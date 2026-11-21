// Exercise the actual browser load functions against an unfinished gzip body.
// No browser globals or external packages are required.
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {createServer} from 'node:http';
import vm from 'node:vm';
import {gzipSync} from 'node:zlib';

for (const file of ['dashboard.html', 'views.html']) {
    const page = await readFile(new URL('../src/tcren/complementarity/' + file, import.meta.url), 'utf8');
    assert.match(page, /<link rel="stylesheet" href="dashboard\.css">/);
    assert.doesNotMatch(page, /<style[ >]/);
    assert.match(page, /<meta name="viewport" content="width=device-width">/);
}
const projection = await readFile(new URL('../src/tcren/complementarity/dashboard.html', import.meta.url), 'utf8');
assert.equal((projection.match(/<svg /g) ?? []).length, 1);
assert.match(projection, /id="outlines" type="checkbox">/);
assert.match(projection, /id="depth" type="checkbox">/);
await readFile(new URL('../src/tcren/complementarity/dashboard.css', import.meta.url), 'utf8');
const viewer = await readFile(new URL('../src/tcren/complementarity/viewer.js', import.meta.url), 'utf8');
const tables = await readFile(new URL('../src/tcren/complementarity/views-browser.js', import.meta.url), 'utf8');
const gridHelpers = tables.slice(tables.indexOf('function residueAxis('),
    tables.indexOf('function heatmapGrid('));
const {residueAxis, heatmapState} = vm.runInThisContext(
    `(() => {${gridHelpers}; return {residueAxis,heatmapState};})()`);
assert.deepEqual(residueAxis([
    {tcr_index:3,tcr_chain:'B',tcr_aa:'E',tcr_number:92,tcr_insertion:'A'},
    {tcr_index:0,tcr_chain:'B',tcr_aa:'G',tcr_number:90,tcr_insertion:''},
], 'tcr', 'tcr_index'), [[0,'B:G90'], [3,'B:E92A']]);
assert.equal(heatmapState({count:0,contacting:false},'count',false),'value');
assert.equal(heatmapState({count:0,contacting:false},'count',true),'masked');
assert.equal(heatmapState({count:null,contacting:true},'count',true),'missing');
assert.equal(heatmapState(undefined,'count',true),'missing');
let started = 0, closed = 0;
const server = createServer((req, res) => {
    started++;
    res.writeHead(200);
    // Send the gzip header but leave the stream unfinished until cancellation.
    res.write(gzipSync('{}').subarray(0, 10));
    res.on('close', () => { closed++; });
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const base = `http://127.0.0.1:${server.address().port}/`;
const until = async predicate => {
    for (let i = 0; i < 200; i++) {
        if (predicate()) return;
        await new Promise(resolve => setTimeout(resolve, 5));
    }
    throw Error('Timed out waiting for stream cancellation');
};
const elements = new Map();
const element = id => {
    if (!elements.has(id)) elements.set(id, {selectedOptions: [], value: '1', textContent: '',
        setAttribute() {}, replaceChildren() {}, append() {}});
    return elements.get(id);
};
const context = vm.createContext({fetch: (url, options) => fetch(new URL(url, base), options),
    Response, DecompressionStream, AbortController, document: {documentElement: {dataset: {}},
        dispatchEvent() {}, createElement: () => ({style:{}})}, CustomEvent: class {}, $: element, choices() {}, checked: () => new Set(),
    nodeTable() {}, summary() {}, draw() {}, viewFor: scene => ({scene, nodes:[], shown:() => true}), Event, DEFAULT_REGIONS: new Set(),
    location: base, history: {replaceState() {}}, URL,
    shortLabels: new Map([['old', 'old']]), summaryById: new Map(), render() {}});
const atomDetails = viewer.slice(viewer.indexOf('export function residueAtomDetails('),
    viewer.indexOf('function inspect(')).replace('export ', '');
const inspectResidue = vm.runInContext(atomDetails + '\nresidueAtomDetails', context);
const residueScene = {
    residues: [{id:0, chain:'A'}, {id:1, chain:'A'}, {id:2, chain:'B'}],
    atoms: [{id:0, residue:0}, {id:1, residue:0}, {id:2, residue:1}, {id:3, residue:2}],
    atom_contacts: [{source:3, target:0, distance:3}, {source:1, target:3, distance:2},
        {source:1, target:2, distance:1}],
};
const inspected = inspectResidue(residueScene, residueScene.residues[0]);
assert.deepEqual(Array.from(inspected.atoms, atom => atom.id), [0, 1]);
assert.deepEqual(Array.from(inspected.contacts, edge => [edge.own.id, edge.partner.id, edge.distance]),
    [[1, 3, 2], [0, 3, 3]]);
assert.equal(inspectResidue(residueScene, residueScene.residues[1]).contacts.length, 0);
const outcomes = viewer.slice(viewer.indexOf('export function matchesOutcome('),
    viewer.indexOf('const formatOutcome')).replace('export ', '');
const matchesOutcome = vm.runInContext(outcomes + '\nmatchesOutcome', context);
assert.equal(matchesOutcome(0, 'binary', {discrete:'0'}), true);
assert.equal(matchesOutcome(1, 'binary', {discrete:'0'}), false);
assert.equal(matchesOutcome(undefined, 'binary', {discrete:'missing'}), true);
const range = {minimum:-2.5, maximum:1.25, unlabelled:false};
assert.equal(matchesOutcome(-2.5, 'continuous', range), true);
assert.equal(matchesOutcome(1.25, 'continuous', range), true);
assert.equal(matchesOutcome(1.26, 'continuous', range), false);
assert.equal(matchesOutcome(undefined, 'continuous', range), false);
const loader = viewer.slice(viewer.indexOf('export async function loadScene('),
    viewer.indexOf('// Coordinate ownership')).replace('export ', '');
const render = viewer.slice(viewer.indexOf('function redraw()'),
    viewer.indexOf("choices('edges'"));
vm.runInContext(`const byId=new Map(),outcomeColor=()=>"#ccc",formatOutcome=()=>"Unlabelled";let pinned=null,cameraExtent=null,focusedContact=null,vertical=2,generation=0,renderController=null,groupFiles=null,currentViews=[],nodePage=0;${loader}${render}`, context);
try {
    element('structures').selectedOptions = [{value: 'old-scene'}];
    const oldRender = vm.runInContext('render()', context);
    await until(() => started === 1);
    element('structures').selectedOptions = [];
    await vm.runInContext('render()', context);
    await oldRender;
    await until(() => closed === 1);
    assert.equal(element('status').textContent, 'Select a structure to inspect.');

    const calls = [];
    const modelIndex = viewer.slice(viewer.indexOf('async function readModelIndex()'),
        viewer.indexOf('async function loadDownloads()'));
    for (const [query, status, expected] of [
        ['', 200, ['groups.json']], ['?group=second&id=model-2', 200, ['groups.json']],
        ['?id=model-2', 200, ['index.json']], ['', 404, ['groups.json', 'index.json']],
    ]) {
        const requested = new URL(base + query).searchParams;
        const fetch = async path => {
            calls.push(path);
            return {status: path === 'groups.json' ? status : 200};
        };
        await vm.runInNewContext(`${modelIndex};readModelIndex()`, {requested, fetch});
        assert.deepEqual(calls.splice(0), expected);
    }

    // Successful rendering of a legacy selection must publish its resolved group
    // before the navigation controls synchronize their cross-view links.
    const observed = [];
    context.location = base + '?id=model-2';
    context.history.replaceState = (_state, _title, url) => observed.push(url.search);
    context.fetch = async () => new Response(gzipSync(JSON.stringify({
        schema:'tcren.complementarity/2', id:'model-2', components:[], residues:[], warnings:[],
        frame:{reference:'second', origin_component:'peptide', origin:'apex', anchors:3, rmsd:0},
    })));
    element('group').value = 'second';
    element('structures').selectedOptions = [{value:'model-2.gz'}];
    await vm.runInContext('render()', context);
    assert.deepEqual(observed, ['?id=model-2&group=second']);
    context.fetch = (url, options) => fetch(new URL(url, base), options);

    const tableLoad = tables.slice(tables.indexOf('async function load()'),
        tables.indexOf('const choices = grouped'));
    vm.runInContext(`generation=0;let loadController=null,active=null;
        let candidateScene=null,motifFocus=null;const grouped=false;
        function updateCandidateMotif(){}
        const visibleEntries=[{id:'old',cache:'old-table'}];${tableLoad}`, context);
    element('structure').selectedIndex = 0;
    const oldTable = vm.runInContext('load()', context);
    await until(() => started === 2);
    // The replacement fails promptly; the superseded gzip stream must stop and
    // its AbortError must not replace the latest selection's error status.
    vm.runInContext("visibleEntries.push({id:'latest',cache:'latest-table'});", context);
    context.fetch = (url, options) => url === 'latest-table'
        ? Promise.reject(Error('latest selection failure')) : fetch(new URL(url, base), options);
    element('structure').selectedIndex = 1;
    await vm.runInContext('load()', context);
    await oldTable;
    await until(() => closed === 2);
    assert.equal(element('status').textContent, 'Cannot load latest: latest selection failure');
    context.fetch = async () => new Response(gzipSync(JSON.stringify({
        tables:{contact_heatmap:[]}, sources:{},
    })));
    vm.runInContext("visibleEntries.push({id:'model-2',reference:'second',cache:'model-2.gz'})", context);
    element('structure').selectedIndex = 2;
    await vm.runInContext('load()', context);
    assert.deepEqual(observed, ['?id=model-2&group=second', '?id=model-2&group=second']);
    console.log('Passed: molecular and table selection abort unfinished gzip bodies; stale errors preserve latest status.');
} finally {
    server.closeAllConnections();
    await new Promise(resolve => server.close(resolve));
}
