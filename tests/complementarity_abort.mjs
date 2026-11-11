// Exercise the actual browser load functions against an unfinished gzip body.
// No browser globals or external packages are required.
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {createServer} from 'node:http';
import vm from 'node:vm';
import {gzipSync} from 'node:zlib';

const viewer = await readFile(new URL('../src/tcren/complementarity/viewer.js', import.meta.url), 'utf8');
const tables = await readFile(new URL('../src/tcren/complementarity/views-browser.js', import.meta.url), 'utf8');
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
        setAttribute() {}, replaceChildren() {}});
    return elements.get(id);
};
const context = vm.createContext({fetch: (url, options) => fetch(new URL(url, base), options),
    Response, DecompressionStream, AbortController, document: {documentElement: {dataset: {}},
        dispatchEvent() {}}, CustomEvent: class {}, $: element, choices() {}, checked: () => new Set(),
    nodeTable() {}, summary() {}, draw() {}, DEFAULT_REGIONS: new Set(),
    location: base, history: {replaceState() {}}, URL,
    shortLabels: new Map([['old', 'old']]), summaryById: new Map(), render() {}});
const loader = viewer.slice(viewer.indexOf('export async function loadScene('),
    viewer.indexOf('// Coordinate ownership')).replace('export ', '');
const render = viewer.slice(viewer.indexOf('async function render()'),
    viewer.indexOf("choices('edges'"));
vm.runInContext(`let generation=0,renderController=null,currentViews=[],nodePage=0;${loader}${render}`, context);
try {
    element('structures').selectedOptions = [{value: 'old-scene'}];
    const oldRender = vm.runInContext('render()', context);
    await until(() => started === 1);
    element('structures').selectedOptions = [];
    await vm.runInContext('render()', context);
    await oldRender;
    await until(() => closed === 1);
    assert.equal(element('status').textContent, 'Select a structure to inspect.');

    const tableLoad = tables.slice(tables.indexOf('async function load()'),
        tables.indexOf('const requested='));
    vm.runInContext(`generation=0;let loadController=null,active=null;
        const entries=[{id:'old',cache:'old-table'}];${tableLoad}`, context);
    element('structure').selectedIndex = 0;
    const oldTable = vm.runInContext('load()', context);
    await until(() => started === 2);
    // The replacement fails promptly; the superseded gzip stream must stop and
    // its AbortError must not replace the latest selection's error status.
    vm.runInContext("entries.push({id:'latest',cache:'latest-table'});", context);
    context.fetch = (url, options) => url === 'latest-table'
        ? Promise.reject(Error('latest selection failure')) : fetch(new URL(url, base), options);
    element('structure').selectedIndex = 1;
    await vm.runInContext('load()', context);
    await oldTable;
    await until(() => closed === 2);
    assert.equal(element('status').textContent, 'Cannot load latest: latest selection failure');
    console.log('Passed: molecular and table selection abort unfinished gzip bodies; stale errors preserve latest status.');
} finally {
    server.closeAllConnections();
    await new Promise(resolve => server.close(resolve));
}
