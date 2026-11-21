// Both projections use the same normalized molecular scene.
const $ = id => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const COLORS = {alpha:'#0072b2', beta:'#56b4e9', peptide:'#d55e00', mhc:'#009e73', b2m:'#999999'};
const DEFAULT_REGIONS = new Set(['CDR1','CDR2','CDR3','PEPTIDE','HELIX_A1','HELIX_A2','HELIX_B1','GROOVE_FLOOR']);
const edgeName = (a,b) => [a,b].sort().join(':');
let entries = [], generation = 0, currentViews = [], nodePage = 0;
let renderController = null;
let vertical = 2;
let focusedContact = null;
let pinned = null, cameraExtent = null, outcomeDomain = null;
let byId = new Map();
const SHAPES = {
    alpha:'M -1 0 A 1 1 0 1 0 1 0 A 1 1 0 1 0 -1 0 Z',
    beta:'M -1 -1 H 1 V 1 H -1 Z',
    peptide:'M 0 -1.2 L 1.15 1 H -1.15 Z',
    mhc:'M 0 -1.3 L 1.2 0 0 1.3 -1.2 0 Z',
    b2m:'M -.33 -1 H .33 V -.33 H 1 V .33 H .33 V 1 H -.33 V .33 H -1 V -.33 H -.33 Z',
};
let groupController = null;
const requested = new URL(location).searchParams;
let groupFiles = null;
async function selectGroup(id) {
    groupController?.abort();
    const controller = groupController = new AbortController();
    document.dispatchEvent(new CustomEvent('tcren:cancel-detail'));
    entries = [];
    pinned = cameraExtent = null;
    $('structures').replaceChildren();
    $('projection').replaceChildren();
    $('status').textContent = 'Loading group…';
    try {
        const response = await fetch(groupFiles.get($('group').value), {signal:controller.signal});
        if(!response.ok)throw Error(`Group HTTP ${response.status}`);
        const payload = await response.json();
        if(controller.signal.aborted)return;
        entries = payload.entries;
        setupOutcomes();
        updateStructureOptions();
        const options = [...$('structures').options];
        const selected = options.find(option=>option.title===id) ?? options[0];
        if(selected)selected.selected=true;
        await render();
    } catch(error) {
        if(!controller.signal.aborted)$('status').textContent=String(error);
    }
}

export function choices(id, values, initiallyVisible) {
    const existing = new Map([...$(id).querySelectorAll('input')].map(n => [n.value,n.checked]));
    $(id).replaceChildren();
    for (const value of [...new Set(values)].sort()) {
        const label = document.createElement('label'), input = document.createElement('input');
        input.type = 'checkbox'; input.value = value;
        input.checked = existing.has(value) ? existing.get(value) : initiallyVisible(value);
        label.append(input,document.createTextNode(' '+value)); $(id).append(label);
    }
}
const checked = id => new Set([...$(id).querySelectorAll('input:checked')].map(n => n.value));

function svgElement(root, tag, attributes = {}, title) {
    const node = document.createElementNS(NS,tag);
    for (const [key,value] of Object.entries(attributes)) node.setAttribute(key,value);
    if (title) {
        const label = document.createElementNS(NS,'title'); label.textContent = title; node.append(label);
    }
    root.append(node); return node;
}

export function matchesOutcome(value, kind, filter) {
    if (kind === 'binary')
        return !filter.discrete || (filter.discrete === 'missing' ? !Number.isFinite(value)
            : value === Number(filter.discrete));
    if (kind === 'continuous')
        return Number.isFinite(value) ? value >= filter.minimum && value <= filter.maximum : filter.unlabelled;
    return true;
}

const formatOutcome = value => Number.isFinite(value) ? String(value) : 'Unlabelled';

function setupOutcomes() {
    byId = new Map(entries.map(entry => [entry.id,entry]));
    const kind = entries.find(entry => entry.outcome_kind)?.outcome_kind;
    $('outcome-controls').hidden = !kind;
    $('discrete-label').hidden = kind !== 'binary';
    $('continuous-label').hidden = kind !== 'continuous';
    $('discrete').value = '';
    const domain = outcomeDomain ?? [0,1];
    for (const [id,value] of [['minimum',domain[0]],['maximum',domain[1]]]) {
        const input = $(id);
        input.min = domain[0]; input.max = domain[1]; input.value = value;
        $(id+'-value').textContent = value;
    }
}

function outcomeColor(entry) {
    if (!Number.isFinite(entry?.outcome)) return '#979da6';
    if (entry.outcome_kind === 'binary') return entry.outcome === 0 ? '#287bd4' : '#ce6418';
    const [low,high] = outcomeDomain ?? [0,1];
    const fraction = high === low ? .5 : (entry.outcome-low)/(high-low);
    return `hsl(${240-213*fraction} 74% 54%)`;
}

export function updateStructureOptions() {
    const select = $('structures');
    const selected = select.value;
    select.replaceChildren();
    const query = $('search').value.toLowerCase(), group = $('group').value;
    const kind = entries.find(entry => entry.outcome_kind)?.outcome_kind;
    const filter = {discrete:$('discrete').value, minimum:Number($('minimum').value),
        maximum:Number($('maximum').value), unlabelled:$('unlabelled').checked};
    for (const entry of entries.filter(e => e.reference === group &&
            (e.id+' '+(e.display_label||'')+' '+(e.search||'')).toLowerCase().includes(query) &&
            matchesOutcome(e.outcome,kind,filter))) {
        const label = entry.display_label || entry.id;
        const option = new Option(kind ? label+' | Label '+formatOutcome(entry.outcome) : label,entry.file);
        option.title = entry.id;
        select.add(option);
    }
    if ([...select.options].some(option => option.value === selected)) select.value = selected;
    else if (select.options.length) select.selectedIndex = 0;
}

export async function loadScene(file, {signal} = {}) {
    const response = await fetch(file, {signal});
    if (!response.ok) throw Error(`${file}: ${response.status}`);
    const text = await new Response(response.body.pipeThrough(new DecompressionStream('gzip'), {signal})).text();
    signal?.throwIfAborted();
    const scene = JSON.parse(text);
    if (scene.schema !== 'tcren.complementarity/2') throw Error('This prototype uses an older schema; regenerate it with the current exporter.');
    return scene;
}

// Coordinate ownership is exclusively atomic; residue positions resolve their CA.
export function nodesFor(scene, atomic) {
    return atomic ? scene.atoms : scene.residues.map(r => ({...r,
        xyz:r.ca_atom === null ? null : scene.atoms[r.ca_atom].xyz}));
}

function viewFor(scene, atomic, visible, groups) {
    const nodes = nodesFor(scene,atomic);
    const residue = node => atomic ? scene.residues[node.residue] : node;
    const shown = node => node.xyz && visible(residue(node));
    const contacts = (atomic ? scene.atom_contacts : scene.residue_contacts).filter(edge =>
        shown(nodes[edge.source]) && shown(nodes[edge.target]) &&
        groups.has(edgeName(residue(nodes[edge.source]).role,residue(nodes[edge.target]).role)));
    const visibleCounts = new Map(), visiblePartners = new Map();
    for (const edge of contacts) {
        for (const [id,partner] of [[edge.source,edge.target],[edge.target,edge.source]]) {
            visibleCounts.set(id,(visibleCounts.get(id)||0)+(atomic?1:edge.contact_count));
            if (!visiblePartners.has(id)) visiblePartners.set(id,new Set());
            visiblePartners.get(id).add(partner);
        }
    }
    return {scene,nodes,residue,shown,contacts,visibleCounts,visiblePartners,atomic};
}

export function residueAtomDetails(scene, residue) {
    const atoms = scene.atoms.filter(atom => atom.residue === residue.id);
    const members = new Set(atoms.map(atom => atom.id));
    const contacts = [];
    for (const edge of scene.atom_contacts) {
        if (!members.has(edge.source) && !members.has(edge.target)) continue;
        const [source, target] = members.has(edge.source)
            ? [edge.source, edge.target] : [edge.target, edge.source];
        const own = scene.atoms[source], partner = scene.atoms[target];
        if (scene.residues[partner.residue].chain !== residue.chain)
            contacts.push({own, partner, distance: edge.distance});
    }
    contacts.sort((a, b) => a.distance - b.distance);
    return {atoms, contacts};
}

function inspect(view,node) {
    const residue = view.residue(node), atom = view.atomic;
    const values={Structure:view.scene.id,Type:atom?'Atom':'Residue',ID:node.id,
        'Source identity':node.key,Chain:residue.chain,Region:residue.region,
        Residue:residue.resname+' ('+residue.aa+')','Author number':residue.number+residue.insertion,
        'Sequence index (0-based)':residue.seq_index,
        'Covalent degree (global)':node.covalent_degree,'Atom-pair contacts (global)':node.contact_count,
        'Contact partner degree (global)':node.contact_degree,'Partner chains (global)':node.partner_chain_count,
        'Atom-pair contacts (after filters)':view.visibleCounts.get(node.id)||0,
        'Contact partner degree (after filters)':view.visiblePartners.get(node.id)?.size||0};
    if(atom)Object.assign(values,{'Atom name':node.name,Element:node.element,
        'Alternate location':node.altloc||'shared','Atom ordinal (0-based)':node.ordinal,
        Backbone:node.backbone,'Partner residues (global)':node.partner_residue_count});
    else Object.assign(values,{'Member atoms':node.atom_count,'Contacting atoms':node.contacting_atom_count,
        'Cα atom ID':node.ca_atom});
    if (focusedContact) {
        focusedContact = null;
        $('projection').querySelector('.focused-contact')?.remove();
    }
    const panel = $('inspector');
    panel.replaceChildren();
    const title = document.createElement('h3');
    title.textContent = view.scene.id + ' | ' + residue.chain + ' ' + residue.resname + ' ' + residue.number + residue.insertion;
    panel.append(title);
    const {atoms, contacts} = residueAtomDetails(view.scene, residue);
    const atomName = item => item.name + (item.altloc ? ' (' + item.altloc + ')' : '');
    const members = document.createElement('p');
    members.textContent = atoms.map(atomName).join(', ') || 'No exported heavy atoms.';
    panel.append(members);
    const rows = document.createElement('div');
    rows.className = 'atom-contacts';
    rows.textContent = contacts.length ? '' : 'No recorded interchain atom contacts.';
    for (const contact of contacts) {
        const partner = view.scene.residues[contact.partner.residue];
        const button = document.createElement('button');
        button.type = 'button';
        button.setAttribute('aria-pressed', 'false');
        button.textContent = atomName(contact.own) + ' → ' + partner.chain + ' ' +
            partner.aa + partner.number + partner.insertion + ' ' + atomName(contact.partner) +
            ' · ' + contact.distance.toFixed(2) + ' Å';
        button.onclick = () => {
            rows.querySelector('[aria-pressed="true"]')?.setAttribute('aria-pressed', 'false');
            focusedContact = focusedContact === contact ? null : contact;
            button.setAttribute('aria-pressed', String(Boolean(focusedContact)));
            redraw();
        };
        rows.append(button);
    }
    const details = document.createElement('details');
    const caption = document.createElement('summary');
    caption.textContent = 'Identity and contact counts';
    const table = document.createElement('table');
    for(const [key,value]of Object.entries(values)){
        const row=table.insertRow();row.insertCell().textContent=key;row.insertCell().textContent=String(value);
    }
    details.append(caption, table);
    panel.append(rows, details);
}

function summary(views,atomic) {
    $('summary').replaceChildren();
    const table = document.createElement('table');
    const header = table.insertRow();
    for (const label of ['Structure','Residues','Atoms','Covalent bonds','Atom contacts','Residue contacts','TCR–peptide atom contacts','Nodes after filters','Atom contacts after filters']) {
        const th = document.createElement('th'); th.textContent=label; header.append(th);
    }
    for (const view of views) {
        const s=view.scene, row=table.insertRow();
        const values=[s.id,s.summary.residues,s.summary.atoms,s.summary.covalent_bonds,
            s.summary.atom_contacts,s.summary.residue_contacts,s.summary.tcr_peptide_atom_contacts,view.nodes.filter(view.shown).length,
            view.contacts.reduce((sum,e)=>sum+(atomic?1:e.contact_count),0)];
        for (const [index,value] of values.entries()){
            const cell=row.insertCell();cell.textContent=index===0&&String(value).length>32?String(value).slice(0,12)+'…':String(value);
            if(index===0)cell.title=s.id;
        }
        if(s.summary.no_tcr_peptide_interface)row.title='No TCR–peptide atom contacts at this cutoff.';
    }
    $('summary').append(table);
}

function nodeTable() {
    const query=$('nodesearch').value.toLowerCase();
    const rows=[];
    for(const view of currentViews)for(const node of view.nodes){
        const r=view.residue(node);
        if($('visibleonly').checked&&!view.shown(node))continue;
        if(!`${view.scene.id} ${node.id} ${node.name||''} ${r.chain} ${r.region} ${r.resname} ${r.number}${r.insertion}`.toLowerCase().includes(query))continue;
        rows.push({view,node,residue:r});
    }
    nodePage=Math.min(nodePage,Math.max(0,Math.ceil(rows.length/50)-1));
    $('nodetable').replaceChildren();const table=document.createElement('table'),header=table.insertRow();
    for(const name of ['Structure','ID','Chain / region','Name','Author number','Atoms','Global contacts','Partner degree','Covalent degree']){
        const th=document.createElement('th');th.textContent=name;header.append(th);
    }
    for(const {view,node,residue:r}of rows.slice(nodePage*50,(nodePage+1)*50)){
        const row=table.insertRow();row.tabIndex=0;
        for(const value of [view.scene.id,node.id,r.chain+' / '+r.region,node.name||r.resname,
                r.number+r.insertion,node.atom_count??1,node.contact_count,node.contact_degree,node.covalent_degree])row.insertCell().textContent=String(value);
        row.onclick=()=>inspect(view,node);row.onkeydown=e=>{if(e.key==='Enter')inspect(view,node);};
    }
    $('nodetable').append(table);$('pageinfo').textContent=`${rows.length} nodes · page ${nodePage+1} / ${Math.max(1,Math.ceil(rows.length/50))}`;
    $('previous').disabled=nodePage===0;$('next').disabled=(nodePage+1)*50>=rows.length;
}

function draw(svg,vertical,views,atomic,scale) {
    svg.replaceChildren();
    const labelBoxes=[],labelBuckets=[[],[],[]],labelLimit=250;let requestedLabels=0;
    function queueLabel(priority,...label){
        const [cx,cy]=label;if(cx<5||cx>695||cy<5||cy>510)return;
        requestedLabels++;
        if(labelBuckets[priority].length<labelLimit)labelBuckets[priority].push(label);
    }
    function labelAt(cx,cy,text,color,title) {
        if(cx<5||cx>695||cy<5||cy>510)return;
        const width=text.length*8+4,height=16;
        let box;
        for(const dy of [-6,-22,12,-38,28,-54,44,-70,60,-86,76]) {
            for(const dx of [6,-width-6,22,-width-22]) {
                const candidate={x:cx+dx,y:cy+dy-height,w:width,h:height};
                if(candidate.x<5||candidate.x+width>695||candidate.y<5||candidate.y+height>515)continue;
                if(labelBoxes.some(b=>candidate.x<b.x+b.w+2&&candidate.x+width+2>b.x&&candidate.y<b.y+b.h+2&&candidate.y+height+2>b.y))continue;
                box=candidate;break;
            }
            if(box)break;
        }
        box??={x:cx+6,y:cy-22,w:width,h:height};
        labelBoxes.push(box);
        if(Math.abs(box.y+height-cy)>12||Math.abs(box.x-cx)>width+8)
            svgElement(svg,'line',{x1:cx,y1:cy,x2:Math.max(box.x,Math.min(cx,box.x+width)),y2:box.y+height/2,stroke:color,'stroke-width':.5,opacity:.65});
        svgElement(svg,'text',{x:box.x,y:box.y+height-3,fill:'var(--map-ink)','font-size':13,
            'paint-order':'stroke',stroke:'var(--map-paper)','stroke-width':2},title).textContent=text;
    }
    const position=p=>[350+(p[0]-Number($('panx').value))*scale,
        280-(p[vertical]-Number($(vertical===2?'panz':'pany').value))*scale];
    const [zeroX,zeroY]=position([0,0,0]);
    for (const axis of [0,1]) svgElement(svg,'line',{x1:axis?zeroX:25,y1:axis?25:zeroY,
        x2:axis?zeroX:675,y2:axis?535:zeroY,stroke:'var(--map-grid)'});
    svgElement(svg,'text',{x:660,y:Math.min(540,Math.max(20,zeroY+18)),'font-size':13}).textContent='+x';
    svgElement(svg,'text',{x:Math.min(650,Math.max(10,zeroX+8)),y:30,'font-size':13}).textContent=vertical===2?'+z':'+y';
    const maximum=130/scale,power=10**Math.floor(Math.log10(maximum));
    const bar=[5,2,1].find(value=>value*power<=maximum)*power;
    svgElement(svg,'line',{x1:30,y1:525,x2:30+bar*scale,y2:525,stroke:'var(--map-ink)','stroke-width':2});
    svgElement(svg,'text',{x:30,y:545,'font-size':12}).textContent=`${Number(bar.toPrecision(2))} Å`;
    const hiddenAxis=3-vertical;
    let near=-Infinity,far=Infinity;
    if($('depth').checked)for(const view of views)for(const node of view.nodes){
        if(!view.shown(node))continue;
        near=Math.max(near,node.xyz[hiddenAxis]);
        far=Math.min(far,node.xyz[hiddenAxis]);
    }
    for (const view of [...views].reverse()) {
        const {scene,nodes,residue,shown,contacts}=view;
        const opacity=pinned?.id===scene.id && views.length>1?.46:1;
        const fill=outcomeColor(byId.get(scene.id));
        function line(a,b,attributes,title) {
            if (!shown(a)||!shown(b)) return;
            const [x1,y1]=position(a.xyz),[x2,y2]=position(b.xyz);
            svgElement(svg,'line',{x1,y1,x2,y2,...attributes},title);
        }
        for (const edge of contacts) {
            const count=atomic?1:edge.contact_count;
            const copies=!atomic&&$('weight').value==='parallel'?count:1;
            const title=`${scene.id}: ${count} atom-pair contacts${atomic?' · '+edge.distance.toFixed(2)+' Å':''}`;
            const attributes={stroke:'#666','stroke-width':!atomic&&$('weight').value==='weighted'?count*.5:.7,
                'stroke-dasharray':'3 3',opacity,fill:'none'};
            for (let k=0;k<copies;k++) {
                if(copies===1) {line(nodes[edge.source],nodes[edge.target],attributes,title);continue;}
                const [x1,y1]=position(nodes[edge.source].xyz),[x2,y2]=position(nodes[edge.target].xyz);
                const length=Math.hypot(x2-x1,y2-y1)||1,offset=(k-(copies-1)/2)*2;
                const cx=(x1+x2)/2-(y2-y1)/length*offset,cy=(y1+y2)/2+(x2-x1)/length*offset;
                svgElement(svg,'path',{d:`M ${x1} ${y1} Q ${cx} ${cy} ${x2} ${y2}`,...attributes},title);
            }
        }
        if (atomic&&$('covalent').checked) for (const edge of scene.bonds) line(nodes[edge.source],nodes[edge.target],{
            stroke:'#aaa','stroke-width':.8,opacity},`${edge.kind}; ${edge.method}; ${edge.distance.toFixed(2)} Å`);
        if ($('backbone').checked) {
            const links=atomic?scene.bonds.filter(e=>e.kind==='backbone').map(e=>[e.source,e.target]):scene.trace;
            for (const [a,b] of links) line(nodes[a],nodes[b],{
                stroke:'var(--map-trace)','stroke-width':2,opacity},'Backbone');
        }
        for (const node of nodes.filter(shown)) {
            const r=residue(node),[cx,cy]=position(node.xyz);
            const name=`${scene.id} | ${r.chain} ${r.region} ${r.aa}${r.number}${r.insertion}${atomic?' '+node.name+' '+node.altloc:''}\nID ${node.id}; ${node.contact_count} global atom contacts; degree ${node.contact_degree}; covalent degree ${node.covalent_degree}\n${view.visibleCounts.get(node.id)||0} atom contacts after filters`;
            const radius=atomic?2.5:4.7;
            const point=svgElement(svg,'path',{d:SHAPES[r.role]??SHAPES.alpha,
                transform:`translate(${cx} ${cy}) scale(${radius})`,fill,
                'vector-effect':'non-scaling-stroke',
                stroke:$('outlines').checked?'var(--map-outline, #333)':'none','stroke-width':1,
                opacity:opacity*(near>far?.22+.78*(node.xyz[hiddenAxis]-far)/(near-far):1),tabindex:0,
                'data-structure-id':scene.id,'data-residue-id':r.id},name);
            point.onclick=()=>inspect(view,node);
            point.onkeydown=e=>{if(e.key==='Enter'||e.key===' ')inspect(view,node);};
            const scope=$('labelscope').value;
            const peptideLabel=r.role==='peptide'&&(!atomic||node.name==='CA');
            const label=scope==='peptide'?peptideLabel&&view===views[0]:scope==='allpeptides'?peptideLabel:
                scope==='all'||(scope==='contacts'&&node.contact_count>0&&(!atomic||node.name==='CA'))||peptideLabel;
            if ($('labels').checked&&label) queueLabel(r.role==='peptide'?0:view.visibleCounts.get(node.id)>0?1:2,cx,cy,
                atomic&&scope==='all'?node.name:`${r.aa}${r.number}${r.insertion}${atomic?' CA':''}`,COLORS[r.role],name);
        }
    }
    for(const bucket of labelBuckets)for(const label of bucket){
        if(labelBoxes.length===labelLimit)break;
        labelAt(...label);
    }
    if (focusedContact) {
        const group = svgElement(svg, 'g', {class:'focused-contact', 'pointer-events':'none'});
        const [x1,y1] = position(focusedContact.own.xyz);
        const [x2,y2] = position(focusedContact.partner.xyz);
        svgElement(group,'line',{x1,y1,x2,y2,stroke:'var(--contact-focus)', 'stroke-width':3});
        for (const [cx,cy] of [[x1,y1],[x2,y2]])
            svgElement(group,'circle',{cx,cy,r:6,fill:'var(--contact-focus)',stroke:'var(--map-outline, #333)'});
    }
    if(requestedLabels>labelBoxes.length)svgElement(svg,'text',{
        x:350,y:556,'font-size':11,'text-anchor':'middle',class:'label-limit',fill:'var(--map-ink)'
    }).textContent=`${requestedLabels-labelBoxes.length} labels omitted (limit ${labelLimit}); filter chains/regions. All points remain inspectable.`;
}

export function fitCameraCoordinates(positions) {
    const lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];let extent=10,count=0;
    for(const p of positions){count++;for(let axis=0;axis<3;axis++){
        lo[axis]=Math.min(lo[axis],p[axis]);hi[axis]=Math.max(hi[axis],p[axis]);extent=Math.max(extent,Math.abs(p[axis]));
    }}
    if(!count)return false;
    for(const [axis,id]of [[0,'panx'],[1,'pany'],[2,'panz']])$(id).value=String(Math.max(-100,Math.min(100,Math.round((lo[axis]+hi[axis])/2))));
    const radius=Math.max(1,...lo.map((value,i)=>(hi[i]-value)/2))*1.1;
    $('zoom').value=String(Math.min(5,Math.max(1,Math.round((cameraExtent??extent)/radius*4)/4)));return true;
}

export function resetCameraCoordinates(){for(const id of ['panx','pany','panz'])$(id).value='0';$('zoom').value='1';}

function redraw() {
    let extent=10;
    for(const view of currentViews)for(const node of view.nodes){
        if(view.shown(node))extent=Math.max(extent,...node.xyz.map(Math.abs));
    }
    if (currentViews.length) cameraExtent ??= extent;
    const scale=230/(cameraExtent??extent)*Number($('zoom').value);
    draw($('projection'),vertical,currentViews,currentViews[0]?.atomic??false,scale);
}

async function render() {
    const token=++generation;
    renderController?.abort();
    const controller = renderController = new AbortController();
    $('projection').setAttribute('aria-busy','true');
    $('status').textContent='Loading selection…';
    try {
        const selected=[...$('structures').selectedOptions];
        const candidate=selected[0];
        if (pinned && pinned.file !== candidate?.value) selected.push({value:pinned.file,title:pinned.id});
        const loaded = new Map(currentViews.map(view => [view.scene.id,view.scene]));
        const scenes=await Promise.all(selected.map(o=>loaded.get(o.title) ??
            loadScene(o.value,{signal:controller.signal})));
        if(token!==generation)return;
        choices('chains',scenes.flatMap(s=>s.components.map(c=>c.role+':'+c.chain)),()=>true);
        choices('regions',scenes.flatMap(s=>s.residues.map(r=>r.role+':'+r.region)),v=>DEFAULT_REGIONS.has(v.split(':')[1]));
        const regions=checked('regions'),chains=checked('chains'),groups=checked('edges');
        const visible=r=>regions.has(r.role+':'+r.region)&&chains.has(r.role+':'+r.chain);
        const atomic=$('detail').value==='atom';
        const views=scenes.map(s=>viewFor(s,atomic,visible,groups));
        currentViews=views;nodePage=0;nodeTable();
        const outcomeKey = $('outcome-key');
        outcomeKey.replaceChildren();
        for (const [index,scene] of scenes.entries()) {
            const entry=byId.get(scene.id);
            const item=document.createElement('span');
            item.style.borderInlineStartColor=outcomeColor(entry);
            item.textContent=(index===0 && candidate?'Active: ':'Reference: ')+scene.id+' | Label '+formatOutcome(entry?.outcome);
            outcomeKey.append(item);
        }
        focusedContact = null;
        $('inspector').textContent='Select a residue to see its atoms and interchain contacts.';
        summary(views,atomic);
        redraw();
        const url=new URL(location);
        if(candidate || url.searchParams.has('id')) {
            if($('group').value)url.searchParams.set('group',$('group').value);
            if(candidate)url.searchParams.set('id',scenes[0].id);
            else url.searchParams.delete('id');
            history.replaceState(null,'',url);
            document.dispatchEvent(new Event('tcren:model'));
        }
        document.dispatchEvent(new CustomEvent('tcren:render',{detail:{scenes,atomic}}));
        $('status').textContent=scenes.length?scenes.map(s=>`${s.id}\nReference ${s.frame.reference}; fixed origin ${s.frame.origin_component}: ${s.frame.origin}; ${s.frame.anchors} anchors; RMSD ${s.frame.rmsd.toFixed(3)} Å\n${s.warnings.length?s.warnings.join('\n'):'No export warnings'}`).join('\n\n')+`\nCamera offsets x/y/z: ${$('panx').value}/${$('pany').value}/${$('panz').value} Å; zoom ${$('zoom').value}. Molecular reference coordinates remain fixed.`:'Select a structure to inspect.';
    } catch(error) {
        controller.abort(); // Also stop sibling loads when one selection member fails.
        if(token!==generation)return;
        currentViews=[];nodeTable();summary([],false);
        $('projection').replaceChildren();
        $('inspector').textContent='Selection unavailable.';$('status').textContent=String(error);
    } finally {
        if(renderController===controller)renderController=null;
        if(token===generation)$('projection').setAttribute('aria-busy','false');
    }
}

choices('edges',['alpha:peptide','beta:peptide','alpha:mhc','beta:mhc','mhc:peptide','alpha:beta','b2m:mhc','mhc:mhc'],v=>!['alpha:beta','b2m:mhc','mhc:mhc'].includes(v));
document.addEventListener('tcren:cancel-detail',()=>{
    generation++;
    renderController?.abort();renderController=null;
    $('projection').setAttribute('aria-busy','false');
});
async function readModelIndex() {
    // An ID-only bookmark needs the legacy index to discover its group.
    if(requested.has('group') || !requested.has('id')){
        const response=await fetch('groups.json');
        if(response.status!==404)return response;
    }
    return fetch('index.json');
}

async function loadDownloads() {
    const exported=await fetch('tables/manifest.json');
    if(!exported.ok)throw Error(`Table index HTTP ${exported.status}`);
    for(const item of (await exported.json()).files){
        if(!/^[a-z_]+\/part-[0-9]+\.parquet$/.test(item.file))continue;
        const row=document.createElement('div');
        for(const filename of [item.file,item.tsv_file].filter(Boolean)){
            if(!/^[a-z_]+\/part-[0-9]+\.(parquet|tsv\.gz)$/.test(filename))continue;
            const link=document.createElement('a');link.href='tables/'+filename;link.download='';
            link.textContent=`${filename} (${item.rows} rows) `;row.append(link);
        }
        $('downloads').append(row);
    }
}

const downloads=$('downloads');
downloads?.closest('details')?.addEventListener('toggle',async function loadWhenOpened() {
    if(!this.open)return;
    this.removeEventListener('toggle',loadWhenOpened);
    try {
        await loadDownloads();
    } catch(error) {
        downloads.textContent=String(error);
        this.addEventListener('toggle',loadWhenOpened);
    }
});

try {
    const response=await readModelIndex();
    if(!response.ok)throw Error(`Index ${response.status}`);
    const index=await response.json();
    const groups=Array.isArray(index)
        ? [...new Map(index.map(entry=>[entry.reference,{id:entry.reference,label:entry.label}])).values()]
        : index.groups;
    if(Array.isArray(index)){
        entries=index;
        for (const entry of entries) {
            if (!Number.isFinite(entry.outcome)) continue;
            outcomeDomain ??= [entry.outcome,entry.outcome];
            outcomeDomain[0]=Math.min(outcomeDomain[0],entry.outcome);
            outcomeDomain[1]=Math.max(outcomeDomain[1],entry.outcome);
        }
        setupOutcomes();
    }
    else {
        outcomeDomain=index.outcome_domain??null;
        groupFiles=new Map(groups.map(group=>[group.id,group.file]));
    }
    $('population').textContent=`${index.count ?? entries.length} structures · ${groups.length} compatible reference groups · independent local dashboard`;
    const diagnostics=await fetch('errors.json');
    $('errors').textContent=diagnostics.ok?JSON.stringify(await diagnostics.json(),null,2):'No command-level diagnostics supplied.';
    for(const {id:group,label:display} of groups){
        const label=display||group;
        const text=/^[a-f0-9]{32,}$/i.test(group)?label.replace(group,group.slice(0,8)):label;
        const option=new Option(text,group);option.title=`${label}\nReference ${group}`;
        $('group').add(option);
    }
    if(groups.some(group=>group.id===requested.get('group')))$('group').value=requested.get('group');
    if(groupFiles && groups.length)await selectGroup(requested.get('id'));
    else {
        const match=entries.find(entry=>entry.id===requested.get('id'));
        if(match)$('group').value=match.reference;
        updateStructureOptions();
        if(match)$('structures').value=match.file;
        await render();
    }
} catch(error) {$('status').textContent=`Serve this directory over localhost. ${error}`;}
function updatePin() {
    $('pin-model').textContent=pinned?'Clear reference':'Pin reference';
    $('reference-model').textContent=pinned?'Pinned reference: '+pinned.id:'';
}
function applyFilters() {
    updateStructureOptions();
    render();
}
$('group').onchange=()=>{
    pinned = cameraExtent = null;
    updatePin();
    if(groupFiles)return selectGroup();
    setupOutcomes();
    applyFilters();
};
$('search').oninput=applyFilters;
$('discrete').onchange=$('unlabelled').onchange=applyFilters;
for(const id of ['minimum','maximum']){
    $(id).oninput=()=>{
        const min=$('minimum'),max=$('maximum');
        if(Number(min.value)>Number(max.value)) {
            if(id==='minimum')min.value=max.value;
            else max.value=min.value;
        }
        $(id+'-value').textContent=$(id).value;
    };
    $(id).onchange=applyFilters;
}
$('pin-model').onclick=()=>{
    const option=$('structures').selectedOptions[0];
    if(!pinned && !option)return;
    pinned=pinned?null:{file:option.value,id:option.title};
    updatePin();
    render();
};
function stepModel(direction) {
    const select=$('structures');
    const next=select.selectedIndex+direction;
    if(next<0 || next>=select.options.length)return;
    select.selectedIndex=next;
    render();
}
$('previous-model').onclick=()=>stepModel(-1);
$('next-model').onclick=()=>stepModel(1);
document.addEventListener('keydown',event=>{
    if(event.defaultPrevented || event.altKey || event.metaKey || event.ctrlKey ||
        event.target.closest?.('input, textarea, select, [contenteditable]'))return;
    if(event.key==='[' || event.key===']'){
        event.preventDefault();
        stepModel(event.key==='['?-1:1);
    }
});
for(const id of ['structures','detail','weight','labels','backbone','covalent','regions','chains','edges','labelscope'])$(id).addEventListener('change',render);
$('fit-camera').onclick=()=>{
    function* positions(){for(const view of currentViews)for(const node of view.nodes)if(view.shown(node))yield node.xyz;}
    if(fitCameraCoordinates(positions()))redraw();
};
$('reset-camera').onclick=()=>{resetCameraCoordinates();redraw();};
$('nodesearch').oninput=()=>{nodePage=0;nodeTable();};
$('visibleonly').onchange=()=>{nodePage=0;nodeTable();};
$('previous').onclick=()=>{nodePage--;nodeTable();};$('next').onclick=()=>{nodePage++;nodeTable();};

for(const id of ['zoom','panx','pany','panz','outlines','depth'])$(id).addEventListener('input',redraw);
for(const [id,axis] of [['side-view',2],['top-view',1]])$(id).onclick=()=>{
    vertical=axis;
    $('side-view').setAttribute('aria-pressed',String(axis===2));
    $('top-view').setAttribute('aria-pressed',String(axis===1));
    $('projection').setAttribute('aria-label',axis===2?'Side contact view':'Top contact view');
    $('depth-label').textContent=axis===2?'Fade by y depth':'Fade by z depth';
    redraw();
};
function stepZoom(delta) {
    const control=$('zoom');
    const previous=control.value;
    control.value=String(Number(previous)+delta);
    if(control.value===previous)return false;
    redraw();
    return true;
}
$('zoom-out').onclick=()=>stepZoom(-.25);
$('zoom-in').onclick=()=>stepZoom(.25);
$('projection').addEventListener('wheel',event=>{
    if(!event.deltaY||event.ctrlKey||event.metaKey)return;
    if(stepZoom(event.deltaY<0?.25:-.25))event.preventDefault();
},{passive:false});
