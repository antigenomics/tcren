// Rendering derives both views from the single normalized molecular scene.
const $ = id => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const COLORS = {alpha:'#0072b2', beta:'#56b4e9', peptide:'#d55e00', mhc:'#009e73', b2m:'#999999'};
const DEFAULT_REGIONS = new Set(['CDR1','CDR2','CDR3','PEPTIDE','HELIX_A1','HELIX_A2','HELIX_B1','GROOVE_FLOOR']);
const edgeName = (a,b) => [a,b].sort().join(':');
let entries = [], generation = 0, currentViews = [], nodePage = 0;

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

export function updateStructureOptions() {
    const selected = new Set([...$('structures').selectedOptions].map(o => o.value));
    $('structures').replaceChildren();
    const query = $('search').value.toLowerCase();
    for (const entry of entries.filter(e => e.reference === $('group').value &&
            (e.id+' '+(e.search||'')).toLowerCase().includes(query))) {
        const option = new Option(entry.display_label||entry.id,entry.file);
        option.title=entry.id;option.selected = selected.has(entry.file);
        $('structures').add(option);
    }
}

export async function loadScene(file) {
    const response = await fetch(file);
    if (!response.ok) throw Error(`${file}: ${response.status}`);
    const scene = JSON.parse(await new Response(response.body.pipeThrough(new DecompressionStream('gzip'))).text());
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
    $('inspector').replaceChildren();const table=document.createElement('table');
    for(const [key,value]of Object.entries(values)){
        const row=table.insertRow();row.insertCell().textContent=key;row.insertCell().textContent=String(value);
    }
    $('inspector').append(table);
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
        svgElement(svg,'text',{x:box.x,y:box.y+height-3,fill:color,'font-size':13,
            'paint-order':'stroke',stroke:'white','stroke-width':2},title).textContent=text;
    }
    const position=p=>[350+(p[0]-Number($('panx').value))*scale,
        280-(p[vertical]-Number($(vertical===2?'panz':'pany').value))*scale];
    const [zeroX,zeroY]=position([0,0,0]);
    for (const axis of [0,1]) svgElement(svg,'line',{x1:axis?zeroX:25,y1:axis?25:zeroY,
        x2:axis?zeroX:675,y2:axis?535:zeroY,stroke:'#d6dce3'});
    svgElement(svg,'text',{x:660,y:Math.min(540,Math.max(20,zeroY+18)),'font-size':13}).textContent='+x';
    svgElement(svg,'text',{x:Math.min(650,Math.max(10,zeroX+8)),y:30,'font-size':13}).textContent=vertical===2?'+z':'+y';
    const maximum=130/scale,power=10**Math.floor(Math.log10(maximum));
    const bar=[5,2,1].find(value=>value*power<=maximum)*power;
    svgElement(svg,'line',{x1:30,y1:525,x2:30+bar*scale,y2:525,stroke:'#333','stroke-width':2});
    svgElement(svg,'text',{x:30,y:545,'font-size':12}).textContent=`${Number(bar.toPrecision(2))} Å`;
    for (const view of views) {
        const {scene,nodes,residue,shown,contacts}=view, opacity=views.length>1?.55:1;
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
                stroke:COLORS[residue(nodes[a]).role],'stroke-width':2,opacity},'Backbone');
        }
        for (const node of nodes.filter(shown)) {
            const r=residue(node),[cx,cy]=position(node.xyz);
            const name=`${scene.id} | ${r.chain} ${r.region} ${r.aa}${r.number}${r.insertion}${atomic?' '+node.name+' '+node.altloc:''}\nID ${node.id}; ${node.contact_count} global atom contacts; degree ${node.contact_degree}; covalent degree ${node.covalent_degree}\n${view.visibleCounts.get(node.id)||0} atom contacts after filters`;
            const point=svgElement(svg,'circle',{cx,cy,r:atomic?2.5:5,fill:COLORS[r.role],
                stroke:atomic&&node.backbone?'#111':'white','stroke-width':1,opacity,tabindex:0,
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
    if(requestedLabels>labelBoxes.length)svgElement(svg,'text',{
        x:350,y:556,'font-size':11,'text-anchor':'middle',class:'label-limit',fill:'#6c3c12'
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
    $('zoom').value=String(Math.min(5,Math.max(1,Math.round(extent/radius*4)/4)));return true;
}

export function resetCameraCoordinates(){for(const id of ['panx','pany','panz'])$(id).value='0';$('zoom').value='1';}

async function render() {
    const token=++generation;
    for(const id of ['side','top'])$(id).setAttribute('aria-busy','true');
    $('status').textContent='Loading selection…';
    try {
        const selected=[...$('structures').selectedOptions];
        const limit=Math.min(16,Math.max(1,parseInt(document.documentElement.dataset.maxOverlay)||16));
        if(selected.length>limit)throw Error(`Select at most ${limit} structures per detailed overlay.`);
        const scenes=await Promise.all(selected.map(o=>loadScene(o.value)));
        if(token!==generation)return;
        choices('chains',scenes.flatMap(s=>s.components.map(c=>c.role+':'+c.chain)),()=>true);
        choices('regions',scenes.flatMap(s=>s.residues.map(r=>r.role+':'+r.region)),v=>DEFAULT_REGIONS.has(v.split(':')[1]));
        const regions=checked('regions'),chains=checked('chains'),groups=checked('edges');
        const visible=r=>regions.has(r.role+':'+r.region)&&chains.has(r.role+':'+r.chain);
        const atomic=$('detail').value==='atom';
        const views=scenes.map(s=>viewFor(s,atomic,visible,groups));
        currentViews=views;nodePage=0;nodeTable();
        $('inspector').textContent='Click a point to inspect its identity and global counts and counts after filters.';
        let extent=10;
        for(const view of views)for(const node of view.nodes.filter(view.shown))extent=Math.max(extent,...node.xyz.map(Math.abs));
        const scale=230/extent*Number($('zoom').value);
        summary(views,atomic);
        draw($('side'),2,views,atomic,scale);draw($('top'),1,views,atomic,scale);
        document.dispatchEvent(new CustomEvent('tcren:render',{detail:{scenes,atomic}}));
        $('status').textContent=scenes.length?scenes.map(s=>`${s.id}\nReference ${s.frame.reference}; fixed origin ${s.frame.origin_component}: ${s.frame.origin}; ${s.frame.anchors} anchors; RMSD ${s.frame.rmsd.toFixed(3)} Å\n${s.warnings.length?s.warnings.join('\n'):'No export warnings'}`).join('\n\n')+`\nCamera offsets x/y/z: ${$('panx').value}/${$('pany').value}/${$('panz').value} Å; zoom ${$('zoom').value}. Molecular reference coordinates remain fixed.`:'Select a structure to inspect.';
    } catch(error) {
        if(token!==generation)return;
        currentViews=[];nodeTable();summary([],false);
        $('side').replaceChildren();$('top').replaceChildren();
        $('inspector').textContent='Selection unavailable.';$('status').textContent=String(error);
    } finally {
        if(token===generation)for(const id of ['side','top'])$(id).setAttribute('aria-busy','false');
    }
}

choices('edges',['alpha:peptide','beta:peptide','alpha:mhc','beta:mhc','mhc:peptide','alpha:beta','b2m:mhc','mhc:mhc'],v=>!['alpha:beta','b2m:mhc','mhc:mhc'].includes(v));
try {
    const response=await fetch('index.json');if(!response.ok)throw Error(`Index ${response.status}`);
    entries=await response.json();
    $('population').textContent=`${entries.length} structures · ${new Set(entries.map(e=>e.reference)).size} compatible reference groups · independent local dashboard`;
    const diagnostics=await fetch('errors.json');
    $('errors').textContent=diagnostics.ok?JSON.stringify(await diagnostics.json(),null,2):'No command-level diagnostics supplied.';
    const exported=await fetch('tables/manifest.json');
    if(exported.ok)for(const item of (await exported.json()).files){
        if(!/^[a-z_]+\/part-[0-9]+\.parquet$/.test(item.file))continue;
        const row=document.createElement('div');
        for(const filename of [item.file,item.tsv_file].filter(Boolean)){
            if(!/^[a-z_]+\/part-[0-9]+\.(parquet|tsv\.gz)$/.test(filename))continue;
            const link=document.createElement('a');link.href='tables/'+filename;link.download='';
            link.textContent=`${filename} (${item.rows} rows) `;row.append(link);
        }
        $('downloads').append(row);
    }
    for(const group of new Set(entries.map(e=>e.reference)))$('group').add(new Option(entries.find(e=>e.reference===group).label||group,group));
    updateStructureOptions();await render();
} catch(error) {$('status').textContent=`Serve this directory over localhost. ${error}`;}
$('group').onchange=()=>{updateStructureOptions();render();};$('search').oninput=()=>{updateStructureOptions();render();};
for(const id of ['structures','detail','weight','labels','backbone','covalent','regions','chains','edges','zoom','panx','pany','panz','labelscope'])$(id).addEventListener('change',render);
$('fit-camera').onclick=()=>{
    function* positions(){for(const view of currentViews)for(const node of view.nodes)if(view.shown(node))yield node.xyz;}
    if(fitCameraCoordinates(positions()))render();
};
$('reset-camera').onclick=()=>{resetCameraCoordinates();render();};
$('nodesearch').oninput=()=>{nodePage=0;nodeTable();};
$('visibleonly').onchange=()=>{nodePage=0;nodeTable();};
$('previous').onclick=()=>{nodePage--;nodeTable();};$('next').onclick=()=>{nodePage++;nodeTable();};
