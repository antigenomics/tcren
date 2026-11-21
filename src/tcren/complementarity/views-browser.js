// All scientific values are serialized table rows. This module only filters and colors them.
const $ = id => document.getElementById(id);
async function readJSON(path, signal) {
  const response = await fetch(path, {cache:'no-cache', signal});
  if (!response.ok) throw Error(`${path}: HTTP ${response.status}`);
  return response.json();
}
const requested = new URL(location).searchParams;
const groupIndex = requested.has('group') || !requested.has('id')
  ? await fetch('groups.json', {cache:'no-cache'}) : null;
if (groupIndex && !groupIndex.ok && groupIndex.status !== 404) {
  throw Error(`groups.json: HTTP ${groupIndex.status}`);
}
const index = groupIndex?.ok ? await groupIndex.json() : await readJSON('index.json');
const grouped = !Array.isArray(index);
const manifest = await readJSON(grouped ? index.manifest : 'manifest.json');
const entries = grouped ? [] : index;
const summaries = grouped ? [] : await readJSON('summary.json');
let ranking = grouped ? index.ranking : entries[0]?.task;
const fixedSurfaces = grouped ? Object.create(null) : entries[0]?.fixed_surface
  ? await readJSON(entries[0].fixed_surface) : {};
const groups = new Map();
for (const entry of entries) {
  if (!groups.has(entry.reference)) groups.set(entry.reference, []);
  groups.get(entry.reference).push(entry);
}
let visibleEntries = [], groupController;
let summaryById = new Map(summaries.map(row=>[row.structure_id,row]));
const squareSize = new Map(Object.entries(manifest.sizes ?? {}));
for (const entry of entries) for (const facet of entry.facets ?? []) {
  const partner = facet.partner_facet;
  squareSize.set(partner, Math.max(squareSize.get(partner) ?? 1, facet.tcr_residues, facet.partner_residues));
}
const shortLabels = new Map(entries.map((entry,i)=>[entry.id, modelLabel(entry,i)]));
function modelLabel(entry,i) {
  return entry.shortLabel || (entry.display_label!==entry.id && entry.display_label) ||
    `${entry.task==='tcr'?'Receptor':'Peptide model'} ${i+1} · ${entry.id.length>16?entry.id.slice(0,8):entry.id}`;
}
let active = null, generation = 0;
let candidateScene = null, motifFocus = null;
let loadController = null;
const palette = [[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
function color(value, lo, hi, reverse=false) {
  let t=Math.min(1,Math.max(0,(value-lo)/(hi-lo||1)));
  if(reverse)t=1-t;
  const q=t*(palette.length-1), i=Math.min(palette.length-2,Math.floor(q)), f=q-i;
  return `rgb(${palette[i].map((v,k)=>Math.round(v*(1-f)+palette[i+1][k]*f)).join(',')})`;
}
function canvasPanel(parent,title,nx,ny,records,coordinates,value,lo,hi,reverse=false,flip=false,highlight=()=>false) {
  const figure=document.createElement('figure'), caption=document.createElement('figcaption'), canvas=document.createElement('canvas');
  caption.textContent=title; figure.append(caption,canvas); parent.append(figure);
  canvas.width=500;canvas.height=canvas.width;
  const context=canvas.getContext('2d'), lookup=new Map();
  const outlines=[];
  for(const r of records){const [x,y]=coordinates(r);lookup.set(`${x},${y}`,r);
    const measurement=value(r);if(!Number.isFinite(measurement))continue;
    context.fillStyle=color(measurement,lo,hi,reverse);
    const plotY = flip ? ny-1-y : y;
    const left = Math.ceil(x*canvas.width/nx), right = Math.ceil((x+1)*canvas.width/nx);
    const top = Math.ceil(plotY*canvas.height/ny), bottom = Math.ceil((plotY+1)*canvas.height/ny);
    context.fillRect(left,top,right-left,bottom-top);
    if(highlight(r))outlines.push([left+1,top+1,right-left-2,bottom-top-2]);
  }
  for (const [stroke,width] of [['#172535',4],['#ffae46',2]]) {
    context.strokeStyle=stroke;context.lineWidth=width;
    for (const rect of outlines) context.strokeRect(...rect);
  }
  canvas.onmousemove=event=>{const bounds=canvas.getBoundingClientRect();
    const x=Math.floor((event.clientX-bounds.left)*nx/bounds.width), yp=Math.floor((event.clientY-bounds.top)*ny/bounds.height), y=flip?ny-1-yp:yp;
    $('details').textContent=JSON.stringify(lookup.get(`${x},${y}`)??{missing:true},null,2);};
}
function selected(group){return [...$(group).querySelectorAll('input:checked')].map(e=>e.value);}
function renderSurface(parent, cells) {
  const heights=cells.map(r=>r.surface_z).filter(v=>v!==null);
  const [ny,nx]=manifest.grid_ny_nx, zlo=heights.length?Math.min(...heights):0,zhi=heights.length?Math.max(...heights):1;
  canvasPanel(parent,'Hydropathy KD: −4.5 (purple) → +4.5 (yellow)',nx,ny,cells,r=>[r.ix,r.iy],r=>r.hydropathy,-4.5,4.5,false,true);
  canvasPanel(parent,`Canonical surface z: ${zlo.toFixed(2)} → ${zhi.toFixed(2)} Å`,nx,ny,cells,r=>[r.ix,r.iy],r=>r.surface_z,zlo,zhi,false,true);
}
function render(){
  if(!active)return;
  $('topography').replaceChildren();$('heatmaps').replaceChildren();
  renderSurface($('topography'), active.tables.topography_cells);
  $('definitions').textContent=`${manifest.definitions[ranking]}. x left→right, y bottom→top; extent ${manifest.extent_xy.join(', ')} Å. ${manifest.definitions.height}. Height colors span this structure; KD colors are fixed.`;
  const rows=selected('rows'),columns=selected('columns'),metric=$('metric').value,masked=$('contacts').checked,pairs=active.tables.contact_heatmap;
  const hi=metric==='ca_distance'?25:pairs.reduce((m,r)=>Math.max(m,r.atom_contact_count),1);
  $('scale').textContent=`${metric}: 0 → ${hi}${metric==='ca_distance'?' Å (saturated at 25; yellow=near, purple=far)':' witnesses (purple=0, yellow=max)'}. Each panel: TCR N→C left→right; partner N→C top→bottom. Empty space pads each partner row to a fixed square across models. Hover for sequence, author numbering and chain IDs.`;
  $('heatmaps').style.setProperty('--columns',Math.max(1,columns.length));
  const focused = new Set(focusedResidues().map(residue => residue.id));
  for(const row of rows)for(const column of columns){
    const panel=document.createElement('div');panel.className='facet';$('heatmaps').append(panel);
    const source=pairs.filter(r=>r.partner_facet===row&&r.loop_facet===column);
    const visible=source.filter(r=>!masked||r.contacting);
    const note=document.createElement('p');
    const missing=visible.filter(r=>r[metric]===null).length;
    note.textContent=!source.length?'Absent facet':!visible.length?`${row} · ${column}: No contacting pairs in this facet`:`${visible.length} pairs shown · ${source.length-visible.length} masked · ${missing} missing ${metric==='ca_distance'?'Cα distances':'values'}`;
    panel.append(note);
    const size = squareSize.get(row) ?? 1;
    canvasPanel(panel,`${row} · ${column}`,size,size,visible,r=>[r.tcr_index,r.partner_index],r=>r[metric],0,hi,metric==='ca_distance',false,
      pair=>pair.contacting && focused.has(ranking === 'tcr' ? pair.tcr_residue_id : pair.partner_residue_id));
  }
}
for(const [group,values] of [['rows',manifest.partner_order],['columns',manifest.loop_order]])for(const value of values){
  const label=document.createElement('label'),input=document.createElement('input');input.type='checkbox';input.value=value;input.checked=true;
  input.onchange=render;label.append(input,document.createTextNode(value));$(group).append(label);
}
function showTrace(parent, trace, caption, contactLinks=[], markers=[]) {
  parent.replaceChildren();
  if (!trace.length && !markers.length) return;
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  const colors = {alpha:'#0072b2', beta:'#56b4e9', peptide:'#d55e00', mhc:'#009e73', b2m:'#666'};
  const points = [...trace, ...contactLinks].flatMap(([, first, second]) => [first, second]).concat(markers);
  const xmin = Math.min(...points.map(p => p[0])) - 2, xmax = Math.max(...points.map(p => p[0])) + 2;
  const zmin = Math.min(...points.map(p => -p[2])) - 2, zmax = Math.max(...points.map(p => -p[2])) + 2;
  svg.setAttribute('viewBox', `${xmin} ${zmin} ${xmax-xmin} ${zmax-zmin}`);
  svg.setAttribute('role', 'img');
  svg.setAttribute('aria-label', caption);
  for (const [role, first, second] of [...trace, ...contactLinks]) {
    const path = document.createElementNS(svg.namespaceURI, 'line');
    for (const [key, value] of Object.entries({
      x1:first[0], y1:-first[2], x2:second[0], y2:-second[2], stroke:colors[role] ?? '#9b460e',
      'stroke-dasharray':role === 'contact' ? '4 2' : 'none',
      'stroke-width':role === 'contact' ? 2.5 : 1.5, 'vector-effect':'non-scaling-stroke',
    })) path.setAttribute(key, value);
    svg.append(path);
  }
  for (const point of markers) {
    const marker = document.createElementNS(svg.namespaceURI, 'circle');
    for (const [key,value] of Object.entries({
      cx:point[0],cy:-point[2],r:1.6,fill:'#ffae46',stroke:'#172535',
      'stroke-width':1.5,'vector-effect':'non-scaling-stroke','data-motif-hit':'',
    })) marker.setAttribute(key,value);
    svg.append(marker);
  }
  const label = document.createElement('small');
  label.textContent = caption;
  parent.append(label, svg);
}
const format = (value, digits, suffix = '') =>
  Number.isFinite(value) ? value.toFixed(digits) + suffix : 'unavailable';

const aminoAcids = 'ACDEFGHIKLMNPQRSTVWY';

function motifText(parent, tag, text) {
  const node = document.createElement(tag);
  node.textContent = text;
  parent.append(node);
  return node;
}

function setMotifFocus(focus, button) {
  motifFocus = focus;
  for (const element of document.querySelectorAll('[data-motif-control][aria-pressed="true"]')) {
    element.setAttribute('aria-pressed','false');
  }
  updateCandidateMotif();
  if (button?.isConnected) button.setAttribute('aria-pressed','true');
  render();
}

function motifControl(parent, text, title, role, columns, index, aa) {
  const button = motifText(parent, 'button', text);
  button.type = 'button';
  button.dataset.motifControl = '';
  button.title = title;
  button.setAttribute('aria-label', title);
  button.setAttribute('aria-pressed','false');
  button.onclick = () => setMotifFocus(
    {role, length:columns.length, positions:[{index, aa}]}, button);
}

function motifStrip(parent, columns, marks, role) {
  const strip = motifText(parent, 'div', '');
  strip.className = 'motif-strip';
  strip.style.setProperty('--positions', columns.length);
  for (const [index, column] of columns.entries()) {
    const peak = Math.max(...Object.values(column));
    const leaders = Object.keys(column).filter(aa => column[aa] === peak).sort();
    const cell = motifText(strip, 'span', '');
    const title = `Position ${index+1} · ${leaders.join('/')} · ${peak.toLocaleString()}`;
    motifControl(cell, leaders.length === 1 ? leaders[0] : '·', title, role, columns, index, leaders.join(''));
    if (leaders.some(aa => marks[index+1]?.includes(aa))) cell.className = 'motif-mark';
  }
}

function motifMatrix(parent, columns, marks, label, range, role) {
  const cells = columns.flatMap(column => Object.values(column));
  const lo = range ? Math.min(...cells, 0) : 0, hi = range ? Math.max(...cells, 0) : 1;
  const table = document.createElement('table');
  table.className = 'motif-matrix';
  table.setAttribute('aria-label', label);
  const heading = table.createTHead().insertRow();
  motifText(heading, 'th', 'AA');
  for (const [index] of columns.entries()) motifText(heading, 'th', index+1).scope = 'col';
  const body = table.createTBody();
  for (const aa of aminoAcids) {
    const row = body.insertRow();
    motifText(row, 'th', aa).scope = 'row';
    for (const [index, column] of columns.entries()) {
      const value = column[aa] ?? 0;
      const cell = row.insertCell();
      cell.style.background = value < 0
        ? `rgba(158,72,37,${value/lo})` : `rgba(9,98,108,${value/(hi||1)})`;
      const exact = range ? value.toString() : `${(100*value).toFixed(2)}%`;
      cell.title = `${aa} · position ${index+1} · ${exact}`;
      motifControl(cell, '', cell.title, role, columns, index, aa);
      if (marks[index+1]?.includes(aa)) cell.className = 'motif-mark';
    }
  }
  parent.append(table);
  if (range) motifText(parent, 'small', `Response scale within this scan · ${lo.toLocaleString()} to ${hi.toLocaleString()}`);
}

function motifSource(parent, note, url, linkLabel) {
  if (!note && !url) return;
  const p = motifText(parent, 'p', note ?? '');
  p.className = 'motif-source';
  if (url && /^https?:\/\//i.test(url)) {
    p.append(' ');
    motifText(p, 'a', linkLabel).href = url;
  }
}

function renderMotif(fixed) {
  const panel = $('fixed-motif');
  panel.replaceChildren();
  const records = fixed?.motif;
  panel.hidden = !records?.length;
  if (panel.hidden) return;
  if (ranking === 'epitope') {
    const {note, url, marks = {}} = fixed.publication ?? {};
    motifText(panel, 'h3', 'Positional-scanning CPL');
    motifText(panel, 'small', 'Peak residue per position · signed MIP-1β response below');
    motifStrip(panel, records, marks, "peptide");
    motifMatrix(panel, records, marks, 'CPL response by amino acid and position', true, 'peptide');
    motifSource(panel, note, url, 'Publication');
    return;
  }

  motifText(panel, 'h3', 'TCR paratope motifs');
  motifText(panel, 'small', 'Repertoire motif per chain · CDR3 excluding anchors');
  for (const {chain, id, members, v, v_count, pwm} of records) {
    const card = document.createElement('article');
    panel.append(card);
    motifText(card, 'strong', chain === 'TRA' ? 'CDR3α' : 'CDR3β');
    const known = chain === 'TRB' ? fixed.known_motif?.slice(0,3) : null;
    const marks = known?.[1] ?? {};
    if (known) motifText(card, 'small', known[0] + ' highlighted');
    const role = chain === 'TRA' ? 'alpha' : 'beta';
    motifStrip(card, pwm, marks, role);
    motifText(card, 'small', `Top V · ${v} · ${(100*v_count/members).toFixed(1)}% (${v_count}/${members})`);
    const details = document.createElement('details');
    card.append(details);
    motifText(details, 'summary', `Full PWM · ${id}`);
    motifMatrix(details, pwm, marks, `${chain} CDR3 amino-acid frequency PWM`, false, role);
    if (known) {
      motifSource(card, fixed.known_motif[3], known[2], 'Publication');
    }
  }
}

async function selectGroup(reference, requested) {
  groupController?.abort();
  loadController?.abort();
  generation++;
  const controller = groupController = new AbortController();
  active = null;
  candidateScene = null;
  $('structure').replaceChildren();
  for(const panel of ['fixed-topography','fixed-trace','candidate-trace','topography','heatmaps'])$(panel).replaceChildren();
  if(grouped) {
    $('status').textContent = 'Loading reference group…';
    try {
      const group = index.groups.find(item=>item.id===reference);
      const payload = await readJSON(group.file, controller.signal);
      if(controller.signal.aborted)return;
      visibleEntries = payload.entries;
      fixedSurfaces[reference] = payload.fixed;
      summaryById = new Map((payload.summaries ?? []).map(row=>[row.structure_id,row]));
      shortLabels.clear();
      visibleEntries.forEach((entry,i)=>shortLabels.set(entry.id,modelLabel(entry,i)));
      ranking = visibleEntries[0]?.task;
    } catch(error) {
      if(!controller.signal.aborted)$('status').textContent = `Cannot load reference group: ${error.message}`;
      return;
    }
  } else visibleEntries = groups.get(reference);
  if(!visibleEntries.length)return;
  $('fixed-partner').value = reference;
  $('structure').replaceChildren(...visibleEntries.map(entry => {
    const label = [shortLabels.get(entry.id), entry.outcome].filter(Boolean).join(' — ');
    const option = new Option(label, entry.id);
    option.title = entry.id;
    return option;
  }));
  if (visibleEntries.some(entry => entry.id === requested)) $('structure').value = requested;
  const fixed = fixedSurfaces[reference];
  motifFocus = null;
  renderMotif(fixed);
  $('fixed-topography').replaceChildren();
  $('fixed-trace').replaceChildren();
  $('reference-angles').hidden = ranking !== 'epitope' || !fixed;
  $('fixed-surface-heading').hidden = !fixed;
  if (fixed) {
    $('fixed-surface-heading').textContent = ranking === 'tcr' ? 'Peptide upper face' : 'TCR variable-domain lower face';
    renderSurface($('fixed-topography'), fixed.cells);
    showTrace($('fixed-trace'), fixed.trace, 'Reference side view (x/z, fitted)');
    const {crossing, incident} = fixed.model_metrics ?? {};
    $('reference-angles').textContent =
      `Reference complex · Crossing ${format(crossing, 1, '°')} · Incidence ${format(incident, 1, '°')}`;
  }
  return load();
}

function candidatePositions(role) {
  if (!candidateScene) return [];
  return candidateScene.residues.filter(residue => residue.role === role &&
    (role === 'peptide' || residue.region === 'CDR3'))
    .sort((a,b) => a.seq_index-b.seq_index);
}

function focusedResidues() {
  if (!motifFocus) return [];
  const residues = candidatePositions(motifFocus.role);
  if (residues.length !== motifFocus.length || new Set(residues.map(r=>r.chain)).size !== 1) return [];
  return motifFocus.positions.map(({index, aa}) => residues[index]?.aa && aa.includes(residues[index].aa)
    ? residues[index] : null).filter(Boolean);
}

function updateCandidateMotif() {
  const known = $('known-motif'), status = $('motif-selection'), clear = $('clear-motif');
  const jump = $('jump-motif-heatmaps');
  jump.hidden = true;
  known.replaceChildren();
  known.hidden = true;
  clear.hidden = !motifFocus;
  status.textContent = motifFocus && !candidateScene ? 'Waiting for the candidate structure.' : '';
  if (!candidateScene) return;
  const entry = visibleEntries[$('structure').selectedIndex];
  const reference = fixedSurfaces[entry.reference]?.known_motif;
  if (ranking === 'tcr' && reference) {
    const residues = candidatePositions('beta');
    const sequence = residues.map(r=>r.aa).join('');
    const cdr3 = fixedSurfaces[entry.reference].motif?.find(r=>r.chain==='TRB')?.pwm ?? [];
    const positions = Object.entries(reference[1]).map(([index,aa])=>({index:Number(index)-1,aa}));
    const motif = positions.map(p=>p.aa).join('');
    const sameLength = residues.length === cdr3.length;
    const oneChain = new Set(residues.map(r=>r.chain)).size === 1;
    const aligned = oneChain && sameLength && positions.every(p=>residues[p.index]?.aa===p.aa);
    const occurrence = oneChain ? sequence.indexOf(motif) : -1;
    const message = !oneChain ? 'CDR3β unavailable' : aligned ? 'present at reference positions'
      : occurrence >= 0 ? `observed at CDR3 position ${occurrence+1}${sameLength ? ' (shifted)' : ' (different length)'}`
      : 'absent';
    known.hidden = false;
    const button = motifText(known,'button',`${reference[0]} · ${message}`);
    button.type = 'button';
    button.dataset.motifControl = '';
    button.setAttribute('aria-pressed', String(Boolean(motifFocus?.known)));
    button.disabled = !aligned && occurrence < 0;
    button.onclick = () => setMotifFocus({role:'beta', length:residues.length, known:true,
      positions:aligned ? positions : [...motif].map((aa,i)=>({aa,index:occurrence+i}))}, button);
  }

  const selected = focusedResidues();
  const ids = new Set(selected.map(residue=>residue.id));
  const coords = residue => candidateScene.atoms[residue.ca_atom]?.xyz;
  const contacts = candidateScene.residue_contacts.filter(edge=>ids.has(edge.source)||ids.has(edge.target));
  if (motifFocus) {
    const residues = candidatePositions(motifFocus.role);
    if (residues.length !== motifFocus.length || new Set(residues.map(r=>r.chain)).size !== 1) {
      status.textContent = `Position mapping unavailable · candidate length ${residues.length}, selected length ${motifFocus.length}`;
    } else if (!selected.length) {
      status.textContent = motifFocus.positions.map(({index,aa})=>
        `Position ${index+1}: candidate ${residues[index]?.aa ?? 'unavailable'}, selected ${aa.split('').join('/')}`).join(' · ');
    } else {
      const hits = selected.map(r=>`${r.aa} ${r.chain}:${r.number}${r.insertion ?? ''}`).join(', ');
      const witnesses = contacts.reduce((sum,edge)=>sum+edge.contact_count,0);
      status.textContent = `${hits} · ${contacts.length} residue contacts · ${witnesses} atom-pair witnesses`;
      jump.hidden = !contacts.length;
    }
  }
  const roles = ranking === 'tcr' ? ['alpha','beta'] : ['peptide','mhc','b2m'];
  const trace = candidateScene.trace.flatMap(([a,b]) => {
    const first = candidateScene.residues[a], second = candidateScene.residues[b];
    return roles.includes(first.role) ? [[first.role, coords(first), coords(second)]] : [];
  });
  const links = contacts.flatMap(({source,target})=>{
    const first=coords(candidateScene.residues[source]), second=coords(candidateScene.residues[target]);
    return first && second ? [['contact',first,second]] : [];
  });
  const preview = $('candidate-trace');
  preview.hidden = false;
  showTrace(preview, trace, 'Candidate side view (x/z, fitted)', links, selected.map(coords).filter(Boolean));
}

async function load(){const entry=visibleEntries[$('structure').selectedIndex];if(!entry)return;
  const current=++generation;
  const outcome = $('outcome');
  outcome.hidden = !entry.outcome;
  outcome.textContent = entry.outcome ?? '';
  if ('library_response' in entry) {
    outcome.textContent += ` · Library response ${format(entry.library_response, 1)}`;
  }
  const molecular = $('molecular-view');
  molecular.hidden = !new URL(location).pathname.endsWith('/table-views/');
  if (!molecular.hidden) molecular.href = `../?${new URLSearchParams({id:entry.id, ...(grouped ? {group:entry.reference} : {})})}`;
  const metrics = $('model-metrics');
  metrics.hidden = !entry.model_metrics;
  if (entry.model_metrics) {
    const {crossing, crossing_signed, incident, iptm, plddt} = entry.model_metrics;
    metrics.textContent = [
      `Crossing ${format(crossing, 1, '°')}`,
      `Signed crossing ${format(crossing_signed, 1, '°')}`,
      `Incidence ${format(incident, 1, '°')}`,
      `ipTM ${format(iptm, 3)}`,
      `Mean pLDDT ${format(plddt, 1)}`,
    ].join(' · ');
  }
  loadController?.abort();const controller=loadController=new AbortController();
  const preview = $('candidate-trace');
  candidateScene = null;
  if (motifFocus?.known) motifFocus = null;
  preview.replaceChildren();
  preview.hidden = true;
  updateCandidateMotif();
  if (!molecular.hidden && entry.molecular_file) {
    fetch(`../${entry.molecular_file}`, {signal:controller.signal}).then(async response => {
      if (!response.ok) throw Error(response.status);
      const scene = await new Response(response.body.pipeThrough(
        new DecompressionStream('gzip'), {signal:controller.signal})).json();
      if (current !== generation) return;
      candidateScene = scene;
      updateCandidateMotif();
      render();
    }).catch(error => {
      if (current !== generation || controller.signal.aborted) return;
      preview.hidden = false;
      preview.textContent = `Side view unavailable: ${error.message}`;
    });
  }
  active=null;$('topography').replaceChildren();$('heatmaps').replaceChildren();
  $('status').textContent=`Loading ${shortLabels.get(entry.id)}…`;
  try{const response=await fetch(entry.cache,{signal:controller.signal});if(!response.ok)throw Error(response.status);
    const stream=response.body.pipeThrough(new DecompressionStream('gzip'),{signal:controller.signal});
    const payload=await new Response(stream).json();controller.signal.throwIfAborted();if(current!==generation)return;
    active=payload;render();
    $('structure').title=entry.id;
    const warning=(entry.clipped ?? summaryById.get(entry.id)?.surface_may_be_clipped)?' · Occupied grid boundary: inspect extent':'';
    $('status').textContent=`${shortLabels.get(entry.id)} · ${entry.task} ranking · reference ${shortLabels.get(entry.reference)||entry.reference.slice(0,8)} · ${active.tables.contact_heatmap.length} residue pairs${warning} · source: ${Object.values(payload.sources).join(', ')}`;
    const url=new URL(location);url.searchParams.set('id',entry.id);url.searchParams.set('group',entry.reference);history.replaceState(null,'',url);
    document.dispatchEvent(new Event('tcren:model'));
  }catch(error){controller.abort();if(current===generation){active=null;$('topography').replaceChildren();$('heatmaps').replaceChildren();$('status').textContent=`Cannot load ${entry.id}: ${error.message}`;}}
}
const choices = grouped ? index.groups : [...groups].map(([id,models])=>
  ({id,label:models[0].group_label ?? shortLabels.get(id) ?? id}));
if (choices.length) {
  const tcrRanking = ranking === 'tcr';
  $('fixed-label').textContent = tcrRanking ? 'pMHC' : 'TCR';
  $('candidate-label').textContent = tcrRanking ? 'TCR model' : 'pMHC model';
  $('fixed-heading').textContent = tcrRanking ? 'Fixed pMHC' : 'Fixed TCR';
  $('candidate-heading').textContent = tcrRanking ? 'Candidate TCR' : 'Candidate pMHC';
  for (const {id,label} of choices) $('fixed-partner').add(new Option(label,id));
  $('fixed-partner').onchange = () => selectGroup($('fixed-partner').value);
  $('clear-motif').onclick = () => setMotifFocus(null);
  $('structure').onchange=load;$('metric').onchange=render;$('contacts').onchange=render;
  const id = requested.get('id');
  const initial = grouped ? requested.get('group') : entries.find(entry=>entry.id===id)?.reference;
  await selectGroup(choices.find(group=>group.id===initial)?.id ?? choices[0].id,id);
}
