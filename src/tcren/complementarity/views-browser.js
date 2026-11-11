// All scientific values are serialized table rows. This module only filters and colors them.
const $ = id => document.getElementById(id);
const [entries, manifest, summaries] = await Promise.all(['index.json','manifest.json','summary.json'].map(p=>fetch(p,{cache:'no-cache'}).then(r=>r.json())));
const summaryById=new Map(summaries.map(r=>[r.structure_id,r]));
const shortLabels=new Map(entries.map((e,i)=>[e.id,(e.display_label&&e.display_label!==e.id?e.display_label:null)||`${e.task==='tcr'?'Receptor':'Peptide model'} ${i+1} · ${e.id.length>16?e.id.slice(0,8):e.id}`]));
let active = null, generation = 0;
let loadController = null;
const palette = [[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
function color(value, lo, hi, reverse=false) {
  if(value===null || !Number.isFinite(value)) return '#eee';
  let t=Math.min(1,Math.max(0,(value-lo)/(hi-lo||1)));
  if(reverse)t=1-t;
  const q=t*(palette.length-1), i=Math.min(palette.length-2,Math.floor(q)), f=q-i;
  return `rgb(${palette[i].map((v,k)=>Math.round(v*(1-f)+palette[i+1][k]*f)).join(',')})`;
}
function canvasPanel(parent,title,nx,ny,records,coordinates,value,lo,hi,reverse=false,flip=false) {
  const figure=document.createElement('figure'), caption=document.createElement('figcaption'), canvas=document.createElement('canvas');
  caption.textContent=title; figure.append(caption,canvas); parent.append(figure);
  canvas.width=500;canvas.height=parent.id==='topography'?500:230;
  const context=canvas.getContext('2d'), lookup=new Map();
  context.fillStyle='#eee';context.fillRect(0,0,canvas.width,canvas.height);
  for(const r of records){const [x,y]=coordinates(r);lookup.set(`${x},${y}`,r);
    context.fillStyle=color(value(r),lo,hi,reverse);
    context.fillRect(x*canvas.width/nx,(flip?ny-1-y:y)*canvas.height/ny,Math.ceil(canvas.width/nx),Math.ceil(canvas.height/ny));}
  canvas.onmousemove=event=>{const bounds=canvas.getBoundingClientRect();
    const x=Math.floor((event.clientX-bounds.left)*nx/bounds.width), yp=Math.floor((event.clientY-bounds.top)*ny/bounds.height), y=flip?ny-1-yp:yp;
    $('details').textContent=JSON.stringify(lookup.get(`${x},${y}`)??{missing:true},null,2);};
}
function selected(group){return [...$(group).querySelectorAll('input:checked')].map(e=>e.value);}
function render(){
  if(!active)return;
  $('topography').replaceChildren();$('heatmaps').replaceChildren();
  const cells=active.tables.topography_cells, heights=cells.map(r=>r.surface_z).filter(v=>v!==null);
  const [ny,nx]=manifest.grid_ny_nx, zlo=heights.length?Math.min(...heights):0,zhi=heights.length?Math.max(...heights):1;
  canvasPanel($('topography'),'Hydropathy KD: −4.5 (purple) → +4.5 (yellow)',nx,ny,cells,r=>[r.ix,r.iy],r=>r.hydropathy,-4.5,4.5,false,true);
  canvasPanel($('topography'),`Canonical surface z: ${zlo.toFixed(2)} → ${zhi.toFixed(2)} Å`,nx,ny,cells,r=>[r.ix,r.iy],r=>r.surface_z,zlo,zhi,false,true);
  $('definitions').textContent=`${manifest.definitions[entries[$('structure').selectedIndex].task]}. x left→right, y bottom→top; extent ${manifest.extent_xy.join(', ')} Å. ${manifest.definitions.height}. Height colors span this structure; KD colors are fixed.`;
  const rows=selected('rows'),columns=selected('columns'),metric=$('metric').value,masked=$('contacts').checked,pairs=active.tables.contact_heatmap;
  const hi=metric==='ca_distance'?25:pairs.reduce((m,r)=>Math.max(m,r.atom_contact_count),1);
  $('scale').textContent=`${metric}: 0 → ${hi}${metric==='ca_distance'?' Å (saturated at 25; yellow=near, purple=far)':' witnesses (purple=0, yellow=max)'}. Each panel: TCR N→C left→right; partner N→C top→bottom. Hover for sequence, author numbering and chain IDs.`;
  $('heatmaps').style.setProperty('--columns',Math.max(1,columns.length));
  for(const row of rows)for(const column of columns){
    const panel=document.createElement('div');panel.className='facet';$('heatmaps').append(panel);
    const source=pairs.filter(r=>r.partner_facet===row&&r.loop_facet===column);
    const visible=source.filter(r=>!masked||r.contacting);
    if(!source.length){panel.textContent=`${row} · ${column}: absent facet`;continue;}
    const note=document.createElement('p');
    const missing=visible.filter(r=>r[metric]===null).length;
    note.textContent=!visible.length?`${row} · ${column}: No contacting pairs in this facet`:`${visible.length} pairs shown · ${source.length-visible.length} masked · ${missing} missing ${metric==='ca_distance'?'Cα distances':'values'}`;
    panel.append(note);
    if(!visible.length)continue;
    const nx=1+source.reduce((m,r)=>Math.max(m,r.tcr_index),0),ny=1+source.reduce((m,r)=>Math.max(m,r.partner_index),0);
    canvasPanel(panel,`${row} · ${column}`,nx,ny,visible,r=>[r.tcr_index,r.partner_index],r=>r[metric],0,hi,metric==='ca_distance');
  }
}
for(const [group,values] of [['rows',manifest.partner_order],['columns',manifest.loop_order]])for(const value of values){
  const label=document.createElement('label'),input=document.createElement('input');input.type='checkbox';input.value=value;input.checked=true;
  input.onchange=render;label.append(input,document.createTextNode(value));$(group).append(label);
}
for(const entry of entries){const option=document.createElement('option');option.value=entry.id;option.textContent=shortLabels.get(entry.id);option.title=entry.id;$('structure').append(option);}
async function load(){const current=++generation,entry=entries[$('structure').selectedIndex];
  loadController?.abort();const controller=loadController=new AbortController();
  active=null;$('topography').replaceChildren();$('heatmaps').replaceChildren();
  $('status').textContent=`Loading ${shortLabels.get(entry.id)}…`;
  try{const response=await fetch(entry.cache,{signal:controller.signal});if(!response.ok)throw Error(response.status);
    const stream=response.body.pipeThrough(new DecompressionStream('gzip'),{signal:controller.signal});
    const payload=await new Response(stream).json();controller.signal.throwIfAborted();if(current!==generation)return;
    active=payload;render();
    $('structure').title=entry.id;
    const warning=summaryById.get(entry.id)?.surface_may_be_clipped?' · Occupied grid boundary: inspect extent':'';
    $('status').textContent=`${shortLabels.get(entry.id)} · ${entry.task} ranking · reference ${shortLabels.get(entry.reference)||entry.reference.slice(0,8)} · ${active.tables.contact_heatmap.length} residue pairs${warning} · source: ${Object.values(payload.sources).join(', ')}`;
    const url=new URL(location);url.searchParams.set('id',entry.id);history.replaceState(null,'',url);
  }catch(error){controller.abort();if(current===generation){active=null;$('topography').replaceChildren();$('heatmaps').replaceChildren();$('status').textContent=`Cannot load ${entry.id}: ${error.message}`;}}
  finally{if(loadController===controller)loadController=null;}
}
const requested=new URL(location).searchParams.get('id');if(entries.some(e=>e.id===requested))$('structure').value=requested;
$('structure').onchange=load;$('metric').onchange=render;$('contacts').onchange=render;
await load();
