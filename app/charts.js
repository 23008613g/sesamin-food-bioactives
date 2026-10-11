// Dependency-free charts. Every mark is generated from the supplied study data.
const INK = '#183b38', MUTED = '#768579', GRID = '#e6e9df';
export const colors = ['#14796f','#c7793b','#819775','#687f9b','#bb797c','#a99a5c','#7a68a1','#5d999a','#bd9776','#a5aa82'];
export const escapeHtml = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
export const finite = value => typeof value === 'number' && Number.isFinite(value);
export function expressionColor(value,max){if(!finite(value))return '#c4c4bd';if(value===0)return '#e0e5d9';const f=value/(max||1),a=[163,202,166],b=[20,87,74];return `rgb(${a.map((v,j)=>Math.round(v+(b[j]-v)*f)).join(',')})`;}
export function format(value, digits = 2) { return finite(value) ? value.toLocaleString('en-US', {minimumFractionDigits:digits, maximumFractionDigits:digits}) : 'NA'; }
export function signed(value, digits = 2) { return finite(value) ? `${value > 0 ? '+' : value < 0 ? '−' : ''}${format(Math.abs(value),digits)}` : 'NA'; }
export function qformat(value) { if(!finite(value))return 'NA'; if(value === 0)return '0'; if(value<.0001)return value.toExponential(2);return value.toFixed(value<.01?5:4).replace(/0+$/,'').replace(/\.$/,''); }
function svgOpen(w,h,label) {return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" role="img" aria-label="${escapeHtml(label)}"><title>${escapeHtml(label)}</title>`;}
function txt(x,y,text,attrs='') {return `<text x="${x}" y="${y}" font-family="Segoe UI,Arial,sans-serif" font-size="10" fill="${MUTED}" ${attrs}>${escapeHtml(text)}</text>`;}
function line(x1,y1,x2,y2,attrs='') {return `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${GRID}" ${attrs}/>`;}
export function emptyChart(message) {return `<div class="empty-chart"><strong>NA</strong><span>${escapeHtml(message)}</span></div>`;}
function tickFormat(v,range) {return Math.abs(v)<1e-10?'0':Math.abs(range)<.1?v.toFixed(3):Math.abs(range)<1?v.toFixed(2):Math.abs(range)<10?v.toFixed(1):v.toFixed(0);}
export function dotChart({samples,valueKey,groups,labels,unit,title,startAtZero=true}) {
  if(!samples.length)return '';
  const valid = samples.filter(s=>finite(s[valueKey])); if(!valid.length)return '';
  const w=470,h=235,left=49,right=22,top=28,bottom=45,plotH=h-top-bottom;
  const vals=valid.map(s=>s[valueKey]);let min=startAtZero?0:Math.min(...vals),max=Math.max(...vals);let range=max-min;if(range===0)range=Math.max(Math.abs(max)*.3,1);max+=range*.18;if(!startAtZero)min-=range*.18;
  const y=v=>h-bottom-(v-min)/(max-min)*plotH;
  const x=i=>left+(i+.5)*(w-left-right)/groups.length;
  let out=svgOpen(w,h,title)+txt(left,13,unit,'font-size="9"');
  for(let i=0;i<5;i++){const val=min+(max-min)*i/4,yy=y(val);out+=line(left,yy,w-right,yy)+txt(left-9,yy+3,tickFormat(val,max-min),'text-anchor="end" font-size="9"');}
  const palette = groups.length===3?['#94a389','#c79969','#14796f']:['#94a389','#c7793b'];
  groups.forEach((g,i)=>{const points=valid.filter(s=>s.group===g);if(!points.length)return;const mean=points.reduce((sum,s)=>sum+s[valueKey],0)/points.length;out+=line(x(i)-28,y(mean),x(i)+28,y(mean),`style="stroke:${palette[i]};stroke-width:2"`);points.forEach((s,j)=>{const jitter=(j-(points.length-1)/2)*12;out+=`<circle cx="${x(i)+jitter}" cy="${y(s[valueKey])}" r="4.7" fill="${palette[i]}" stroke="#fffefa" stroke-width="1.5"><title>${escapeHtml(`${s.sample}: ${format(s[valueKey],4)} ${unit}`)}</title></circle>`;});out+=txt(x(i),h-24,labels[i],'text-anchor="middle" fill="#183b38" font-size="10"')+txt(x(i),h-9,`n = ${points.length}`,'text-anchor="middle" font-size="8"');});
  return out+'</svg>';
}
export function localizationChart(rows,lang) {
  const zh=lang==='zh', libs=[...new Set(rows.map(r=>r.library))];
  const pairs=libs.map(lib=>({lib,h:rows.find(r=>r.library===lib&&r.category==='Hepatocytes'),c:rows.find(r=>r.library===lib&&r.category==='Cholangiocytes')})).filter(p=>p.h&&p.c&&finite(p.h.detection_pct)&&finite(p.c.detection_pct));
  if(!pairs.length)return emptyChart(zh?'沒有此候選的單細胞檢出資料。':'No measured single-cell detection data for this candidate.');
  const w=470,h=228,left=93,right=30,top=37,bottom=40;const observedMax=Math.max(...pairs.flatMap(p=>[p.h.detection_pct,p.c.detection_pct]));const xmax=observedMax<=2?2:observedMax<=5?5:observedMax<=10?10:observedMax<=25?25:observedMax<=50?50:100;
  const x=v=>left+v/xmax*(w-left-right);let out=svgOpen(w,h,zh?'各技術文庫內膽管細胞與肝細胞的檢出比例':'Detection in cholangiocytes and hepatocytes within each technical library');
  out+=`<circle cx="106" cy="10" r="4" fill="#14796f"/>`+txt(117,13,zh?'肝細胞':'Hepatocytes','font-size="9"')+`<circle cx="252" cy="10" r="4" fill="#c7793b"/>`+txt(263,13,zh?'膽管細胞':'Cholangiocytes','font-size="9"');
  for(let i=0;i<=4;i++){const val=xmax*i/4;out+=line(x(val),top,x(val),h-bottom)+txt(x(val),h-23,`${tickFormat(val,xmax)}%`,'text-anchor="middle" font-size="9"');}
  const conditionCounts={};
  pairs.forEach((p,i)=>{const condition=p.h.condition,rep=conditionCounts[condition]=(conditionCounts[condition]||0)+1;const label=condition==='CTRL'?`${zh?'對照':'Control'} ${rep}`:`${condition} ${rep}`;const yy=top+17+i*34;out+=txt(left-11,yy-1,label,'text-anchor="end" fill="#183b38" font-size="9"')+txt(left-11,yy+11,p.lib,'text-anchor="end" font-size="7"');out+=line(x(p.h.detection_pct),yy,x(p.c.detection_pct),yy,'style="stroke:#becabc;stroke-width:2.5"');for(const [r,c]of[[p.h,'#14796f'],[p.c,'#c7793b']])out+=`<circle cx="${x(r.detection_pct)}" cy="${yy}" r="5" fill="${c}" stroke="#fffefa" stroke-width="1"><title>${escapeHtml(`${p.lib} · ${r.category}: ${format(r.detection_pct,3)}% (${r.n_detected}/${r.n_cells})`)}</title></circle>`;});
  return out+txt((left+w-right)/2,h-6,zh?'表達被檢出的細胞比例':'Cells with detected expression','text-anchor="middle" font-size="9"')+'</svg>';
}
export function compositionChart(data,lang) {
  const zh=lang==='zh';if(!data||!finite(data.total))return emptyChart(zh?'沒有可分解的單細胞資料。':'No measured single-cell decomposition is available.');
  const w=470,h=228,left=136,right=66,top=42,bottom=41,vals=[data.composition,data.within,data.total];
  let min=Math.min(0,...vals),max=Math.max(0,...vals),span=max-min;if(span===0)span=.001;min-=span*.12;max+=span*.18;const x=v=>left+(v-min)/(max-min)*(w-left-right);
  let out=svgOpen(w,h,zh?'Kitagawa描述性分解：組成項加類內項等於總差值':'Descriptive Kitagawa decomposition: composition plus within-category term equals total');
  out+=txt(left,17,'Δ mean log1p(CP10k)','font-size="9"');
  for(let i=0;i<=3;i++){const val=min+(max-min)*i/3;out+=line(x(val),top-8,x(val),h-bottom)+txt(x(val),h-23,tickFormat(val,max-min),'text-anchor="middle" font-size="8"');}
  out+=line(x(0),top-8,x(0),h-bottom,'style="stroke:#97aa99;stroke-width:1.2"');
  const labels=zh?['細胞大類組成','大類內表達','總差值']:['Category composition','Within categories','Total difference'];
  vals.forEach((v,i)=>{const yy=top+20+i*42,c=i===0?'#bf9c6c':i===1?'#7b9f82':'#14796f';const a=x(Math.min(0,v)),b=x(Math.max(0,v));out+=txt(left-11,yy+3,labels[i],'text-anchor="end" font-size="10"');out+=`<rect x="${a}" y="${yy-9}" width="${Math.max(.8,b-a)}" height="18" rx="2" fill="${c}"><title>${escapeHtml(`${labels[i]}: ${signed(v,7)}`)}</title></rect>`;out+=txt(w-right+10,yy+3,signed(v,4),'font-size="10" fill="#183b38"');});
  return out+txt(left,h-6,zh?'描述性差值；不是機制貢獻率':'Descriptive differences; not causal contribution rates','font-size="8"')+'</svg>';
}
export class UmapChart {
  constructor(canvas,tooltip,onInfo){this.canvas=canvas;this.ctx=canvas.getContext('2d');this.tooltip=tooltip;this.onInfo=onInfo;this.data=null;this.zoom=1;this.dx=0;this.dy=0;this.drag=null;this.screen=[];this.opts={gene:'VNN1',mode:'expression',library:'all',lang:'zh'};this.raf=0;this.install();this.observer=new ResizeObserver(()=>this.draw());this.observer.observe(canvas.parentElement);}
  install(){const c=this.canvas;c.addEventListener('wheel',e=>{e.preventDefault();this.zoomBy(e.deltaY<0?1.16:1/1.16,e.offsetX,e.offsetY);},{passive:false});c.addEventListener('pointerdown',e=>{this.drag={x:e.clientX,y:e.clientY,dx:this.dx,dy:this.dy,moved:false};c.setPointerCapture(e.pointerId);this.tooltip.classList.add('hidden');});c.addEventListener('pointermove',e=>{if(this.drag){const mx=e.clientX-this.drag.x,my=e.clientY-this.drag.y;this.dx=this.drag.dx+mx;this.dy=this.drag.dy+my;this.drag.moved=Math.abs(mx)+Math.abs(my)>3;this.draw();}else this.hover(e.offsetX,e.offsetY);});c.addEventListener('pointerup',e=>{if(this.drag&&!this.drag.moved)this.hover(e.offsetX,e.offsetY);this.drag=null;});c.addEventListener('pointercancel',()=>this.drag=null);c.addEventListener('pointerleave',()=>{if(!this.drag)this.tooltip.classList.add('hidden');});c.addEventListener('keydown',e=>{if(['+','=','-','0','ArrowLeft','ArrowRight','ArrowUp','ArrowDown'].includes(e.key)){e.preventDefault();if(e.key==='+'||e.key==='=')this.zoomBy(1.3);else if(e.key==='-')this.zoomBy(1/1.3);else if(e.key==='0')this.reset();else{this.dx+=e.key==='ArrowLeft'?30:e.key==='ArrowRight'?-30:0;this.dy+=e.key==='ArrowUp'?30:e.key==='ArrowDown'?-30:0;this.draw();}}});}
  setData(data){this.data=data;const xs=data.coordinates.x,ys=data.coordinates.y;this.bounds={xmin:Math.min(...xs),xmax:Math.max(...xs),ymin:Math.min(...ys),ymax:Math.max(...ys)};this.draw();}
  setOptions(opts){this.opts={...this.opts,...opts};this.tooltip.classList.add('hidden');this.draw();}
  reset(){this.zoom=1;this.dx=0;this.dy=0;this.draw();}
  zoomBy(factor,cx,cy){const w=this.canvas.clientWidth,h=this.canvas.clientHeight;cx??=w/2;cy??=h/2;const next=Math.min(12,Math.max(.65,this.zoom*factor));const ratio=next/this.zoom;this.dx=cx-w/2-(cx-w/2-this.dx)*ratio;this.dy=cy-h/2-(cy-h/2-this.dy)*ratio;this.zoom=next;this.draw();}
  draw(){if(this.raf)cancelAnimationFrame(this.raf);this.raf=requestAnimationFrame(()=>this.paint());}
  paint(){this.raf=0;if(!this.data)return;const c=this.canvas,w=c.clientWidth,h=c.clientHeight;if(!w||!h)return;const dpr=Math.min(window.devicePixelRatio||1,2);c.width=Math.round(w*dpr);c.height=Math.round(h*dpr);const ctx=this.ctx;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.fillStyle='#fdfdf9';ctx.fillRect(0,0,w,h);const b=this.bounds,baseScale=Math.min((w-82)/(b.xmax-b.xmin),(h-76)/(b.ymax-b.ymin)),scale=baseScale*this.zoom;const tx=x=>w/2+this.dx+(x-(b.xmax+b.xmin)/2)*scale,ty=y=>h/2+this.dy-(y-(b.ymax+b.ymin)/2)*scale;const {data,opts}=this,values=data.expression[opts.gene],mode=opts.mode,libIndex=data.categories.library.indexOf(opts.library);let max=0,detected=0,visible=0,missing=0;const filtered=[];for(let i=0;i<data.n_cells;i++){if(libIndex>=0&&data.codes.library[i]!==libIndex)continue;filtered.push(i);visible++;const val=values?.[i];if(finite(val)){if(val>0)detected++;max=Math.max(max,val);}else missing++;}this.screen=[];const radius=Math.min(3,1.45*Math.sqrt(this.zoom));const records=filtered.map(i=>({i,x:tx(data.coordinates.x[i]),y:ty(data.coordinates.y[i]),v:values?.[i]}));if(mode==='expression')records.sort((a,b)=>(a.v??-1)-(b.v??-1));for(const p of records){if(p.x<0||p.x>w||p.y<0||p.y>h)continue;let color;if(mode==='expression')color=expressionColor(p.v,max);else color=colors[data.codes[mode][p.i]%colors.length];ctx.fillStyle=color;ctx.beginPath();ctx.arc(p.x,p.y,radius,0,Math.PI*2);ctx.fill();this.screen.push(p);}ctx.strokeStyle='#819485';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(28,h-52);ctx.lineTo(73,h-52);ctx.moveTo(28,h-52);ctx.lineTo(28,h-97);ctx.stroke();ctx.fillStyle='#6d8272';ctx.font='9px Segoe UI,Arial';ctx.fillText('UMAP 1',29,h-37);ctx.save();ctx.translate(17,h-51);ctx.rotate(-Math.PI/2);ctx.fillText('UMAP 2',0,0);ctx.restore();this.onInfo?.({visible,detected,missing,max,mode});}
  hover(x,y){let nearest=null,best=49;for(const p of this.screen){const d=(p.x-x)**2+(p.y-y)**2;if(d<best){best=d;nearest=p;}}if(!nearest){this.tooltip.classList.add('hidden');return;}const i=nearest.i,d=this.data,zh=this.opts.lang==='zh',label=key=>d.categories[key][d.codes[key][i]];this.tooltip.innerHTML=`<strong>${escapeHtml(label('subtype'))}</strong><br>${escapeHtml(label('library'))} · ${escapeHtml(label('condition'))}<br>${escapeHtml(this.opts.gene)}: ${finite(nearest.v)?format(nearest.v,4):'NA'}<br><span style="opacity:.65">${zh?'細胞索引':'Cell index'} ${i+1} · log1p(CP10k)</span>`;this.tooltip.classList.remove('hidden');const tw=this.tooltip.offsetWidth,th=this.tooltip.offsetHeight;this.tooltip.style.left=`${Math.max(5,Math.min(x+13,this.canvas.clientWidth-tw-8))}px`;this.tooltip.style.top=`${Math.max(5,Math.min(y+13,this.canvas.clientHeight-th-8))}px`;}
  download(){
    if(!this.data||!this.canvas.clientWidth||!this.canvas.clientHeight)return;
    // Flush the pending frame so selected labels cannot accompany the preceding gene's plot.
    if(this.raf)cancelAnimationFrame(this.raf);this.raf=0;this.paint();
    const out=document.createElement('canvas'),w=1600,cssW=this.canvas.clientWidth,cssH=this.canvas.clientHeight,k=w/cssW,plotH=Math.round(cssH*k),header=154;
    out.width=w;out.height=header+plotH+600;const ctx=out.getContext('2d');ctx.fillStyle='#fffefa';ctx.fillRect(0,0,out.width,out.height);
    ctx.fillStyle=INK;ctx.font='600 30px Segoe UI,Arial';ctx.fillText(`Sesamin Atlas · ${this.opts.gene}`,42,47);
    const modeNames={expression:'Gene expression',subtype:'Author cell subtype',broad:'Broad cell category',condition:'Condition'},library=this.opts.library==='all'?'All four technical libraries':this.opts.library;
    ctx.font='20px Segoe UI,Arial';ctx.fillText(`${modeNames[this.opts.mode]} | ${library} | Current viewport`,42,85);
    ctx.fillStyle=MUTED;ctx.font='18px Segoe UI,Arial';ctx.fillText('GSE207889 · Author-saved UMAP · Human liver organoid reference',42,120);
    const values=this.data.expression[this.opts.gene],libIndex=this.data.categories.library.indexOf(this.opts.library);let max=0,count=0,detected=0,missing=0;
    for(let i=0;i<this.data.n_cells;i++){if(libIndex>=0&&this.data.codes.library[i]!==libIndex)continue;count++;const value=values?.[i];if(!finite(value))missing++;else{max=Math.max(max,value);if(value>0)detected++;}}
    // Draw every visible point again at export resolution. Do not enlarge a phone screenshot.
    ctx.save();ctx.translate(0,header);ctx.scale(k,k);ctx.beginPath();ctx.rect(0,0,cssW,cssH);ctx.clip();ctx.fillStyle='#fdfdf9';ctx.fillRect(0,0,cssW,cssH);
    const radius=Math.min(3,1.45*Math.sqrt(this.zoom));
    for(const p of this.screen){ctx.fillStyle=this.opts.mode==='expression'?expressionColor(p.v,max):colors[this.data.codes[this.opts.mode][p.i]%colors.length];ctx.beginPath();ctx.arc(p.x,p.y,radius,0,Math.PI*2);ctx.fill();}
    ctx.strokeStyle='#819485';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(28,cssH-52);ctx.lineTo(73,cssH-52);ctx.moveTo(28,cssH-52);ctx.lineTo(28,cssH-97);ctx.stroke();ctx.fillStyle=MUTED;ctx.font='9px Segoe UI,Arial';ctx.fillText('UMAP 1',29,cssH-37);ctx.save();ctx.translate(17,cssH-51);ctx.rotate(-Math.PI/2);ctx.fillText('UMAP 2',0,0);ctx.restore();ctx.restore();
    let yy=header+plotH+35;ctx.font='18px Segoe UI,Arial';
    const note=text=>{let part='';for(const char of text){if(ctx.measureText(part+char).width>w-84){ctx.fillText(part,42,yy);yy+=25;part=char;}else part+=char;}if(part){ctx.fillText(part,42,yy);yy+=25;}};
    if(this.opts.mode==='expression'){
      if(missing===count){ctx.fillStyle='#c4c4bd';ctx.fillRect(42,yy-15,18,18);ctx.fillStyle=INK;ctx.fillText('NA: this gene feature is unavailable, not zero expression.',75,yy);}
      else if(detected===0){ctx.fillStyle='#e0e5d9';ctx.fillRect(42,yy-15,18,18);ctx.fillStyle=INK;ctx.fillText('Measured zero in all selected cells; no positive expression values.',75,yy);}
      else{ctx.fillStyle='#e0e5d9';ctx.fillRect(42,yy-15,18,18);ctx.fillStyle=INK;ctx.fillText('Measured zero',75,yy);const grad=ctx.createLinearGradient(290,0,590,0);grad.addColorStop(0,'#a3caa6');grad.addColorStop(1,'#14574a');ctx.fillStyle=grad;ctx.fillRect(290,yy-15,300,18);ctx.fillStyle=INK;ctx.fillText(`Detected > 0 to ${format(max,2)} · log1p(CP10k) · linear color scale`,610,yy);}
      yy+=31;ctx.fillStyle=MUTED;note(missing===count?`${count.toLocaleString('en-US')} cells selected; feature unavailable in all selected cells.`:`${count.toLocaleString('en-US')} cells selected; ${detected.toLocaleString('en-US')} with detected expression. Counts cover selected libraries, not the zoomed viewport.`);
    }else{
      let xx=42;for(const [i,label]of this.data.categories[this.opts.mode].entries()){const sw=ctx.measureText(label).width+50;if(xx+sw>w-42){xx=42;yy+=28;}ctx.fillStyle=colors[i%colors.length];ctx.beginPath();ctx.arc(xx+6,yy-6,6,0,Math.PI*2);ctx.fill();ctx.fillStyle=INK;ctx.fillText(label,xx+23,yy);xx+=sw;}yy+=31;ctx.fillStyle=MUTED;note(`${count.toLocaleString('en-US')} cells selected. Counts cover selected libraries, not the zoomed viewport. Source category labels retained.`);
    }
    yy+=13;ctx.fillStyle=MUTED;ctx.font='17px Segoe UI,Arial';
    note('No sesamin treatment. Source: four technical libraries from one differentiation; no independent biological replicates.');
    note('UMAP is not spatial transcriptomics. Distances are not effect sizes or statistical significance.');
    note(`App ${this.opts.appVersion||'1.1.0'} | Data ${this.opts.dataVersion||'2026-10-09-analysis-strengthening-15'} | Exported ${new Date().toISOString().slice(0,10)} UTC`);
    note('https://fc28-sesamin-atlas.static.hf.space/'+location.hash);
    note('The link restores selected options; this image records the current pan and zoom separately.');
    const final=document.createElement('canvas');final.width=w;final.height=Math.min(out.height,Math.ceil(yy+20));final.getContext('2d').drawImage(out,0,0);
    const link=document.createElement('a');link.download=`sesamin-atlas_${this.opts.gene}_${this.opts.mode}_${this.opts.library}_v${this.opts.appVersion||'1.1.0'}.png`;link.href=final.toDataURL('image/png');link.click();
  }
}
