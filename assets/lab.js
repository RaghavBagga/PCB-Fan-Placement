/* PCB Fan Placement Lab: finite-volume heat-equation solver running in the browser.
   Same discretisation as simulation.py: 2 mm cells, 5-point stencil, insulated edges, PCG solve. */
(function(){
"use strict";
const NX=50, NY=40, DX=0.002, N=NX*NY;
const DEF={dh:60,sig:12,lim:60,p:[2,1,0.8],fanOn:true,fx:38,fy:50};
const S={k:15,t:0.0016,hnat:10,Ta:25,dh:60,sig:12,lim:60,fanOn:true,fx:38,fy:50,model:"gauss",showFlux:false,
  chips:[{n:"CPU",x:30,y:50,w:12,d:12,P:2},{n:"Regulator",x:70,y:55,w:8,d:8,P:1},{n:"Driver",x:65,y:20,w:10,d:6,P:0.8}]};
const $=id=>document.getElementById(id);
const nnCount=new Uint8Array(N);
for(let i=0;i<NX;i++)for(let j=0;j<NY;j++){nnCount[i*NY+j]=(i>0)+(i<NX-1)+(j>0)+(j<NY-1);}

/* ---------- fan profile: Gaussian or truncated Taylor polynomial of exp(-u) ---------- */
const ucut={};
function taylor(u,m){let s=0,term=1;for(let k=0;k<=m;k++){s+=term;term*=-u/(k+1);}return s;}
function cutoffU(m){ if(ucut[m]!==undefined)return ucut[m];
  let u=0;while(u<30&&Math.abs(taylor(u,m)-Math.exp(-u))<=0.02)u+=0.0005;return ucut[m]=u;}
function profile(r2){
  const u=r2/(2*S.sig*S.sig);
  if(S.model==="gauss")return Math.exp(-u);
  const m=+S.model;return u<=cutoffU(m)?taylor(u,m):0;
}

/* ---------- solver ---------- */
let Q=new Float64Array(N);
function buildQ(){
  Q.fill(0);
  for(const c of S.chips){
    const cells=[];
    for(let i=0;i<NX;i++)for(let j=0;j<NY;j++){
      // same arithmetic as simulation.py (metres, IEEE doubles) so cells exactly on a chip edge are included/excluded identically
      const x=(i+.5)*DX,y=(j+.5)*DX;
      if(Math.abs(x-c.x/1000)<=(c.w/1000)/2&&Math.abs(y-c.y/1000)<=(c.d/1000)/2)cells.push(i*NY+j);
    }
    if(!cells.length)continue;
    const per=c.P/(cells.length*DX*DX);
    for(const k of cells)Q[k]+=per;
  }
}
function buildH(on,fx,fy){
  fx=fx===undefined?S.fx:fx;fy=fy===undefined?S.fy:fy;
  const h=new Float64Array(N);
  for(let i=0;i<NX;i++){const x=(i+.5)*2;
    for(let j=0;j<NY;j++){const y=(j+.5)*2;
      h[i*NY+j]=S.hnat+(on?S.dh*profile((x-fx)**2+(y-fy)**2):0);}}
  return h;
}
const bv=new Float64Array(N), rr=new Float64Array(N), zz=new Float64Array(N), pp=new Float64Array(N), Ap=new Float64Array(N), dg=new Float64Array(N);
function solveCG(h,x){           // theta = T - Ta, warm-started from x
  const kt=S.k*S.t, d2=DX*DX;
  for(let m=0;m<N;m++){dg[m]=kt*nnCount[m]+2*h[m]*d2;bv[m]=Q[m]*d2;}
  const apply=(v,out)=>{
    for(let i=0;i<NX;i++){const base=i*NY;
      for(let j=0;j<NY;j++){const m=base+j,c=v[m];let s=0;
        if(i>0)s+=c-v[m-NY]; if(i<NX-1)s+=c-v[m+NY]; if(j>0)s+=c-v[m-1]; if(j<NY-1)s+=c-v[m+1];
        out[m]=kt*s+2*h[m]*d2*c;}}
  };
  apply(x,Ap);
  let bn=0,rz=0;
  for(let m=0;m<N;m++){rr[m]=bv[m]-Ap[m];bn+=bv[m]*bv[m];zz[m]=rr[m]/dg[m];pp[m]=zz[m];rz+=rr[m]*zz[m];}
  bn=Math.sqrt(bn)||1;
  for(let it=0;it<800;it++){
    let rn=0;for(let m=0;m<N;m++)rn+=rr[m]*rr[m];
    if(Math.sqrt(rn)<1e-9*bn)break;
    apply(pp,Ap);
    let pAp=0;for(let m=0;m<N;m++)pAp+=pp[m]*Ap[m];
    const a=rz/pAp;let rzn=0;
    for(let m=0;m<N;m++){x[m]+=a*pp[m];rr[m]-=a*Ap[m];zz[m]=rr[m]/dg[m];rzn+=rr[m]*zz[m];}
    const b=rzn/rz;rz=rzn;
    for(let m=0;m<N;m++)pp[m]=zz[m]+b*pp[m];
  }
  return x;
}
function measure(th,h){
  let mx=-1,mi=0,sum=0,cnt=0,out=0;
  for(let m=0;m<N;m++){const v=th[m];if(v>mx){mx=v;mi=m;}sum+=v;if(v+S.Ta>S.lim)cnt++;out+=2*h[m]*v*DX*DX;}
  let inn=0;for(let m=0;m<N;m++)inn+=Q[m]*DX*DX;
  return {Tmax:mx+S.Ta,mi,hx:(Math.floor(mi/NY)+.5)*2,hy:((mi%NY)+.5)*2,area:cnt*4,mean:sum/N+S.Ta,heatIn:inn,heatOut:out};
}
const scratch=new Float64Array(N);
function tmaxAt(fx,fy){const h=buildH(true,fx,fy);solveCG(h,scratch);let mx=-1;for(let m=0;m<N;m++)if(scratch[m]>mx)mx=scratch[m];return mx+S.Ta;}

/* ---------- colour ---------- */
const INF=[[0,0,4],[87,16,110],[188,55,84],[249,142,9],[252,255,164]];
const VIR=[[68,1,84],[59,82,139],[33,145,140],[94,201,98],[253,231,37]];
function cmap(st,u){u=Math.max(0,Math.min(1,u))*(st.length-1);const i=Math.min(st.length-2,Math.floor(u)),f=u-i,a=st[i],b=st[i+1];
  return [a[0]+(b[0]-a[0])*f,a[1]+(b[1]-a[1])*f,a[2]+(b[2]-a[2])*f];}
$("bar").style.background="linear-gradient(90deg,"+INF.map((c,i)=>`rgb(${c}) ${i*25}%`).join(",")+")";

/* ---------- drawing ---------- */
const heat=$("heat"),hctx=heat.getContext("2d");
const off=document.createElement("canvas");off.width=NX;off.height=NY;const octx=off.getContext("2d");
function fitCanvas(cv){const r=cv.getBoundingClientRect(),d=Math.min(2,window.devicePixelRatio||1);
  const w=Math.max(1,Math.round(r.width*d)),h=Math.max(1,Math.round(r.height*d));
  if(cv.width!==w||cv.height!==h){cv.width=w;cv.height=h;}return [w,h];}
function drawContour(th,W,H){
  hctx.strokeStyle="#6fe3ff";hctx.lineWidth=Math.max(1.5,W/450);hctx.beginPath();
  const L=S.lim-S.Ta,px=i=>((i+.5)*2/100)*W,py=j=>H-((j+.5)*2/80)*H;
  const v=(i,j)=>th[i*NY+j];
  for(let i=0;i<NX-1;i++)for(let j=0;j<NY-1;j++){
    const c=[v(i,j),v(i+1,j),v(i+1,j+1),v(i,j+1)];
    const idx=(c[0]>L?1:0)|(c[1]>L?2:0)|(c[2]>L?4:0)|(c[3]>L?8:0);
    if(idx===0||idx===15)continue;
    const X=[i,i+1,i+1,i],Y=[j,j,j+1,j+1];
    const edge=e=>{const a=e,b=(e+1)%4;
      const t=(L-c[a])/(c[b]-c[a]);return [px(X[a]+(X[b]-X[a])*t),py(Y[a]+(Y[b]-Y[a])*t)];};
    const T={1:[[3,0]],2:[[0,1]],3:[[3,1]],4:[[1,2]],5:[[3,2],[0,1]],6:[[0,2]],7:[[3,2]],8:[[2,3]],9:[[0,2]],10:[[0,3],[1,2]],11:[[1,2]],12:[[3,1]],13:[[0,1]],14:[[3,0]]}[idx];
    for(const [e1,e2] of T){const p=edge(e1),q=edge(e2);hctx.moveTo(p[0],p[1]);hctx.lineTo(q[0],q[1]);}
  }
  hctx.stroke();
}
function arrow(x0,y0,dx,dy,len,head,col,lw){
  const m=Math.hypot(dx,dy);if(m<1e-12)return;dx/=m;dy/=m;
  const x1=x0+dx*len,y1=y0+dy*len;hctx.strokeStyle=col;hctx.fillStyle=col;hctx.lineWidth=lw;
  hctx.beginPath();hctx.moveTo(x0,y0);hctx.lineTo(x1,y1);hctx.stroke();
  hctx.beginPath();hctx.moveTo(x1,y1);hctx.lineTo(x1-head*(dx*0.9-dy*0.5),y1-head*(dy*0.9+dx*0.5));
  hctx.lineTo(x1-head*(dx*0.9+dy*0.5),y1-head*(dy*0.9-dx*0.5));hctx.closePath();hctx.fill();
}
function drawFlux(th,W,H){            // F = -kt grad T : direction of heat flow, brightness ~ |F|
  const g=(i,j)=>th[Math.max(0,Math.min(NX-1,i))*NY+Math.max(0,Math.min(NY-1,j))];
  let mx=0;const pts=[];
  for(let i=2;i<NX;i+=4)for(let j=2;j<NY;j+=4){
    const gx=(g(i+1,j)-g(i-1,j))/(2*DX*1000),gy=(g(i,j+1)-g(i,j-1))/(2*DX*1000);   // K/mm
    const m=Math.hypot(gx,gy);pts.push([i,j,-gx,-gy,m]);if(m>mx)mx=m;}
  const L=W/42;
  for(const [i,j,fx,fy,m] of pts){
    if(m<0.02*mx)continue;
    const a=0.25+0.75*Math.sqrt(m/mx);
    arrow(((i+.5)*2/100)*W,H-((j+.5)*2/80)*H,fx,-fy,L,L*0.38,`rgba(255,226,120,${a.toFixed(2)})`,Math.max(1,W/700));
  }
}
let trail=[],gradVec=null;
function drawBoard(th,m,hi){
  const [W,H]=fitCanvas(heat);
  const img=octx.createImageData(NX,NY);
  for(let i=0;i<NX;i++)for(let j=0;j<NY;j++){
    const u=(th[i*NY+j])/(hi-S.Ta),c=cmap(INF,u),o=((NY-1-j)*NX+i)*4;
    img.data[o]=c[0];img.data[o+1]=c[1];img.data[o+2]=c[2];img.data[o+3]=255;}
  octx.putImageData(img,0,0);
  hctx.imageSmoothingEnabled=true;hctx.imageSmoothingQuality="high";
  hctx.drawImage(off,0,0,W,H);
  drawContour(th,W,H);
  if(S.showFlux)drawFlux(th,W,H);
  const x=(m.hx/100)*W,y=H-(m.hy/80)*H,r=Math.max(6,W/70);
  hctx.strokeStyle="#fff";hctx.lineWidth=Math.max(1.5,W/450);
  hctx.beginPath();hctx.arc(x,y,r,0,7);hctx.moveTo(x-r*1.5,y);hctx.lineTo(x+r*1.5,y);hctx.moveTo(x,y-r*1.5);hctx.lineTo(x,y+r*1.5);hctx.stroke();
  if(trail.length>1){hctx.strokeStyle="rgba(111,227,255,.9)";hctx.setLineDash([4,3]);hctx.lineWidth=Math.max(1.5,W/450);hctx.beginPath();
    trail.forEach((p,k)=>{const px=p[0]/100*W,py=H-p[1]/80*H;k?hctx.lineTo(px,py):hctx.moveTo(px,py);});hctx.stroke();hctx.setLineDash([]);}
  if(gradVec&&S.fanOn){const fxp=S.fx/100*W,fyp=H-S.fy/80*H;       // steepest-descent direction of Tmax
    arrow(fxp,fyp,-gradVec[0],gradVec[1],W/16,W/55,"#ffd24a",Math.max(2,W/260));}
}
function layoutOverlay(){
  const box=$("chips");
  if(!box.children.length){
    for(const c of S.chips){const d=document.createElement("div");d.className="chip";d.innerHTML="<span></span>";box.appendChild(d);}
  }
  S.chips.forEach((c,k)=>{const d=box.children[k];
    d.style.left=c.x+"%";d.style.top=(100-c.y/80*100)+"%";d.style.width=c.w+"%";d.style.height=(c.d/80*100)+"%";
    d.firstChild.textContent=c.n+" "+c.P.toFixed(1)+" W";});
  const f=$("fan");
  f.style.left=S.fx+"%";f.style.top=(100-S.fy/80*100)+"%";f.style.width=(2*S.sig)+"%";
  f.classList.toggle("off",!S.fanOn);
  f.setAttribute("aria-valuetext","x "+S.fx.toFixed(1)+" mm, y "+S.fy.toFixed(1)+" mm");
}

/* ---------- main update ---------- */
let base=null,baseTh=new Float64Array(N),fanTh=new Float64Array(N),baseDirty=true,pending=false;
function fmt(v,d=1){return v.toFixed(d)}
function delta(el,v,unit,lowerIsGood=true,d=1){
  if(Math.abs(v)<0.05&&d===1){el.textContent="no change vs no fan";el.className="d";return;}
  el.textContent=(v>0?"+":"−")+Math.abs(v).toFixed(d)+" "+unit+" vs no fan";
  el.className="d "+((v<0)===lowerIsGood?"good":"bad");
}
function update(){
  pending=false;
  if(baseDirty){buildQ();const h0=buildH(false);solveCG(h0,baseTh);base=measure(baseTh,h0);baseDirty=false;}
  const h=buildH(S.fanOn);
  solveCG(h,fanTh);
  const m=measure(fanTh,h);
  // gradient of Tmax w.r.t. fan position by central differences (1 mm)
  if(S.fanOn){const e=1;
    gradVec=[(tmaxAt(S.fx+e,S.fy)-tmaxAt(S.fx-e,S.fy))/(2*e),(tmaxAt(S.fx,S.fy+e)-tmaxAt(S.fx,S.fy-e))/(2*e)];}
  else gradVec=null;
  const hi=Math.ceil((base.Tmax+1)/5)*5;
  $("hi").textContent=hi+" °C";$("lo").textContent=S.Ta+" °C";
  drawBoard(fanTh,m,hi);
  layoutOverlay();
  $("mTmax").innerHTML=fmt(m.Tmax)+"<small>°C</small>";
  delta($("dTmax"),m.Tmax-base.Tmax,"°C");
  $("mArea").innerHTML=m.area+"<small>mm²</small>";
  delta($("dArea"),m.area-base.area,"mm²",true,0);
  $("mHot").textContent="("+fmt(m.hx,0)+", "+fmt(m.hy,0)+") mm";
  let who="",bd=1e9;S.chips.forEach(c=>{const d=Math.hypot(c.x-m.hx,c.y-m.hy);if(d<bd){bd=d;who=c.n;}});
  $("mHotWho").textContent="at the "+who;
  $("mMean").innerHTML=fmt(m.mean)+"<small>°C</small>";
  delta($("dMean"),m.mean-base.mean,"°C");
  $("balance").innerHTML="<span>Heat in "+fmt(m.heatIn,3)+" W · heat removed "+fmt(m.heatOut,3)+" W · mismatch "+Math.abs(m.heatOut-m.heatIn).toExponential(1)+" W</span>"+
    (gradVec?"<span>∇T<sub>max</sub> w.r.t. fan position = ("+gradVec[0].toFixed(3)+", "+gradVec[1].toFixed(3)+") K/mm, yellow arrow = steepest descent</span>":"");
  drawScan();
}
function request(dirty){if(dirty)baseDirty=true;if(!pending){pending=true;requestAnimationFrame(update);}}

/* ---------- inputs ---------- */
const board=$("board");let drag=false;
function stopOpt(){optRunning=false;}
function placeFan(e){stopOpt();trail=[];const r=board.getBoundingClientRect();
  S.fx=Math.max(0,Math.min(100,(e.clientX-r.left)/r.width*100));
  S.fy=Math.max(0,Math.min(80,(1-(e.clientY-r.top)/r.height)*80));
  S.fx=Math.round(S.fx*2)/2;S.fy=Math.round(S.fy*2)/2;request(false);}
board.addEventListener("pointerdown",e=>{drag=true;try{board.setPointerCapture(e.pointerId);}catch(_){}placeFan(e);});
board.addEventListener("pointermove",e=>{if(drag)placeFan(e);});
board.addEventListener("pointerup",()=>{drag=false;});
board.addEventListener("pointercancel",()=>{drag=false;});
$("fan").addEventListener("keydown",e=>{
  const s=e.shiftKey?10:2;let used=true;stopOpt();trail=[];
  if(e.key==="ArrowLeft")S.fx=Math.max(0,S.fx-s);else if(e.key==="ArrowRight")S.fx=Math.min(100,S.fx+s);
  else if(e.key==="ArrowUp")S.fy=Math.min(80,S.fy+s);else if(e.key==="ArrowDown")S.fy=Math.max(0,S.fy-s);else used=false;
  if(used){e.preventDefault();request(false);}
});
function bindRange(id,outId,unit,apply,dec=0){
  const el=$(id),o=$(outId);
  const sync=()=>{o.textContent=(+el.value).toFixed(dec)+" "+unit;};
  el.addEventListener("input",()=>{apply(+el.value);sync();onParam();});sync();
}
function onParam(){scan=null;trail=[];$("scanBtn").disabled=false;updateScanUI();request(true);}
bindRange("dh","dhO","W/m²K",v=>S.dh=v);
bindRange("sig","sigO","mm",v=>S.sig=v);
bindRange("lim","limO","°C",v=>S.lim=v);
S.chips.forEach((c,k)=>bindRange("p"+k,"p"+k+"O","W",v=>c.P=v,1));
$("fanOn").addEventListener("change",e=>{S.fanOn=e.target.checked;request(false);});
$("fluxOn").addEventListener("change",e=>{S.showFlux=e.target.checked;request(false);});
$("model").addEventListener("change",e=>{S.model=e.target.value;onParam();});
$("reset").addEventListener("click",()=>{
  stopOpt();S.dh=DEF.dh;S.sig=DEF.sig;S.lim=DEF.lim;S.fanOn=true;S.fx=DEF.fx;S.fy=DEF.fy;S.model="gauss";$("model").value="gauss";S.chips.forEach((c,k)=>c.P=DEF.p[k]);
  $("dh").value=S.dh;$("sig").value=S.sig;$("lim").value=S.lim;$("fanOn").checked=true;
  S.chips.forEach((c,k)=>{$("p"+k).value=c.P;});
  ["dh","sig","lim","p0","p1","p2"].forEach(id=>$(id).dispatchEvent(new Event("input")));
});

/* ---------- live optimiser: compass (pattern) search on Tmax(x_f, y_f) ---------- */
let optRunning=false;
$("optBtn").addEventListener("click",()=>{
  if(optRunning){optRunning=false;return;}
  if(!S.fanOn){$("optOut").textContent="Turn the fan on first.";return;}
  optRunning=true;$("optBtn").textContent="Stop";trail=[[S.fx,S.fy]];
  let step=10,f=tmaxAt(S.fx,S.fy),evals=1,iters=0;
  const dirs=[[1,0],[-1,0],[0,1],[0,-1],[.7071,.7071],[-.7071,.7071],[.7071,-.7071],[-.7071,-.7071]];
  const go=()=>{
    if(!optRunning){$("optBtn").textContent="Optimise fan position";return;}
    let moved=false;
    for(const d of dirs){
      const nx=S.fx+step*d[0],ny=S.fy+step*d[1];
      if(nx<0||nx>100||ny<0||ny>80)continue;
      const g=tmaxAt(nx,ny);evals++;
      if(g<f-1e-6){f=g;S.fx=nx;S.fy=ny;moved=true;trail.push([nx,ny]);break;}
    }
    iters++;
    if(!moved){step/=2;}
    $("optOut").textContent="Compass search: step "+step.toFixed(2)+" mm, fan at ("+S.fx.toFixed(2)+", "+S.fy.toFixed(2)+") mm, Tmax "+f.toFixed(3)+" °C, "+evals+" solves";
    request(false);
    if(step<0.1||iters>400){optRunning=false;$("optBtn").textContent="Optimise fan position";
      $("optOut").textContent+=" (converged). Compare with the scan result and the report: optimum ≈ (37.7, 50.25) mm for the default fan.";return;}
    setTimeout(go,30);
  };
  go();
});

/* ---------- scan ---------- */
let scan=null,scanning=false;
const scv=$("scanCv"),sctx=scv.getContext("2d");
const XS=[],YS=[];for(let x=0;x<=100;x+=4)XS.push(x);for(let y=0;y<=80;y+=4)YS.push(y);
const soff=document.createElement("canvas");soff.width=XS.length;soff.height=YS.length;
function objective(){return document.querySelector('input[name=obj]:checked').value;}
function bestOf(sc){let bi=0;
  for(let n=1;n<sc.tmax.length;n++){
    const a=objective()==="tmax"?[sc.tmax[n],sc.area[n]]:[sc.area[n],sc.tmax[n]];
    const b=objective()==="tmax"?[sc.tmax[bi],sc.area[bi]]:[sc.area[bi],sc.tmax[bi]];
    if(a[0]<b[0]-1e-9||(Math.abs(a[0]-b[0])<=1e-9&&a[1]<b[1]))bi=n;}
  return bi;}
function updateScanUI(){
  $("scanEmpty").hidden=!!scan;
  $("goBest").disabled=$("csvBtn").disabled=!scan;
  if(!scan){$("bestOut").textContent=scanning?"Scanning…":"No scan yet. Changing a model control clears the scan.";$("scanRead").textContent="";}
  else{
    const b=bestOf(scan),xi=Math.floor(b/YS.length),yi=b%YS.length;
    $("bestOut").innerHTML='<div class="big">Fan at ('+XS[xi]+', '+YS[yi]+') mm</div><div>Max temperature '+scan.tmax[b].toFixed(1)+' °C, area above '+S.lim+' °C: '+scan.area[b]+' mm². With no fan: '+scan.baseT.toFixed(1)+' °C and '+scan.baseA+' mm².</div>';
  }
  drawScan();
}
function drawScan(){
  const [W,H]=fitCanvas(scv);sctx.clearRect(0,0,W,H);
  if(!scan)return;
  const nxs=XS.length,nys=YS.length,img=soff.getContext("2d").createImageData(nxs,nys);
  let lo=1e9,hi=-1e9;for(const v of scan.tmax){if(v<lo)lo=v;if(v>hi)hi=v;}
  for(let xi=0;xi<nxs;xi++)for(let yi=0;yi<nys;yi++){
    const c=cmap(VIR,(scan.tmax[xi*nys+yi]-lo)/(hi-lo||1)),o=((nys-1-yi)*nxs+xi)*4;
    img.data[o]=c[0];img.data[o+1]=c[1];img.data[o+2]=c[2];img.data[o+3]=255;}
  soff.getContext("2d").putImageData(img,0,0);
  sctx.imageSmoothingEnabled=true;sctx.drawImage(soff,0,0,W,H);
  const pos=(x,y)=>[x/100*W,H-y/80*H],u=Math.max(5,W/90);
  for(const c of S.chips){const [x,y]=pos(c.x,c.y);sctx.strokeStyle="rgba(255,255,255,.85)";sctx.lineWidth=1.5;sctx.strokeRect(x-c.w/200*W,y-c.d/160*H,c.w/100*W,c.d/80*H);}
  const b=bestOf(scan),[bx,by]=pos(XS[Math.floor(b/YS.length)],YS[b%YS.length]);
  sctx.fillStyle="#ff4d4d";sctx.strokeStyle="#fff";sctx.lineWidth=1.5;sctx.beginPath();
  for(let k=0;k<10;k++){const r=k%2?u*.5:u*1.3,a=-Math.PI/2+k*Math.PI/5;sctx.lineTo(bx+r*Math.cos(a),by+r*Math.sin(a));}
  sctx.closePath();sctx.fill();sctx.stroke();
  const [fx,fy]=pos(S.fx,S.fy);sctx.strokeStyle="#6fe3ff";sctx.lineWidth=2;sctx.beginPath();sctx.arc(fx,fy,u,0,7);sctx.stroke();
  $("scanLegend").textContent="Colour: max T from "+lo.toFixed(1)+" °C (dark) to "+hi.toFixed(1)+" °C (light). Red star: best. Cyan ring: current fan.";
}
function startScan(){
  if(scanning)return;stopOpt();scanning=true;scan=null;$("scanBtn").disabled=true;$("prog").hidden=false;$("prog").value=0;updateScanUI();
  buildQ();
  const nys=YS.length,res={tmax:new Float64Array(XS.length*nys),area:new Float64Array(XS.length*nys),hx:new Float64Array(XS.length*nys),hy:new Float64Array(XS.length*nys),baseT:base.Tmax,baseA:base.area};
  const saveX=S.fx,saveY=S.fy,th=new Float64Array(N);let xi=0;
  const step=()=>{
    const t0=performance.now();
    while(xi<XS.length&&performance.now()-t0<25){
      for(let yi=0;yi<nys;yi++){S.fx=XS[xi];S.fy=YS[yi];const h=buildH(true);solveCG(h,th);const m=measure(th,h),n=xi*nys+yi;
        res.tmax[n]=m.Tmax;res.area[n]=m.area;res.hx[n]=m.hx;res.hy[n]=m.hy;}
      xi++;
    }
    $("prog").value=xi/XS.length;
    if(xi<XS.length){setTimeout(step,0);return;}
    S.fx=saveX;S.fy=saveY;scanning=false;scan=res;$("prog").hidden=true;$("scanBtn").disabled=false;updateScanUI();
  };
  setTimeout(step,0);
}
$("scanBtn").addEventListener("click",startScan);
document.querySelectorAll('input[name=obj]').forEach(r=>r.addEventListener("change",updateScanUI));
$("goBest").addEventListener("click",()=>{if(!scan)return;const b=bestOf(scan);S.fx=XS[Math.floor(b/YS.length)];S.fy=YS[b%YS.length];trail=[];request(false);});
function scanPick(e){const r=scv.getBoundingClientRect();
  const x=Math.max(0,Math.min(100,(e.clientX-r.left)/r.width*100)),y=Math.max(0,Math.min(80,(1-(e.clientY-r.top)/r.height)*80));
  return [x,y];}
scv.addEventListener("click",e=>{if(!scan)return;const [x,y]=scanPick(e);S.fx=Math.round(x*2)/2;S.fy=Math.round(y*2)/2;trail=[];request(false);});
scv.addEventListener("pointermove",e=>{if(!scan)return;const [x,y]=scanPick(e);
  const xi=Math.min(XS.length-1,Math.round(x/4)),yi=Math.min(YS.length-1,Math.round(y/4)),n=xi*YS.length+yi;
  $("scanRead").textContent="Fan at ("+XS[xi]+", "+YS[yi]+") mm → max "+scan.tmax[n].toFixed(1)+" °C, "+scan.area[n]+" mm² above limit";});
$("csvBtn").addEventListener("click",()=>{
  if(!scan)return;const rows=["x_fan_mm,y_fan_mm,Tmax_C,x_hot_mm,y_hot_mm,area_over_limit_mm2"];
  for(let xi=0;xi<XS.length;xi++)for(let yi=0;yi<YS.length;yi++){const n=xi*YS.length+yi;
    rows.push([XS[xi],YS[yi],scan.tmax[n].toFixed(3),scan.hx[n],scan.hy[n],scan.area[n]].join(","));}
  const blob=new Blob([rows.join("\n")],{type:"text/csv"}),a=document.createElement("a");
  a.href=URL.createObjectURL(blob);a.download="fan_scan_limit_"+S.lim+"C.csv";document.body.appendChild(a);a.click();a.remove();
  setTimeout(()=>URL.revokeObjectURL(a.href),2000);
});

/* ---------- start ---------- */
window.addEventListener("resize",()=>request(false));
updateScanUI();
request(true);
})();
