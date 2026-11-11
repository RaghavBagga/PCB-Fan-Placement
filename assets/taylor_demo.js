/* Live Taylor demo: fan profile exp(-u), u = r^2/(2 sigma^2), against its degree-m Taylor polynomial in u. */
(function(){
"use strict";
const cv=document.getElementById("tplot");if(!cv)return;
const ctx=cv.getContext("2d"),$=id=>document.getElementById(id);
function P(u,m){let s=0,t=1;for(let k=0;k<=m;k++){s+=t;t*=-u/(k+1);}return s;}
function css(n){return getComputedStyle(document.documentElement).getPropertyValue(n).trim();}
function draw(){
  const d=Math.min(2,window.devicePixelRatio||1),r=cv.getBoundingClientRect();
  const W=Math.round(r.width*d),H=Math.round(r.height*d);if(cv.width!==W||cv.height!==H){cv.width=W;cv.height=H;}
  const m=+$("tm").value,rs=+$("tr").value,X0=46*d,X1=W-14*d,Y0=H-34*d,Y1=12*d,XM=4,YLO=-0.4,YHI=1.4;
  const px=x=>X0+(x/XM)*(X1-X0),py=y=>Y0-((y-YLO)/(YHI-YLO))*(Y0-Y1);
  ctx.clearRect(0,0,W,H);ctx.font=`${11*d}px ${css("--mono")}`;ctx.textBaseline="middle";
  // grid + axes
  ctx.strokeStyle=css("--line");ctx.fillStyle=css("--muted");ctx.lineWidth=1*d;
  for(const y of [0,0.5,1]){ctx.beginPath();ctx.moveTo(X0,py(y));ctx.lineTo(X1,py(y));ctx.stroke();ctx.textAlign="right";ctx.fillText(y.toFixed(1),X0-6*d,py(y));}
  ctx.textAlign="center";for(let x=0;x<=4;x++){ctx.fillText(x+"σ",px(x),Y0+14*d);}
  ctx.fillText("r (distance from fan centre)",(X0+X1)/2,H-6*d);
  // 2% band: region where polynomial is within 0.02 of exact, scanning outward from 0
  let rc=0;while(rc<XM&&Math.abs(P(rc*rc/2,m)-Math.exp(-rc*rc/2))<=0.02)rc+=0.002;
  ctx.fillStyle=css("--panel");ctx.globalAlpha=1;ctx.fillStyle="rgba(11,127,155,.14)";ctx.fillRect(px(0),Y1,px(Math.min(rc,XM))-px(0),Y0-Y1);
  // curves
  const line=(f,col,w)=>{ctx.beginPath();ctx.strokeStyle=col;ctx.lineWidth=w*d;let on=false;
    for(let i=0;i<=400;i++){const x=XM*i/400,y=f(x);if(y<YLO||y>YHI){on=false;continue;}on?ctx.lineTo(px(x),py(y)):ctx.moveTo(px(x),py(y));on=true;}ctx.stroke();};
  line(x=>Math.exp(-x*x/2),css("--ink"),2.4);
  line(x=>P(x*x/2,m),css("--accent"),2);
  // probe
  const u=rs*rs/2,ex=Math.exp(-u),ap=P(u,m);
  ctx.strokeStyle=css("--muted");ctx.setLineDash([4*d,4*d]);ctx.beginPath();ctx.moveTo(px(rs),Y1);ctx.lineTo(px(rs),Y0);ctx.stroke();ctx.setLineDash([]);
  ctx.fillStyle=css("--ink");ctx.beginPath();ctx.arc(px(rs),py(ex),4*d,0,7);ctx.fill();
  if(ap>=YLO&&ap<=YHI){ctx.fillStyle=css("--accent");ctx.beginPath();ctx.arc(px(rs),py(ap),4*d,0,7);ctx.fill();}
  // readout
  const err=Math.abs(ap-ex);let bound=1;for(let k=1;k<=m+1;k++)bound*=u/k;
  let need=0;{let b=1;for(let n=1;n<200;n++){b*=u/n;if(b<1e-3){need=n;break;}}}   // alternating-series bound u^n/n! after n terms
  $("tmO").textContent="m = "+m+" (order "+2*m+" in r)";$("trO").textContent=rs.toFixed(2)+" σ";
  $("tread").innerHTML=
    "At r = "+rs.toFixed(2)+"σ (u = "+u.toFixed(3)+"): exact "+ex.toFixed(4)+", polynomial "+ap.toFixed(4)+
    "<br>actual error "+err.toExponential(2)+" ≤ bound u<sup>m+1</sup>/(m+1)! = "+bound.toExponential(2)+
    "<br>Polynomial stays within 0.02 of the exact fan out to r = "+rc.toFixed(2)+"σ (shaded)."+
    (rs>0?"<br>Terms needed here for error below 10<sup>−3</sup>: about "+need+".":"");
}
["tm","tr"].forEach(i=>$(i).addEventListener("input",draw));
window.addEventListener("resize",draw);
const mq=window.matchMedia&&window.matchMedia("(prefers-color-scheme: dark)");if(mq&&mq.addEventListener)mq.addEventListener("change",draw);
document.addEventListener("themechange",draw);
draw();
})();
