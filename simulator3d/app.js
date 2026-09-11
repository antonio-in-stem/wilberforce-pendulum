"use strict";

// ---- presets (mirror src/wilberforce_presets.f90) --------------------------
const PRESETS = WilberforcePresets.map(p => ({...p,n:`Set ${p.id} · ${p.name}`,th0:p.theta0}));

// ---- model -----------------------------------------------------------------
const P = {delta:1.6, eps:0.075, m:0.07, I:0.5};
let z0 = 1.0, th0 = 0.0, speed = 1.0;
let S = {z:0, v:0, th:0, om:0, t:0};
let E0 = 0, running = true;
const hist = {t:[], z:[], th:[], E:[]};
const HMAX = 1400;

const kOf = () => P.m*P.delta/P.I;
const w0 = () => Wilberforce.normalModes(P).natural;
const w1 = () => Wilberforce.normalModes(P).plus;
const w2 = () => Wilberforce.normalModes(P).minus;
const Tbeat = () => P.eps === 0 ? Infinity : 2*Math.PI/(w1()-w2());
let accumulatedTime = 0;
function rk4(s,h){
  const u=Wilberforce.rk4Step(P,[s.z,s.v,s.th,s.om],h);
  return {z:u[0],v:u[1],th:u[2],om:u[3],t:s.t+h};
}
function energy(s){
  const e=Wilberforce.energyOf(P,[s.z,s.v,s.th,s.om]);
  return {Ez:e.verticalEnergy,Et:e.twistingEnergy,Ec:e.couplingEnergy,tot:e.totalEnergy};
}
function resetSim(){
  accumulatedTime=0;
  S={z:z0,v:0,th:th0,om:0,t:0}; E0=energy(S).tot;
  hist.t.length=hist.z.length=hist.th.length=hist.E.length=0;
}

// ============================================================================
//  Three.js scene
// ============================================================================
const CEIL_Y = 3.0, BOB_REST = -0.3, GROUND = -3.2, VZ = 1.35;
const BOB_R = 0.82, BOB_H = 1.0, SPR_R = 0.44, COILS = 9;

let renderer, scene, camera, controls, springMesh, bobGroup, stripe, rodGroup;
const stage = document.getElementById("stage");

function gradientTexture(){
  const c = document.createElement("canvas"); c.width = c.height = 512;
  const g = c.getContext("2d").createRadialGradient(256,180,60, 256,300,520);
  g.addColorStop(0, "#1b2740"); g.addColorStop(0.55, "#111a2c"); g.addColorStop(1, "#070a12");
  const cx = c.getContext("2d"); cx.fillStyle = g; cx.fillRect(0,0,512,512);
  const tex = new THREE.CanvasTexture(c); tex.encoding = THREE.sRGBEncoding; return tex;
}

function initThree(){
  const w = stage.clientWidth, h = stage.clientHeight;
  renderer = new THREE.WebGLRenderer({antialias:true, preserveDrawingBuffer:true});
  renderer.setPixelRatio(Math.min(window.devicePixelRatio||1, 2));
  renderer.setSize(w, h);
  renderer.outputEncoding = THREE.sRGBEncoding;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.12;
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.domElement.id = "gl";
  stage.appendChild(renderer.domElement);

  scene = new THREE.Scene();
  scene.background = gradientTexture();
  scene.fog = new THREE.Fog(0x0a0e18, 12, 26);

  camera = new THREE.PerspectiveCamera(40, w/h, 0.1, 100);
  camera.position.set(6.4, 2.9, 8.2);

  controls = new THREE.OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true; controls.dampingFactor = 0.08;
  controls.target.set(0, 0.05, 0);
  controls.minDistance = 4; controls.maxDistance = 18;
  controls.maxPolarAngle = Math.PI*0.92;

  // ---- lighting (studio) ----
  scene.add(new THREE.HemisphereLight(0x8fb0ff, 0x0e1220, 0.55));
  const key = new THREE.DirectionalLight(0xffffff, 2.0);
  key.position.set(6, 10, 6); key.castShadow = true;
  key.shadow.mapSize.set(2048, 2048);
  key.shadow.camera.near = 1; key.shadow.camera.far = 40;
  key.shadow.camera.left = -8; key.shadow.camera.right = 8;
  key.shadow.camera.top = 8;  key.shadow.camera.bottom = -8;
  key.shadow.bias = -0.0004; key.shadow.radius = 6;
  scene.add(key);
  const fill = new THREE.DirectionalLight(0x5f8bff, 0.5); fill.position.set(-7, 3, 2); scene.add(fill);
  const rim  = new THREE.PointLight(0xff7ab0, 0.6, 30);   rim.position.set(-3, 1, -5); scene.add(rim);

  // ---- ground with soft contact shadow ----
  const ground = new THREE.Mesh(
    new THREE.CircleGeometry(11, 64),
    new THREE.MeshStandardMaterial({color:0x0c1120, roughness:1.0, metalness:0.0}));
  ground.rotation.x = -Math.PI/2; ground.position.y = GROUND; ground.receiveShadow = true;
  scene.add(ground);
  const grid = new THREE.GridHelper(16, 32, 0x33507f, 0x1e2a44);
  grid.position.y = GROUND + 0.002; grid.material.opacity = 0.28; grid.material.transparent = true;
  scene.add(grid);

  // ---- ceiling mount ----
  const ceil = new THREE.Mesh(
    new THREE.BoxGeometry(3.0, 0.28, 3.0),
    new THREE.MeshStandardMaterial({color:0x2b3350, roughness:0.5, metalness:0.4}));
  ceil.position.y = CEIL_Y + 0.14; ceil.castShadow = true; ceil.receiveShadow = true;
  scene.add(ceil);

  // ---- spring (rebuilt each frame) ----
  const springMat = new THREE.MeshStandardMaterial({
    color:0x5b9bff, roughness:0.28, metalness:0.75, emissive:0x0c1a33, emissiveIntensity:0.6});
  springMesh = new THREE.Mesh(makeSpringGeometry(CEIL_Y, BOB_REST+BOB_H/2, 0), springMat);
  springMesh.castShadow = true; scene.add(springMesh);

  // ---- bob group (cylinder + stripe + rod + sliders) ----
  bobGroup = new THREE.Group(); scene.add(bobGroup);

  const bob = new THREE.Mesh(
    new THREE.CylinderGeometry(BOB_R, BOB_R, BOB_H, 64, 1, false),
    new THREE.MeshStandardMaterial({color:0xcdd3ea, roughness:0.32, metalness:0.9}));
  bob.castShadow = true; bob.receiveShadow = true; bobGroup.add(bob);
  // bevel rings on the rims for a machined look
  for (const sy of [BOB_H/2, -BOB_H/2]) {
    const ring = new THREE.Mesh(new THREE.TorusGeometry(BOB_R, 0.045, 12, 64),
      new THREE.MeshStandardMaterial({color:0x9aa2c4, roughness:0.4, metalness:0.9}));
    ring.rotation.x = Math.PI/2; ring.position.y = sy; ring.castShadow = true; bobGroup.add(ring);
  }
  // reference stripe on the rim (emissive) to reveal the twist
  stripe = new THREE.Mesh(new THREE.BoxGeometry(0.10, BOB_H*0.96, 0.06),
    new THREE.MeshStandardMaterial({color:0xff5c8a, emissive:0xff2d6a, emissiveIntensity:1.4, roughness:0.5}));
  stripe.position.set(BOB_R, 0, 0); bobGroup.add(stripe);

  // rod + sliders
  rodGroup = new THREE.Group(); bobGroup.add(rodGroup);
  const rod = new THREE.Mesh(new THREE.CylinderGeometry(0.045, 0.045, 3.5, 20),
    new THREE.MeshStandardMaterial({color:0x4a5170, roughness:0.4, metalness:0.85}));
  rod.rotation.z = Math.PI/2; rod.castShadow = true; rodGroup.add(rod);
  const slMat = new THREE.MeshStandardMaterial({color:0x9d6bff, roughness:0.35, metalness:0.6,
    emissive:0x2a145c, emissiveIntensity:0.5});
  for (const sx of [1.55, -1.55]) {
    const sl = new THREE.Mesh(new THREE.BoxGeometry(0.34, 0.42, 0.42), slMat);
    sl.position.x = sx; sl.castShadow = true; rodGroup.add(sl);
  }

  window.addEventListener("resize", onResize);
}

function makeSpringGeometry(topY, botY, twist){
  const pts = [], N = COILS*16;
  for (let i=0;i<=N;i++){
    const u=i/N, a=u*COILS*2*Math.PI + twist*0.4;
    pts.push(new THREE.Vector3(SPR_R*Math.cos(a), topY - u*(topY-botY), SPR_R*Math.sin(a)));
  }
  const curve = new THREE.CatmullRomCurve3(pts);
  return new THREE.TubeGeometry(curve, N, 0.058, 9, false);
}

function zAmp(){ let a=Math.abs(z0); for(const zz of hist.z) if(Math.abs(zz)>a) a=Math.abs(zz); return Math.max(a,0.02); }

function updateScene(){
  const drop = Math.max(-VZ, Math.min(VZ, (S.z/zAmp())*VZ));
  const bobY = BOB_REST - drop, th = S.th;
  bobGroup.position.y = bobY;
  bobGroup.rotation.y = th;
  // rebuild the spring between the ceiling and the bob top
  const old = springMesh.geometry;
  springMesh.geometry = makeSpringGeometry(CEIL_Y - 0.05, bobY + BOB_H/2 - 0.02, th);
  old.dispose();
  controls.update();
  renderer.render(scene, camera);
}

function onResize(){
  const w = stage.clientWidth, h = stage.clientHeight;
  renderer.setSize(w, h); camera.aspect = w/h; camera.updateProjectionMatrix();
  setupPlots();
}

// ============================================================================
//  2D plots
// ============================================================================
const pTime=document.getElementById("pTime"), pPhase=document.getElementById("pPhase"),
      pEnergy=document.getElementById("pEnergy");
let ctxT, ctxP, ctxE;
function setupCanvas(cv){ const r=cv.getBoundingClientRect(), dpr=window.devicePixelRatio||1;
  cv.width=r.width*dpr; cv.height=r.height*dpr; const c=cv.getContext("2d"); c.setTransform(dpr,0,0,dpr,0,0); return c; }
function setupPlots(){ ctxT=setupCanvas(pTime); ctxP=setupCanvas(pPhase); ctxE=setupCanvas(pEnergy); }
function grid(c,w,h){ c.fillStyle="#0d1526"; c.fillRect(0,0,w,h);
  c.strokeStyle="rgba(90,120,170,0.12)"; c.lineWidth=1;
  for(let i=1;i<6;i++){ c.beginPath(); c.moveTo(0,h*i/6); c.lineTo(w,h*i/6); c.stroke(); } }
function series(c,w,h,arr,scale,col){ c.strokeStyle=col; c.lineWidth=1.6; c.beginPath();
  const N=arr.length; for(let i=0;i<N;i++){ const X=i/(HMAX-1)*w, Y=h/2-arr[i]/scale*h*0.42; i?c.lineTo(X,Y):c.moveTo(X,Y);} c.stroke(); }
function drawPlots(){
  const N=hist.t.length;
  let cv=pTime,c=ctxT,w=cv.clientWidth,h=cv.clientHeight; grid(c,w,h);
  if(N>1){ const zm=Math.max(1e-9,...hist.z.map(Math.abs)), tm=Math.max(1e-9,...hist.th.map(Math.abs));
    series(c,w,h,hist.z,zm,"#4f9bff"); series(c,w,h,hist.th,tm,"#b57cff"); }
  cv=pPhase;c=ctxP;w=cv.clientWidth;h=cv.clientHeight; grid(c,w,h);
  if(N>1){ const zm=Math.max(1e-9,...hist.z.map(Math.abs)), tm=Math.max(1e-9,...hist.th.map(Math.abs));
    c.strokeStyle="#33d6a6"; c.lineWidth=1.5; c.beginPath();
    for(let i=0;i<N;i++){ const X=w/2+hist.th[i]/tm*w*0.42, Y=h/2-hist.z[i]/zm*h*0.42; i?c.lineTo(X,Y):c.moveTo(X,Y);} c.stroke();
    const X=w/2+hist.th[N-1]/tm*w*0.42, Y=h/2-hist.z[N-1]/zm*h*0.42;
    c.fillStyle="#ff5c8a"; c.beginPath(); c.arc(X,Y,3.6,0,7); c.fill(); }
  cv=pEnergy;c=ctxE;w=cv.clientWidth;h=cv.clientHeight; grid(c,w,h);
  const e=energy(S), tot=Math.max(1e-12,Math.abs(E0));
  const bars=[["Ez",e.Ez,"#4f9bff"],["Eθ",e.Et,"#b57cff"],["Ec",e.Ec,"#ffb44a"],["E",e.tot,"#33d6a6"]];
  const bw=w*0.84/bars.length*0.6;
  bars.forEach((b,i)=>{ const x=w*0.08+i*(w*0.84/bars.length); const bh=Math.max(-h*0.42,Math.min(h*0.42,b[1]/tot*h*0.42));
    c.fillStyle=b[2]; c.fillRect(x,h*0.55-Math.max(bh,0),bw,Math.abs(bh));
    c.fillStyle="#93a0bf"; c.font="10px sans-serif"; c.textAlign="center"; c.fillText(b[0],x+bw/2,h*0.72); });
  c.strokeStyle="rgba(150,170,210,0.25)"; c.beginPath(); c.moveTo(0,h*0.55); c.lineTo(w,h*0.55); c.stroke();
}

// ============================================================================
//  UI
// ============================================================================
const $=id=>document.getElementById(id);
function fmt(x){ if(x===Infinity)return "∞"; const a=Math.abs(x);
  if(a!==0&&(a<1e-3||a>=1e4))return x.toExponential(2); return x.toFixed(a<1?4:3); }
function syncLabels(){
  $("vdelta").textContent=fmt(P.delta); $("veps").textContent=fmt(P.eps);
  $("vm").textContent=fmt(P.m); $("vI").textContent=fmt(P.I);
  $("vz0").textContent=fmt(z0); $("vth0").textContent=fmt(th0); $("vspeed").textContent=speed.toFixed(2);
  $("sw0").textContent=fmt(w0()); $("sw1").textContent=fmt(w1()); $("sw2").textContent=fmt(w2());
  $("stb").textContent=fmt(Tbeat())+" s"; $("sk").textContent=fmt(kOf());
}
function pushSliders(){
  const values={...P,z0,th0,speed};
  for(const id of ["delta","eps","m","I","z0","th0","speed"]){
    const el=$(id),value=values[id];
    // HTML range constraints must not silently rewrite a scientific preset.
    el.min=Math.min(Number(el.min),value);el.max=Math.max(Number(el.max),value);
    el.step="any";el.value=String(value);
  }
}
function loadPreset(i){
  const p=PRESETS[i]; if(!p)throw new RangeError("Unknown preset");
  Wilberforce.validate(p);
  P.delta=p.delta;P.eps=p.eps;P.m=p.m;P.I=p.I;
  z0=p.z0;th0=p.th0;
  $("model-error").textContent="";
  pushSliders();resetSim();syncLabels();
}
function bindSlider(id,obj,key){
  $(id).addEventListener("input",()=>{
    const value=Number($(id).value);
    if(!Number.isFinite(value))return;
    if(obj==="P"){
      try{Wilberforce.validate({...P,[key]:value});}
      catch(error){$("model-error").textContent=error.message;pushSliders();return;}
      P[key]=value;resetSim();
    } else if(key==="z0"){z0=value;resetSim();}
      else if(key==="th0"){th0=value;resetSim();}
      else if(key==="speed")speed=value;
    $("model-error").textContent="";syncLabels();
  });
}

function initUI(){
  PRESETS.forEach((p,i)=>{ const o=document.createElement("option"); o.value=i; o.textContent=p.n; $("preset").appendChild(o); });
  $("preset").value=11;
  $("preset").addEventListener("change",e=>loadPreset(+e.target.value));
  ["delta","eps","m","I"].forEach(k=>bindSlider(k,"P",k));
  bindSlider("z0",null,"z0"); bindSlider("th0",null,"th0"); bindSlider("speed",null,"speed");
  $("play").addEventListener("click",()=>{running=!running;$("play").textContent=running?"Pause":"Play";});
  $("reset").addEventListener("click",resetSim);
  $("view").addEventListener("click",()=>{ camera.position.set(6.4,2.9,8.2); controls.target.set(0,0.05,0); });
  $("periodic").addEventListener("click",()=>{ th0=z0*Math.sqrt(P.m/P.I);
    pushSliders(); resetSim(); syncLabels(); });
}

// ============================================================================
//  main loop
// ============================================================================
let last = performance.now();
function frame(now){
  const dt=Math.min(0.05,(now-last)/1000); last=now;
  if(running && !document.hidden){
    const h=Math.min(0.004,0.05/w1());
    const advance=Wilberforce.advanceClock(accumulatedTime,dt,speed,h);
    accumulatedTime=advance.remainder;const steps=advance.steps;
    for(let i=0;i<steps;i++){ S=rk4(S,h);
      hist.t.push(S.t); hist.z.push(S.z); hist.th.push(S.th); hist.E.push(energy(S).tot);
      if(hist.t.length>HMAX){ hist.t.shift(); hist.z.shift(); hist.th.shift(); hist.E.shift(); } }
  }
  updateScene(); drawPlots();
  const e=energy(S);
  $("st").textContent=S.t.toFixed(2)+" s";
  $("sE").textContent=fmt(e.tot)+"  ("+((e.tot-E0)).toExponential(1)+")";
  requestAnimationFrame(frame);
}

function init(){ initThree(); setupPlots(); initUI(); loadPreset(11); syncLabels(); requestAnimationFrame(frame); }
window.addEventListener("DOMContentLoaded", init);
