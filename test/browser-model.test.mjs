import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import assert from 'node:assert/strict';
import test from 'node:test';
import { analyticState, rk4Step, energyOf, normalModes, validate, advanceClock } from '../web/model.mjs';
import { PRESETS } from '../web/presets.mjs';
const fixtures=JSON.parse(readFileSync(new URL('./fixtures/fortran.json',import.meta.url),'utf8'));
const close=(a,b,tolerance=2e-10)=>assert.ok(Math.abs(a-b)<=tolerance,`${a} != ${b}, delta=${a-b}`);
for (const {parameters:p,checkpoints} of fixtures) test(`Fortran/JS parity: preset ${p.id}, full horizon`,()=>{
  let u=[p.z0,0,p.theta0,0],previous=0;
  for(const q of checkpoints){
    for(let i=previous;i<q.step;i++)u=rk4Step(p,u,p.h);
    previous=q.step;
    u.forEach((v,i)=>close(v,q.numeric[i]));
    analyticState(p,q.time).forEach((v,i)=>close(v,q.analytic[i],2e-11));
    close(energyOf(p,u).totalEnergy,q.energy,2e-11);
  }
});
test('all fourteen presets are present, including the unmodified small inertia',()=>{
  assert.equal(PRESETS.length,14);assert.equal(PRESETS[13].I,Math.fround(1.39e-4));
  assert.equal(PRESETS[12].h,.001);
});
test('zero/tiny/negative coupling and arbitrary initial velocities',()=>{
  for(const eps of [0,1e-20,-.8]){
    const p={...PRESETS[6],eps};let u=[.2,.3,-.4,.5],initial=[...u];
    for(let i=0;i<100;i++)u=rk4Step(p,u,.001);
    const a=analyticState(p,.1,initial);u.forEach((v,i)=>close(v,a[i],1e-11));
  }
});
test('invalid states, unstable systems, off-resonance and invalid timesteps fail explicitly',()=>{
  const p=PRESETS[6];
  for(const patch of [{I:0},{m:-1},{delta:NaN},{eps:10},{k:2},{I:Infinity}])assert.throws(()=>validate({...p,...patch}),RangeError);
  for(const h of [0,-1,Infinity,5])assert.throws(()=>rk4Step(p,[0,0,0,0],h),RangeError);
  assert.throws(()=>analyticState(p,NaN),RangeError);
  assert.throws(()=>rk4Step(p,[0,NaN,0,0],.01),RangeError);
});
test('single mode retains constant total energy including coupling',()=>{
  const p=PRESETS[6],u0=[p.z0,0,p.z0*Math.sqrt(p.m/p.I),0],e0=energyOf(p,u0);
  for(let t=0;t<=20;t+=.1){const e=energyOf(p,analyticState(p,t,u0));
    // Individual coordinate energies are not generally constant in a coupled eigenmode.
    // The physical invariant is the complete energy, including coupling.
    close(e.totalEnergy,e0.totalEnergy,1e-16);
  }
});
test('clock does not advance at zero speed and is refresh-rate independent',()=>{
  assert.equal(advanceClock(.003,.02,0,.004).steps,0);
  for(const fps of [30,60,120]){
    let remainder=0,n=0;
    for(let i=0;i<fps;i++){const a=advanceClock(remainder,1/fps,1,.004);remainder=a.remainder;n+=a.steps;}
    assert.equal(n,250);close(remainder,0,1e-12);
  }
  assert.equal(advanceClock(0,100,1,.004).steps,400);
});
