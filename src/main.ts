import './style.css'
import { initUeff, drawUeffGraph, resizeUeff } from './ueff';

let canvas: HTMLCanvasElement;
let ctx: CanvasRenderingContext2D;
const dpr = window.devicePixelRatio || 1;

// ==============================================================================================================================
// Simulator UI
// ==============================================================================================================================
let hudStatus: string = ''; // HUD status text (rendered on canvas). Kept in module scope so draw() can access it.
let center: [number, number];
let scaleKM: number = 60; // km per pixel
let running = false;
let history1: Vec[] = [];
let history2: Vec[] = [];

// ==============================================================================================================================
// Orbital parameters
// ==============================================================================================================================
type SimState = {
  mu: number;
  r1: Vec; // Earth position in COM frame
  v1: Vec; // Earth velocity in COM frame
  r2: Vec; // Spacecraft position in COM frame
  v2: Vec; // Spacecraft velocity in COM frame
  t: number;
};
let simState: SimState = {
  mu: 398600.4418, // Earth's gravitational parameter, km^3/s^2
  r1: vec(0, 0),
  v1: vec(0, 0),
  r2: vec(12000, 0),
  v2: vec(0, 0), // chosen velocity (chose depending on theme. ie. solar system)
  t: 0
};

function resize() {
  const rect = canvas.getBoundingClientRect();
  canvas.width = rect.width * dpr;
  canvas.height = rect.height * dpr;
  ctx.setTransform(dpr, 0, 0, -dpr, canvas.width / 2, canvas.height / 2);
  center = [0,0];
}

function resizeAll() {
  resize();         // main canvas sizing
  resizeUeff(dpr);  // guarded single-resize for ueff
  draw();
}

// ==============================================================================================================================
// Utilities
// ==============================================================================================================================
type Vec = [number, number];
function vec(x:number,y:number): Vec{return[x,y];}
//component-wise operations
function vadd(a:Vec,b:Vec):Vec{return[a[0]+b[0],a[1]+b[1]];}
function vsub(a:Vec,b:Vec):Vec{return[a[0]-b[0],a[1]-b[1]];}
function vscale(s:number,v:Vec):Vec{return[s*v[0],s*v[1]];}
function vdot(a:Vec,b:Vec):number{return a[0]*b[0]+a[1]*b[1];}
function vmag(v:Vec):number{return Math.sqrt(vdot(v,v));}

// ==============================================================================================================================
// Stumpff functions C(z) and S(z)
// ==============================================================================================================================
function stumpffC(z: number): number {
  if (z > 0) {
    return (1 - Math.cos(Math.sqrt(z))) / z;
  } else if (z < 0) {
    return (Math.cosh(Math.sqrt(-z)) - 1) / -z;
  } else {
    return 1 / 2;
  }
}
function stumpffS(z: number): number {
  if (z > 0) {
    return (Math.sqrt(z) - Math.sin(Math.sqrt(z))) / (Math.sqrt(z ** 3));
  } else if (z < 0) {
    return (Math.sinh(Math.sqrt(-z)) - Math.sqrt(-z)) / (Math.sqrt((-z) ** 3));
  } else {
    return 1 / 6;
  }
}

// ==============================================================================================================================
// Universal kepler equation & its derivative
// ==============================================================================================================================
function uniKeplerF(chi:number, r0:number, vr0:number, alpha:number, dt:number, mu:number):number{
    const z = alpha * chi**2;
    return (r0 * vr0 / Math.sqrt(mu)) * chi**2 * stumpffC(z) + (1 - alpha * r0) * chi**3 * stumpffS(z) + r0 * chi - Math.sqrt(mu) * dt;
}
function uniKeplerD(chi:number, r0:number, vr0:number, alpha:number, mu:number):number{
    const z = alpha * chi**2;
    return (r0 * vr0 / Math.sqrt(mu)) * chi * (1 - alpha * chi**2 * stumpffS(z)) + (1 - alpha * r0) * chi**2 * stumpffC(z) + r0;
}

// ==============================================================================================================================
// Solving universal kepler equation for chi using newton-raphson method
// ==============================================================================================================================
type StateVector = {r:Vec, v:Vec};
function solveUniversalKepler(r0vec:Vec, vr0vec:Vec, dt:number, mu:number):StateVector{
  const r0 = vmag(r0vec);
  const v0 = vmag(vr0vec);
  const vr0 = vdot(r0vec, vr0vec)/r0;
  const alpha = 2/r0 - (v0**2)/mu;

  let chi = Math.sqrt(mu) * Math.abs(alpha) * dt; // initial guess
  if (!isFinite(chi) || chi === 0) chi = Math.cbrt(mu)*dt; // fallback initial guess - if chi is not finite or is zero (parabolic case)
  const maxIter = 1000;
  const tol = 1e-8;

  for(let i=0; i<maxIter; i++){
    const f = uniKeplerF(chi, r0, vr0, alpha, dt, mu);
    const df = uniKeplerD(chi, r0, vr0, alpha, mu);
    const dchi = -f/df;
    chi += dchi;
    if(Math.abs(dchi) < tol){
      break;
    }
  }
  // computing f,g,fdot,gdot
  const z = alpha * chi**2;
  const C = stumpffC(z);
  const S = stumpffS(z);
  const fcoef = 1 - (chi**2 / r0) * C;
  const gcoef = dt - (chi**3 / Math.sqrt(mu)) * S;
  const rvec = vadd(vscale(fcoef, r0vec), vscale(gcoef, vr0vec));
  const rmag = vmag(rvec)
  if (rmag === 0) {
    throw new Error("Resulting position vector has zero magnitude.");
  }
  const fdotcoef = (chi * Math.sqrt(mu) / (r0 * rmag)) * (alpha * chi**2 * S - 1);
  const gdotcoef = 1 - (chi**2 / rmag) * C;
  const vvec = vadd(vscale(fdotcoef, r0vec), vscale(gdotcoef, vr0vec));
  return {r:rvec, v:vvec};
}

function parseVec(str: string): number[] {
  const parts: number[] = str.split(',').map(s => parseFloat(s.trim()));
  if (parts.length < 2) return [0, 0];
  return [parts[0], parts[1]];
}

// ==============================================================================================================================
// Reset from UI
// ==============================================================================================================================
function resetFromUI(){
  // Read G and masses from UI and compute gravitational parameter mu = G*(m1 + m2)
  const G = parseFloat(((document.getElementById('G') as HTMLInputElement)?.value) || '0');
  const m1 = parseFloat(((document.getElementById('m1') as HTMLInputElement)?.value) || '0');
  const m2 = parseFloat(((document.getElementById('m2') as HTMLInputElement)?.value) || '0');
  const computedMu = G * (m1 + m2);
  (document.getElementById('computedMu') as HTMLOutputElement).textContent = computedMu.toFixed(4)
  simState.mu = computedMu;

  // reference frame choice
  const refFrame = (document.getElementById('refFrame') as HTMLSelectElement).value;

  // Initial relative position and velocity (as user specified)
  let rRel = vec(...(parseVec((document.getElementById('r0') as HTMLInputElement).value) as [number, number]));
  let vRel: [number, number];
  /*
  let vRel = vec(...(parseVec((document.getElementById('v0') as HTMLInputElement).value) as [number, number]));
  use if you want user controlled initial velocity
  */

  const rmag = vmag(rRel);
  const mu = simState.mu;
  const e = parseFloat((document.getElementById('eccInput') as HTMLInputElement).value);
  const vmagnitude = Math.sqrt(mu * (1 + e) / rmag); // vis-viva equation for desired eccentricity at periapsis
  const tangential = vec(-rRel[1]/rmag, rRel[0]/rmag); // unit vector perpendicular to r
  vRel = vscale(vmagnitude, tangential);
  
  if (refFrame === 'COM') {
    // center of mass inertial frame
    simState.r1 = vscale(-m2 / (m1 + m2), rRel);
    simState.r2 = vscale(m1 / (m1 + m2), rRel);
    simState.v1 = vscale(-m2 / (m1 + m2), vRel);
    simState.v2 = vscale(m1 / (m1 + m2), vRel);
  } else {
    // relative non inertial frame
    simState.r1 = vec(0,0); // earth at origin
    simState.r2 = rRel; // spacecraft at relative distance
    simState.v1 = vec(0,0);
    simState.v2 = vRel;
  }
  
  simState.t = 0;
  history1 = [simState.r1.slice() as Vec];
  history2 = [simState.r2.slice() as Vec];

  // Ueff graph
  const M = parseFloat(((document.getElementById('m1') as HTMLInputElement)?.value) || '0');
  const m = parseFloat(((document.getElementById('m2') as HTMLInputElement)?.value) || '0');

  drawUeffGraph(M, m, G, rRel, vRel);

  fitScale();
  draw();
  updateStatus('reset');
}

// ==============================================================================================================================
// Scaling (fitting full orbit in view)
// ==============================================================================================================================
function fitScale(){
  const rmag = vmag(vsub(simState.r2, simState.r1));
  const mu = simState.mu;

  // for relative orbit
  const vRel = vsub(simState.v2, simState.v1);
  const eps = (vmag(vRel)**2)/2 - mu/rmag; // specific energy
  const h = (simState.r2[0] - simState.r1[0]) * (simState.v2[1] - simState.v1[1]) - (simState.r2[1] - simState.r1[1]) * (simState.v2[0] - simState.v1[0]); // specific angular momentum

  // eccentricity (clamp to zero for near circular orbits to avoid numerical issues)
  const eRaw = Math.sqrt(Math.max(0, 1 + (2*eps*h*h)/(mu*mu)));
  const e = Math.min(eRaw, 1); // clamping zoom at parabolic limit

  let rmax = rmag;
  if (eps < 0) {
    // elliptical bound orbit
    const a = -mu/(2*eps); // semi-major axis
    rmax = a * (1 + e); // apoapsis distance
  } else {
    // parabolic/hyperbolic - zoom out cap
    rmax = rmag * 19999; // zoom out cap - (1 + 0.9999)/(1 - 0.9999) = 19999
  }

  const margin = 0.65 // orbit takes up 65% of canvas size
  const halfCanvas = Math.min(canvas.width, canvas.height) / 2 * margin;
  scaleKM = Math.max(5, rmax / halfCanvas);
}

// ==============================================================================================================================
// World to canvas (coordinate system mapping)
// ==============================================================================================================================
function worldToCanvas(pos:Vec):[number, number]{
  return [center[0] + pos[0]/scaleKM, center[1] + pos[1]/scaleKM];
}

// ==============================================================================================================================
// Drawing
// ==============================================================================================================================
function draw() {
  ctx.setTransform(dpr, 0, 0, -dpr, canvas.width / 2, canvas.height / 2);
  ctx.clearRect(-canvas.width / 2, -canvas.height / 2, canvas.width, canvas.height);

  // get reference frame choice
  const refFrame = (document.getElementById('refFrame') as HTMLSelectElement).value;

  // draw Earth (r1 - position depends on frame of reference)
  const drawEarthEl = document.getElementById('drawEarth') as HTMLInputElement | null;
  const drawEarth = drawEarthEl?.checked ?? false;
  if (drawEarth) {
    const earthRadiusKM = 3000; // chosen Earth's radius in km
    const earthRadiusPixels = earthRadiusKM / scaleKM;

    // if condition: COM frame is true use simState.r1, if condition: Relative frame is true use [center[0], center[1]]
    const earthDisplayPos = refFrame === 'COM' ? worldToCanvas(simState.r1) : [center[0], center[1]];

    ctx.beginPath();
    ctx.arc(earthDisplayPos[0], earthDisplayPos[1], earthRadiusPixels, 0, 2 * Math.PI);
    
    // earth colour gradient
    const gradient = ctx.createRadialGradient(earthDisplayPos[0], earthDisplayPos[1], earthRadiusPixels * 0.1, earthDisplayPos[0], earthDisplayPos[1], earthRadiusPixels);
    gradient.addColorStop(0, 'rgb(37, 71, 74)');
    gradient.addColorStop(1, 'rgb(14, 27, 27)');
    ctx.fillStyle = gradient;
    ctx.fill();
  }

  // draw history
  const showTrailEl = document.getElementById('showTrails') as HTMLInputElement | null;
    // earth trail
    if (showTrailEl?.checked && history1.length > 1) {
    ctx.beginPath();
    const startCanvasPos = worldToCanvas(history1[0]);
    ctx.moveTo(startCanvasPos[0], startCanvasPos[1]);
    for (let i = 1; i < history1.length; i++) {
      const canvasPos = worldToCanvas(history1[i]);
      ctx.lineTo(canvasPos[0], canvasPos[1]);
    }
    ctx.strokeStyle = 'rgba(255, 252, 248, 0.3)';
    ctx.lineWidth = 1.5;
    ctx.stroke();
  }
  // spacecraft trail
  if (showTrailEl?.checked && history2.length > 1) {
    ctx.beginPath();
    const startCanvasPos = worldToCanvas(history2[0]);
    ctx.moveTo(startCanvasPos[0], startCanvasPos[1]);
    for (let i = 1; i < history2.length; i++) {
      const canvasPos = worldToCanvas(history2[i]);
      ctx.lineTo(canvasPos[0], canvasPos[1]);
    }
    ctx.strokeStyle = 'rgba(255, 252, 248, 0.3)';
    ctx.lineWidth = 1.5;
    ctx.stroke();
  }

  if (refFrame === 'COM' && showTrailEl?.checked) {
    // Size in pixels, gently shrinking as you zoom out (scaleKM increases).
    const notchLen = Math.max(3, Math.min(4, 4 * Math.sqrt(60 / scaleKM)));
    const gap = Math.max(3, Math.min(5, 3 * Math.sqrt(60 / scaleKM)));

    // Opacity also drops a bit when zoomed out to avoid visual dominance.
    const opacity = Math.max(0.25, Math.min(0.6, 0.6 * Math.sqrt(60 / scaleKM)));

    ctx.save();
    ctx.strokeStyle = `rgba(160, 160, 160, ${opacity})`;
    ctx.fillStyle = `rgba(160, 160, 160, ${Math.min(0.9, opacity + 0.2)})`;
    ctx.lineWidth = 1.5;

    // Notches (gapped plus)
    ctx.beginPath();
    ctx.moveTo(- gap - notchLen, 0); ctx.lineTo(- gap, 0);
    ctx.moveTo(gap, 0);           ctx.lineTo(gap + notchLen, 0);
    ctx.moveTo(0, - gap - notchLen); ctx.lineTo(0, - gap);
    ctx.moveTo(0, gap);            ctx.lineTo(0, gap + notchLen);
    ctx.stroke();

    // Center dot
    ctx.beginPath();
    ctx.arc(0, 0, 2, 0, 2 * Math.PI);
    ctx.fill();

    ctx.restore();
  }

  // draw spacecraft (r2)
  const spacecraftPos = worldToCanvas(simState.r2);
  ctx.beginPath();
  // scale spacecraft size with scaleKM so it shrinks when zoomed out
  ctx.arc(spacecraftPos[0], spacecraftPos[1], 400 / scaleKM, 0, 2 * Math.PI);
  const gradient = ctx.createRadialGradient(spacecraftPos[0]-2, spacecraftPos[1]-2, 2, spacecraftPos[0],spacecraftPos[1],10);
  gradient.addColorStop(0, '#a12719');
  gradient.addColorStop(1, 'rgb(116, 30, 15)');
  ctx.fillStyle = gradient;
  ctx.fill();

  // HUD
  ctx.setTransform(dpr, 0, 0, dpr, canvas.width / 2, canvas.height / 2);
  ctx.fillStyle = 'white';
  ctx.font = '14px Arial';

  ctx.fillText(`t: ${simState.t.toFixed(1)} s`, -canvas.width / 2 + 300, canvas.height / 2 -990);
  ctx.fillText(`r: ${vmag(simState.r2).toFixed(1)} km`, -canvas.width / 2 + 300, canvas.height / 2 -970);
  ctx.fillText(`v: ${vmag(simState.v2).toFixed(3)} km/s`, -canvas.width / 2 + 300, canvas.height / 2 -950);

  // status (HUD)
  ctx.fillText(`Spacecraft status: ${hudStatus}`, - canvas.width / 2 + 300, -canvas.height / 2 + 1000);

  ctx.restore();
}

// real time textual data (HUD) update
function updateStatus(s: string){
  // update HUD string
  hudStatus = s;
}

// ==============================================================================================================================
// Step the simulation
// ==============================================================================================================================
function stepSimulation(dt: number) {
  // reference frame choice
  const refFrame = (document.getElementById('refFrame') as HTMLSelectElement).value;
  // previously computed constants and variables
  const m1 = parseFloat(((document.getElementById('m1') as HTMLInputElement)?.value) || '0');
  const m2 = parseFloat(((document.getElementById('m2') as HTMLInputElement)?.value) || '0');
  const rRel = vsub(simState.r2, simState.r1);
  const vRel = vsub(simState.v2, simState.v1);

  // propagate the relative motion using the universal Kepler solver
  const relState = solveUniversalKepler(rRel, vRel, dt, simState.mu);
  
  // update relative state
  const rRelnew = relState.r;
  const vRelnew = relState.v;

  /*
  // for drifting center of mass feature
  // const COM
  const rCOM = vscale(1 / (m1 + m2), vadd(vscale(m1, simState.r1), vscale(m2, simState.r2)));
  const vCOM = vscale(1 / (m1 + m2), vadd(vscale(m1, simState.v1), vscale(m2, simState.v2)));
  const rCOMnew = vadd(rCOM, vscale(dt, vCOM));
  */

  if (refFrame === 'COM') {
    // center of mass frame - both bodies orbiting common barycenter (inertial frame)
    simState.r1 = vscale(-m2 / (m1 + m2), rRelnew);
    simState.r2 = vscale(m1 / (m1 + m2), rRelnew);
    simState.v1 = vscale(-m2 / (m1 + m2), vRelnew);
    simState.v2 = vscale(m1 / (m1 + m2), vRelnew);
  } else {
    // perifocal frame - Earth fixed at origin, spacecraft at seperation distance (relative non inertial frame)
    simState.r1 = vec(0,0); // earth at origin
    simState.r2 = rRelnew; // spacecraft at relative distance
    simState.v1 = vec(0,0);
    simState.v2 = vRelnew;
  }

  simState.t += dt;
  history1.push(simState.r1.slice() as Vec);
  history2.push(simState.r2.slice() as Vec);

  // prevent history from growing indefinitely
  if (history1.length > 10000) {
    history1.shift();
  }
  if (history2.length > 10000) {
    history2.shift();
  }

  draw();
  updateStatus('advanced by ' + dt.toFixed(2) + ' s to t=' + simState.t.toFixed(2) + ' s');
}

// ==============================================================================================================================
// Animation loop
// ==============================================================================================================================
function startAnimationLoop(): void {
  let lastTimestamp: number | null = null;

  const loop = (ts: DOMHighResTimeStamp): void => {
    if (lastTimestamp === null) {
      lastTimestamp = ts;
    }
    // calculating change in time - converting from milliseconds (ts) to seconds
    const elapsed = (ts - lastTimestamp) / 1000;
    lastTimestamp = ts;

    if (running) {
      // get time scale from UI
      const timeScaleInput = document.getElementById('dt') as HTMLInputElement | null;
      if (timeScaleInput && timeScaleInput.value) {
        const dt = parseFloat(timeScaleInput.value);
        if (!isNaN(dt) && dt > 0) {
          // Decouple physics from frame rate:
          // Run multiple steps if elapsed time was long
          // Formula: steps = max(1, floor(1 + elapsed * 30))
          // This ensures ~30 steps per second target
          const steps = Math.max(1, Math.floor(1 + elapsed * 30));
          for (let i = 0; i < steps; i++) {
            stepSimulation(dt);
          }
        }
      }
    }

    draw();

    requestAnimationFrame(loop);
  };
  
  requestAnimationFrame(loop);
}

// ==============================================================================================================================
// UI wiring
// ==============================================================================================================================
function wireUI() {
  // wiring start button
  const startBtn = document.getElementById('start') as HTMLButtonElement | null;
  if (!startBtn) throw new Error('Missing #start button in the DOM');

  const updateStartLabel = () => { startBtn.innerText = running ? 'Stop' : 'Start'; };

  startBtn.onclick = () => {
    running = !running;
    updateStartLabel();
  };

  // ensure the label matches initial state
  updateStartLabel();


  (document.getElementById('step') as HTMLButtonElement).onclick = () => {
    const timeScaleInput = document.getElementById('dt') as HTMLInputElement | null;
    if (timeScaleInput && timeScaleInput.value) {
      const dt = parseFloat(timeScaleInput.value);
      if (!isNaN(dt) && dt > 0) {
        stepSimulation(dt);
      }
    }
  };
  (document.getElementById('reset') as HTMLButtonElement).onclick = () => {
    resetFromUI();
  };

  // eccentricity slider and input box synchronization -----------------------------------------------------------------------

  const eccSlider = document.getElementById('ecc') as HTMLInputElement | null;
  const eccInput  = document.getElementById('eccInput') as HTMLInputElement | null;

  if (!eccSlider || !eccInput) {
    throw new Error("Missing #ecc or #eccInput in the DOM");
  }
  // Read bounds from HTML attributes
  const maxEcc = parseFloat(eccInput.max);
  const minEcc = parseFloat(eccInput.min);
  const power = 3; // exponent for non-linear power scaling

  // Slider [0,1] → Eccentricity [minEcc, maxEcc]
  const sliderToEcc = (sliderVal: number): number => {
    const normalized = parseFloat(eccSlider.max) > 0 
      ? sliderVal / parseFloat(eccSlider.max) 
      : 0; // normalize slider to [0,1]
    return minEcc + (maxEcc - minEcc) * Math.pow(normalized, power);
  };

  // Eccentricity → Slider [0, eccSlider.max]
  const eccToSlider = (ecc: number): number => {
    const clamped = Math.max(minEcc, Math.min(maxEcc, ecc));
    const normalized = (clamped - minEcc) / (maxEcc - minEcc); // [0,1]
    const sliderRange = parseFloat(eccSlider.max);
    return sliderRange * Math.pow(normalized, 1 / power);
  };

  eccSlider.addEventListener('input', () => {
    const sliderVal = parseFloat(eccSlider.value);
    const ecc = sliderToEcc(sliderVal);
    eccInput.value = ecc.toFixed(4);
    resetFromUI();
  });

  eccInput.addEventListener('input', () => {
    const min = parseFloat(eccInput.min);
    const max = parseFloat(eccInput.max);

    let v = parseFloat(eccInput.value);
    if (!Number.isFinite(v)) return; // ignore invalid typing

    if (Number.isFinite(min) && v < min) v = min;
    if (Number.isFinite(max) && v > max) v = max;

    // round to the input step (4 decimal places)
    const step = parseFloat(eccInput.step) || 0.0001;
    v = Math.round(v / step) * step;

    eccInput.value = v.toFixed(4);
    eccSlider.value = String(eccToSlider(v));

    resetFromUI();
  });

  // initialize input to slider value on setup
  const initSliderVal = parseFloat(eccSlider.value);
  const initEcc = sliderToEcc(initSliderVal);
  eccInput.value = initEcc.toFixed(4);
  // ecc on refresh
  const INIT_ECC = 0.3;
  eccInput.value = INIT_ECC.toFixed(4);
  eccSlider.value = String(eccToSlider(INIT_ECC));
  // Wire G/mass inputs to update computed mu and reset simulation
  const gEl = document.getElementById('G') as HTMLInputElement | null;
  const m1El = document.getElementById('m1') as HTMLInputElement | null;
  const m2El = document.getElementById('m2') as HTMLInputElement | null;
  [gEl, m1El, m2El].forEach(el => {
    if (el) el.addEventListener('input', () => {
      resetFromUI();
    });
  });
  resetFromUI();
  // ------------------------------------------------------------------------------------------------------------------------
};

// ==============================================================================================================================
// DOM Ready - Initialize everything after HTML loads
// ==============================================================================================================================
document.addEventListener('DOMContentLoaded', () => {
  canvas = document.getElementById('orbit') as HTMLCanvasElement;
  ctx = canvas.getContext('2d')!;

  const ueff = document.getElementById('ueff') as HTMLCanvasElement;
  initUeff(ueff);

  resize();
  wireUI();

  // ONE initial sizing pass, after layout is stable
  requestAnimationFrame(() => {
    resizeAll();  // ueff canvas
    resetFromUI();
    draw();            // draw once with correct sizes
    startAnimationLoop();
  });

  window.addEventListener('resize', () => {
    resizeAll();
    draw();
  });
});
