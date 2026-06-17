let ueffCanvas!: HTMLCanvasElement;
let ueffCtx!: CanvasRenderingContext2D;

const ueffPadding = 50;
const ueffPaddingRight = 10;

export function initUeff(canvas: HTMLCanvasElement): void {
  ueffCanvas = canvas;
  ueffCtx = canvas.getContext('2d')!;
}
export function resizeUeff(dpr: number): void {
  const rect = ueffCanvas.getBoundingClientRect();
  ueffCanvas.width = Math.round(rect.width * dpr);
  ueffCanvas.height = Math.round(rect.height * dpr);
  ueffCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
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
// Equation for the effective potential
// ==============================================================================================================================
function calcUeff(r: number, M: number, m: number, G: number, rRel: Vec, vRel: Vec): number {
  const h = Math.abs(rRel[0]*vRel[1] - rRel[1]*vRel[0]);
  const term1 = -(G * M * m) / r;
  const term2 = (M * m * (h**2)) / (2 * (M + m) * (r**2));
  return term1 + term2;
}

// ==============================================================================================================================
// Algebraicly finding the radius corresponding to the minimum effective potential (done by diffrentiating Ueff)
// ==============================================================================================================================
function findUeffMinimum(M: number, m: number, G: number, rRel: Vec, vRel: Vec): { r: number; U: number } {
  const h = Math.abs(rRel[0] * vRel[1] - rRel[1] * vRel[0]);
  const mu = G * (M + m);
  const rMin = h**2 / mu;  // where dU_eff/dr = 0
  return { r: rMin, U: calcUeff(rMin, M, m, G, rRel, vRel) };
}

// ==============================================================================================================================
// Draw the effective potential graph
// ==============================================================================================================================
export function drawUeffGraph(M: number, m: number, G: number, rRel: Vec, vRel: Vec): void {
  // compute min + plotting bounds
  const { r: rMin, U: UMin } = findUeffMinimum(M, m, G, rRel, vRel);
  const rMax = rMin * 20;
  const rMin_plot = rMin * 0.1;
  const UMax = 0;
  const UMin_plot = UMin * 1.5;

  // logical CSS size (drawing coordinates use logical units)
  const rect = ueffCanvas.getBoundingClientRect();
  const logicalWidth = rect.width;
  const logicalHeight = rect.height;

  const worldToCanvasX = (r: number): number =>
  ueffPadding + ((r - rMin_plot) / (rMax - rMin_plot)) * (logicalWidth - ueffPadding - ueffPaddingRight);

  const worldToCanvasY = (U: number): number =>
  logicalHeight - ueffPadding - ((U - UMin_plot) / (UMax - UMin_plot)) * (logicalHeight - 2 * ueffPadding);

  // background (logical coords)
  ueffCtx.clearRect(0, 0, logicalWidth, logicalHeight);


  // draw U_eff curve — sample in logical pixels (no DPR math here)
  // If you want curve thickness to be physically N pixels, set curveLineWidth = N / dpr
  const dpr = window.devicePixelRatio || 1;
  const curveLineWidth = 2; // leave as logical units (change to 2 / dpr for 2 physical px)
  ueffCtx.strokeStyle = '#fefefe';
  ueffCtx.lineWidth = curveLineWidth;
  ueffCtx.beginPath();

  const samples = Math.max(Math.floor(logicalWidth - 2 * ueffPadding), 2);
  for (let i = 0; i <= samples; i++) {
    const t = i / samples;
    const r = rMin_plot + t * (rMax - rMin_plot);
    const U = calcUeff(r, M, m, G, rRel, vRel);
    const x = worldToCanvasX(r);
    const y = worldToCanvasY(U);
    if (i === 0) ueffCtx.moveTo(x, y);
    else ueffCtx.lineTo(x, y);
  }
  ueffCtx.stroke();

  // mark r_min
  const rMinX = worldToCanvasX(rMin);
  const rMinY = worldToCanvasY(UMin);
  ueffCtx.fillStyle = '#f4faf4';
  ueffCtx.fillRect(rMinX - 4, rMinY - 4, 8, 8);

  // Optional: crisp 1px axes (aligned to physical pixels)
  // translate and lineWidth must be in logical units because the ctx was previously scaled by dpr in resizeUeff()
  const axisTranslate = 0.5 / dpr;   // logical units -> 0.5 physical px
  const axisLineWidth = 1 / dpr;     // logical units -> 1 physical px

  ueffCtx.save();
  ueffCtx.translate(axisTranslate, axisTranslate);
  ueffCtx.strokeStyle = '#666';
  ueffCtx.lineWidth = axisLineWidth;

  ueffCtx.beginPath();
  // round coordinates if you prefer exact integer locations in logical space
  ueffCtx.moveTo(ueffPadding, worldToCanvasY(0));
  ueffCtx.lineTo(logicalWidth - ueffPaddingRight, worldToCanvasY(0));
  ueffCtx.stroke();

  ueffCtx.beginPath();
  ueffCtx.moveTo(ueffPadding, ueffPadding);
  ueffCtx.lineTo(ueffPadding, logicalHeight - ueffPadding);
  ueffCtx.stroke();
  ueffCtx.restore();


  // labels and text (logical coordinates)
  
  // Generate sensible tick marks: aim for ~5 ticks on each axis
  // Adaptive axis labels

  function formatNum(v: number): string {
    if (Math.abs(v) >= 1e12) return (v / 1e12).toFixed(2) + 'T';
    if (Math.abs(v) >= 1e9) return (v / 1e9).toFixed(2) + 'G';
    if (Math.abs(v) >= 1e6) return (v / 1e6).toFixed(1) + 'M';
    if (Math.abs(v) >= 1e3) return (v / 1e3).toFixed(1) + 'k';
    return v.toFixed(0);
  }

  function makeTicks(min: number, max: number): number[] {
    const range = max - min;
    const mag = Math.pow(10, Math.floor(Math.log10(range)));
    const norm = range / mag;
    const step = norm < 1.5 ? mag * 0.2 : norm < 3 ? mag * 0.5 : norm < 7 ? mag : mag * 2;
    const ticks = [];
    for (let t = Math.ceil(min / step) * step; t <= max + step / 2; t += step) ticks.push(t);
    return ticks;
  }

  // Draw ticks
  ueffCtx.fillStyle = '#666';
  ueffCtx.font = '10px system-ui';
  ueffCtx.textAlign = 'center';
  makeTicks(rMin_plot, rMax).forEach(r => {
    const x = worldToCanvasX(r);
    ueffCtx.fillText(formatNum(r) + ' km', x, logicalHeight - ueffPadding + 12);
  });

  ueffCtx.textAlign = 'right';
  makeTicks(UMin_plot, UMax).forEach(u => {
    const y = worldToCanvasY(u);
    ueffCtx.fillText(formatNum(u) + ' J', ueffPadding - 6, y + 3);
  });

  // Labels
  ueffCtx.fillStyle = '#cacaca';
  ueffCtx.font = '14px system-ui';
  ueffCtx.textAlign = 'center';
  ueffCtx.fillText('r (km)', logicalWidth - ueffPaddingRight - 20, worldToCanvasY(0) - 15);
  ueffCtx.save();
  ueffCtx.translate(10, logicalHeight / 2);
  ueffCtx.rotate(-Math.PI / 2);
  ueffCtx.fillText('Ueff ( J )', 50, 10);
  ueffCtx.restore();

  //====================================================================================
  // autoscaling grid lines
  ueffCtx.save();
  ueffCtx.strokeStyle = 'rgba(102, 102, 102, 0.6)';
  ueffCtx.lineWidth = 0.5 / dpr;

  // Vertical grid lines (at r ticks)
  makeTicks(rMin_plot, rMax).forEach(r => {
    const x = worldToCanvasX(r);
    ueffCtx.beginPath();
    ueffCtx.moveTo(x, ueffPadding);
    ueffCtx.lineTo(x, logicalHeight - ueffPadding);
    ueffCtx.stroke();
  });

  // Horizontal grid lines (at U ticks)
  makeTicks(UMin_plot, UMax).forEach(u => {
    const y = worldToCanvasY(u);
    ueffCtx.beginPath();
    ueffCtx.moveTo(ueffPadding, y);
    ueffCtx.lineTo(logicalWidth - ueffPadding, y);
    ueffCtx.stroke();
  });

  ueffCtx.restore();
  //=====================================================================================

  // r_min marker
  ueffCtx.fillStyle = '#f4faf4';
  ueffCtx.fillRect(rMinX - 4, rMinY - 4, 8, 8);

}
