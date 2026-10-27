/**
 * Dependency-free, interactive molecular renderer.
 * Coordinates and bonds are displayed as supplied; the viewer does not alter
 * molecular geometry, connectivity, or the exported model.
 */

const ELEMENTS = {
  C: { color: '#50616b', light: '#a8b5b9', dark: '#263943', radius: 0.31 },
  N: { color: '#528dc2', light: '#c2e0f7', dark: '#255578', radius: 0.33 },
  O: { color: '#d47567', light: '#ffd0bc', dark: '#883d35', radius: 0.33 },
  H: { color: '#d9e0df', light: '#ffffff', dark: '#a1b0ae', radius: 0.19 },
  S: { color: '#d6b457', light: '#fff0b9', dark: '#917627', radius: 0.38 },
  P: { color: '#c290b2', light: '#f3d4e8', dark: '#825576', radius: 0.38 },
};

const clamp = (value, min, max) => Math.max(min, Math.min(max, value));
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const normalize = (v) => {
  const length = Math.hypot(...v);
  return length > 1e-12 ? v.map((value) => value / length) : [1, 0, 0];
};
const elementStyle = (element) => ELEMENTS[element] || ELEMENTS.C;

function segmentIntersectsBox(a, b, box) {
  let lower = 0, upper = 1;
  for (const [axis, minimum, maximum] of [['x', box.x, box.x + box.width], ['y', box.y, box.y + box.height]]) {
    const delta = b[axis] - a[axis];
    if (Math.abs(delta) < 1e-9) {
      if (a[axis] < minimum || a[axis] > maximum) return false;
    } else {
      const first = (minimum - a[axis]) / delta;
      const last = (maximum - a[axis]) / delta;
      lower = Math.max(lower, Math.min(first, last));
      upper = Math.min(upper, Math.max(first, last));
      if (lower > upper) return false;
    }
  }
  return true;
}

function principalFrame(points) {
  // Orient the longest dimension horizontally without changing the model.
  const center = [0, 0, 0];
  for (const point of points) for (let axis = 0; axis < 3; axis++) center[axis] += point[axis] / points.length;
  const centered = points.map((point) => point.map((value, axis) => value - center[axis]));
  const matrix = Array.from({ length: 3 }, () => [0, 0, 0]);
  for (const point of centered) {
    for (let row = 0; row < 3; row++) {
      for (let column = 0; column < 3; column++) matrix[row][column] += point[row] * point[column];
    }
  }
  const multiply = (vector) => matrix.map((row) => dot(row, vector));
  let primary = normalize([1, 0.73, 0.37]);
  for (let iteration = 0; iteration < 32; iteration++) primary = normalize(multiply(primary));
  if (primary[0] < 0) primary = primary.map((value) => -value);
  let secondary = normalize(cross(primary, Math.abs(primary[2]) < 0.8 ? [0, 0, 1] : [0, 1, 0]));
  for (let iteration = 0; iteration < 32; iteration++) {
    const next = multiply(secondary);
    const parallel = dot(next, primary);
    const perpendicular = next.map((value, axis) => value - parallel * primary[axis]);
    if (Math.hypot(...perpendicular) < 1e-10) break;
    secondary = normalize(perpendicular);
  }
  const tertiary = normalize(cross(primary, secondary));
  return centered.map((point) => [dot(point, primary), dot(point, secondary), dot(point, tertiary)]);
}

export class MoleculeViewer {
  constructor(canvas, { onInteraction } = {}) {
    if (!canvas || typeof canvas.getContext !== 'function') throw new TypeError('A canvas element is required.');
    this.canvas = canvas;
    this.context = canvas.getContext('2d', { alpha: true });
    if (!this.context) throw new Error('A 2D canvas context is unavailable.');
    this.onInteraction = typeof onInteraction === 'function' ? onInteraction : null;
    this.options = { showHydrogens: false, style: 'ball-stick', autoRotate: false };
    this.model = null;
    this.points = [];
    this.projected = [];
    this.endCapIndices = [];
    this.sprites = new Map();
    this.pointers = new Map();
    this.listeners = [];
    this.yaw = -0.16;
    this.pitch = -0.28;
    this.zoom = 1;
    this.pan = [0, 0];
    this.frame = 0;
    this.previousTime = 0;
    this.destroyed = false;
    this.fitScale = 1;
    this.radius = 1;
    this.selected = -1;
    this.dragDistance = 0;
    this.width = 1;
    this.height = 1;
    this.canvas.style.touchAction = 'none';
    this.canvas.style.cursor = 'grab';
    if (!canvas.hasAttribute('tabindex')) canvas.tabIndex = 0;
    if (!canvas.hasAttribute('aria-label')) canvas.setAttribute('aria-label', 'Interactive molecular structure. Drag to rotate, scroll to zoom, and Shift-drag to pan.');
    this.listen(canvas, 'pointerdown', (event) => this.pointerDown(event));
    this.listen(canvas, 'pointermove', (event) => this.pointerMove(event));
    this.listen(canvas, 'pointerup', (event) => this.pointerUp(event));
    this.listen(canvas, 'pointercancel', (event) => this.pointerUp(event, true));
    this.listen(canvas, 'lostpointercapture', (event) => this.pointerUp(event, true));
    this.listen(canvas, 'wheel', (event) => {
      event.preventDefault();
      const multiplier = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? this.height : 1;
      this.zoom = clamp(this.zoom * Math.exp(-event.deltaY * multiplier * 0.0015), 0.08, 15);
      this.interact('zoom');
      this.requestRender();
    }, { passive: false });
    this.listen(canvas, 'dblclick', () => this.resetView());
    this.listen(canvas, 'contextmenu', (event) => event.preventDefault());
    this.listen(canvas, 'keydown', (event) => this.keyDown(event));
    this.listen(document, 'visibilitychange', () => {
      this.previousTime = 0;
      if (!document.hidden) this.requestRender();
    });
    if (typeof ResizeObserver !== 'undefined') {
      this.resizeObserver = new ResizeObserver(() => this.resize());
      this.resizeObserver.observe(canvas);
    } else {
      this.listen(window, 'resize', () => this.resize());
    }
    this.resize();
  }

  listen(target, type, listener, options) {
    target.addEventListener(type, listener, options);
    this.listeners.push(() => target.removeEventListener(type, listener, options));
  }

  interact(type, extra = {}) {
    if (this.onInteraction) this.onInteraction({ type, ...extra });
  }

  setModel(model) {
    if (!model || !Array.isArray(model.atoms) || !Array.isArray(model.bonds)) throw new TypeError('A model must contain atoms and bonds.');
    for (const atom of model.atoms) {
      if (![atom.x, atom.y, atom.z].every(Number.isFinite)) throw new TypeError('Atom coordinates must be finite numbers.');
    }
    for (const bond of model.bonds) {
      if (!Number.isInteger(bond.a) || !Number.isInteger(bond.b) || !model.atoms[bond.a] || !model.atoms[bond.b]) {
        throw new TypeError('Bond endpoints must be valid zero-based atom indices.');
      }
    }
    this.model = model;
    this.endCapIndices = model.atoms.flatMap((atom, index) => atom.role === 'endCap' && atom.element === 'C' ? [index] : []);
    this.points = model.atoms.length ? principalFrame(model.atoms.map(({ x, y, z }) => [x, y, z])) : [];
    this.radius = Math.max(1, ...this.points.map((point) => Math.hypot(...point)));
    this.selected = -1;
    this.resetView(false);
  }

  setOptions(options = {}) {
    if (typeof options.showHydrogens === 'boolean') this.options.showHydrogens = options.showHydrogens;
    if (typeof options.autoRotate === 'boolean') this.options.autoRotate = options.autoRotate;
    if (options.style === 'wireframe' || options.style === 'ball-stick') this.options.style = options.style;
    this.previousTime = 0;
    this.requestRender();
  }

  resetView(notify = true) {
    this.yaw = -0.16;
    this.pitch = -0.28;
    this.zoom = 1;
    this.pan = [0, 0];
    this.computeFit();
    this.requestRender();
    if (notify) this.interact('reset');
  }

  resize() {
    if (this.destroyed) return;
    const bounds = this.canvas.getBoundingClientRect();
    this.visible = bounds.width > 0 && bounds.height > 0;
    if (!this.visible) {
      if (this.frame) cancelAnimationFrame(this.frame);
      this.frame = 0;
      this.previousTime = 0;
      return;
    }
    this.width = Math.max(1, bounds.width || this.canvas.clientWidth || 640);
    this.height = Math.max(1, bounds.height || this.canvas.clientHeight || 440);
    this.dpr = Math.min(window.devicePixelRatio || 1, 2);
    const width = Math.round(this.width * this.dpr);
    const height = Math.round(this.height * this.dpr);
    if (this.canvas.width !== width || this.canvas.height !== height) {
      this.canvas.width = width;
      this.canvas.height = height;
    }
    this.computeFit();
    this.requestRender();
  }

  rotatePoint(point) {
    const cy = Math.cos(this.yaw), sy = Math.sin(this.yaw);
    const cp = Math.cos(this.pitch), sp = Math.sin(this.pitch);
    const x = cy * point[0] + sy * point[2];
    const z = -sy * point[0] + cy * point[2];
    return [x, cp * point[1] - sp * z, sp * point[1] + cp * z];
  }

  computeFit() {
    if (!this.points.length) return;
    let extentX = 0, extentY = 0;
    const cameraDistance = this.radius * 5 + 10;
    for (const point of this.points) {
      const rotated = this.rotatePoint(point);
      const perspective = cameraDistance / (cameraDistance - rotated[2]);
      extentX = Math.max(extentX, Math.abs(rotated[0] * perspective));
      extentY = Math.max(extentY, Math.abs(rotated[1] * perspective));
    }
    this.fitScale = Math.min(this.width * 0.41 / Math.max(extentX + 0.6, 1), this.height * 0.36 / Math.max(extentY + 0.6, 1));
  }

  requestRender() {
    if (this.destroyed || this.frame || document.hidden || this.visible === false) return;
    this.frame = requestAnimationFrame((time) => {
      this.frame = 0;
      if (this.destroyed || this.visible === false) return;
      if (this.options.autoRotate && this.pointers.size === 0) {
        const elapsed = this.previousTime ? Math.min(time - this.previousTime, 50) : 0;
        this.yaw += elapsed * 0.00016;
      }
      this.previousTime = time;
      this.render();
      if (this.options.autoRotate && this.pointers.size === 0) this.requestRender();
    });
  }

  sprite(element) {
    if (this.sprites.has(element)) return this.sprites.get(element);
    const canvas = document.createElement('canvas');
    canvas.width = 64;
    canvas.height = 64;
    const context = canvas.getContext('2d');
    const material = elementStyle(element);
    const gradient = context.createRadialGradient(22, 18, 1, 33, 32, 30);
    gradient.addColorStop(0, material.light);
    gradient.addColorStop(0.33, material.color);
    gradient.addColorStop(1, material.dark);
    context.fillStyle = gradient;
    context.beginPath();
    context.arc(32, 32, 30, 0, Math.PI * 2);
    context.fill();
    context.strokeStyle = 'rgba(255,255,255,0.28)';
    context.lineWidth = 1;
    context.stroke();
    this.sprites.set(element, canvas);
    return canvas;
  }

  render() {
    const context = this.context;
    context.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    context.clearRect(0, 0, this.width, this.height);
    if (!this.model || !this.points.length) return;
    const scale = this.fitScale * this.zoom;
    const cameraDistance = this.radius * 5 + 10;
    const wireframe = this.options.style === 'wireframe';
    const showHydrogens = this.options.showHydrogens;
    this.projected = this.points.map((point, index) => {
      const rotated = this.rotatePoint(point);
      const perspective = cameraDistance / (cameraDistance - rotated[2]);
      const atom = this.model.atoms[index];
      const material = elementStyle(atom.element);
      return {
        x: this.width / 2 + this.pan[0] + rotated[0] * scale * perspective,
        y: this.height / 2 + this.pan[1] - rotated[1] * scale * perspective,
        z: rotated[2],
        radius: clamp(material.radius * scale * perspective, wireframe ? 0.6 : 1.05, 20),
        visible: showHydrogens || atom.element !== 'H',
        perspective,
      };
    });
    const layers = [];
    for (let index = 0; index < this.model.bonds.length; index++) {
      const bond = this.model.bonds[index];
      const a = this.projected[bond.a], b = this.projected[bond.b];
      if (a.visible && b.visible) layers.push({ bond: index, depth: (a.z + b.z) / 2 });
    }
    for (let index = 0; index < this.model.atoms.length; index++) {
      if (this.projected[index].visible) layers.push({ atom: index, depth: this.projected[index].z + 0.01 });
    }
    layers.sort((a, b) => a.depth - b.depth);
    context.lineCap = 'round';
    for (const layer of layers) {
      if (layer.bond !== undefined) {
        this.drawBond(this.model.bonds[layer.bond], scale, wireframe);
      } else {
        const index = layer.atom;
        const projected = this.projected[index];
        const atom = this.model.atoms[index];
        const radius = wireframe ? Math.max(1.15, Math.min(projected.radius * 0.28, 2.5)) : projected.radius;
        if (projected.x < -radius || projected.x > this.width + radius || projected.y < -radius || projected.y > this.height + radius) continue;
        if (wireframe || radius < 1.8) {
          context.fillStyle = elementStyle(atom.element).color;
          context.beginPath();
          context.arc(projected.x, projected.y, radius, 0, Math.PI * 2);
          context.fill();
        } else {
          context.drawImage(this.sprite(atom.element), projected.x - radius, projected.y - radius, radius * 2, radius * 2);
        }
        if (this.selected === index) {
          context.strokeStyle = '#648f76';
          context.lineWidth = 1.5;
          context.beginPath();
          context.arc(projected.x, projected.y, radius + 5, 0, Math.PI * 2);
          context.stroke();
        }
      }
    }
    this.drawEndCapLabels();
    if (this.selected >= 0 && this.projected[this.selected]?.visible) this.drawSelectedLabel();
  }

  drawEndCapLabels() {
    if (!this.endCapIndices.length) return;
    const context = this.context;
    const label = 'CH₃';
    context.font = '500 10px system-ui, sans-serif';
    const width = context.measureText(label).width + 12;
    const height = 18;
    const boxes = [];
    const angles = [0, -Math.PI / 4, Math.PI / 4, -Math.PI / 2, Math.PI / 2, Math.PI];
    for (const index of this.endCapIndices) {
      const point = this.projected[index];
      if (!point?.visible || index === this.selected || point.x < 0 || point.y < 0 || point.x > this.width || point.y > this.height) continue;
      const bond = this.model.bonds.find((entry) => (entry.a === index && this.model.atoms[entry.b].element !== 'H') || (entry.b === index && this.model.atoms[entry.a].element !== 'H'));
      const neighbor = bond ? this.projected[bond.a === index ? bond.b : bond.a] : null;
      // Prefer the space beyond the terminal bond, then try nearby directions.
      const direction = neighbor && Math.hypot(point.x - neighbor.x, point.y - neighbor.y) > 1
        ? Math.atan2(point.y - neighbor.y, point.x - neighbor.x)
        : Math.atan2(point.y - this.height / 2, point.x - this.width / 2);
      let chosen;
      for (const distance of [point.radius + 24, point.radius + 38, point.radius + 52]) {
        for (const angle of angles) {
          const x = point.x + Math.cos(direction + angle) * distance - width / 2;
          const y = point.y + Math.sin(direction + angle) * distance - height / 2;
          const box = { x, y, width, height };
          if (x < 6 || y < 6 || x + width > this.width - 6 || y + height > this.height - 6) continue;
          if (boxes.some((other) => x < other.x + other.width + 4 && x + width + 4 > other.x && y < other.y + other.height + 4 && y + height + 4 > other.y)) continue;
          const occupied = this.projected.some((atom) => {
            if (!atom.visible) return false;
            const radius = this.options.style === 'wireframe' ? 3 : atom.radius;
            return Math.hypot(atom.x - clamp(atom.x, x, x + width), atom.y - clamp(atom.y, y, y + height)) < radius + 3;
          });
          if (occupied) continue;
          if (this.model.bonds.some((entry) => {
            const a = this.projected[entry.a], b = this.projected[entry.b];
            return a.visible && b.visible && segmentIntersectsBox(a, b, box);
          })) continue;
          chosen = box;
          break;
        }
        if (chosen) break;
      }
      // Dense or strongly overlapping projections remain uncluttered.
      if (!chosen) continue;
      boxes.push(chosen);
      const centerX = chosen.x + width / 2, centerY = chosen.y + height / 2;
      const distance = Math.hypot(centerX - point.x, centerY - point.y);
      const dx = (centerX - point.x) / distance, dy = (centerY - point.y) / distance;
      const inset = Math.min(Math.abs(width / (2 * dx)), Math.abs(height / (2 * dy)));
      const radius = this.options.style === 'wireframe' ? 3 : point.radius;
      context.strokeStyle = 'rgba(105, 132, 110, 0.48)';
      context.lineWidth = 0.8;
      context.beginPath();
      context.moveTo(point.x + dx * (radius + 3), point.y + dy * (radius + 3));
      context.lineTo(centerX - dx * (inset + 2), centerY - dy * (inset + 2));
      context.stroke();
      context.fillStyle = 'rgba(250, 253, 247, 0.94)';
      context.beginPath();
      if (context.roundRect) context.roundRect(chosen.x, chosen.y, width, height, 4);
      else context.rect(chosen.x, chosen.y, width, height);
      context.fill();
      context.fillStyle = '#69836b';
      context.textAlign = 'center';
      context.textBaseline = 'middle';
      context.fillText(label, centerX, centerY + 0.5);
    }
  }

  drawBond(bond, scale, wireframe) {
    const context = this.context;
    const a = this.projected[bond.a], b = this.projected[bond.b];
    if ((a.x < 0 && b.x < 0) || (a.x > this.width && b.x > this.width) || (a.y < 0 && b.y < 0) || (a.y > this.height && b.y > this.height)) return;
    const atomA = this.model.atoms[bond.a], atomB = this.model.atoms[bond.b];
    const isHydrogen = atomA.element === 'H' || atomB.element === 'H';
    const perspective = (a.perspective + b.perspective) / 2;
    const width = wireframe ? clamp(scale * 0.045 * perspective, 0.8, 1.7) : clamp(scale * (isHydrogen ? 0.075 : 0.13) * perspective, 0.75, 5);
    const dx = b.x - a.x, dy = b.y - a.y;
    const distance = Math.hypot(dx, dy) || 1;
    const nx = -dy / distance, ny = dx / distance;
    const count = bond.order === 2 ? 2 : bond.order === 3 ? 3 : 1;
    const separation = Math.max(1.8, width * 1.5);
    for (let line = 0; line < count; line++) {
      const offset = (line - (count - 1) / 2) * separation;
      const x1 = a.x + nx * offset, y1 = a.y + ny * offset;
      const x2 = b.x + nx * offset, y2 = b.y + ny * offset;
      const midX = (x1 + x2) / 2, midY = (y1 + y2) / 2;
      context.lineWidth = width;
      context.strokeStyle = elementStyle(atomA.element).color;
      context.beginPath();
      context.moveTo(x1, y1);
      context.lineTo(midX, midY);
      context.stroke();
      context.strokeStyle = elementStyle(atomB.element).color;
      context.beginPath();
      context.moveTo(midX, midY);
      context.lineTo(x2, y2);
      context.stroke();
    }
  }

  drawSelectedLabel() {
    const context = this.context;
    const atom = this.model.atoms[this.selected];
    const point = this.projected[this.selected];
    const label = atom.role === 'endCap' && atom.element === 'C'
      ? `CH₃ ${atom.endCap === 'tail' ? 'tail' : 'head'} cap · ${atom.name || atom.element}`
      : `${atom.name || atom.element} · Unit ${atom.residue || 1}${atom.charge ? ` · ${atom.charge > 0 ? '+' : ''}${atom.charge}` : ''}`;
    context.font = '500 11px system-ui, sans-serif';
    const width = context.measureText(label).width + 20;
    const x = clamp(point.x - width / 2, 10, Math.max(10, this.width - width - 10));
    const y = clamp(point.y - point.radius - 36, 10, this.height - 32);
    context.fillStyle = 'rgba(255,255,255,0.96)';
    context.strokeStyle = '#dbe4dd';
    context.lineWidth = 1;
    context.beginPath();
    if (context.roundRect) context.roundRect(x, y, width, 26, 6);
    else context.rect(x, y, width, 26);
    context.fill();
    context.stroke();
    context.fillStyle = '#385447';
    context.textAlign = 'left';
    context.textBaseline = 'middle';
    context.fillText(label, x + 10, y + 13);
  }

  pointerDown(event) {
    this.canvas.focus({ preventScroll: true });
    this.pointers.set(event.pointerId, { x: event.clientX, y: event.clientY, button: event.button });
    this.canvas.setPointerCapture(event.pointerId);
    this.dragDistance = 0;
    this.canvas.style.cursor = 'grabbing';
  }

  pointerMove(event) {
    const previous = this.pointers.get(event.pointerId);
    if (!previous) return;
    const deltaX = event.clientX - previous.x, deltaY = event.clientY - previous.y;
    this.dragDistance += Math.hypot(deltaX, deltaY);
    if (this.pointers.size > 1) {
      const other = [...this.pointers.entries()].find(([id]) => id !== event.pointerId)?.[1];
      if (other) {
        const before = Math.hypot(previous.x - other.x, previous.y - other.y);
        const after = Math.hypot(event.clientX - other.x, event.clientY - other.y);
        if (before > 5 && after > 5) this.zoom = clamp(this.zoom * after / before, 0.08, 15);
        this.pan[0] += deltaX / 2;
        this.pan[1] += deltaY / 2;
      }
    } else if (event.shiftKey || previous.button === 1 || previous.button === 2) {
      this.pan[0] += deltaX;
      this.pan[1] += deltaY;
    } else {
      this.yaw += deltaX * 0.008;
      this.pitch = clamp(this.pitch + deltaY * 0.008, -Math.PI / 2, Math.PI / 2);
    }
    this.pointers.set(event.pointerId, { ...previous, x: event.clientX, y: event.clientY });
    this.interact('drag');
    this.requestRender();
  }

  pointerUp(event, cancelled = false) {
    if (!this.pointers.has(event.pointerId)) return;
    const previous = this.pointers.get(event.pointerId);
    if (!cancelled && this.pointers.size === 1 && this.dragDistance < 5 && previous.button === 0) {
      const bounds = this.canvas.getBoundingClientRect();
      this.selectAt(event.clientX - bounds.left, event.clientY - bounds.top);
    }
    this.pointers.delete(event.pointerId);
    if (this.canvas.hasPointerCapture(event.pointerId)) this.canvas.releasePointerCapture(event.pointerId);
    if (!this.pointers.size) {
      this.canvas.style.cursor = 'grab';
      this.previousTime = 0;
    }
    this.requestRender();
  }

  selectAt(x, y) {
    let chosen = -1;
    let depth = -Infinity;
    for (let index = 0; index < this.projected.length; index++) {
      const point = this.projected[index];
      if (!point.visible || point.z < depth) continue;
      if (Math.hypot(x - point.x, y - point.y) < Math.max(5, point.radius + 3)) {
        chosen = index;
        depth = point.z;
      }
    }
    this.selected = chosen;
    this.interact('select', { index: chosen, atom: chosen < 0 ? null : this.model.atoms[chosen] });
    this.requestRender();
  }

  keyDown(event) {
    const step = event.shiftKey ? 0.18 : 0.08;
    if (event.key === 'ArrowLeft') this.yaw -= step;
    else if (event.key === 'ArrowRight') this.yaw += step;
    else if (event.key === 'ArrowUp') this.pitch = clamp(this.pitch - step, -Math.PI / 2, Math.PI / 2);
    else if (event.key === 'ArrowDown') this.pitch = clamp(this.pitch + step, -Math.PI / 2, Math.PI / 2);
    else if (event.key === '+' || event.key === '=') this.zoom = clamp(this.zoom * 1.12, 0.08, 15);
    else if (event.key === '-' || event.key === '_') this.zoom = clamp(this.zoom / 1.12, 0.08, 15);
    else if (event.key === '0' || event.key === 'Home') this.resetView();
    else if (event.key === 'Escape') this.selected = -1;
    else return;
    event.preventDefault();
    this.interact('keyboard');
    this.requestRender();
  }

  getPNG() {
    this.render();
    return this.canvas.toDataURL('image/png');
  }

  destroy() {
    this.destroyed = true;
    if (this.frame) cancelAnimationFrame(this.frame);
    this.resizeObserver?.disconnect();
    for (const remove of this.listeners) remove();
    this.listeners = [];
    this.pointers.clear();
    this.sprites.clear();
  }
}
