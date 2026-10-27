/** Real, explicit-hydrogen molecular graphs assembled from audited templates.
 * Coordinates are in angstroms. Bond endpoint indices are zero-based.
 * This is a deterministic initial-geometry builder, not a force-field engine.
 */
const GROUP_NAMES = {
  PA: 'Primary ammonium', QA: 'Quaternary ammonium',
  GD: 'Guanidinium', IL: 'Imidazolium',
};
const GROUPS = new Set(Object.keys(GROUP_NAMES));
const LENGTHS = new Set([4, 6, 8, 10]);
export const MAX_REPEATS = 30;
const CC_LENGTH = 1.53;
const CH_LENGTH = 1.09;
const subtract = (a, b) => a.map((v, i) => v - b[i]);
const add = (a, b) => a.map((v, i) => v + b[i]);
const scale = (a, s) => a.map(v => v * s);
const dot = (a, b) => a.reduce((sum, v, i) => sum + v * b[i], 0);
const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
const norm = a => Math.hypot(...a);
const unit = a => scale(a, 1 / norm(a));
const xyz = a => [a.x, a.y, a.z];

function rotate(point, axis, angle) {
  const c = Math.cos(angle), s = Math.sin(angle);
  return add(add(scale(point, c), scale(cross(axis, point), s)), scale(axis, dot(axis, point)*(1-c)));
}

function align(point, from, to) {
  const cosine = Math.max(-1, Math.min(1, dot(from, to)));
  if (cosine > 1 - 1e-12) return point.slice();
  let axis = cross(from, to);
  if (norm(axis) < 1e-10) axis = cross(from, Math.abs(from[0]) < 0.8 ? [1,0,0] : [0,1,0]);
  return rotate(point, unit(axis), Math.acos(cosine));
}

function dihedral(a, b, c, d) {
  const axis = unit(subtract(c, b));
  const left = subtract(a, b), right = subtract(d, c);
  const v = subtract(left, scale(axis, dot(left, axis)));
  const w = subtract(right, scale(axis, dot(right, axis)));
  return Math.atan2(dot(cross(axis, v), w), dot(v, w));
}

export function getTemplate(templates, group, carbons) {
  if (!GROUPS.has(group) || !LENGTHS.has(carbons)) throw new RangeError('Choose primary ammonium, quaternary ammonium, guanidinium or imidazolium and 4, 6, 8 or 10 alkyl carbons.');
  const template = templates[`${group}_${carbons}`];
  if (!template) throw new Error(`Missing chemical template for ${GROUP_NAMES[group]}, C${carbons}.`);
  return template;
}

function formulaFor(atoms, charge) {
  const counts = {};
  for (const atom of atoms) counts[atom.element] = (counts[atom.element] || 0) + 1;
  const keys = ['C', 'H', ...Object.keys(counts).filter(x => x !== 'C' && x !== 'H').sort()];
  const formula = keys.filter(k => counts[k]).map(k => k + (counts[k] === 1 ? '' : counts[k])).join('');
  return formula + (charge ? (Math.abs(charge) === 1 ? '' : `^${Math.abs(charge)}`) + (charge > 0 ? '+' : '-') : '');
}

/** Build one linear, homotypic chain with two pendent groups per repeat.
 * All 2n source attachment hydrogens are removed. Add n-1 C4--C1 bonds
 * and one neutral CH3 group at each main-chain endpoint.
 */
export function buildPolymer(templates, {group, carbons, repeats}) {
  const template = getTemplate(templates, group, carbons);
  if (!Number.isInteger(repeats) || repeats < 1 || repeats > MAX_REPEATS) throw new RangeError(`Repeat count must be an integer from 1 to ${MAX_REPEATS}.`);
  const {atoms: source, bonds: sourceBonds, attachment} = template;
  const {head, tail, headCap, tailCap} = attachment;
  const indexByName = Object.fromEntries(source.map((a, i) => [a.name, i]));
  const neighborSets = source.map(() => new Set());
  for (const bond of sourceBonds) { neighborSets[bond.a].add(bond.b); neighborSets[bond.b].add(bond.a); }
  for (const [atom, cap] of [[head, headCap], [tail, tailCap]]) {
    if (source[atom]?.element !== 'C' || source[cap]?.element !== 'H' || !neighborSets[atom].has(cap)) {
      throw new Error('Invalid attachment definition: a carbon-bound cap hydrogen is required.');
    }
  }
  const original = source.map(xyz), origin = original[head];
  const firstDirection = unit(subtract(original[tail], origin));
  let previous = original.map(p => align(subtract(p, origin), firstDirection, [1,0,0]));
  const firstPlaced = previous;
  const atoms = [], bonds = [], units = [];
  // A spatial hash makes candidate checking scale with nearby atoms, not the
  // entire growing chain. Hydrogens are retained in the graph and exports.
  const grid = new Map(), binSize = 3.5;
  const bin = p => p.map(v => Math.floor(v/binSize));
  const gridKey = b => b.join(',');
  const insert = entry => {
    const key = gridKey(bin(entry.p));
    if (!grid.has(key)) grid.set(key, []);
    grid.get(key).push(entry);
  };
  const heavy = source.map((a, i) => a.element === 'H' ? -1 : i).filter(i => i >= 0);
  let assemblyClashPenalty = 0;
  for (let repeat = 0; repeat < repeats; repeat++) {
    let placed = previous;
    if (repeat > 0) {
      const outgoing = unit(subtract(previous[tailCap], previous[tail]));
      const nextHead = add(previous[tail], scale(outgoing, CC_LENGTH));
      const incoming = unit(subtract(original[headCap], origin));
      const oriented = original.map(p => align(subtract(p, origin), incoming, scale(outgoing, -1)));
      let best = null;
      // Include an exact all-trans candidate and nearby rotamers. Prefer an
      // extended backbone while avoiding severe inter-repeat heavy-atom clashes.
      const starting = oriented.map(p => add(p, nextHead));
      const phi = dihedral(previous[indexByName.C3], previous[tail], starting[head], starting[indexByName.C2]);
      const transAngle = Math.PI - phi;
      for (const delta of [0,30,-30,60,-60,90,-90,120,-120,150,-150,180]) {
        const angle = transAngle + delta*Math.PI/180;
        const candidate = oriented.map(p => add(rotate(p, outgoing, angle), nextHead));
        let clash = 0;
        for (const i of heavy) {
          const p = candidate[i], bucket = bin(p);
          for (let x=-1; x<=1; x++) for (let y=-1; y<=1; y++) for (let z=-1; z<=1; z++) {
            const nearby = grid.get(gridKey([bucket[0]+x,bucket[1]+y,bucket[2]+z]));
            if (!nearby) continue;
            for (const entry of nearby) {
              if (entry.repeat === repeat-1 && (
                (i === head && (entry.local === tail || neighborSets[tail].has(entry.local))) ||
                (entry.local === tail && neighborSets[head].has(i)))) continue;
              const distance = norm(subtract(p, entry.p));
              // 2.6 A is a conservative heavy-atom starting-geometry clearance.
              if (distance < 2.6) clash += (2.6-distance)**2;
            }
          }
        }
        const forward = candidate[tail][0]-previous[tail][0];
        const score = clash*100 + (delta/180)**2*1.5 + Math.max(0, 3-forward)**2*2;
        if (!best || score < best.score) best = {candidate, score, clash};
      }
      placed = best.candidate;
      assemblyClashPenalty += best.clash;
    }
    const mapping = new Map();
    for (let i=0; i<source.length; i++) {
      if (i === headCap || i === tailCap) continue;
      const atom = source[i], point = placed[i];
      mapping.set(i, atoms.length);
      atoms.push({id: atoms.length+1, name: atom.name, element: atom.element,
        x: point[0], y: point[1], z: point[2], charge: atom.charge,
        residue: repeat+1, residueName: atom.residueName});
      if (atom.element !== 'H') insert({p: point, repeat, local: i});
    }
    for (const b of sourceBonds) if (mapping.has(b.a) && mapping.has(b.b)) {
      bonds.push({a: mapping.get(b.a), b: mapping.get(b.b), order: b.order});
    }
    if (repeat > 0) bonds.push({a: units[repeat-1].tail, b: mapping.get(head), order: 1});
    units.push({head: mapping.get(head), tail: mapping.get(tail)});
    previous = placed;
  }
  const endCaps = {};
  let endCapClashPenalty = 0;
  function addMethylCap(end, attachmentIndex, removedCapPosition, capNumber) {
    const boundary = atoms[attachmentIndex], attachmentPosition = xyz(boundary);
    const outward = unit(subtract(removedCapPosition, attachmentPosition));
    const carbonPosition = add(attachmentPosition, scale(outward, CC_LENGTH));
    const perpendicular = unit(cross(outward, Math.abs(outward[0]) < 0.8 ? [1,0,0] : [0,1,0]));
    const second = cross(outward, perpendicular);
    let best = null;
    // A tetrahedral methyl has three equivalent 120-degree-spaced hydrogens.
    // Sample its unique 120-degree torsional range to reduce cap contacts.
    for (let degrees=0; degrees<120; degrees+=10) {
      const hydrogens = [0,1,2].map(i => {
        const theta = (degrees + i*120)*Math.PI/180;
        const radial = add(scale(perpendicular, Math.cos(theta)), scale(second, Math.sin(theta)));
        const direction = add(scale(outward, 1/3), scale(radial, Math.sqrt(8/9)));
        return add(carbonPosition, scale(direction, CH_LENGTH));
      });
      let penalty = 0;
      for (const point of hydrogens) for (let i=0; i<atoms.length; i++) {
        if (i === attachmentIndex) continue; // This is a bonded 1-3 contact.
        const clearance = atoms[i].element === 'H' ? 1.6 : 2.0;
        const distance = norm(subtract(point, xyz(atoms[i])));
        if (distance < clearance) penalty += (clearance-distance)**2;
      }
      if (!best || penalty < best.penalty) best = {hydrogens, penalty};
    }
    const append = (name, element, point) => {
      const index = atoms.length;
      atoms.push({id:index+1, name, element, x:point[0], y:point[1], z:point[2], charge:0,
        residue:boundary.residue, residueName:boundary.residueName, role:'endCap', endCap:end});
      return index;
    };
    const carbon = append(`CM${capNumber}`, 'C', carbonPosition);
    const hydrogens = best.hydrogens.map((point,i) => append(`HM${capNumber}${i+1}`, 'H', point));
    bonds.push({a:attachmentIndex,b:carbon,order:1}, ...hydrogens.map(h=>({a:carbon,b:h,order:1})));
    endCaps[end] = {attachment:attachmentIndex,carbon,hydrogens};
    endCapClashPenalty += best.penalty;
  }
  addMethylCap('head', units[0].head, firstPlaced[headCap], 1);
  addMethylCap('tail', units[repeats-1].tail, previous[tailCap], 2);
  // Keep each PDB residue contiguous: the head cap belongs beside repeat 1,
  // rather than appearing after the last repeat with a restarted residue ID.
  const ordered = atoms.map((_,i)=>i).sort((a,b)=>atoms[a].residue-atoms[b].residue || a-b);
  const remap = new Map(ordered.map((old,index)=>[old,index]));
  const sortedAtoms = ordered.map((old,index)=>({...atoms[old],id:index+1}));
  atoms.splice(0,atoms.length,...sortedAtoms);
  for (const bond of bonds) {bond.a=remap.get(bond.a);bond.b=remap.get(bond.b);}
  for (const item of units) {item.head=remap.get(item.head);item.tail=remap.get(item.tail);}
  for (const cap of Object.values(endCaps)) {
    cap.attachment=remap.get(cap.attachment);
    cap.carbon=remap.get(cap.carbon);
    cap.hydrogens=cap.hydrogens.map(i=>remap.get(i));
  }
  const charge = atoms.reduce((sum, a) => sum+a.charge, 0);
  return {atoms, bonds,
    stats: {atoms: atoms.length, bonds: bonds.length, formula: formulaFor(atoms, charge), charge, repeats},
    metadata: {name: `${template.groupName} C${carbons} DP${repeats} CH3`, group, carbons, repeats,
      template: template.name, source: template.source, templateGeometry: template.geometry,
      connection: 'C4(i)-C1(i+1); C5 is retained as the pendent methyl group',
      removedHydrogens: 2*repeats, addedCapAtoms:8, addedCapBonds:8, endCaps,
      endCapChemistry:'Neutral CH3 at both main-chain endpoints; C5 pendent methyl retained',
      units, bondLengthAngstrom: CC_LENGTH, capCHBondLengthAngstrom:CH_LENGTH,
      geometry: 'Rigid conformer assembly with tetrahedral attachment vectors, inter-repeat rotamer screening and methyl-cap torsion screening; no full-polymer minimization',
      stereochemistry: 'One representative template geometry repeated; tacticity is not user-controlled',
      assemblyClashPenalty, endCapClashPenalty, counterions: 'Not included'}};
}

function validateExport(model) {
  if (!model?.atoms?.length || !Array.isArray(model.bonds)) throw new TypeError('A nonempty molecular model is required.');
  if (model.atoms.some(a => ![a.x,a.y,a.z].every(Number.isFinite))) throw new Error('Atom coordinates must be finite.');
  for (const b of model.bonds) {
    if (![b.a,b.b].every(i => Number.isInteger(i) && i >= 0 && i < model.atoms.length) || b.a === b.b || ![1,2,3].includes(b.order)) {
      throw new Error('Invalid bond endpoints or order.');
    }
  }
}

/** PDB includes formal charges and bond orders as repeated CONECT targets.
 * This convention is supported by RDKit and Open Babel; SDF also expresses
 * bond orders directly, for consumers that ignore repeated PDB targets.
 */
export function toPDB(model) {
  validateExport(model);
  if (model.atoms.length > 99998 || model.atoms.some(a => a.residue > 9999)) throw new RangeError('Model exceeds PDB serial/residue field limits; use SDF.');
  const lines = ['HEADER    POLYMER PROGRAMMING INITIAL STRUCTURE',
    `TITLE     ${model.metadata?.name || 'POLYMER'}`,
    'REMARK 200 INITIAL CONFORMER ASSEMBLY; NO FULL-POLYMER MINIMIZATION',
    'REMARK 200 C4(I)-C1(I+1); C5 IS A PENDENT METHYL; ANGSTROM UNITS',
    'REMARK 200 BOTH MAIN-CHAIN ENDS CAPPED BY NEUTRAL CH3 GROUPS',
    'REMARK 200 HEAD CAP CM1/HM11/HM12/HM13; TAIL CAP CM2/HM21/HM22/HM23',
    'REMARK 200 FORMAL CHARGES AND REPEATED CONECT TARGET BOND ORDERS',
    'REMARK 200 SDF ALSO AVAILABLE FOR EXPLICIT BOND ORDER INTERCHANGE'];
  for (let i=0; i<model.atoms.length; i++) {
    const a = model.atoms[i];
    const coords = [a.x,a.y,a.z].map(v => v.toFixed(3));
    if (coords.some(v => v.length > 8)) throw new RangeError('Coordinates exceed PDB field width; use SDF.');
    if (a.name.length > 4 || a.residueName.length > 3 || Math.abs(a.charge) > 9) throw new RangeError('Atom labels or charges exceed PDB fields; use SDF.');
    const charge = a.charge ? `${Math.abs(a.charge)}${a.charge>0?'+':'-'}` : '  ';
    lines.push(`HETATM${String(i+1).padStart(5)} ${a.name.padStart(4)} ${a.residueName.padStart(3)} A${String(a.residue).padStart(4)}    ${coords.map(v=>v.padStart(8)).join('')}  1.00  0.00          ${a.element.padStart(2)}${charge}`);
  }
  lines.push(`TER   ${String(model.atoms.length+1).padStart(5)}`);
  const neighbors = model.atoms.map(() => []);
  for (const b of model.bonds) for (let order=0; order<b.order; order++) {
    neighbors[b.a].push(b.b+1);
    neighbors[b.b].push(b.a+1);
  }
  for (let i=0; i<neighbors.length; i++) {
    const row = neighbors[i].sort((a,b)=>a-b);
    for (let j=0; j<row.length; j+=4) lines.push(`CONECT${String(i+1).padStart(5)}${row.slice(j,j+4).map(n=>String(n).padStart(5)).join('')}`);
  }
  lines.push('END');
  return lines.join('\n')+'\n';
}

/** V3000 is used at every size: it supports >999 atoms and bonds at DP30. */
export function toSDF(model) {
  validateExport(model);
  const lines = [model.metadata?.name || 'Polymer', '  PolyProg          3D',
    'Initial structure; no full-polymer minimization',
    '  0  0  0     0  0            999 V3000', 'M  V30 BEGIN CTAB',
    `M  V30 COUNTS ${model.atoms.length} ${model.bonds.length} 0 0 0`, 'M  V30 BEGIN ATOM'];
  model.atoms.forEach((a,i) => lines.push(`M  V30 ${i+1} ${a.element} ${a.x.toFixed(6)} ${a.y.toFixed(6)} ${a.z.toFixed(6)} 0${a.charge ? ` CHG=${a.charge}` : ''}`));
  lines.push('M  V30 END ATOM', 'M  V30 BEGIN BOND');
  model.bonds.forEach((b,i) => lines.push(`M  V30 ${i+1} ${b.order} ${b.a+1} ${b.b+1}`));
  lines.push('M  V30 END BOND', 'M  V30 END CTAB', 'M  END',
    '>  <REPEAT_COUNT>', String(model.stats.repeats), '',
    '>  <FORMAL_CHARGE>', String(model.stats.charge), '',
    '>  <FORMULA>', model.stats.formula, '',
    '>  <ATOM_LABELS>', model.atoms.map(a=>`${a.id}:${a.residue}:${a.name}`).join(' '), '',
    '>  <END_CAPS>', model.metadata?.endCapChemistry || '', '',
    '>  <END_CAP_ATOM_INDICES_ZERO_BASED>', JSON.stringify(model.metadata?.endCaps || {}), '',
    '>  <GEOMETRY_METHOD>', model.metadata?.geometry || 'Initial structure', '',
    '>  <SOURCE_SHA256>', JSON.stringify(model.metadata?.source || {}), '', '$$$$');
  return lines.join('\n')+'\n';
}
