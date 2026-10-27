import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { buildPolymer, getTemplate, toPDB, toSDF } from '../src/polymer.js';

// Independent chemical expectations, counted from the source PDB CONECT records.
const templates = JSON.parse(readFileSync(new URL('../data/templates.json', import.meta.url), 'utf8'));
const families = {
  PA: { carbonOffset: 10, hydrogenOffset: 21, nitrogens: 2, rings: 0 },
  QA: { carbonOffset: 13, hydrogenOffset: 27, nitrogens: 2, rings: 0 },
  GD: { carbonOffset: 11, hydrogenOffset: 23, nitrogens: 4, rings: 0 },
  IL: { carbonOffset: 13, hydrogenOffset: 22, nitrogens: 3, rings: 1 },
};

function graph(model) {
  const adjacency = model.atoms.map(() => []);
  const unique = new Set();
  for (const bond of model.bonds) {
    assert.ok(Number.isInteger(bond.a) && Number.isInteger(bond.b));
    assert.ok(bond.a >= 0 && bond.b >= 0 && bond.a < model.atoms.length && bond.b < model.atoms.length);
    assert.notEqual(bond.a, bond.b, 'No atom may bond to itself');
    assert.ok([1, 2, 3].includes(bond.order));
    const key = [bond.a, bond.b].sort((a, b) => a - b).join(':');
    assert.ok(!unique.has(key), `Duplicate covalent bond ${key}`);
    unique.add(key);
    adjacency[bond.a].push({ index: bond.b, order: bond.order });
    adjacency[bond.b].push({ index: bond.a, order: bond.order });
  }
  const visited = new Set();
  const pending = [0];
  while (pending.length) {
    const index = pending.pop();
    if (visited.has(index)) continue;
    visited.add(index);
    pending.push(...adjacency[index].map(({ index: neighbor }) => neighbor));
  }
  assert.equal(visited.size, model.atoms.length, 'The polymer must be one connected molecule');
  return adjacency;
}

function atomIndex(model, residue, name) {
  const matches = model.atoms.flatMap((atom, index) => atom.residue === residue && atom.name === name ? [index] : []);
  assert.equal(matches.length, 1, `Expected one ${name} in residue ${residue}`);
  return matches[0];
}

function hydrogenCount(model, adjacency, index) {
  return adjacency[index].filter(({ index: neighbor }) => model.atoms[neighbor].element === 'H').length;
}

function angle(a, center, b) {
  const u = [a.x - center.x, a.y - center.y, a.z - center.z];
  const v = [b.x - center.x, b.y - center.y, b.z - center.z];
  const cosine = u.reduce((sum, value, i) => sum + value * v[i], 0) / (Math.hypot(...u) * Math.hypot(...v));
  return Math.acos(Math.max(-1, Math.min(1, cosine))) * 180 / Math.PI;
}

test('All source-derived 3D templates have ordinary covalent distances and tetrahedral saturated carbons', () => {
  for (const [name, template] of Object.entries(templates)) {
    const adjacency = graph(template);
    for (const { a, b } of template.bonds) {
      const first = template.atoms[a];
      const second = template.atoms[b];
      const distance = Math.hypot(first.x - second.x, first.y - second.y, first.z - second.z);
      const hydrogenBond = first.element === 'H' || second.element === 'H';
      assert.ok(distance >= (hydrogenBond ? 0.85 : 1.1) && distance <= (hydrogenBond ? 1.25 : 1.8), `${name} ${first.name}-${second.name}: unexpected covalent distance ${distance}`);
    }
    for (const [index, atom] of template.atoms.entries()) {
      if (atom.element !== 'C' || adjacency[index].length !== 4) continue;
      const neighbors = adjacency[index].map(({ index: neighbor }) => template.atoms[neighbor]);
      for (let a = 0; a < neighbors.length; a += 1) {
        for (let b = a + 1; b < neighbors.length; b += 1) {
          const value = angle(neighbors[a], atom, neighbors[b]);
          assert.ok(value >= 90 && value <= 130, `${name} ${neighbors[a].name}-${atom.name}-${neighbors[b].name}: saturated carbon angle ${value} is not approximately tetrahedral`);
        }
      }
    }
  }
});

for (const [group, expected] of Object.entries(families)) {
  for (const carbons of [4, 6, 8, 10]) {
    for (const repeats of [1, 2, 7]) {
      test(`${group}_${carbons}, ${repeats} repeats: formula, connectivity, charge, caps and junction geometry`, () => {
        const template = getTemplate(templates, group, carbons);
        assert.ok(template);
        const model = buildPolymer(templates, { group, carbons, repeats });
        const adjacency = graph(model);
        const counts = model.atoms.reduce((result, atom) => {
          result[atom.element] = (result[atom.element] || 0) + 1;
          for (const coordinate of [atom.x, atom.y, atom.z]) assert.ok(Number.isFinite(coordinate), 'Coordinates must be finite');
          return result;
        }, {});
        assert.deepEqual(counts, {
          C: repeats * (carbons + expected.carbonOffset),
          H: repeats * (2 * carbons + expected.hydrogenOffset) - 2 * (repeats - 1),
          N: repeats * expected.nitrogens,
          O: repeats * 3,
        });
        assert.equal(model.bonds.length - model.atoms.length + 1, expected.rings * repeats, 'Cycle count must preserve each imidazolium ring');
        assert.equal(model.atoms.reduce((sum, atom) => sum + atom.charge, 0), repeats);
        assert.equal(model.stats.charge, repeats);
        assert.equal(model.stats.atoms, model.atoms.length);
        assert.equal(model.stats.bonds, model.bonds.length);
        assert.equal(model.stats.repeats, repeats);

        for (const [index, atom] of model.atoms.entries()) {
          const valence = adjacency[index].reduce((sum, bond) => sum + bond.order, 0);
          const expectedValence = { H: 1, C: 4, O: 2, N: atom.charge === 1 ? 4 : 3 }[atom.element];
          assert.equal(valence, expectedValence, `Bad valence at ${atom.residue}:${atom.name}`);
        }

        for (let residue = 1; residue <= repeats; residue += 1) {
          const indices = Object.fromEntries(['C1', 'C2', 'C3', 'C4', 'C5', 'C11', 'C12', 'C13', 'C14', 'N1', 'N2', 'O2'].map(name => [name, atomIndex(model, residue, name)]));
          assert.equal(hydrogenCount(model, adjacency, indices.C1), residue === 1 ? 3 : 2, 'Only the head cap loses a hydrogen');
          assert.equal(hydrogenCount(model, adjacency, indices.C4), residue === repeats ? 1 : 0, 'Only the tail cap loses a hydrogen');
          assert.equal(hydrogenCount(model, adjacency, indices.C5), 3, 'The pendant methyl must be retained');
          assert.equal(hydrogenCount(model, adjacency, indices.N1), 1, 'The amide retains its hydrogen');
          for (const carbon of ['C12', 'C13', 'C14']) assert.equal(hydrogenCount(model, adjacency, indices[carbon]), 2, 'The cation linker has three CH2 groups');

          let previous = indices.O2;
          let current = adjacency[previous].find(({ index }) => model.atoms[index].element === 'C' && hydrogenCount(model, adjacency, index) > 0)?.index;
          assert.notEqual(current, undefined, 'The ester oxygen must connect to the alkyl tail');
          const tail = [];
          while (current !== undefined) {
            assert.equal(model.atoms[current].element, 'C');
            tail.push(current);
            const onward = adjacency[current].filter(({ index }) => index !== previous && model.atoms[index].element !== 'H');
            assert.ok(onward.length <= 1, 'The alkyl tail must be unbranched');
            previous = current;
            current = onward[0]?.index;
            assert.ok(tail.length <= carbons, 'Unexpected tail cycle or length');
          }
          assert.equal(tail.length, carbons, 'The suffix counts all alkyl carbons, including terminal CH3');
          for (let index = 0; index < tail.length; index += 1) assert.equal(hydrogenCount(model, adjacency, tail[index]), index === tail.length - 1 ? 3 : 2);

          if (group === 'PA') assert.equal(hydrogenCount(model, adjacency, indices.N2), 3);
          if (group === 'QA') {
            assert.equal(hydrogenCount(model, adjacency, indices.N2), 0);
            const methyls = adjacency[indices.N2].filter(({ index }) => model.atoms[index].element === 'C' && hydrogenCount(model, adjacency, index) === 3);
            assert.equal(methyls.length, 3, 'Quaternary ammonium needs three methyl substituents');
          }
          if (group === 'GD') {
            assert.equal(hydrogenCount(model, adjacency, indices.N2), 1);
            for (const nitrogen of ['N3', 'N4']) assert.equal(hydrogenCount(model, adjacency, atomIndex(model, residue, nitrogen)), 2);
          }
          if (group === 'IL') {
            assert.equal(hydrogenCount(model, adjacency, indices.N2), 0);
            assert.equal(hydrogenCount(model, adjacency, atomIndex(model, residue, 'N3')), 1);
            for (const carbon of ['C15', 'C16', 'C17']) assert.equal(hydrogenCount(model, adjacency, atomIndex(model, residue, carbon)), 1);
          }

          if (residue < repeats) {
            const nextHead = atomIndex(model, residue + 1, 'C1');
            const nextCarbon = atomIndex(model, residue + 1, 'C2');
            const junction = adjacency[indices.C4].find(({ index }) => index === nextHead);
            assert.equal(junction?.order, 1, 'Polymerization must connect C4 to the next C1');
            const a = model.atoms[indices.C4];
            const b = model.atoms[nextHead];
            const distance = Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
            assert.ok(distance >= 1.4 && distance <= 1.65, `Unphysical junction bond length ${distance}`);
            for (const neighbor of [indices.C3, indices.C5, indices.C11]) {
              const value = angle(model.atoms[neighbor], a, b);
              assert.ok(value >= 95 && value <= 125, `Junction C4 bond angle ${value} must be approximately tetrahedral`);
            }
            const headAngle = angle(a, b, model.atoms[nextCarbon]);
            assert.ok(headAngle >= 95 && headAngle <= 125, `Junction C1 bond angle ${headAngle} must be approximately tetrahedral`);
          }
        }

        const crossResidue = model.bonds.filter(({ a, b }) => model.atoms[a].residue !== model.atoms[b].residue);
        assert.equal(crossResidue.length, repeats - 1, 'A linear polymer has one covalent junction per adjacent pair');
      });
    }
  }
}

test('PDB and SDF exports retain all atoms for a multi-repeat charged polymer', () => {
  const model = buildPolymer(templates, { group: 'IL', carbons: 10, repeats: 7 });
  const pdb = toPDB(model);
  const atomLines = pdb.split(/\r?\n/).filter(line => /^(ATOM  |HETATM)/.test(line));
  assert.equal(atomLines.length, model.atoms.length);
  assert.equal(atomLines.filter(line => line.slice(78, 80).trim() === '1+').length, 7);
  assert.ok(pdb.includes('CONECT'), 'PDB needs explicit polymer connectivity');
  const exportedEdges = new Set();
  const exportedOrder = new Map();
  for (const line of pdb.split(/\r?\n/).filter(line => line.startsWith('CONECT'))) {
    const indices = [];
    for (let offset = 6; offset < line.length; offset += 5) {
      const value = line.slice(offset, offset + 5).trim();
      if (value) indices.push(Number(value) - 1);
    }
    for (const neighbor of indices.slice(1)) {
      exportedEdges.add([indices[0], neighbor].sort((a, b) => a - b).join(':'));
      const direction = `${indices[0]}:${neighbor}`;
      exportedOrder.set(direction, (exportedOrder.get(direction) || 0) + 1);
    }
  }
  const expectedEdges = new Set(model.bonds.map(({ a, b }) => [a, b].sort((first, second) => first - second).join(':')));
  assert.deepEqual(exportedEdges, expectedEdges, 'PDB must retain every covalent edge, including all repeat junctions');
  for (const { a, b, order } of model.bonds) {
    assert.equal(exportedOrder.get(`${a}:${b}`), order, 'Repeated PDB CONECT targets must preserve double-bond order');
    assert.equal(exportedOrder.get(`${b}:${a}`), order);
  }
  const sdf = toSDF(model);
  assert.ok(sdf.includes('M  END'));
  assert.ok(sdf.includes('$$$$'));
  assert.ok(sdf.includes('CHG'), 'SDF must preserve formal cation charges');
});

test('Long polymers use an SDF format that supports more than 999 atoms', () => {
  const model = buildPolymer(templates, { group: 'QA', carbons: 10, repeats: 20 });
  assert.ok(model.atoms.length > 999);
  const sdf = toSDF(model);
  assert.ok(sdf.includes('V3000'), 'V2000 count fields cannot represent this polymer');
  assert.match(sdf, new RegExp(`M  V30 COUNTS ${model.atoms.length} ${model.bonds.length} `));
  assert.ok(sdf.includes('M  END'));
});
