#!/usr/bin/env python3
"""Import verified structures without modifying the source files.

Run with RDKit 2026.09.1, for example:
  PYTHONPATH=/home/watermiu/copolymer/.depict_packages python3 scripts/import-templates.py \
    --source /home/watermiu/copolymer
The saved conformers are starting geometries, not minimized polymer chains.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
from itertools import combinations
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdDepictor, rdMolDescriptors, rdMolTransforms
from rdkit.Chem.Draw import rdMolDraw2D

GROUPS = {'PA': 'Primary ammonium', 'QA': 'Quaternary ammonium',
          'GD': 'Guanidinium', 'IL': 'Imidazolium'}


def load_source(pdb, mol2):
    atoms, edges = [], set()
    for line in pdb.read_text().splitlines():
        if line.startswith(('ATOM  ', 'HETATM')):
            raw = line[78:80].strip()
            charge = (int(raw[:-1]) * (1 if raw[-1] == '+' else -1)) if raw else 0
            atoms.append(dict(id=len(atoms)+1, name=line[12:16].strip(),
                              sourceSerial=int(line[6:11]), element=line[76:78].strip(),
                              charge=charge, residue=1, residueName=line[17:20].strip()))
        elif line.startswith('CONECT'):
            ids = [int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip()]
            edges.update(tuple(sorted((ids[0], j))) for j in ids[1:])
    index = {a['sourceSerial']: i for i, a in enumerate(atoms)}
    bonds, elements, section = [], {}, ''
    for line in mol2.read_text().splitlines():
        if line.startswith('@<TRIPOS>'):
            section = line.strip()
            continue
        fields = line.split()
        if not fields:
            continue
        if section == '@<TRIPOS>ATOM':
            elements[int(fields[0])] = fields[5].split('.')[0]
        elif section == '@<TRIPOS>BOND':
            bonds.append(dict(a=index[int(fields[1])], b=index[int(fields[2])],
                              order=1 if fields[3] == 'am' else int(fields[3])))
    assert all(a['element'] == elements[a['sourceSerial']] for a in atoms)
    assert edges == {tuple(sorted((atoms[b['a']]['sourceSerial'], atoms[b['b']]['sourceSerial']))) for b in bonds}
    return atoms, bonds


def make_template(root, out, group, carbons):
    key = f'{group}_{carbons}'
    pdb, mol2 = (root/key/(key+ext) for ext in ('.pdb', '.mol2'))
    atoms, bonds = load_source(pdb, mol2)
    editable = Chem.RWMol()
    for a in atoms:
        atom = Chem.Atom(a['element'])
        atom.SetNoImplicit(True)
        atom.SetFormalCharge(a['charge'])
        editable.AddAtom(atom)
    for b in bonds:
        editable.AddBond(b['a'], b['b'], {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
                                         3: Chem.BondType.TRIPLE}[b['order']])
    mol = editable.GetMol()
    Chem.SanitizeMol(mol)
    by_name = {a['name']: i for i, a in enumerate(atoms)}
    attach = dict(head=by_name['C1'], tail=by_name['C4'], headCap=by_name['H1'],
                  tailCap=by_name['H7'], methyl=by_name['C5'])
    params = AllChem.ETKDGv3()
    params.randomSeed = 20261027
    params.useRandomCoords = True
    conformers = AllChem.EmbedMultipleConfs(mol, numConfs=24, params=params)
    assert conformers, key
    candidates = []
    for conf_id in conformers:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
        status = ff.Minimize(maxIts=3000, energyTol=1e-7, forceTol=1e-5)
        if status:
            continue
        conf = mol.GetConformer(conf_id)
        valid = True
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.GetDegree() == 4:
                for a, b in combinations([n.GetIdx() for n in atom.GetNeighbors()], 2):
                    angle = rdMolTransforms.GetAngleDeg(conf, a, atom.GetIdx(), b)
                    if not 95 < angle < 125:
                        valid = False
        for bond in mol.GetBonds():
            length = rdMolTransforms.GetBondLength(conf, bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            is_h = any(a.GetSymbol() == 'H' for a in (bond.GetBeginAtom(), bond.GetEndAtom()))
            if not ((0.9 < length < 1.25) if is_h else (1.15 < length < 1.8)):
                valid = False
        if not valid:
            continue
        # C1 hydrogens are equivalent. Select its most outward-facing H as the
        # named H1 cap without changing chemistry or stereochemistry.
        cap_choices = [n.GetIdx() for n in mol.GetAtomWithIdx(attach['head']).GetNeighbors() if n.GetSymbol() == 'H']
        cap = max(cap_choices, key=lambda h: (conf.GetAtomPosition(h)-conf.GetAtomPosition(attach['tailCap'])).Length())
        span = (conf.GetAtomPosition(cap)-conf.GetAtomPosition(attach['tailCap'])).Length()
        candidates.append((span, conf_id, cap, ff.CalcEnergy()))
    assert candidates, (key, 'No converged conformer passed geometry validation')
    span, selected_id, selected_cap, selected_energy = max(candidates)
    conf = mol.GetConformer(selected_id)
    p1, p2 = conf.GetAtomPosition(attach['headCap']), conf.GetAtomPosition(selected_cap)
    conf.SetAtomPosition(attach['headCap'], p2)
    conf.SetAtomPosition(selected_cap, p1)
    for i, atom in enumerate(atoms):
        p = conf.GetAtomPosition(i)
        atom.update(x=round(p.x, 7), y=round(p.y, 7), z=round(p.z, 7))
    # Illustrations use the same verified chemical graph, without explicit H.
    drawing = Chem.RemoveHs(mol)
    rdDepictor.Compute2DCoords(drawing)
    drawer = rdMolDraw2D.MolDraw2DSVG(1000, 520)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = 0.08
    options.bondLineWidth = 2.1
    options.baseFontSize = 0.8
    drawer.DrawMolecule(drawing)
    drawer.FinishDrawing()
    (out/'assets/templates'/f'{key}.svg').write_text(drawer.GetDrawingText())
    sources = {}
    for file in (pdb, mol2):
        shutil.copyfile(file, out/'data/source'/file.name)
        sources[file.suffix[1:]] = dict(path=f'data/source/{file.name}',
            sha256=hashlib.sha256(file.read_bytes()).hexdigest())
    return dict(name=key, group=group, groupName=GROUPS[group], carbons=carbons,
                atoms=atoms, bonds=bonds, attachment=attach,
                formula=rdMolDescriptors.CalcMolFormula(mol), charge=sum(a['charge'] for a in atoms),
                source=sources, image=f'assets/templates/{key}.svg',
                geometry=dict(method='RDKit ETKDGv3; unconstrained UFF; extended validated conformer selected from 24 candidates',
                              rdkitVersion=rdBase.rdkitVersion, randomSeed=20261027,
                              optimizationConverged=True, energyKcalMol=round(selected_energy, 6),
                              attachmentSpanAngstrom=round(span, 6), selectedConformer=selected_id,
                              stereochemistry='Representative conformer; tacticity is not specified'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    for path in ('data/source', 'assets/templates'):
        (args.output/path).mkdir(parents=True, exist_ok=True)
    templates = {}
    for group in GROUPS:
        for carbons in (4, 6, 8, 10):
            key = f'{group}_{carbons}'
            templates[key] = make_template(args.source, args.output, group, carbons)
            print(key, len(templates[key]['atoms']), templates[key]['formula'], flush=True)
    (args.output/'data/templates.json').write_text(json.dumps(templates, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    main()
