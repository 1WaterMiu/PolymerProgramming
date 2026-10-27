#!/usr/bin/env python3
"""Render source-derived charged headgroups as font-independent SVG paths.

Example (RDKit 2026.09.1):
  PYTHONPATH=/home/watermiu/copolymer/.depict_packages \
    python3 scripts/generate-headgroups.py

The source templates are read only. Each group's original linker carbon C14
is replaced by an R attachment marker; no linker methylenes are drawn.
"""
from pathlib import Path
import argparse
import io
import json
import math
import re
from xml.etree import ElementTree

from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D

WIDTH, HEIGHT = 360, 180
INK = (32/255, 53/255, 46/255)  # #20352e, matching the website's primary ink.
GROUPS = [
    ('PA', 'primary-ammonium', '*[NH3+]'),
    ('QA', 'quaternary-ammonium', '*[N+](C)(C)C'),
    ('GD', 'guanidinium', '*NC(=[NH2+])N'),
    ('IL', 'imidazolium', '*[n+]1cc[nH]c1'),
]


def extract_headgroup(template):
    """Copy only the graph beyond N2, retaining source H counts and charges."""
    atoms, bonds = template['atoms'], template['bonds']
    names = {a['name']: i for i, a in enumerate(atoms)}
    nitrogen, linker = names['N2'], names['C14']
    neighbors = {i: set() for i in range(len(atoms))}
    for b in bonds:
        neighbors[b['a']].add(b['b'])
        neighbors[b['b']].add(b['a'])
    assert linker in neighbors[nitrogen]
    selected, pending = set(), [nitrogen]
    while pending:
        i = pending.pop()
        if i in selected or i == linker:
            continue
        selected.add(i)
        pending.extend(neighbors[i] - selected - {linker})
    graph = Chem.RWMol()
    wildcard = Chem.Atom(0)
    wildcard.SetNoImplicit(True)
    graph.AddAtom(wildcard)
    indices = {}
    for i in sorted(selected):
        original = atoms[i]
        atom = Chem.Atom(original['element'])
        atom.SetFormalCharge(original['charge'])
        atom.SetNoImplicit(True)
        indices[i] = graph.AddAtom(atom)
    for b in bonds:
        if b['a'] in selected and b['b'] in selected:
            graph.AddBond(indices[b['a']], indices[b['b']],
                          {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
                           3: Chem.BondType.TRIPLE}[b['order']])
    graph.AddBond(0, indices[nitrogen], Chem.BondType.SINGLE)
    molecule = graph.GetMol()
    Chem.SanitizeMol(molecule)
    molecule = Chem.RemoveHs(molecule)
    assert Chem.GetFormalCharge(molecule) == 1
    assert len(Chem.GetMolFrags(molecule)) == 1
    return molecule


def orient(molecule):
    rdDepictor.Compute2DCoords(molecule)
    conf = molecule.GetConformer()
    r = conf.GetAtomPosition(0)
    n_index = molecule.GetAtomWithIdx(0).GetNeighbors()[0].GetIdx()
    n = conf.GetAtomPosition(n_index)
    angle = -math.atan2(n.y-r.y, n.x-r.x)
    for i in range(molecule.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        x, y = p.x-r.x, p.y-r.y
        conf.SetAtomPosition(i, Point3D(x*math.cos(angle)-y*math.sin(angle),
                                       x*math.sin(angle)+y*math.cos(angle), 0))


def draw(molecule, kind):
    drawer = (rdMolDraw2D.MolDraw2DSVG if kind == 'svg' else rdMolDraw2D.MolDraw2DCairo)(WIDTH, HEIGHT)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.setAtomPalette({i: INK for i in range(119)})
    options.setSymbolColour(INK)
    options.setQueryColour(INK)
    options.explicitMethyl = True
    options.fixedFontSize = 37
    options.bondLineWidth = 2.8
    options.padding = 0.06
    options.additionalAtomLabelPadding = 0.04
    options.atomLabels[0] = 'R'
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    templates = json.loads((args.root/'data/templates.json').read_text())
    out = args.root/'assets/headgroups'
    out.mkdir(parents=True, exist_ok=True)
    provenance, previews = {}, []
    for code, filename, expected in GROUPS:
        source = templates[f'{code}_4']
        molecule = extract_headgroup(source)
        canonical = Chem.MolToSmiles(molecule)
        assert canonical == Chem.MolToSmiles(Chem.MolFromSmiles(expected)), (code, canonical)
        # Every carbon-length variant must contain the same charged headgroup.
        for carbons in (6, 8, 10):
            assert Chem.MolToSmiles(extract_headgroup(templates[f'{code}_{carbons}'])) == canonical
        orient(molecule)
        svg = draw(molecule, 'svg')
        tree = ElementTree.fromstring(svg)
        assert tree.attrib['viewBox'] == f'0 0 {WIDTH} {HEIGHT}'
        assert not any(e.tag.rsplit('}', 1)[-1] in ('text', 'image', 'foreignObject') for e in tree.iter())
        assert '<path' in svg
        assert not re.search(r'url\(\s*["\']?https?://', svg)
        (out/f'{filename}.svg').write_text(svg)
        png = draw(molecule, 'png')
        (out/f'{filename}.png').write_bytes(png)
        previews.append(png)
        provenance[filename] = dict(sourceTemplate=source['name'], source=source['source'],
            headgroupSmiles=canonical, formalCharge=1,
            formulaWithWildcard=rdMolDescriptors.CalcMolFormula(molecule),
            attachment='R replaces C14, directly bonded to original N2',
            viewBox=f'0 0 {WIDTH} {HEIGHT}', glyphs='SVG paths; no text nodes or external fonts',
            rdkitVersion=rdBase.rdkitVersion)
        print(filename, canonical, 'charge=+1; SVG glyphs=paths')
    (out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    # A small scientific drawing contact sheet is convenient for visual QA.
    from PIL import Image
    sheet = Image.new('RGBA', (WIDTH*2, HEIGHT*2), 'white')
    for i, data in enumerate(previews):
        png = Image.open(io.BytesIO(data)).convert('RGBA')
        sheet.alpha_composite(png, ((i % 2)*WIDTH, (i // 2)*HEIGHT))
    sheet.convert('RGB').save(out/'preview.png')


if __name__ == '__main__':
    main()
