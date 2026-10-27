# Polymer Programming

A local-first, English-language molecular workspace for assembling cationic polymers from 16 verified PDB/MOL2 building blocks.

Choose **primary ammonium, quaternary ammonium, guanidinium, or imidazolium**, choose **4, 6, 8, or 10 alkyl carbons**, and connect **1–30 repeat units** with **CH3 caps at both ends**. The selected monomer's 2D structure is always shown above the interactive 3D polymer, so both structures can be inspected together. Download the polymer as PDB, SDF, or a JSON build recipe.

Headgroup selectors use fixed-aspect-ratio SVG images generated from the source structures. Every chemical glyph is an outlined path, preventing font substitution or wrapping from deforming the formulas. The R label indicates attachment to the side chain. The interface, shareable URLs, and download filenames use full chemical names; original source IDs remain unchanged for provenance. Older abbreviation-based URLs remain compatible, and repeat counts above 30 in saved URLs are clamped to 30.

## Run locally

On Windows, double-click **`Start-Website.cmd`**. The launcher finds Node.js on PATH or the existing Codex Node runtime, starts the local server, and opens your browser. Keep its terminal window open while using the website. Press Enter in that window to stop the server. If this website is already running on port 4173, the launcher simply opens it.

Alternatively, with Node.js 20 or newer:

```sh
node server.mjs
```

Open **http://127.0.0.1:4173**. No package installation, API key, account, internet connection, Python, or WSL is needed to run the website. The server binds only to the local loopback address. Set the `PORT` environment variable to use a different port.

Use the HTTP server rather than opening `index.html` directly: browser ES modules, template fetching, and the background builder require HTTP. All model generation happens inside the browser. Generated structures are not uploaded.

## Available chemistry

| Headgroup | Terminal group | Alkyl carbon options |
| --- | --- | --- |
| Primary ammonium | NH3+ | 4, 6, 8, 10 |
| Quaternary ammonium | N+(CH3)3 | 4, 6, 8, 10 |
| Guanidinium | NH–C(=NH2+)–NH2 | 4, 6, 8, 10 |
| Imidazolium | Cationic five-membered aromatic ring | 4, 6, 8, 10 |

The carbon count is the **total number of carbons after the ester oxygen**: C4 means `O–(CH2)3–CH3`. Each building block has an ester-linked alkyl side chain and an amide-linked cationic side chain with a three-carbon spacer. Each repeat contributes one positive formal charge.

The current builder produces a **single linear chain made from one selected template**. Mixed sequences, multiple chains, crosslinks, counterions, solvent, tacticity controls, and molecular dynamics topology generation are not implemented.

## How assembly works

The original PDB files define atom identities, hydrogens, formal charges, and connectivity. Matching MOL2 files supply bond orders. Their atom identities and complete bond graphs are checked against one another during template import. The original files are copied unchanged to `data/source/`; SHA-256 hashes are recorded in `data/templates.json`.

The supplied coordinates are molecular sketches. The importer regenerates 3D conformers with RDKit ETKDG, relaxes individual templates using UFF, and rejects templates with abnormal bond lengths or saturated-carbon angles. This preprocessing is already done and the validated coordinates are included in the repository.

For each successive repeat, the browser builder:

1. Removes the designated cap hydrogen on C4 of the previous unit and on C1 of the new unit.
2. Preserves the C5 pendant methyl and every headgroup, alkyl side chain, and other atom.
3. Rigidly places the new template using the cap-bond directions and adds a 1.53 Å C4–C1 bond.
4. Searches candidate torsions around the new bond to reduce close contacts and favor an extended starting chain.

After the repeat units are connected, the two remaining outer cap hydrogens are replaced by **new, neutral CH3 groups**: one bonds to the first unit's C1 and one to the last unit's C4. The original C5 pendant methyl remains intact. This also applies to a single-repeat build. The original source files and monomer reference diagrams stay unchanged and H-capped.

For a source template with `A` atoms and `B` bonds, the methyl-capped chain contains `n*A − 2*n + 8` atoms, `n*B − n + 7` bonds, and formal charge `+n`. Relative to the earlier H-capped chain, the difference is `C2H4` (+6 atoms and +6 bonds). The builder replaces 2*n source cap hydrogens in total and adds two carbons plus six methyl hydrogens. Within every repeat, C1 is CH2 and C4 has no hydrogen.

The head methyl atoms are named `CM1/HM11/HM12/HM13` in the first residue; the tail methyl atoms are `CM2/HM21/HM22/HM23` in the last residue. Methyl-cap bonds and hydrogen geometry are explicitly constructed. End-cap atom indices are recorded in the build metadata and the exported PDB/SDF includes all cap atoms and bonds.

**These are initial assembled coordinates, not energy-minimized polymer conformations.** A valid molecule graph does not establish an equilibrium conformation or a simulation-ready force field. Tacticity is not controlled, no counterions or solvent are added, and source monomer partial charges/parameters are not reused as a polymer topology. Minimize and parameterize the final chain before molecular dynamics.

## Viewer and exports

- Drag to rotate; Shift-drag or right-drag to pan; scroll or pinch to zoom.
- Use the hydrogen toggle, ball-and-stick/wireframe selector, automatic rotation, reset, and expanded view.
- Keyboard: arrow keys rotate, `+`/`-` zoom, and `0` resets the focused viewer.
- The **Monomer structure** panel continuously displays the selected original H-capped reference above the **CH3-capped polymer**. Select its diagram to open a larger view, or download the source PDB.
- **PDB** exports all atoms, coordinates in angstroms, residue identities, formal charge fields, and CONECT records. Repeated CONECT targets encode multiple bonds using the convention supported by RDKit and other molecular tools; SDF provides an explicit bond-order format.
- **SDF V3000** preserves explicit bond orders and supports models with more than 999 atoms.
- **JSON recipe** records settings, counts, methyl end caps, assembly notes, and source provenance. It is a record of the build; importing recipes is not currently supported. Polymer download filenames include `_CH3` to distinguish these structures from earlier H-capped exports.

Changing a selection does not silently replace the built polymer. Click **Build polymer** to apply changes. Downloads always correspond to the last successfully built polymer, and the preview status indicates pending changes. The URL records the last successful configuration for easy reopening.

## Validate

```sh
node --test tests/*.test.mjs
```

Tests cover all 16 templates, polymer connectivity and formula, headgroup protonation, charge, exact end-cap handling, pendant methyl retention, alkyl lengths, ordinary covalent geometry, inter-unit bond lengths/angles, determinism, invalid inputs, source hashes, and fixed-column PDB / V3000 SDF exports. Boundary cases verify that 30 repeat units succeed and 31 are rejected. The GitHub workflow runs these checks on pushes and pull requests.

## Repository structure

```text
index.html                    English user interface
styles.css                    Responsive workspace styling
src/app.js                    Configuration, UI state, and downloads
src/builder-worker.js         Background molecular assembly
src/polymer.js                Molecular graph assembly and export
src/viewer.js                 Dependency-free interactive 3D canvas viewer
data/templates.json          Validated templates and provenance
data/source/                  The 16 original PDB/MOL2 pairs
assets/templates/             Source-derived 2D structure SVGs
assets/headgroups/             Fixed-aspect headgroup structure images
scripts/import-templates.py   Reproducible template preparation
scripts/generate-headgroups.py Source-derived headgroup image generation
tests/                        Chemical and export validation
server.mjs                    Small local-only static HTTP server
Start-Website.cmd             Windows launcher
```

The static site can also be served under a repository subpath, such as a GitHub Pages project site. All application asset URLs are relative. No hosted deployment is required for local use.
