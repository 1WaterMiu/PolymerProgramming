import { MoleculeViewer } from './viewer.js';
import { buildPolymer, toPDB, toSDF, MAX_REPEATS } from './polymer.js';

const $ = id => document.getElementById(id);
const groups = {
  PA: { name: 'Primary ammonium', slug: 'primary-ammonium', description: 'A compact, protonated amine.' },
  QA: { name: 'Quaternary ammonium', slug: 'quaternary-ammonium', description: 'Three methyls. A permanent charge.' },
  GD: { name: 'Guanidinium', slug: 'guanidinium', description: 'A resonance-stabilized cation.' },
  IL: { name: 'Imidazolium', slug: 'imidazolium', description: 'An aromatic, five-membered ring.' },
};
const lengths = [4, 6, 8, 10];
const search = new URLSearchParams(location.search);
const requestedGroup = search.get('group');
const initialGroup = Object.keys(groups).find(group => groups[group].slug === requestedGroup) || (Object.hasOwn(groups, requestedGroup) ? requestedGroup : 'PA');
const requestedRepeats = Number(search.get('repeats'));
const state = {
  group: initialGroup,
  carbons: lengths.includes(Number(search.get('carbons'))) ? Number(search.get('carbons')) : 4,
  repeats: Number.isInteger(requestedRepeats) && requestedRepeats >= 1 ? Math.min(requestedRepeats, MAX_REPEATS) : 12,
  templates: null, model: null, builtConfig: null, busy: false, error: null,
};
let viewer, worker, job = 0, toastTimer;
const pending = new Map();
const format = new Intl.NumberFormat('en-US');
const config = () => ({ group: state.group, carbons: state.carbons, repeats: state.repeats });
const key = () => `${state.group}_${state.carbons}`;

function notify(message) {
  clearTimeout(toastTimer);
  $('toast').textContent = message;
  $('toast').hidden = false;
  toastTimer = setTimeout(() => { $('toast').hidden = true; }, 3500);
}

function setStatus(message, dirty = false) {
  $('model-status').replaceChildren();
  const dot = document.createElement('span'); dot.className = 'green-dot';
  $('model-status').append(dot, document.createTextNode(message));
  $('model-status').classList.toggle('dirty', dirty);
}

function updateTitle() {
  const current = state.builtConfig || config();
  $('model-title').replaceChildren(document.createTextNode(groups[current.group].name));
  const detail = document.createElement('span');
  detail.textContent = `C${current.carbons} · ${current.repeats} ${current.repeats === 1 ? 'repeat' : 'repeats'}`;
  $('model-title').append(detail);
}

function updateSelection() {
  document.querySelectorAll('[data-group]').forEach(button => {
    const active = button.dataset.group === state.group;
    button.classList.toggle('selected', active); button.setAttribute('aria-pressed', active);
  });
  document.querySelectorAll('[data-carbons]').forEach(button => {
    const active = Number(button.dataset.carbons) === state.carbons;
    button.classList.toggle('active', active); button.setAttribute('aria-pressed', active);
  });
  document.querySelectorAll('[data-library-group]').forEach(card => card.classList.toggle('current', card.dataset.libraryGroup === state.group));
  document.querySelectorAll('[data-template]').forEach(button => {
    const active = button.dataset.template === key();
    button.classList.toggle('active', active); button.setAttribute('aria-pressed', active);
  });
  $('selected-template').textContent = `${groups[state.group].name} · C${state.carbons}`;
  const imagePath = `./assets/templates/${key()}.svg`;
  if ($('template-image').getAttribute('src') !== imagePath) $('template-image').src = imagePath;
  $('template-image').alt = `${groups[state.group].name}, C${state.carbons}: original hydrogen-capped monomer reference`;
  $('template-caption').textContent = 'Original H-capped reference';
  $('monomer-title').textContent = groups[state.group].name;
  $('monomer-name').textContent = `C${state.carbons} alkyl side chain`;
  $('monomer-image-link').href = imagePath;
  $('monomer-source-link').href = `./data/source/${key()}.pdb`;
  $('monomer-source-link').download = `${groups[state.group].slug}_C${state.carbons}_source.pdb`;
  $('source-pdb-link').href = `./data/source/${key()}.pdb`;
  $('source-pdb-link').download = `${groups[state.group].slug}_C${state.carbons}_source.pdb`;
  $('source-pdb-link').textContent = 'Download original monomer PDB ↗';
  const dirty = state.builtConfig && JSON.stringify(config()) !== JSON.stringify(state.builtConfig);
  if (!state.busy) {
    if (state.error) setStatus('Build failed', true);
    else if (dirty) setStatus('Build to apply changes', true);
    else if (state.model) setStatus('Ready');
  }
  updateTitle();
}

function fillLibrary() {
  for (const [group, info] of Object.entries(groups)) {
    const card = document.createElement('article');
    card.className = 'library-card'; card.dataset.libraryGroup = group;
    card.innerHTML = `<div class="library-card-top"><img class="library-headgroup-image" src="./assets/headgroups/${info.slug}.svg" alt="${info.name} headgroup"><span class="library-charge">+1 e / unit</span></div><h3>${info.name}</h3><p class="library-description">${info.description}</p><div class="library-options" role="group" aria-label="${info.name} templates"></div>`;
    for (const carbons of lengths) {
      const button = document.createElement('button'); button.type = 'button';
      button.textContent = `C${carbons}`; button.dataset.template = `${group}_${carbons}`;
      button.setAttribute('aria-label', `Select ${info.name}, ${carbons} alkyl carbons`);
      button.addEventListener('click', () => {
        state.group = group; state.carbons = carbons; updateSelection();
        $('workspace').scrollIntoView({ behavior: matchMedia('(prefers-reduced-motion: reduce)').matches ? 'instant' : 'smooth' });
        notify(`${info.name} · C${carbons} selected. Build to update the polymer.`);
      });
      card.querySelector('.library-options').append(button);
    }
    $('library-grid').append(card);
  }
}

function requestBuild(settings) {
  if (!worker) return Promise.resolve().then(() => buildPolymer(state.templates, settings));
  const id = ++job;
  return new Promise((resolve, reject) => {
    pending.set(id, { resolve, reject });
    worker.postMessage({ type: 'build', id, config: settings });
  });
}

async function generate() {
  if (!state.templates || state.busy) return;
  state.repeats = Number($('repeat-count').value);
  if (!Number.isInteger(state.repeats) || state.repeats < 1 || state.repeats > MAX_REPEATS) {
    $('repeat-error').hidden = false;
    $('repeat-count').setAttribute('aria-invalid', 'true');
    $('repeat-count').setAttribute('aria-describedby', 'repeat-error');
    $('repeat-count').focus();
    return;
  }
  $('repeat-error').hidden = true;
  $('repeat-count').removeAttribute('aria-invalid');
  $('repeat-count').removeAttribute('aria-describedby');
  $('repeat-slider').value = state.repeats;
  const requested = config();
  state.busy = true; state.error = null;
  $('build-button').disabled = true; $('build-button-label').textContent = 'Building your polymer…';
  $('stage-loading').hidden = false;
  $('stage-loading').querySelector('p').textContent = `Connecting ${requested.repeats} ${groups[requested.group].name.toLowerCase()} · C${requested.carbons} repeat units`;
  $('stage-error').hidden = true;
  setStatus('Building');
  try {
    const model = await requestBuild(requested);
    state.model = model; state.builtConfig = requested;
    viewer?.setModel(model);
    $('stat-repeats').textContent = format.format(requested.repeats);
    $('stat-atoms').textContent = format.format(model.atoms.length);
    $('stat-bonds').textContent = format.format(model.bonds.length);
    $('stat-charge').textContent = `${model.stats.charge > 0 ? '+' : ''}${model.stats.charge}`;
    const formula = model.stats.formula.replace(/\^\d+[+-]$|[+-]$/, '').replace(/\d/g, digit => '₀₁₂₃₄₅₆₇₈₉'[Number(digit)]);
    $('formula-display').textContent = `${formula} · ${requested.repeats} positive ${requested.repeats === 1 ? 'charge' : 'charges'}`;
    ['download-pdb', 'download-sdf', 'download-recipe'].forEach(id => { $(id).disabled = false; });
    const address = new URL(location.href);
    for (const [name, value] of Object.entries(requested)) address.searchParams.set(name, name === 'group' ? groups[value].slug : value);
    history.replaceState(null, '', address);
    viewer?.resize();
    setStatus('Ready');
  } catch (error) {
    state.error = error.message;
    $('stage-error').textContent = `The structure could not be built. ${error.message}`;
    $('stage-error').hidden = false;
    setStatus('Build failed', true);
  } finally {
    state.busy = false;
    $('stage-loading').hidden = true;
    $('build-button').disabled = false; $('build-button-label').textContent = 'Build polymer';
    updateSelection();
  }
}

function download(content, extension, mime) {
  if (!state.model) return;
  const built = state.builtConfig;
  const filename = `${groups[built.group].slug}_C${built.carbons}_DP${built.repeats}_CH3.${extension}`;
  const url = URL.createObjectURL(new Blob([content], { type: mime }));
  const link = document.createElement('a'); link.href = url; link.download = filename;
  document.body.append(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 30000);
  notify(`${filename} downloaded`);
}

function exportSafely(extension, mime, create) {
  if (!state.model) return;
  try { download(create(state.model), extension, mime); }
  catch (error) { notify(`Export failed: ${error.message}`); }
}

$('builder-form').addEventListener('submit', event => { event.preventDefault(); generate(); });
$('group-options').addEventListener('click', event => {
  const button = event.target.closest('[data-group]'); if (!button) return;
  state.group = button.dataset.group; updateSelection();
});
$('carbon-options').addEventListener('click', event => {
  const button = event.target.closest('[data-carbons]'); if (!button) return;
  state.carbons = Number(button.dataset.carbons); updateSelection();
});
$('repeat-count').addEventListener('input', () => {
  const value = Number($('repeat-count').value);
  if (Number.isInteger(value) && value >= 1 && value <= MAX_REPEATS) { state.repeats = value; $('repeat-slider').value = value; $('repeat-error').hidden = true; $('repeat-count').removeAttribute('aria-invalid'); $('repeat-count').removeAttribute('aria-describedby'); updateSelection(); }
});
$('repeat-slider').addEventListener('input', () => {
  state.repeats = Number($('repeat-slider').value); $('repeat-count').value = state.repeats; $('repeat-error').hidden = true; $('repeat-count').removeAttribute('aria-invalid'); $('repeat-count').removeAttribute('aria-describedby'); updateSelection();
});
$('show-hydrogens').addEventListener('change', () => viewer?.setOptions({ showHydrogens: $('show-hydrogens').checked }));
$('render-style').addEventListener('change', () => viewer?.setOptions({ style: $('render-style').value }));
$('rotate-toggle').addEventListener('click', () => {
  const rotate = $('rotate-toggle').getAttribute('aria-pressed') !== 'true';
  $('rotate-toggle').setAttribute('aria-pressed', rotate); viewer?.setOptions({ autoRotate: rotate });
});
$('reset-view').addEventListener('click', () => viewer?.resetView());
$('fullscreen').addEventListener('click', async () => {
  const card = $('viewer-card');
  try {
    if (document.fullscreenElement) await document.exitFullscreen();
    else if (card.requestFullscreen && !card.classList.contains('expanded')) await card.requestFullscreen();
    else card.classList.toggle('expanded');
  } catch { card.classList.toggle('expanded'); }
  viewer?.resize();
});
document.addEventListener('keydown', event => { if (event.key === 'Escape') { $('viewer-card').classList.remove('expanded'); viewer?.resize(); } });
$('guide-open').addEventListener('click', () => $('guide-dialog').showModal());
for (const id of ['guide-close', 'guide-done']) $(id).addEventListener('click', () => $('guide-dialog').close());
$('guide-dialog').addEventListener('click', event => { if (event.target === $('guide-dialog')) { const bounds = $('guide-dialog').getBoundingClientRect(); if (event.clientX < bounds.left || event.clientX > bounds.right || event.clientY < bounds.top || event.clientY > bounds.bottom) $('guide-dialog').close(); } });
$('download-pdb').addEventListener('click', () => exportSafely('pdb', 'chemical/x-pdb', toPDB));
$('download-sdf').addEventListener('click', () => exportSafely('sdf', 'chemical/x-mdl-sdfile', toSDF));
$('download-recipe').addEventListener('click', () => exportSafely('json', 'application/json', model => JSON.stringify({
  application: 'Polymer Programming', version: '1.2.0', settings: { headgroup: groups[state.builtConfig.group].name, carbons: state.builtConfig.carbons, repeats: state.builtConfig.repeats, endCaps: 'CH3' },
  template: `${groups[state.builtConfig.group].name} · C${state.builtConfig.carbons}`,
  stats: model.stats, metadata: model.metadata,
  source: state.templates[`${state.builtConfig.group}_${state.builtConfig.carbons}`].source,
  assembly: 'C4 of repeat i is bonded to C1 of repeat i+1. Both outer cap hydrogens are replaced by new neutral CH3 groups, bonded to the first C1 and last C4. The original C5 pendant methyl is retained.',
  geometry: 'Initial assembled coordinates; no whole-polymer energy minimization. No counterions or solvent. Tacticity is not controlled.',
  createdAt: new Date().toISOString(),
}, null, 2) + '\n'));

fillLibrary();
$('repeat-count').max = MAX_REPEATS;
$('repeat-slider').max = MAX_REPEATS;
$('repeat-error').textContent = `Enter a whole number between 1 and ${MAX_REPEATS}.`;
$('repeat-count').value = state.repeats; $('repeat-slider').value = state.repeats;
updateSelection();

async function initialize() {
  try {
    const response = await fetch(new URL('../data/templates.json', import.meta.url));
    if (!response.ok) throw new Error(`Templates could not be loaded (HTTP ${response.status}).`);
    state.templates = await response.json();
    if (Object.keys(state.templates).length !== 16) throw new Error('The complete 16-template library is required.');
    try { viewer = new MoleculeViewer($('molecule-canvas')); }
    catch (error) { notify(`3D preview unavailable: ${error.message}. Structure downloads are still available.`); }
    try {
      worker = new Worker(new URL('./builder-worker.js', import.meta.url), { type: 'module' });
      worker.onmessage = ({ data }) => {
        const task = pending.get(data.id); if (!task) return;
        pending.delete(data.id); data.error ? task.reject(new Error(data.error)) : task.resolve(data.model);
      };
      worker.onerror = () => {
        worker?.terminate(); worker = null;
        for (const { reject } of pending.values()) reject(new Error('The background builder stopped. Please build again.'));
        pending.clear();
      };
      worker.postMessage({ type: 'initialize', templates: state.templates });
    } catch { worker = null; }
    $('build-button').disabled = false;
    await generate();
  } catch (error) {
    $('stage-loading').hidden = true;
    $('stage-error').hidden = false;
    $('stage-error').textContent = `${error.message} Start the local server with Start-Website.cmd, then open http://127.0.0.1:4173.`;
    $('build-button-label').textContent = 'Templates unavailable';
    setStatus('Unable to load', true);
  }
}
initialize();
