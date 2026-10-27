import { buildPolymer } from './polymer.js';
let templates;
self.onmessage = ({ data }) => {
  if (data.type === 'initialize') { templates = data.templates; return; }
  if (data.type === 'build') {
    try { self.postMessage({ id: data.id, model: buildPolymer(templates, data.config) }); }
    catch (error) { self.postMessage({ id: data.id, error: error.message }); }
  }
};
