import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import vm from 'node:vm';

test('opening the HTML directly shows recovery instructions instead of an endless loader', async () => {
  const html = await readFile(new URL('../index.html', import.meta.url), 'utf8');
  const script = html.match(/<script>([\s\S]*?)<\/script>/)[1];
  const elements = Object.fromEntries(['stage-loading', 'stage-error', 'model-status', 'build-button', 'build-button-label'].map(id => [id, {}]));
  let onReady;
  vm.runInNewContext(script, {
    location: { protocol: 'file:' },
    document: {
      addEventListener: (event, callback) => { assert.equal(event, 'DOMContentLoaded'); onReady = callback; },
      getElementById: id => elements[id],
    },
  });
  await onReady();
  assert.equal(elements['stage-loading'].hidden, true);
  assert.equal(elements['stage-error'].hidden, false);
  assert.match(elements['stage-error'].textContent, /Start-Website.cmd/);
  assert.match(elements['stage-error'].textContent, /http:\/\/127.0.0.1:4173/);
  assert.equal(elements['build-button'].disabled, true);
  assert.equal(elements['model-status'].textContent, 'Unable to load');
});
