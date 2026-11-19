// Run with jsdom 24 on NODE_PATH. Loads the published pages via file://.
const {JSDOM} = require('jsdom');
const assert = require('node:assert/strict');
const path = require('node:path');

async function check(stem) {
  const file = path.resolve(__dirname, '../src/docs/images/benchmarks/trees/preview', stem + '_observable.html');
  const dom = await JSDOM.fromFile(file, {resources:'usable', runScripts:'dangerously'});
  const {document, Event, KeyboardEvent} = dom.window;
  await new Promise(ok => dom.window.addEventListener('load', ok));
  const el = id => document.getElementById(id);
  const data = JSON.parse(el('points').textContent);
  const dots = () => document.querySelectorAll('#chart circle').length;
  assert.equal(dots(), data.length);
  const dot = document.querySelector('#chart circle'), tip = document.querySelector('.point-tooltip');
  assert.ok(!dot.querySelector('title'));
  dot.dispatchEvent(new dom.window.MouseEvent('pointermove', {bubbles:true, clientX:80, clientY:120}));
  assert.equal(tip.hidden, false);
  assert.match(tip.textContent, /nodes\n[\d.]+ ns\/(node|path)/);
  el('chart').dispatchEvent(new Event('pointerleave'));
  assert.equal(tip.hidden, true);
  assert.equal(document.querySelectorAll('#chart svg').length, 1);
  const svg = () => el('chart').querySelector('svg');
  const wheel = deltaY => {
    svg().getBoundingClientRect = () => ({left:0, top:0, width:1000, height:450});
    el('chart').dispatchEvent(new dom.window.WheelEvent('wheel', {clientX:400, clientY:200, deltaY, cancelable:true}));
  };
  const full = svg().scale('x').domain.join();
  wheel(-400);
  assert.notEqual(svg().scale('x').domain.join(), full);
  wheel(400);
  assert.equal(svg().scale('x').domain.join(), full);
  wheel(-400);
  const reset = el('picker').nextElementSibling;
  assert.equal(reset.textContent, 'Reset');
  const beforePan = svg().scale('x').domain.join();
  svg().getBoundingClientRect = () => ({left:0, top:0, width:1000, height:450});
  el('chart').setPointerCapture = () => {};
  svg().dispatchEvent(new dom.window.MouseEvent('pointerdown',{bubbles:true,clientX:400,clientY:200}));
  el('chart').dispatchEvent(new dom.window.MouseEvent('pointermove',{clientX:480,clientY:240}));
  assert.notEqual(svg().scale('x').domain.join(), beforePan);
  el('chart').dispatchEvent(new Event('lostpointercapture'));
  reset.click();
  assert.equal(svg().scale('x').domain.join(), full);
  const label = [...document.querySelectorAll('svg text')].find(t => t.textContent ===
    `CPU time (ns/${stem === 'inspected_dfs' ? 'node' : 'path'})`);
  assert.ok(label.closest('[transform*="rotate(-90"]').tagName);
  el('trigger').click();
  await new Promise(ok => setTimeout(ok, 5));
  assert.equal(document.activeElement, el('search'));
  el('search').value = '32 BP';
  el('search').dispatchEvent(new Event('input'));
  const matches = document.querySelectorAll('label:not([hidden]) input');
  assert.equal(matches.length, 1);
  assert.equal(matches[0].value, 'bp:32768');
  el('none').click();
  assert.equal(el('chart').textContent, 'Select a curve to display.');
  matches[0].click();
  assert.equal(dots(), data.filter(p => p.impl === 'bp' && p.b === 32768).length);
  assert.match(el('trigger').textContent, /^1 of /);
  el('all').click();
  assert.equal(dots(), data.length);
  document.dispatchEvent(new KeyboardEvent('keydown', {key:'Escape'}));
  assert.equal(el('picker').open, false);
  assert.equal(document.activeElement, el('trigger'));
  el('trigger').click();
  el('done').click();
  assert.equal(document.activeElement, el('trigger'));
  el('trigger').click();
  document.body.click();
  assert.equal(el('picker').open, false);
  dom.window.close();
  console.log(`${stem}: ${data.length} dots, selector checks passed`);
}
(async () => {
  for (const stem of ['inspected_dfs','inspected_root_path','tree_path']) await check(stem);
})().catch(error => {console.error(error); process.exitCode = 1;});
