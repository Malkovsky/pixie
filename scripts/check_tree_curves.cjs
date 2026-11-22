// With jsdom 24 on NODE_PATH: node scripts/check_tree_curves.cjs [built-page.html]
const {JSDOM, ResourceLoader} = require('jsdom');
const assert = require('node:assert/strict');
const path = require('node:path');
const {pathToFileURL} = require('node:url');

async function check() {
  const assets = path.resolve(__dirname, '../src/docs/images/benchmarks/trees/preview');
  const resources = new class extends ResourceLoader {
    fetch(url, options) { return url.startsWith('file:') ? super.fetch(url, options) : null; }
  }();
  const options = {resources, runScripts:'dangerously'};
  const data = pathToFileURL(path.resolve(__dirname, '../src/docs/benchmarks/data')).href;
  const dom = process.argv[2] ? await JSDOM.fromFile(path.resolve(process.argv[2]), options) :
    new JSDOM(`<link rel="stylesheet" href="tree-curves.css"><script src="d3-7.9.0.min.js"></script>
      <script src="plot-0.6.17.min.js"></script><script src="${data}/tree-matched.js"></script>
      <script src="${data}/tree-history.js"></script><body><script>
      for (const key of [...Object.keys(window.pixieTreeCharts), 'matched-inspected_dfs', 'matched-inspected_root_path', 'matched-tree_path']) {
        const root = document.createElement('div'); root.className = 'pixie-tree-chart';
        root.dataset.chart = key; document.body.append(root);
      }</script><script src="tree-curves.js"></script>`, {...options, url:pathToFileURL(assets + '/fixture.html').href});
  const {document, Event, KeyboardEvent, MouseEvent, WheelEvent} = dom.window;
  if (document.readyState !== 'complete') await new Promise(ok => dom.window.addEventListener('load', ok));
  try {
    const roots = [...document.querySelectorAll('.pixie-tree-chart')];
    const specs = dom.window.pixieTreeCharts;
    assert.equal(roots.length, 9);
    assert.equal(document.querySelector('[data-chart="matched-tree_path_sibling"]'), null);
    if (process.argv[2]) {
      const history = document.querySelector('details#tree-history');
      assert.equal(history.open, false);
      assert.equal(history.querySelectorAll('.pixie-tree-chart').length, 4);
      for (const id of ['history-build-prepared','history-memory']) {
        assert.equal(document.querySelector(`[data-chart="${id}"]`).closest('details'), null);
      }
      history.querySelector('summary').click();
      assert.equal(history.open, true);
    }
    assert.deepEqual(roots.map(r => r.dataset.chart).sort(), Object.keys(specs).sort());
    for (const [i, root] of roots.entries()) {
      const spec = specs[root.dataset.chart], el = s => root.querySelector(s);
      const other = roots[(i + 1) % roots.length], otherSvg = other.querySelector('svg');
      const svg = () => el('.chart-plot svg'), chart = el('.chart-plot');
      const domains = () => ['x','y'].map(k => svg().scale(k).domain.join()).join(';');
      const dots = () => [...root.querySelectorAll('circle')];
      const full = domains();
      assert.equal(dots().length, spec.points.length);
      assert.equal(svg().getAttribute('aria-label'), spec.title);
      assert.equal(svg().scale('y').type, spec.yScale);
      assert.ok([...svg().querySelectorAll('text')].some(t => t.textContent === spec.yLabel && t.closest('[transform*="rotate(-90"]')));
      const dot = dots()[0], tip = el('.point-tooltip');
      dot.dispatchEvent(new MouseEvent('pointermove', {bubbles:true, clientX:80, clientY:120}));
      assert.equal(tip.hidden, false);
      assert.equal(tip.textContent, dot.getAttribute('aria-label'));
      assert.ok(tip.textContent.includes(spec.unit));
      assert.equal(tip.textContent.includes('nodes/path'), spec.points[0].length != null);
      chart.dispatchEvent(new Event('pointerleave'));
      assert.equal(tip.hidden, true);
      const bounds = () => { svg().getBoundingClientRect = () => ({left:0, top:0, width:1000, height:450}); };
      const wheel = deltaY => {
        bounds(); chart.dispatchEvent(new WheelEvent('wheel', {deltaY, clientX:400, clientY:200, cancelable:true}));
      };
      wheel(-400);
      assert.notEqual(domains(), full);
      wheel(400);
      assert.equal(domains(), full);
      wheel(-400);
      const beforePan = domains();
      bounds(); chart.setPointerCapture = () => {};
      svg().dispatchEvent(new MouseEvent('pointerdown', {bubbles:true,clientX:400,clientY:200}));
      chart.dispatchEvent(new MouseEvent('pointermove', {clientX:480,clientY:240}));
      assert.notEqual(domains(), beforePan);
      chart.dispatchEvent(new Event('lostpointercapture'));
      el('.reset').click();
      assert.equal(domains(), full);
      assert.equal(el('.zoom-in, .zoom-out'), null);
      chart.dispatchEvent(new KeyboardEvent('keydown', {key:'+'}));
      assert.notEqual(domains(), full);
      chart.dispatchEvent(new KeyboardEvent('keydown', {key:'-'}));
      assert.equal(domains(), full);
      chart.dispatchEvent(new KeyboardEvent('keydown', {key:'ArrowRight'}));
      assert.notEqual(domains(), full);
      el('.reset').click();
      el('summary').click();
      await new Promise(ok => setTimeout(ok, 5));
      const search = el('[type=search]'), box = el('[type=checkbox]');
      assert.equal(document.activeElement, search);
      search.value = box.getAttribute('aria-label').split(/\s+/).reverse().join(' ');
      search.dispatchEvent(new Event('input'));
      assert.equal(box.closest('label').hidden, false);
      assert.ok(el('label[hidden]'));
      el('.none').click();
      assert.equal(chart.textContent, 'Select a curve to display.');
      box.click();
      assert.match(el('summary').textContent, /^1 of /);
      assert.ok(dots().length > 0 && dots().length < spec.points.length);
      assert.equal(domains(), full);
      el('.all').click();
      assert.equal(dots().length, spec.points.length);
      search.dispatchEvent(new KeyboardEvent('keydown', {key:'Escape', bubbles:true}));
      assert.equal(el('details').open, false);
      assert.equal(document.activeElement, el('summary'));
      el('summary').click(); el('.done').click();
      assert.equal(document.activeElement, el('summary'));
      el('summary').click(); document.body.click();
      assert.equal(el('details').open, false);
      assert.equal(other.querySelector('svg'), otherSvg);
      console.log(`${root.dataset.chart}: ${spec.points.length} points, independent controls passed`);
    }
  } finally { dom.window.close(); }
}
check().catch(error => {console.error(error); process.exitCode = 1;});
