(() => {
function mount(root, {points:data, title, unit, yLabel, yScale}) {
  root.innerHTML = `<div class="chart-toolbar"><details><summary></summary><div class="chart-panel">
    <input type="search" aria-label="Search curves" placeholder="Search curves">
    <fieldset><legend>Visible curves</legend></fieldset>
    <button type="button" class="all">All</button> <button type="button" class="none">None</button>
    <button type="button" class="done">Done</button></div></details>
    <button type="button" class="reset" title="Restore the original axes">Reset</button></div>
    <div class="chart-plot" tabindex="0" role="group" aria-label="Plot: use arrow keys to pan and plus or minus to zoom"></div><div class="point-tooltip" hidden></div>`;
  root.setAttribute('role', 'group');
  root.setAttribute('aria-label', title);
  const el = selector => root.querySelector(selector);
  const sizes = [...new Set(data.map(d => d.n))].sort((a, b) => a - b);
  const names = {bp:'BP', dfuds:'DFUDS', louds:'LOUDS', forest_cover:'ForestCover',
    adj_list:'Adjacency list', pointer:'Pointers', bfs_csr:'BFS CSR'};
  const colors = ['#0072B2','#D55E00','#009E73','#AC548D','#777','#222','#827100'];
  const label = d => d.curve ?? (['adj_list','pointer','bfs_csr'].includes(d.impl) ? 'Reference' :
    d.b === 0 ? 'Baseline' : `${d.b / 1024} KiB`);
  const key = d => `${d.impl}:${d.curve ?? d.b}`;
  const family = d => d.family ?? names[d.impl] ?? d.impl;
  const families = [...new Set([...Object.keys(names), ...data.map(d => d.impl)])];
  const color = d3.scaleOrdinal(families, [...colors, ...d3.schemeTableau10]);
  const curves = d3.group(data, key);
  const controls = el('fieldset'), picker = el('details'), trigger = el('summary');
  const search = el('[type=search]'), chart = el('.chart-plot'), reset = el('.reset');
  const tip = el('.point-tooltip');
  let width = chart.clientWidth || 1000;
  const height = () => Math.min(450, Math.max(290, width * 0.55));
  let zoom = 1, domains, drag;
  reset.onclick = () => { zoom = 1; domains = undefined; draw(); };
  function zoomAt(axes, pivot, factor) {
    const next = Math.max(1, Math.min(32, zoom * factor));
    if (next === zoom) return false;
    domains = next === 1 ? undefined : axes.map((s,i) => s.range.map(edge =>
      s.invert(pivot[i] + (edge - pivot[i]) * zoom / next)));
    zoom = next;
    draw();
    return true;
  }
  chart.onkeydown = e => {
    if (e.ctrlKey || e.metaKey || e.altKey) return;
    const step = {ArrowLeft:[-.1,0], ArrowRight:[.1,0], ArrowDown:[0,.1], ArrowUp:[0,-.1]}[e.key];
    const svg = el('svg');
    if (svg && (e.key === '+' || e.key === '-')) {
      e.preventDefault();
      const axes = [svg.scale('x'), svg.scale('y')];
      zoomAt(axes, axes.map(s => d3.mean(s.range)), e.key === '+' ? 2 : .5);
      return;
    }
    if (!step || !svg) return;
    e.preventDefault();
    domains = ['x','y'].map((k,i) => {
      const s = svg.scale(k);
      return s.range.map(edge => s.invert(edge + step[i] * Math.abs(s.range[1] - s.range[0])));
    });
    draw();
  };
  chart.addEventListener('wheel', e => {
    if (drag) return;
    const svg = chart.querySelector('svg');
    if (!svg || !e.deltaY) return;
    const rect = svg.getBoundingClientRect();
    const x = (e.clientX - rect.left) * width / rect.width;
    const y = (e.clientY - rect.top) * height() / rect.height;
    const axes = [svg.scale('x'), svg.scale('y')];
    if (![x,y].every((p,i) => p >= Math.min(...axes[i].range) && p <= Math.max(...axes[i].range))) return;
    const delta = e.deltaY * (e.deltaMode === 1 ? 16 : e.deltaMode === 2 ? rect.height : 1);
    if (zoomAt(axes, [x,y], 2 ** (-delta / 400))) e.preventDefault();
  }, {passive:false});
  chart.addEventListener('pointerdown', e => {
    const svg = e.target.closest('svg');
    if (!svg || drag || e.button !== 0) return;
    const rect = svg.getBoundingClientRect();
    const axes = [svg.scale('x'), svg.scale('y')];
    const scale = [width / rect.width, height() / rect.height];
    const pos = [(e.clientX - rect.left) * scale[0], (e.clientY - rect.top) * scale[1]];
    if (!pos.every((p,i) => p >= Math.min(...axes[i].range) && p <= Math.max(...axes[i].range))) return;
    drag = {id:e.pointerId, x:e.clientX, y:e.clientY, axes, scale};
    chart.setPointerCapture(e.pointerId);
    chart.classList.add('dragging');
    tip.hidden = true;
    e.preventDefault();
  });
  chart.addEventListener('lostpointercapture', () => {
    drag = undefined;
    chart.classList.remove('dragging');
  });
  chart.addEventListener('pointermove', e => {
    if (drag) {
      if (e.pointerId !== drag.id) return;
      const delta = [e.clientX - drag.x, e.clientY - drag.y];
      domains = drag.axes.map((s,i) => s.range.map(edge => s.invert(edge - delta[i] * drag.scale[i])));
      draw();
      return;
    }
    const text = e.target.closest('circle')?.getAttribute('aria-label');
    tip.hidden = !text;
    if (!text) return;
    if (tip.textContent !== text) tip.textContent = text;
    tip.style.left = `${Math.max(0, Math.min(e.clientX + 12, innerWidth - tip.offsetWidth - 4))}px`;
    tip.style.top = `${e.clientY + 12 + tip.offsetHeight <= innerHeight ? e.clientY + 12 : Math.max(0, e.clientY - tip.offsetHeight - 8)}px`;
  });
  chart.addEventListener('pointerleave', () => { tip.hidden = true; });
  for (const id of families) {
    const points = data.filter(d => d.impl === id);
    if (!points.length) continue;
    const group = document.createElement('div'), heading = document.createElement('strong');
    heading.textContent = family(points[0]);
    heading.style.color = color(id);
    group.append(heading);
    for (const d of d3.group(points, key).values()) {
      const row = document.createElement('label'), box = document.createElement('input');
      box.type = 'checkbox';
      box.checked = true;
      box.value = key(d[0]);
      box.setAttribute('aria-label', `${family(d[0])} ${label(d[0])}`);
      row.append(box, ` ${label(d[0])}`);
      group.append(row);
    }
    controls.append(group);
  }
  function draw() {
    reset.disabled = !domains;
    tip.hidden = true;
    const on = new Set([...controls.querySelectorAll(':checked')].map(c => c.value));
    trigger.textContent = `${on.size} of ${curves.size} curves`;
    if (!on.size) {
      chart.textContent = 'Select a curve to display.';
      return;
    }
    const marks = [];
    for (const k of on) {
      const points = curves.get(k), first = points[0];
      const bySize = new Map(points.map(d => [d.n, d.t]));
      // An undefined y breaks the line instead of bridging a missing curve at a measured size.
      const dash = first.curve === undefined ? ({32768:'7,4',524288:'2,3'})[first.b] :
        ({BFS:'7,4','32 KiB':'2,3'})[first.curve];
      marks.push(Plot.line(sizes, {clip:!!domains, x:n => n, y:n => bySize.get(n), stroke:color(first.impl), strokeDasharray:dash}));
    }
    marks.push(Plot.dot(data.filter(d => on.has(key(d))), {clip:!!domains, x:'n', y:'t', stroke:'impl', r:4,
      pointerEvents:'all', ariaLabel:d => `${family(d)} · ${label(d)}\n${d.n.toLocaleString()} nodes\n${d3.format('.5~g')(d.t)} ${unit}` +
        (d.length == null ? '' : ` · ${d.length.toFixed(1)} nodes/path`)}));
    const stride = Math.max(1, Math.ceil(sizes.length / Math.max(2, (width - 100) / 55)));
    chart.replaceChildren(Plot.plot({width, height:height(), marginLeft:74, marginRight:30,
      ariaLabel:title, color:{domain:families, range:families.map(color)},
      x:{type:'log', domain:domains?.[0] ?? d3.extent(sizes), ticks:domains ? 6 : sizes.filter((_,i) => i % stride === 0),
        tickFormat:domains ? undefined : n => Number.isInteger(Math.log2(n)) ? '2^' + Math.log2(n) : d3.format('~s')(n), label:'Tree size (nodes)'},
      y:{type:yScale, domain:domains?.[1] ?? (yScale === 'linear' ? [0,d3.max(data,d => d.t)] : d3.extent(data,d => d.t)),
        nice:!domains, grid:true, label:yLabel, labelAnchor:'center', labelOffset:62, labelArrow:false}, marks}));
  }
  controls.onchange = draw;
  for (const [id, checked] of [['all',true], ['none',false]]) {
    el('.' + id).onclick = () => { controls.querySelectorAll('input').forEach(c => c.checked = checked); draw(); };
  }
  const close = () => { picker.open = false; trigger.focus(); };
  el('.done').onclick = close;
  search.oninput = () => {
    const terms = search.value.trim().toLowerCase().split(/\s+/);
    for (const group of controls.querySelectorAll('div')) {
      for (const row of group.querySelectorAll('label')) {
        const text = row.firstChild.getAttribute('aria-label').toLowerCase();
        row.hidden = !terms.every(term => text.includes(term));
      }
      group.hidden = !group.querySelector('label:not([hidden])');
    }
  };
  picker.ontoggle = () => { if (picker.open) search.focus(); };
  picker.addEventListener('keydown', e => { if (e.key === 'Escape' && picker.open) { e.preventDefault(); close(); } });
  document.addEventListener('click', e => { if (picker.open && !picker.contains(e.target)) picker.open = false; });
  draw();
  if (window.ResizeObserver) new ResizeObserver(() => {
    if (chart.clientWidth && width !== chart.clientWidth) { width = chart.clientWidth; draw(); }
  }).observe(chart);
}
function init() {
  const charts = window.pixieTreeCharts ??= {};
  // Matched data retains the three source-family medians alongside each plotted mean.
  const titles = {inspected_dfs:'Matched inspected DFS', inspected_root_path:'Matched inspected root paths', tree_path:'Matched random tree paths'};
  for (const [op, points] of Object.entries(window.PIXIE_TREE_MATCHED?.charts ?? {})) {
    if (!Object.hasOwn(titles, op)) continue;
    const dfs = op === 'inspected_dfs', unit = dfs ? 'ns/node' : 'ns/path';
    charts['matched-' + op] = {points, title:titles[op], unit, yLabel:`CPU time (${unit})`, yScale:dfs ? 'linear' : 'log'};
  }
  for (const root of document.querySelectorAll('.pixie-tree-chart[data-chart]')) {
    const spec = charts[root.dataset.chart];
    if (spec) mount(root, spec);
    else root.textContent = `Chart data unavailable: ${root.dataset.chart}`;
  }
}
if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init, {once:true});
else init();
})();
