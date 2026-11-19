(() => {
  const el = id => document.getElementById(id);
  const data = JSON.parse(el('points').textContent);
  const dfs = el('points').dataset.operation === 'inspected_dfs';
  const unit = dfs ? 'ns/node' : 'ns/path';
  const sizes = [...new Set(data.map(d => d.n))].sort((a, b) => a - b);
  const names = {bp:'BP', dfuds:'DFUDS', louds:'LOUDS', forest_cover:'ForestCover',
    adj_list:'Adjacency list', pointer:'Pointers'};
  const colors = ['#0072B2','#D55E00','#009E73','#AC548D','#777','#222'];
  const label = d => ['adj_list','pointer'].includes(d.impl) ? 'Reference' :
    d.b === 0 ? 'Baseline' : `${d.b / 1024} KiB`;
  const key = d => `${d.impl}:${d.b}`;
  const curves = d3.group(data, key);
  const controls = el('controls'), picker = el('picker'), trigger = el('trigger');
  const search = el('search'), chart = el('chart');
  const reset = document.createElement('button');
  reset.textContent = 'Reset';
  reset.title = 'Restore the original axes';
  const toolbar = document.createElement('div');
  toolbar.className = 'chart-toolbar';
  picker.before(toolbar);
  toolbar.append(picker, reset);
  let zoom = 1, domains, drag;
  reset.onclick = () => { zoom = 1; domains = undefined; draw(); };
  chart.addEventListener('wheel', e => {
    if (drag) return;
    const svg = chart.querySelector('svg');
    if (!svg || !e.deltaY) return;
    const rect = svg.getBoundingClientRect();
    const x = (e.clientX - rect.left) * 1000 / rect.width;
    const y = (e.clientY - rect.top) * 450 / rect.height;
    const axes = [svg.scale('x'), svg.scale('y')];
    if (![x,y].every((p,i) => p >= Math.min(...axes[i].range) && p <= Math.max(...axes[i].range))) return;
    const delta = e.deltaY * (e.deltaMode === 1 ? 16 : e.deltaMode === 2 ? rect.height : 1);
    const next = Math.max(1, Math.min(32, zoom * 2 ** (-delta / 400)));
    if (next === zoom) return;
    e.preventDefault();
    domains = next === 1 ? undefined : axes.map((s,i) => {
      const p = i ? y : x;
      return s.range.map(edge => s.invert(p + (edge - p) * zoom / next));
    });
    zoom = next;
    draw();
  }, {passive:false});
  const tip = document.createElement('div');
  tip.className = 'point-tooltip';
  tip.hidden = true;
  document.body.append(tip);
  chart.addEventListener('pointerdown', e => {
    const svg = e.target.closest('svg');
    if (!svg || drag || e.button !== 0) return;
    const rect = svg.getBoundingClientRect();
    const axes = [svg.scale('x'), svg.scale('y')];
    const scale = [1000 / rect.width, 450 / rect.height];
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
  for (const [index, id] of Object.keys(names).entries()) {
    const points = data.filter(d => d.impl === id);
    if (!points.length) continue;
    const group = document.createElement('div'), heading = document.createElement('h2');
    heading.textContent = names[id];
    heading.style.color = colors[index];
    group.append(heading);
    for (const d of d3.group(points, key).values()) {
      const row = document.createElement('label'), box = document.createElement('input');
      box.type = 'checkbox';
      box.checked = true;
      box.value = key(d[0]);
      box.setAttribute('aria-label', `${names[id]} ${label(d[0])}`);
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
      marks.push(Plot.line(sizes, {clip:!!domains, x:n => n, y:n => bySize.get(n), stroke:() => first.impl,
        strokeDasharray:({32768:'7,4',524288:'2,3'})[first.b]}));
    }
    marks.push(Plot.dot(data.filter(d => on.has(key(d))), {clip:!!domains, x:'n', y:'t', stroke:'impl', r:4,
      pointerEvents:'all', ariaLabel:d => `${names[d.impl]} · ${label(d)}\n${d.n.toLocaleString()} nodes\n${d.t.toFixed(1)} ${unit}` +
        (dfs ? '' : ` · ${d.length.toFixed(1)} nodes/path`)}));
    chart.replaceChildren(Plot.plot({width:1000, height:450, marginLeft:74, marginRight:40,
      ariaLabel:document.querySelector('h1').textContent,
      color:{domain:Object.keys(names), range:colors},
      x:{type:'log', domain:domains?.[0] ?? d3.extent(sizes), ticks:domains ? 6 : sizes, tickFormat:domains ? undefined : n => '2^' + Math.log2(n), label:'Tree size (nodes)'},
      y:{type:dfs ? 'linear' : 'log', domain:domains?.[1] ?? (dfs ? [0,d3.max(data,d => d.t)] : d3.extent(data,d => d.t)),
        nice:!domains, grid:true, label:`CPU time (${unit})`, labelAnchor:'center', labelOffset:62, labelArrow:false}, marks}));
  }
  controls.onchange = draw;
  for (const [id, checked] of [['all',true], ['none',false]]) {
    el(id).onclick = () => { controls.querySelectorAll('input').forEach(c => c.checked = checked); draw(); };
  }
  const close = () => { picker.open = false; trigger.focus(); };
  el('done').onclick = close;
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
  document.addEventListener('keydown', e => { if (e.key === 'Escape' && picker.open) { e.preventDefault(); close(); } });
  document.addEventListener('click', e => { if (picker.open && !picker.contains(e.target)) picker.open = false; });
  draw();
})();
