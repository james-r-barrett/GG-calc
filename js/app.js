let PARTS = [], ACCEPT = [], LINKERS = [];

async function loadData() {
  const [parts, acceptors, linkers] = await Promise.all([
    fetch('data/parts.json').then(r => r.json()),
    fetch('data/acceptors.json').then(r => r.json()),
    fetch('data/linkers.json').then(r => r.json()),
  ]);
  PARTS = parts; ACCEPT = acceptors; LINKERS = linkers;
  init();
}

function init() {
  /* ============ Embedded part libraries ============
     PARTS   — L0 parts:      {n:name, s:seq, t:type, s5:5' syntax, s3:3' syntax}
     ACCEPT  — acceptor/resistance vectors: {n, s, enz:enzyme, lvl:level acceptor, pos:position at L2, mark:selectable marker}
     LINKERS — end linkers / dummy fragments: {n, s, t:type, pos:positions occupied, lvl:level}
  */
  
  /* ============ Category colours ============ */
  const CAT_COLORS = {
    'Promoter': '#2F7D4F',
    "Promoter + 5'UTR": '#2F7D4F',
    "5'UTR": '#4E9B6E',
    "3'UTR + Terminator": '#B24545',
    'CTP': '#5C7CBA',
    'Signal Peptide': '#5C7CBA',
    'Epitope': '#9B6EC4',
    'Fluorescent Tag': '#D4568C',
    'Reporter': '#B85C38',
    'Resistance Marker': '#C08A2E',
    'Intron': '#7C8CA6',
    'MicroRNA': '#8A63B0',
    'Multi-STOP': '#A6564B',
    'Leaky Stop': '#C97B6E',
    'Ribosome Reinitiation': '#4E7A8C',
    'End Linker': '#6B7280',
    'Dummy Sequence': '#9CA3AF',
    'Acceptor Vector': '#14213D',
    'Custom': '#8A8F98',
  };
  function catColor(cat){ return CAT_COLORS[cat] || '#8A8F98'; }
  function catDot(cat){ return `<span class="cat-dot" style="background:${catColor(cat)}"></span>`; }
  function catPill(cat){ return `<span class="cat-pill">${catDot(cat)}${escapeHtml(cat||'Custom')}</span>`; }
  
  function escapeHtml(str){
    return String(str==null?'':str).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  }
  
  /* ============ DNA maths (mirrors the workbook's "Vol Req - Sequence" method) ============ */
  function baseCounts(seq){
    seq = (seq||'').toUpperCase();
    let a=0,c=0,g=0,t=0;
    for (let i=0;i<seq.length;i++){
      const ch = seq[i];
      if (ch==='A') a++; else if (ch==='C') c++; else if (ch==='G') g++; else if (ch==='T') t++;
    }
    return {a,c,g,t};
  }
  function molWeight(seq){
    const {a,c,g,t} = baseCounts(seq);
    return a*331.2 + c*307.2 + g*347.2 + t*322.2;
  }
  // targetFmol * MW(g/mol) * 1e-15(fmol->mol) / (conc(ng/uL) * 1e-9(ng->g)) = volume in uL
  function volReqUl(targetFmol, seq, concNgUl){
    if (!seq || !concNgUl || concNgUl<=0 || targetFmol==null) return null;
    return (targetFmol * molWeight(seq) * 1e-6) / concNgUl;
  }
  // Length-only method: average dsDNA MW per bp (617.96/bp + 36.04), halved as in the workbook's "Vol Req - Length" column.
  function molWeightFromLength(lenBp){
    return lenBp*617.96 + 36.04;
  }
  function volReqLengthUl(targetFmol, lenBp, concNgUl){
    if (lenBp==null || !concNgUl || concNgUl<=0 || targetFmol==null) return null;
    return (targetFmol * molWeightFromLength(lenBp) * 1e-6) / (2*concNgUl);
  }
  function cleanSeq(raw){
    return (raw||'').toUpperCase().replace(/[^ACGT]/g, '');
  }
  function fmt(n, dp){
    if (n==null || isNaN(n)) return '&ndash;';
    return n.toFixed(dp==null?2:dp);
  }
  
  /* ============ Combined searchable insert list ============ */
  const INSERTABLE = [];
  PARTS.forEach(p => INSERTABLE.push({ kind:'part', raw:p, cat:(p.t||'').trim() }));
  LINKERS.forEach(l => INSERTABLE.push({ kind:'linker', raw:l, cat:(l.t||'').trim() }));
  
  /* ============ App state ============ */
  let uidCounter = 1;
  const state = {
    accFmol: 25,
    insFmol: 50,
    totalVol: 10,
    calcMode: 'sequence', // 'sequence' | 'length'
    mm: {
      bufferType: 'neb',           // 'neb' | 'homemade'
      bufferVol: null, bufferTouched: false,   // NEB 10x buffer override (default totalVol/10)
      bufferAVol: 1, bufferBVol: 1,            // home-made 2-part buffer, freely editable
      ligaseVol: null, ligaseTouched: false,   // default totalVol/40
      enzymeVol: null, enzymeTouched: false,   // default totalVol/20
    },
    acceptor: { mode:'library', part:null, conc:'', customName:'', customSeq:'', customLen:'' },
    inserts: [ mkInsert(), mkInsert() ],
  };
  function mkInsert(mode){ return { id:'ins'+(uidCounter++), mode: mode||'library', part:null, cat:null, conc:'', customName:'', customSeq:'', customLen:'' }; }
  
  function rowSeq(row){
    if (row.mode==='custom') return cleanSeq(row.customSeq);
    return row.part ? row.part.raw.s : null;
  }
  function rowName(row, fallback){
    if (row.mode==='custom') return (row.customName||'').trim() || fallback;
    return row.part ? row.part.raw.n : fallback;
  }
  function rowCat(row){
    if (row.mode==='custom') return 'Custom';
    return row.part ? row.part.cat : null;
  }
  // Only meaningful for custom-mode rows in length-only calculation mode.
  function rowLenBp(row){
    const n = parseInt(row.customLen, 10);
    return Number.isFinite(n) && n > 0 ? n : null;
  }
  // Fragment length to *display* next to a row, independent of whether it can be used for a volume calc.
  function displayLen(row){
    if (row.mode === 'library') return row.part ? row.part.raw.s.length : null;
    if (state.calcMode === 'length') return rowLenBp(row);
    const seq = rowSeq(row);
    return seq ? seq.length : null;
  }
  
  /* ============ Master mix defaults ============ */
  function mmDefaults(){
    return { buffer: state.totalVol/10, ligase: state.totalVol/40, enzyme: state.totalVol/20 };
  }
  function mmCurrent(){
    const d = mmDefaults();
    return {
      buffer: state.mm.bufferTouched ? state.mm.bufferVol : d.buffer,
      bufferA: state.mm.bufferAVol,
      bufferB: state.mm.bufferBVol,
      ligase: state.mm.ligaseTouched ? state.mm.ligaseVol : d.ligase,
      enzyme: state.mm.enzymeTouched ? state.mm.enzymeVol : d.enzyme,
    };
  }
  function roundStr(n, dp){
    if (n==null || isNaN(n)) return '';
    const f = Math.pow(10, dp==null?3:dp);
    return (Math.round(n*f)/f).toString();
  }
  
  /* ============ Rendering: Acceptor + Insert rows ============ */
  function renderPartRow(container, row, opts){
    // opts: {label, items, onRemove}
    const len = displayLen(row);
    const lengthMode = state.calcMode === 'length';
  
    let html = `<div class="part-row-head">
        <span class="part-row-label">${opts.label}</span>
        <div class="row-actions">
          <div class="mode-toggle">
            <button type="button" class="mode-btn${row.mode==='library'?' active':''}" data-mode="library">Library</button>
            <button type="button" class="mode-btn${row.mode==='custom'?' active':''}" data-mode="custom">Custom</button>
          </div>
          ${opts.onRemove ? `<button type="button" class="remove-btn" data-action="remove" title="Remove">&times;</button>` : ''}
        </div>
      </div>`;
  
    if (row.mode === 'library'){
      const inputVal = row.part ? row.part.raw.n : '';
      html += `<div class="combo-wrap">
        <div class="combo-input-row">
          <input type="text" class="combo-input" placeholder="Search parts by name&hellip;" value="${escapeHtml(inputVal)}" autocomplete="off" role="combobox" aria-expanded="false" aria-autocomplete="list" aria-label="Search parts by name">
          <input type="text" inputmode="decimal" class="conc-input" data-field="conc" placeholder="ng/&micro;L" value="${escapeHtml(row.conc)}" aria-label="Stock concentration in nanograms per microlitre">
        </div>
        <div class="combo-panel" role="listbox" hidden></div>
      </div>`;
      if (row.part){
        const cat = row.part.cat;
        html += `<div class="selected-chip">${catPill(cat)}<span class="chip-len mono">${len} bp</span></div>`;
      }
    } else if (lengthMode){
      html += `<div class="custom-fields">
        <input type="text" data-field="customName" placeholder="Fragment name" value="${escapeHtml(row.customName)}">
        <div class="conc-row">
          <input type="number" min="1" step="1" data-field="customLen" placeholder="Length (bp)" value="${escapeHtml(row.customLen)}">
          <input type="text" inputmode="decimal" data-field="conc" placeholder="Stock conc. (ng/&micro;L)" value="${escapeHtml(row.conc)}">
        </div>
      </div>`;
    } else {
      html += `<div class="custom-fields">
        <input type="text" data-field="customName" placeholder="Fragment name" value="${escapeHtml(row.customName)}">
        <textarea data-field="customSeq" placeholder="Paste DNA sequence (FASTA or plain, whitespace ignored)&hellip;">${escapeHtml(row.customSeq)}</textarea>
        <div class="conc-row">
          <input type="text" inputmode="decimal" data-field="conc" placeholder="Stock conc. (ng/&micro;L)" value="${escapeHtml(row.conc)}">
          <span class="chip-len mono custom-len-badge" style="align-self:center; white-space:nowrap;">${len!=null ? len+' bp' : ''}</span>
        </div>
      </div>`;
    }
  
    container.innerHTML = html;
  
    container.querySelectorAll('[data-mode]').forEach(btn => {
      btn.addEventListener('click', () => {
        if (row.mode !== btn.dataset.mode){ row.mode = btn.dataset.mode; renderAll(); }
      });
    });
    if (opts.onRemove){
      container.querySelector('[data-action="remove"]').addEventListener('click', opts.onRemove);
    }
  
    if (row.mode === 'library'){
      const input = container.querySelector('.combo-input');
      const panel = container.querySelector('.combo-panel');
      wireCombobox(input, panel, opts.items, row, () => renderAll());
      container.querySelector('.conc-input').addEventListener('input', e => { row.conc = e.target.value; renderResults(); });
    } else {
      container.querySelector('[data-field="customName"]').addEventListener('input', e => { row.customName = e.target.value; renderResults(); updateAssemblyName(); });
      const seqEl = container.querySelector('[data-field="customSeq"]');
      if (seqEl){
        seqEl.addEventListener('input', e => {
          row.customSeq = e.target.value;
          const badge = container.querySelector('.custom-len-badge');
          const seq = cleanSeq(row.customSeq);
          if (badge) badge.textContent = seq.length ? seq.length + ' bp' : '';
          renderResults();
        });
      }
      const lenEl = container.querySelector('[data-field="customLen"]');
      if (lenEl){
        lenEl.addEventListener('input', e => { row.customLen = e.target.value; renderResults(); });
      }
      container.querySelector('[data-field="conc"]').addEventListener('input', e => { row.conc = e.target.value; renderResults(); });
    }
  }
  
  // Single delegated listener (set up once, see Init) closes whichever combo panel is open.
  let openCombo = null;
  function closeOpenCombo(){ if (openCombo){ openCombo.panel.hidden = true; openCombo = null; } }
  
  function wireCombobox(input, panel, items, row, onPick){
    let activeIdx = -1;
    let filtered = [];
  
    function renderPanel(){
      const q = input.value.trim().toLowerCase();
      filtered = !q ? items.slice(0, 60) : items.filter(it => it.raw.n.toLowerCase().includes(q)).slice(0, 60);
      if (filtered.length === 0){
        panel.innerHTML = `<div class="combo-empty">No matches</div>`;
      } else {
        panel.innerHTML = filtered.map((it,i) => `
          <div class="combo-item${i===activeIdx?' active':''}" role="option" data-idx="${i}">
            ${catDot(it.cat)}
            <span class="ci-name">${escapeHtml(it.raw.n)}</span>
            <span class="ci-len mono">${it.raw.s.length} bp</span>
          </div>`).join('');
        panel.querySelectorAll('.combo-item').forEach(el => {
          el.addEventListener('mousedown', (e) => {
            e.preventDefault();
            const it = filtered[+el.dataset.idx];
            row.part = it; row.cat = it.cat;
            panel.hidden = true; openCombo = null;
            onPick();
          });
        });
      }
      panel.hidden = false;
      openCombo = { input, panel };
    }
  
    input.addEventListener('focus', () => { activeIdx=-1; renderPanel(); });
    input.addEventListener('input', () => {
      if (row.part && input.value !== row.part.raw.n) row.part = null;
      activeIdx = -1;
      renderPanel();
    });
    input.addEventListener('keydown', (e) => {
      if (panel.hidden) return;
      if (e.key === 'ArrowDown'){ e.preventDefault(); activeIdx = Math.min(activeIdx+1, filtered.length-1); renderPanel(); }
      else if (e.key === 'ArrowUp'){ e.preventDefault(); activeIdx = Math.max(activeIdx-1, 0); renderPanel(); }
      else if (e.key === 'Enter'){ e.preventDefault(); if (activeIdx>=0 && filtered[activeIdx]){ const it=filtered[activeIdx]; row.part=it; row.cat=it.cat; panel.hidden=true; openCombo=null; onPick(); } }
      else if (e.key === 'Escape'){ panel.hidden = true; openCombo = null; }
    });
  }
  
  function renderAcceptorRow(){
    const container = document.getElementById('acceptor-row');
    const items = ACCEPT.map(a => ({ kind:'acceptor', raw:a, cat:'Acceptor Vector' }));
    renderPartRow(container, state.acceptor, { label:'Acceptor', items, onRemove:null });
  }
  
  function renderInsertRows(){
    const wrap = document.getElementById('insert-rows');
    wrap.innerHTML = '';
    state.inserts.forEach((row, i) => {
      const div = document.createElement('div');
      div.className = 'part-row';
      wrap.appendChild(div);
      renderPartRow(div, row, {
        label: 'Insert ' + (i+1),
        items: INSERTABLE,
        onRemove: () => { state.inserts.splice(i,1); renderAll(); }
      });
    });
    document.getElementById('insert-count').textContent = '(' + state.inserts.length + ')';
  }
  
  /* ============ Results ============ */
  function buildRow(row, targetFmol, fallbackName){
    const conc = parseFloat(row.conc);
    const hasConc = conc > 0;
    let len = null, vol = null, hasData = false;
    if (state.calcMode === 'length'){
      len = row.mode === 'custom' ? rowLenBp(row) : (row.part ? row.part.raw.s.length : null);
      hasData = len != null;
      vol = (hasData && hasConc) ? volReqLengthUl(targetFmol, len, conc) : null;
    } else {
      const seq = rowSeq(row);
      hasData = !!seq;
      len = hasData ? seq.length : null;
      vol = (hasData && hasConc) ? volReqUl(targetFmol, seq, conc) : null;
    }
    return {
      name: rowName(row, fallbackName),
      cat: rowCat(row) || 'Custom',
      len, vol, hasData, hasConc, complete: hasData && hasConc
    };
  }
  
  function computeResults(){
    const rows = [];
    const accRow = buildRow(state.acceptor, state.accFmol, 'Acceptor');
    accRow.cat = 'Acceptor Vector';
    rows.push(accRow);
    state.inserts.forEach((ins, i) => rows.push(buildRow(ins, state.insFmol, 'Insert ' + (i+1))));
  
    const cur = mmCurrent();
    const bufferType = state.mm.bufferType;
    const bufferTotal = bufferType === 'neb' ? cur.buffer : (cur.bufferA + cur.bufferB);
  
    const incomplete = rows.filter(r => !r.complete);
    const partVolSum = rows.reduce((s,r) => s + (r.vol||0), 0);
    const water = state.totalVol - (bufferTotal + cur.ligase + cur.enzyme) - partVolSum;
  
    return { rows, bufferType, bufferTotal, bufferA: cur.bufferA, bufferB: cur.bufferB, ligase: cur.ligase, enzyme: cur.enzyme, water, incomplete, partVolSum };
  }
  
  function renderResults(){
    const res = computeResults();
    const tbody = document.getElementById('results-tbody');
    let html = '';
  
    res.rows.forEach(r => {
      html += `<tr>
        <td class="rname">${escapeHtml(r.name)}${r.len!=null?` <span class="dim mono">${r.len}bp</span>`:''}</td>
        <td>${catDot(r.cat)}</td>
        <td class="num">${r.vol!=null ? fmt(r.vol,2) : (r.hasData ? '<span class="dim">need conc.</span>' : '<span class="dim">&ndash;</span>')}</td>
      </tr>`;
    });
  
    html += `<tr class="divider"><td colspan="3"></td></tr>`;
    if (res.bufferType === 'neb'){
      html += `<tr><td>T4 ligase buffer (10&times;)</td><td></td><td class="num">${fmt(res.bufferTotal,2)}</td></tr>`;
    } else {
      html += `<tr><td>Buffer &ndash; Part A</td><td></td><td class="num">${fmt(res.bufferA,2)}</td></tr>`;
      html += `<tr><td>Buffer &ndash; Part B</td><td></td><td class="num">${fmt(res.bufferB,2)}</td></tr>`;
    }
    html += `<tr><td>T4 DNA ligase</td><td></td><td class="num">${fmt(res.ligase,2)}</td></tr>`;
    html += `<tr><td>Restriction enzyme</td><td></td><td class="num">${fmt(res.enzyme,2)}</td></tr>`;
  
    const waterOk = res.incomplete.length === 0;
    html += `<tr><td>Water</td><td></td><td class="num ${waterOk && res.water<0 ? 'neg':''}">${waterOk ? fmt(res.water,2) : '<span class="dim">&ndash;</span>'}</td></tr>`;
    html += `<tr class="total"><td>Total</td><td></td><td class="num">${fmt(state.totalVol,2)}</td></tr>`;
  
    tbody.innerHTML = html;
  
    const warnEl = document.getElementById('warning-msg');
    const noteEl = document.getElementById('note-msg');
    if (waterOk && res.water < 0){
      warnEl.hidden = false;
      warnEl.textContent = `Parts and master mix exceed the reaction volume by ${fmt(Math.abs(res.water),2)} \u00b5L. Increase the reaction volume, or use more concentrated stocks.`;
    } else {
      warnEl.hidden = true;
    }
    if (!waterOk){
      noteEl.hidden = false;
      noteEl.textContent = state.calcMode === 'length'
        ? `Enter a length (bp) and stock concentration for every selected fragment to calculate water volume.`
        : `Enter a stock concentration for every selected fragment to calculate water volume.`;
    } else {
      noteEl.hidden = true;
    }
  
    renderCyclingCaption();
  }
  
  function renderCyclingCaption(){
    const el = document.getElementById('enzyme-caption');
    if (!el) return;
    const enz = (state.acceptor.mode === 'library' && state.acceptor.part) ? state.acceptor.part.raw.enz : null;
    el.hidden = !enz;
    if (enz) el.textContent = `Selected acceptor uses ${enz}.`;
  }
  
  function updateAssemblyName(){
    const parts = [];
    const accName = rowName(state.acceptor, null);
    if (accName) parts.push(accName);
    state.inserts.forEach(ins => { const n = rowName(ins, null); if (n) parts.push(n); });
    document.getElementById('assembly-name').value = parts.join('_');
  }
  
  function renderAll(){
    renderAcceptorRow();
    renderInsertRows();
    renderResults();
    updateAssemblyName();
  }
  
  /* ============ Master mix rendering ============ */
  function mmFieldRow(key, label, value, showReset, defaultHint){
    return `<label class="mm-field">
      <span>${label} <span class="unit">&micro;L</span></span>
      <div class="mm-field-input">
        <input type="text" inputmode="decimal" data-mm="${key}" value="${roundStr(value,3)}">
        ${defaultHint ? `<button type="button" class="reset-link" data-mm-reset="${key}"${showReset?'':' hidden'}>reset (${defaultHint})</button>` : ''}
      </div>
    </label>`;
  }
  function wireMmField(container, key, onInput, onReset){
    const input = container.querySelector(`[data-mm="${key}"]`);
    const resetBtn = container.querySelector(`[data-mm-reset="${key}"]`);
    if (input){
      input.addEventListener('input', e => {
        onInput(parseFloat(e.target.value)||0);
        if (resetBtn) resetBtn.hidden = false;
        renderResults();
      });
    }
    if (resetBtn && onReset){
      resetBtn.addEventListener('click', () => { onReset(); renderMasterMix(); renderResults(); });
    }
  }
  function renderMasterMix(){
    const c = document.getElementById('mastermix-fields');
    const cur = mmCurrent();
    let html = `<div class="radio-row">
        <label class="radio-option"><input type="radio" name="bufferType" value="neb"${state.mm.bufferType==='neb'?' checked':''}> NEB 10&times; T4 ligase buffer</label>
        <label class="radio-option"><input type="radio" name="bufferType" value="homemade"${state.mm.bufferType==='homemade'?' checked':''}> Home-made 2-part buffer</label>
      </div>`;
    if (state.mm.bufferType === 'neb'){
      html += mmFieldRow('buffer', 'T4 ligase buffer (10&times;)', cur.buffer, state.mm.bufferTouched, 'vol&divide;10');
    } else {
      html += mmFieldRow('bufferA', 'Buffer &ndash; Part A', cur.bufferA, false, null);
      html += mmFieldRow('bufferB', 'Buffer &ndash; Part B', cur.bufferB, false, null);
    }
    html += mmFieldRow('ligase', 'T4 DNA ligase', cur.ligase, state.mm.ligaseTouched, 'vol&divide;40');
    html += mmFieldRow('enzyme', 'Restriction enzyme', cur.enzyme, state.mm.enzymeTouched, 'vol&divide;20');
    c.innerHTML = html;
  
    c.querySelectorAll('input[name="bufferType"]').forEach(r => r.addEventListener('change', e => {
      state.mm.bufferType = e.target.value; renderMasterMix(); renderResults();
    }));
    wireMmField(c, 'buffer', v => { state.mm.bufferVol = v; state.mm.bufferTouched = true; }, () => { state.mm.bufferTouched = false; });
    wireMmField(c, 'bufferA', v => { state.mm.bufferAVol = v; });
    wireMmField(c, 'bufferB', v => { state.mm.bufferBVol = v; });
    wireMmField(c, 'ligase', v => { state.mm.ligaseVol = v; state.mm.ligaseTouched = true; }, () => { state.mm.ligaseTouched = false; });
    wireMmField(c, 'enzyme', v => { state.mm.enzymeVol = v; state.mm.enzymeTouched = true; }, () => { state.mm.enzymeTouched = false; });
  }
  
  /* ============ Setup field wiring ============ */
  document.getElementById('in-accFmol').addEventListener('input', e => { state.accFmol = parseFloat(e.target.value)||0; renderResults(); });
  document.getElementById('in-insFmol').addEventListener('input', e => { state.insFmol = parseFloat(e.target.value)||0; renderResults(); });
  document.getElementById('in-totalVol').addEventListener('input', e => {
    state.totalVol = parseFloat(e.target.value)||0;
    renderMasterMix();
    renderResults();
  });
  document.querySelectorAll('input[name="calcMode"]').forEach(r => {
    r.addEventListener('change', e => { state.calcMode = e.target.value; renderAll(); });
  });
  document.getElementById('add-insert-btn').addEventListener('click', () => {
    const lastMode = state.inserts.length ? state.inserts[state.inserts.length-1].mode : 'library';
    state.inserts.push(mkInsert(lastMode));
    renderAll();
  });
  
  /* ============ Tabs ============ */
  document.querySelectorAll('.tab-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('.tab-btn').forEach(b => { b.classList.remove('active'); b.setAttribute('aria-selected','false'); });
      btn.classList.add('active'); btn.setAttribute('aria-selected','true');
      document.querySelectorAll('.tab-panel').forEach(p => p.hidden = true);
      document.getElementById('tab-' + btn.dataset.tab).hidden = false;
      if (btn.dataset.tab === 'db') renderDbTable();
    });
  });
  
  /* ============ Part Library tab ============ */
  const ALL_DB = [];
  ACCEPT.forEach(a => ALL_DB.push({
    kind:'acceptor', raw:a, name:a.n, cat:'Acceptor Vector', len:a.s.length,
    details: [a.enz, a.lvl!=null&&a.lvl!=='' ? 'Level '+a.lvl : null, a.pos&&a.pos!=='-' ? 'Position '+a.pos : null, a.mark].filter(Boolean).join(' \u00b7 ')
  }));
  PARTS.forEach(p => ALL_DB.push({
    kind:'part', raw:p, name:p.n, cat:(p.t||'').trim(), len:p.s.length,
    details: `5\u2032 ${p.s5||''}  \u00b7  3\u2032 ${p.s3||''}`
  }));
  LINKERS.forEach(l => ALL_DB.push({
    kind:'linker', raw:l, name:l.n, cat:(l.t||'').trim(), len:l.s.length,
    details: [l.pos?('Positions '+l.pos):null, l.lvl!=null?('Level '+l.lvl):null].filter(Boolean).join(' \u00b7 ')
  }));
  
  const DB_CATS = Array.from(new Set(ALL_DB.map(d => d.cat))).sort();
  let dbActiveFilter = 'All';
  
  function renderDbFilters(){
    const wrap = document.getElementById('db-filters');
    const cats = ['All', ...DB_CATS];
    wrap.innerHTML = cats.map(c => `<button type="button" class="filter-chip${c===dbActiveFilter?' active':''}" data-cat="${escapeHtml(c)}">${c==='All'?'':catDot(c)}${escapeHtml(c)}</button>`).join('');
    wrap.querySelectorAll('.filter-chip').forEach(chip => {
      chip.addEventListener('click', () => { dbActiveFilter = chip.dataset.cat; renderDbTable(); });
    });
  }
  
  function renderDbTable(){
    const q = document.getElementById('db-search').value.trim().toLowerCase();
    const rows = ALL_DB.filter(d => (dbActiveFilter==='All' || d.cat===dbActiveFilter) && (!q || d.name.toLowerCase().includes(q)));
    const tbody = document.getElementById('db-tbody');
    tbody.innerHTML = rows.map((d, i) => `
      <tr>
        <td class="dname">${escapeHtml(d.name)}</td>
        <td>${catPill(d.cat)}</td>
        <td class="dlen mono">${d.len} bp</td>
        <td class="ddetails">${escapeHtml(d.details)}</td>
        <td class="dactions">
          <button type="button" class="icon-btn" data-copy="${i}">Copy seq</button>
          <button type="button" class="icon-btn" data-use="${i}">Use&nbsp;&rarr;</button>
        </td>
      </tr>`).join('');
  
    tbody.querySelectorAll('[data-copy]').forEach(b => b.addEventListener('click', () => {
      const d = rows[+b.dataset.copy];
      navigator.clipboard.writeText(d.raw.s).then(() => { b.textContent='Copied'; setTimeout(()=>b.textContent='Copy seq',1200); });
    }));
    tbody.querySelectorAll('[data-use]').forEach(b => b.addEventListener('click', () => {
      const d = rows[+b.dataset.use];
      if (d.kind === 'acceptor'){
        state.acceptor = { mode:'library', part:{kind:'acceptor', raw:d.raw, cat:'Acceptor Vector'}, conc:'', customName:'', customSeq:'', customLen:'' };
      } else {
        const ins = mkInsert();
        ins.part = { kind:d.kind, raw:d.raw, cat:d.cat };
        state.inserts.push(ins);
      }
      renderAll();
      document.querySelector('.tab-btn[data-tab="calc"]').click();
    }));
  
    document.getElementById('db-count').textContent = `${rows.length} of ${ALL_DB.length} parts`;
  }
  document.getElementById('db-search').addEventListener('input', renderDbTable);
  
  /* ============ Init ============ */
  document.addEventListener('click', (e) => {
    if (openCombo && !openCombo.panel.contains(e.target) && e.target !== openCombo.input) closeOpenCombo();
  });
  renderDbFilters();
  renderMasterMix();
  renderAll();
}

loadData();
