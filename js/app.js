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
  // Rounds to at most maxDp decimals but drops trailing zeros, so e.g. 1 stays "1"
  // while 0.25 keeps both decimals -- lets a value's own precision show through
  // instead of padding everything to a fixed number of decimal places.
  function fmtSmart(n, maxDp){
    if (n==null || isNaN(n)) return '&ndash;';
    maxDp = maxDp==null ? 3 : maxDp;
    return n.toFixed(maxDp).replace(/(\.\d*?)0+$/,'$1').replace(/\.$/,'');
  }

  /* ============ Combined searchable insert list ============ */
  const INSERTABLE = [];
  PARTS.forEach(p => INSERTABLE.push({ kind:'part', raw:p, cat:(p.t||'').trim() }));
  LINKERS.forEach(l => INSERTABLE.push({ kind:'linker', raw:l, cat:(l.t||'').trim() }));

  /* ============ App state ============ */
  let uidCounter = 1;
  function mkInsert(mode){ return { id:'ins'+(uidCounter++), mode: mode||'custom', part:null, cat:null, conc:'', customName:'', customSeq:'', customLen:'', calcModeOverride:null }; }
  function mkAcceptor(){ return { mode:'custom', part:null, conc:'', customName:'', customSeq:'', customLen:'', calcModeOverride:null }; }
  function mkAssembly(){
    return {
      id: 'asm'+(uidCounter++),
      collapsed: false,
      acceptor: mkAcceptor(),
      inserts: [ mkInsert() ],
      checks: {},
    };
  }
  function duplicateAssembly(asm){
    return {
      id: 'asm'+(uidCounter++),
      collapsed: false,
      acceptor: { ...asm.acceptor },
      inserts: asm.inserts.map(row => ({ ...row, id: 'ins'+(uidCounter++) })),
      checks: {},
    };
  }
  const state = {
    accFmol: 25,
    insFmol: 50,
    totalVol: 10,
    calcMode: 'sequence', // 'sequence' | 'length'
    mm: {
      enabled: false,               // shared master mix batching across assemblies (off = each assembly computed independently)
      bufferType: 'homemade',      // 'neb' | 'homemade'
      bufferVol: null, bufferTouched: false,   // NEB 10x buffer override (default totalVol/10)
      bufferAVol: 1, bufferBVol: 1,            // home-made 2-part buffer, freely editable
      ligaseVol: null, ligaseTouched: false,   // default totalVol/40
      enzymeVol: null, enzymeTouched: false,   // default totalVol/20
    },
    mmOverage: 1,       // extra reactions' worth of common master mix, editable
    mmChecks: {},       // checkbox state for the shared batch rows
    assemblies: [ mkAssembly() ],
    cycling: defaultCycling(),  // rampRate (deg C/sec, rough estimate) + heatupMin + editable step/cycle blocks
  };
  // Recommended digestion / heat-inactivation temperatures per Type IIS enzyme (rough, manufacturer defaults).
  // Ligation always runs at 16°C (T4 ligase) regardless of the restriction enzyme, so it isn't listed here.
  const CYCLING_ENZYMES = [
    { key:'bsai',   label:'BsaI / BsaI-HFv2',      digest:37, heatInact:65 },
    { key:'bbsi',   label:'BbsI-HF (BpiI)',        digest:37, heatInact:65 },
    { key:'bspqi',  label:'BspQI (SapI)',          digest:37, heatInact:65 },
    { key:'bsmbi',  label:'BsmBI-v2 (Esp3I)',      digest:42, heatInact:80 },
    { key:'custom', label:'Custom / other',        digest:null, heatInact:null },
  ];
  function enzymeKeyFromName(name){
    const n = (name||'').toLowerCase();
    if (n.includes('bsai')) return 'bsai';
    if (n.includes('bbsi') || n.includes('bpii')) return 'bbsi';
    if (n.includes('bspqi') || n.includes('sapi')) return 'bspqi';
    if (n.includes('bsmbi') || n.includes('esp3i')) return 'bsmbi';
    return null;
  }
  function applyCyclingEnzyme(key){
    const e = CYCLING_ENZYMES.find(x => x.key === key);
    if (!e || e.digest == null) return;
    state.cycling.blocks[0].steps[0].temp = e.digest;      // Digestion
    state.cycling.blocks[1].steps[0].temp = e.digest;      // Final digestion
    state.cycling.blocks[2].steps[0].temp = e.heatInact;   // Heat inactivation
  }
  function defaultCycling(){
    return {
      rampRate: 3, heatupMin: 5,
      enzyme: 'bsai', enzymeTouched: false,   // enzymeTouched: user picked manually, stop auto-syncing to the selected acceptor
      blocks: [
        { cycles: 35, steps: [ { name:'Digestion', temp:37, time:5 }, { name:'Ligation', temp:16, time:5 } ] },
        { cycles: 1,  steps: [ { name:'Final digestion', temp:37, time:10 } ] },
        { cycles: 1,  steps: [ { name:'Heat inactivation', temp:65, time:20 } ] },
      ],
    };
  }

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
  // A row can override the global sequence/length calculation mode; falls back to the global setting.
  function rowCalcMode(row){
    return row.calcModeOverride || state.calcMode;
  }
  // Fragment length to *display* next to a row, independent of whether it can be used for a volume calc.
  function displayLen(row){
    if (row.mode === 'library') return row.part ? row.part.raw.s.length : null;
    if (rowCalcMode(row) === 'length') return rowLenBp(row);
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

  /* ============ Checkbox cell ============ */
  function checkCell(checksObj, key){
    return `<input type="checkbox" class="row-check" data-check-key="${escapeHtml(key)}"${checksObj[key]?' checked':''}>`;
  }
  function wireCheckCells(container, checksObj, onChange){
    container.querySelectorAll('[data-check-key]').forEach(cb => {
      cb.addEventListener('change', e => {
        checksObj[cb.dataset.checkKey] = e.target.checked;
        if (onChange) onChange();
      });
    });
  }

  /* ============ Rendering: Acceptor + Insert rows ============ */
  function renderPartRow(container, row, opts){
    // opts: {label, items, onRemove}
    const len = displayLen(row);
    const lengthMode = rowCalcMode(row) === 'length';

    let html = `<div class="part-row-head">
        <span class="part-row-label">${opts.label}</span>
        <div class="row-actions">
          <select class="calc-mode-select" title="Volume calculation method for this fragment">
            <option value=""${!row.calcModeOverride?' selected':''}>Global (${state.calcMode==='length'?'Length':'Sequence'})</option>
            <option value="sequence"${row.calcModeOverride==='sequence'?' selected':''}>Sequence composition</option>
            <option value="length"${row.calcModeOverride==='length'?' selected':''}>Length only</option>
          </select>
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
    container.querySelector('.calc-mode-select').addEventListener('change', e => {
      row.calcModeOverride = e.target.value || null;
      renderAll();
    });
    if (opts.onRemove){
      container.querySelector('[data-action="remove"]').addEventListener('click', opts.onRemove);
    }

    if (row.mode === 'library'){
      const input = container.querySelector('.combo-input');
      const panel = container.querySelector('.combo-panel');
      wireCombobox(input, panel, opts.items, row, () => renderAll());
      container.querySelector('.conc-input').addEventListener('input', e => { row.conc = e.target.value; refreshAllResults(); });
    } else {
      container.querySelector('[data-field="customName"]').addEventListener('input', e => { row.customName = e.target.value; renderAssemblyHeaders(); refreshAllResults(); });
      const seqEl = container.querySelector('[data-field="customSeq"]');
      if (seqEl){
        seqEl.addEventListener('input', e => {
          row.customSeq = e.target.value;
          const badge = container.querySelector('.custom-len-badge');
          const seq = cleanSeq(row.customSeq);
          if (badge) badge.textContent = seq.length ? seq.length + ' bp' : '';
          refreshAllResults();
        });
      }
      const lenEl = container.querySelector('[data-field="customLen"]');
      if (lenEl){
        lenEl.addEventListener('input', e => { row.customLen = e.target.value; refreshAllResults(); });
      }
      container.querySelector('[data-field="conc"]').addEventListener('input', e => { row.conc = e.target.value; refreshAllResults(); });
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

  /* ============ Per-assembly computation ============ */
  function buildRow(row, targetFmol, fallbackName){
    const conc = parseFloat(row.conc);
    const hasConc = conc > 0;
    let len = null, vol = null, hasData = false;
    if (rowCalcMode(row) === 'length'){
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

  function computeAssembly(asm){
    const rows = [];
    const accRow = buildRow(asm.acceptor, state.accFmol, 'Acceptor');
    accRow.cat = 'Acceptor Vector';
    accRow.key = 'acceptor';
    rows.push(accRow);
    asm.inserts.forEach((ins, i) => {
      const r = buildRow(ins, state.insFmol, 'Insert ' + (i+1));
      r.key = 'ins' + ins.id;
      rows.push(r);
    });
    const partVolSum = rows.reduce((s,r) => s + (r.vol||0), 0);
    const complete = rows.every(r => r.complete);
    return { rows, partVolSum, complete };
  }

  function assemblyName(asm){
    const parts = [];
    const accName = rowName(asm.acceptor, null);
    if (accName) parts.push(accName);
    asm.inserts.forEach(ins => { const n = rowName(ins, null); if (n) parts.push(n); });
    return parts.join('_');
  }

  /* ============ Global (batched master mix) computation ============ */
  function computeGlobal(){
    const cur = mmCurrent();
    const bufferTotal = state.mm.bufferType === 'neb' ? cur.buffer : (cur.bufferA + cur.bufferB);
    const commonPerRxn = bufferTotal + cur.ligase + cur.enzyme;

    const perAssembly = state.assemblies.map(asm => {
      const res = computeAssembly(asm);
      const water = res.complete ? (state.totalVol - commonPerRxn - res.partVolSum) : null;
      return { asm, res, water };
    });

    const allComplete = perAssembly.every(p => p.res.complete);
    const N = state.assemblies.length;
    const overage = parseFloat(state.mmOverage) || 0;
    const batchFactor = N + overage;
    const enabled = state.mm.enabled;

    let minWater = null;
    if (enabled && allComplete && N > 0){
      minWater = Math.min(...perAssembly.map(p => p.water));
    }

    perAssembly.forEach(p => {
      if (enabled && minWater != null){
        p.aliquot = commonPerRxn + minWater;
        p.topup = p.water - minWater;
        p.total = p.aliquot + p.res.partVolSum + p.topup;
      } else {
        p.aliquot = null; p.topup = null;
        p.total = p.water != null ? (commonPerRxn + p.res.partVolSum + p.water) : null;
      }
    });

    return {
      enabled,
      bufferType: state.mm.bufferType, bufferTotal, bufferA: cur.bufferA, bufferB: cur.bufferB,
      ligase: cur.ligase, enzyme: cur.enzyme, commonPerRxn,
      perAssembly, allComplete, N, overage, batchFactor, minWater,
      bufferBatch: (enabled && minWater!=null) ? bufferTotal*batchFactor : null,
      bufferABatch: (enabled && minWater!=null) ? cur.bufferA*batchFactor : null,
      bufferBBatch: (enabled && minWater!=null) ? cur.bufferB*batchFactor : null,
      ligaseBatch: (enabled && minWater!=null) ? cur.ligase*batchFactor : null,
      enzymeBatch: (enabled && minWater!=null) ? cur.enzyme*batchFactor : null,
      waterBatch: (enabled && minWater!=null) ? minWater*batchFactor : null,
    };
  }

  /* ============ Rendering: Assemblies ============ */
  function renderAssemblies(){
    const wrap = document.getElementById('assembly-list');
    wrap.innerHTML = '';
    const glob = computeGlobal();

    state.assemblies.forEach((asm, ai) => {
      const p = glob.perAssembly[ai];
      const card = document.createElement('div');
      card.className = 'assembly-card';
      wrap.appendChild(card);

      const name = assemblyName(asm) || ('Assembly ' + (ai+1));
      let head = `<div class="assembly-head">
        <button type="button" class="assembly-collapse" data-action="collapse" title="${asm.collapsed?'Expand':'Collapse'}">${asm.collapsed?'&#9656;':'&#9662;'}</button>
        <span class="assembly-name mono">${escapeHtml(name)}</span>
        <button type="button" class="duplicate-btn" data-action="duplicate" title="Duplicate assembly">&#10697;</button>
        ${state.assemblies.length>1 ? `<button type="button" class="remove-btn" data-action="remove" title="Remove assembly">&times;</button>` : ''}
      </div>`;
      card.innerHTML = head;
      card.querySelector('[data-action="collapse"]').addEventListener('click', () => {
        asm.collapsed = !asm.collapsed; renderAssemblies();
      });
      card.querySelector('[data-action="duplicate"]').addEventListener('click', () => {
        state.assemblies.splice(ai+1, 0, duplicateAssembly(asm)); renderAll();
      });
      if (state.assemblies.length>1){
        card.querySelector('[data-action="remove"]').addEventListener('click', () => {
          state.assemblies.splice(ai,1); renderAll();
        });
      }

      if (asm.collapsed) return;

      const body = document.createElement('div');
      body.className = 'assembly-body';
      card.appendChild(body);

      const accWrap = document.createElement('div');
      accWrap.className = 'part-row part-row-acceptor';
      body.appendChild(accWrap);
      renderPartRow(accWrap, asm.acceptor, { label:'Acceptor', items: ACCEPT.map(a => ({ kind:'acceptor', raw:a, cat:'Acceptor Vector' })), onRemove:null });

      const insertsHeader = document.createElement('div');
      insertsHeader.className = 'assembly-subhead';
      insertsHeader.innerHTML = `<span>Insert fragments <span class="count">(${asm.inserts.length})</span></span>`;
      body.appendChild(insertsHeader);

      asm.inserts.forEach((row, i) => {
        const div = document.createElement('div');
        div.className = 'part-row';
        body.appendChild(div);
        renderPartRow(div, row, {
          label: 'Insert ' + (i+1),
          items: INSERTABLE,
          onRemove: () => { asm.inserts.splice(i,1); renderAll(); }
        });
      });

      const addBtn = document.createElement('button');
      addBtn.type = 'button'; addBtn.className = 'btn-add'; addBtn.textContent = '+ Add insert fragment';
      addBtn.addEventListener('click', () => {
        const lastMode = asm.inserts.length ? asm.inserts[asm.inserts.length-1].mode : 'custom';
        asm.inserts.push(mkInsert(lastMode));
        renderAll();
      });
      body.appendChild(addBtn);

      const resultsWrap = document.createElement('div');
      resultsWrap.className = 'assembly-results';
      resultsWrap.innerHTML = renderAssemblyResultsHtml(p, glob);
      body.appendChild(resultsWrap);
      wireCheckCells(resultsWrap, asm.checks, null);
    });

    document.getElementById('assembly-count').textContent = '(' + state.assemblies.length + ')';
  }

  function renderAssemblyHeaders(){
    // cheap re-render of just the name/heading text without losing focus in open inputs
    document.querySelectorAll('.assembly-card').forEach((card, ai) => {
      const asm = state.assemblies[ai];
      if (!asm) return;
      const nameEl = card.querySelector('.assembly-name');
      if (nameEl) nameEl.textContent = assemblyName(asm) || ('Assembly ' + (ai+1));
    });
  }

  function renderAssemblyResultsHtml(p, glob){
    let html = `<table class="results-table">
      <thead><tr><th>Component</th><th></th><th style="text-align:right;text-transform:none">&micro;L</th><th></th></tr></thead>
      <tbody>`;
    p.res.rows.forEach(r => {
      html += `<tr>
        <td class="rname">${escapeHtml(r.name)}${r.len!=null?` <span class="dim mono">${r.len}bp</span>`:''}</td>
        <td>${catDot(r.cat)}</td>
        <td class="num">${r.vol!=null ? fmt(r.vol,2) : (r.hasData ? '<span class="dim">need conc.</span>' : '<span class="dim">&ndash;</span>')}</td>
        <td class="chk">${checkCell(p.asm.checks, r.key)}</td>
      </tr>`;
    });
    html += `<tr class="divider"><td colspan="4"></td></tr>`;

    if (glob.enabled){
      html += `<tr><td>Master mix aliquot</td><td></td><td class="num">${p.aliquot!=null ? fmt(p.aliquot,2) : '<span class="dim">&ndash;</span>'}</td><td class="chk">${checkCell(p.asm.checks,'aliquot')}</td></tr>`;
      const topupOk = p.topup!=null;
      html += `<tr><td>Top-up water</td><td></td><td class="num ${topupOk && p.topup<0 ? 'neg':''}">${topupOk ? fmt(p.topup,2) : '<span class="dim">&ndash;</span>'}</td><td class="chk">${checkCell(p.asm.checks,'topup')}</td></tr>`;
    } else {
      if (glob.bufferType === 'neb'){
        html += `<tr><td>T4 ligase buffer (10&times;)</td><td></td><td class="num">${fmt(glob.bufferTotal,2)}</td><td class="chk">${checkCell(p.asm.checks,'buffer')}</td></tr>`;
      } else {
        html += `<tr><td>Buffer &ndash; Part A</td><td></td><td class="num">${fmt(glob.bufferA,2)}</td><td class="chk">${checkCell(p.asm.checks,'bufferA')}</td></tr>`;
        html += `<tr><td>Buffer &ndash; Part B</td><td></td><td class="num">${fmt(glob.bufferB,2)}</td><td class="chk">${checkCell(p.asm.checks,'bufferB')}</td></tr>`;
      }
      html += `<tr><td>T4 DNA ligase</td><td></td><td class="num">${fmt(glob.ligase,2)}</td><td class="chk">${checkCell(p.asm.checks,'ligase')}</td></tr>`;
      html += `<tr><td>Restriction enzyme</td><td></td><td class="num">${fmt(glob.enzyme,2)}</td><td class="chk">${checkCell(p.asm.checks,'enzyme')}</td></tr>`;
      const waterOk = p.water!=null;
      html += `<tr><td>Water</td><td></td><td class="num ${waterOk && p.water<0 ? 'neg':''}">${waterOk ? fmt(p.water,2) : '<span class="dim">&ndash;</span>'}</td><td class="chk">${checkCell(p.asm.checks,'water')}</td></tr>`;
    }
    html += `<tr class="total"><td>Total</td><td></td><td class="num">${p.total!=null ? fmt(p.total,2) : fmt(state.totalVol,2)}</td><td></td></tr>`;
    html += `</tbody></table>`;
    if (!p.res.complete){
      html += `<p class="note">Enter${state.calcMode==='length'?' a length (bp) and':''} a stock concentration for every selected fragment to calculate water volume.</p>`;
    }
    return html;
  }

  /* ============ Master mix batch box (aside, top) ============ */
  function batchRow(label, batchVal, key, batchFactor, perRxnVal){
    const ok = batchVal!=null;
    const calc = ok ? `<span class="mm-calc-hint mono">${fmtSmart(batchFactor)}&thinsp;&times;&thinsp;${fmtSmart(perRxnVal)}&thinsp;=</span>` : '';
    return `<tr><td>${label}</td><td></td><td class="num ${ok && batchVal<0 ? 'neg':''}">${ok ? calc+fmtSmart(batchVal) : '<span class="dim">&ndash;</span>'}</td><td class="chk">${checkCell(state.mmChecks,key)}</td></tr>`;
  }

  // Builds the skeleton (overage field + results table) once when master mix is switched on;
  // left in place afterwards so typing in the overage field doesn't need to rebuild it.
  function renderMasterMixBatch(){
    const box = document.getElementById('mastermix-batch');
    if (!state.mm.enabled){
      box.innerHTML = '';
      renderCyclingCaption();
      return;
    }
    box.innerHTML = `
      <label class="mm-field">
        <span>Extra reactions <span class="unit">overage</span></span>
        <div class="mm-field-input">
          <input type="text" inputmode="decimal" id="in-mmOverage" value="${roundStr(state.mmOverage,3)}">
        </div>
      </label>
      <p class="note" id="mm-batch-caption"></p>
      <table class="results-table" id="results-table">
        <thead><tr><th>Component</th><th></th><th style="text-align:right;text-transform:none">&micro;L</th><th></th></tr></thead>
        <tbody id="results-tbody"></tbody>
      </table>
      <p id="warning-msg" class="warning" hidden></p>
      <p id="note-msg" class="note" hidden></p>`;
    document.getElementById('in-mmOverage').addEventListener('input', e => {
      state.mmOverage = parseFloat(e.target.value)||0;
      renderBatchResults();
    });
    renderBatchResults();
  }

  // Numbers-only update of the batch results table; assumes renderMasterMixBatch() has already
  // built the skeleton. Safe to call on every keystroke elsewhere without losing input focus.
  function renderBatchResults(){
    if (!state.mm.enabled) return;
    const glob = computeGlobal();
    const tbody = document.getElementById('results-tbody');
    if (!tbody) return;
    const captionEl = document.getElementById('mm-batch-caption');
    const warnEl = document.getElementById('warning-msg');
    const noteEl = document.getElementById('note-msg');
    const bf = glob.batchFactor;

    captionEl.textContent = `Batch for ${glob.N} assembl${glob.N===1?'y':'ies'} + ${fmtSmart(glob.overage)} overage = ${fmtSmart(bf)}× reactions.`;

    let html = '';
    if (glob.bufferType === 'neb'){
      html += batchRow('T4 ligase buffer (10&times;)', glob.bufferBatch, 'buffer', bf, glob.bufferTotal);
    } else {
      html += batchRow('Buffer &ndash; Part A', glob.bufferABatch, 'bufferA', bf, glob.bufferA);
      html += batchRow('Buffer &ndash; Part B', glob.bufferBBatch, 'bufferB', bf, glob.bufferB);
    }
    html += batchRow('T4 DNA ligase', glob.ligaseBatch, 'ligase', bf, glob.ligase);
    html += batchRow('Restriction enzyme', glob.enzymeBatch, 'enzyme', bf, glob.enzyme);
    const waterOk = glob.waterBatch!=null;
    html += batchRow('Water (shared minimum)', glob.waterBatch, 'water', bf, glob.minWater);
    html += `<tr class="total"><td>Aliquot per reaction</td><td></td><td class="num">${waterOk ? fmtSmart(glob.commonPerRxn+glob.minWater) : '<span class="dim">&ndash;</span>'}</td><td></td></tr>`;

    tbody.innerHTML = html;
    wireCheckCells(tbody, state.mmChecks, null);

    const anyNegative = waterOk && (glob.waterBatch < 0 || glob.perAssembly.some(p => p.topup!=null && p.topup<0));
    if (anyNegative){
      warnEl.hidden = false;
      warnEl.textContent = `Parts and master mix exceed the reaction volume for at least one assembly. Increase the reaction volume, or use more concentrated stocks.`;
    } else {
      warnEl.hidden = true;
    }
    if (!glob.allComplete){
      noteEl.hidden = false;
      noteEl.textContent = state.calcMode === 'length'
        ? `Enter a length (bp) and stock concentration for every fragment in every assembly to calculate the shared master mix.`
        : `Enter a stock concentration for every fragment in every assembly to calculate the shared master mix.`;
    } else {
      noteEl.hidden = true;
    }

    renderCyclingCaption();
  }

  function renderCyclingCaption(){
    const el = document.getElementById('enzyme-caption');
    if (!el) return;
    const first = state.assemblies[0];
    const enz = (first && first.acceptor.mode === 'library' && first.acceptor.part) ? first.acceptor.part.raw.enz : null;
    el.hidden = !enz;
    if (enz) el.textContent = `Selected acceptor uses ${enz}.`;

    if (enz && !state.cycling.enzymeTouched){
      const key = enzymeKeyFromName(enz);
      if (key && key !== state.cycling.enzyme){
        state.cycling.enzyme = key;
        applyCyclingEnzyme(key);
        renderCyclingTable();
      }
    }
  }

  /* ============ Cycling table (editable steps/cycles + rough estimated run time) ============ */
  // Sums step times plus ramp time (|Δtemp| / ramp rate) between every consecutive step, including
  // between repeats of a cycled block and between blocks, plus a fixed initial heat-up allowance.
  function cyclingTotalMinutes(){
    const rampPerMin = (state.cycling.rampRate || 0) * 60;
    let total = state.cycling.heatupMin || 0;
    let currentTemp = null;
    state.cycling.blocks.forEach(block => {
      const cycles = Math.max(1, block.cycles || 1);
      for (let c = 0; c < cycles; c++){
        block.steps.forEach(step => {
          if (currentTemp != null && rampPerMin > 0){
            total += Math.abs(step.temp - currentTemp) / rampPerMin;
          }
          total += step.time || 0;
          currentTemp = step.temp;
        });
      }
    });
    return total;
  }
  function fmtDuration(mins){
    if (mins==null || isNaN(mins)) return '&ndash;';
    const totalMin = Math.round(mins);
    const h = Math.floor(totalMin/60), m = totalMin%60;
    return h>0 ? `${h}h ${m}min` : `${m} min`;
  }
  function updateCyclingTotal(){
    const el = document.getElementById('cycling-total');
    if (!el) return;
    el.innerHTML = `<span class="cycling-total-label">Estimated total run time</span>
      <span class="cycling-total-value">~${fmtDuration(cyclingTotalMinutes())}</span>
      <span class="cycling-total-sub">assumes a ${fmtSmart(state.cycling.rampRate)}&deg;C/sec block ramp rate plus ${fmtSmart(state.cycling.heatupMin)} min initial heat-up. Rough estimate only; actual cyclers vary.</span>`;
  }
  function renderCyclingTable(){
    const wrap = document.getElementById('cycling-table-wrap');
    if (!wrap) return;
    let html = `<label class="cycling-enzyme-field">
      <span>Restriction enzyme</span>
      <select id="in-cyclingEnzyme">
        ${CYCLING_ENZYMES.map(e => `<option value="${e.key}"${state.cycling.enzyme===e.key?' selected':''}>${escapeHtml(e.label)}</option>`).join('')}
      </select>
    </label>`;
    html += `<table class="cycling-table">
      <colgroup><col><col class="col-temp"><col class="col-time"><col class="col-cyc"></colgroup>
      <thead><tr><th>Step</th><th>Temp</th><th>Min</th><th>&times;</th></tr></thead><tbody>`;
    state.cycling.blocks.forEach((block, bi) => {
      block.steps.forEach((step, si) => {
        html += `<tr>
          <td>${escapeHtml(step.name)}</td>
          <td class="cyc-temp-static">${fmtSmart(step.temp,1)}&deg;</td>
          <td><input type="text" inputmode="decimal" class="cyc-input cyc-time" data-block="${bi}" data-step="${si}" value="${step.time}"></td>
          ${si===0 ? `<td${block.steps.length>1?` rowspan="${block.steps.length}"`:''}><input type="text" inputmode="numeric" class="cyc-input cyc-cycles" data-block="${bi}" value="${block.cycles}"></td>` : ''}
        </tr>`;
      });
    });
    html += `</tbody></table>`;
    wrap.innerHTML = html;

    wrap.querySelector('#in-cyclingEnzyme').addEventListener('change', e => {
      state.cycling.enzyme = e.target.value;
      state.cycling.enzymeTouched = true;
      applyCyclingEnzyme(state.cycling.enzyme);
      renderCyclingTable();
    });
    wrap.querySelectorAll('.cyc-time').forEach(inp => inp.addEventListener('input', e => {
      state.cycling.blocks[+e.target.dataset.block].steps[+e.target.dataset.step].time = parseFloat(e.target.value) || 0;
      updateCyclingTotal();
    }));
    wrap.querySelectorAll('.cyc-cycles').forEach(inp => inp.addEventListener('input', e => {
      state.cycling.blocks[+e.target.dataset.block].cycles = Math.max(1, parseInt(e.target.value, 10) || 1);
      updateCyclingTotal();
    }));

    updateCyclingTotal();
  }

  function renderAll(){
    renderAssemblies();
    renderMasterMixBatch();
  }

  // Targeted refresh for value-only changes (typing conc/name/length, master mix volumes, fmol/volume
  // targets): updates every assembly's results table + the batch box, without rebuilding the row
  // inputs themselves (which would drop keyboard focus mid-keystroke).
  function refreshAllResults(){
    const glob = computeGlobal();
    document.querySelectorAll('#assembly-list > .assembly-card').forEach((card, ai) => {
      const asm = state.assemblies[ai];
      if (!asm || asm.collapsed) return;
      const p = glob.perAssembly[ai];
      const resultsWrap = card.querySelector('.assembly-results');
      if (resultsWrap){
        resultsWrap.innerHTML = renderAssemblyResultsHtml(p, glob);
        wireCheckCells(resultsWrap, asm.checks, null);
      }
    });
    if (state.mm.enabled) renderBatchResults();
    else renderCyclingCaption();
  }

  /* ============ CSV export ============ */
  function csvEscape(v){
    const s = String(v==null?'':v);
    return /[",\n]/.test(s) ? '"' + s.replace(/"/g,'""') + '"' : s;
  }
  function downloadCsv(){
    const glob = computeGlobal();
    const rows = [['Assembly','Component','Category','Volume (µL)']];

    glob.perAssembly.forEach(p => {
      const name = assemblyName(p.asm) || 'Assembly';
      p.res.rows.forEach(r => rows.push([name, r.name, r.cat, r.vol!=null?r.vol.toFixed(2):'']));
      if (glob.enabled){
        rows.push([name, 'Master mix aliquot', '', p.aliquot!=null?p.aliquot.toFixed(2):'']);
        rows.push([name, 'Top-up water', '', p.topup!=null?p.topup.toFixed(2):'']);
      } else {
        if (glob.bufferType === 'neb'){
          rows.push([name, 'T4 ligase buffer (10x)', '', glob.bufferTotal.toFixed(2)]);
        } else {
          rows.push([name, 'Buffer - Part A', '', glob.bufferA.toFixed(2)]);
          rows.push([name, 'Buffer - Part B', '', glob.bufferB.toFixed(2)]);
        }
        rows.push([name, 'T4 DNA ligase', '', glob.ligase.toFixed(2)]);
        rows.push([name, 'Restriction enzyme', '', glob.enzyme.toFixed(2)]);
        rows.push([name, 'Water', '', p.water!=null?p.water.toFixed(2):'']);
      }
      rows.push([name, 'Total', '', p.total!=null?p.total.toFixed(2):'']);
    });

    if (glob.enabled){
      const batchName = 'Master mix (batch)';
      if (glob.bufferType === 'neb'){
        rows.push([batchName, 'T4 ligase buffer (10x)', '', glob.bufferBatch!=null?glob.bufferBatch.toFixed(2):'']);
      } else {
        rows.push([batchName, 'Buffer - Part A', '', glob.bufferABatch!=null?glob.bufferABatch.toFixed(2):'']);
        rows.push([batchName, 'Buffer - Part B', '', glob.bufferBBatch!=null?glob.bufferBBatch.toFixed(2):'']);
      }
      rows.push([batchName, 'T4 DNA ligase', '', glob.ligaseBatch!=null?glob.ligaseBatch.toFixed(2):'']);
      rows.push([batchName, 'Restriction enzyme', '', glob.enzymeBatch!=null?glob.enzymeBatch.toFixed(2):'']);
      rows.push([batchName, 'Water (shared minimum)', '', glob.waterBatch!=null?glob.waterBatch.toFixed(2):'']);
    }

    const csv = rows.map(r => r.map(csvEscape).join(',')).join('\r\n');
    const blob = new Blob([csv], { type:'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = 'pipetting-volumes.csv';
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
    URL.revokeObjectURL(url);
  }
  document.getElementById('download-csv-btn').addEventListener('click', downloadCsv);

  /* ============ Master mix settings (Global conditions) ============ */
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
        refreshAllResults();
      });
    }
    if (resetBtn && onReset){
      resetBtn.addEventListener('click', () => { onReset(); renderMasterMixSettings(); refreshAllResults(); });
    }
  }
  // These settings (buffer type, per-reaction buffer/ligase/enzyme volumes) are always shown and
  // always in effect — they drive each assembly's own independent calculation, and, when the
  // shared master mix toggle is on, the batch numbers in the Master mix box too.
  function renderMasterMixSettings(){
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
      state.mm.bufferType = e.target.value; renderMasterMixSettings(); refreshAllResults();
    }));
    wireMmField(c, 'buffer', v => { state.mm.bufferVol = v; state.mm.bufferTouched = true; }, () => { state.mm.bufferTouched = false; });
    wireMmField(c, 'bufferA', v => { state.mm.bufferAVol = v; });
    wireMmField(c, 'bufferB', v => { state.mm.bufferBVol = v; });
    wireMmField(c, 'ligase', v => { state.mm.ligaseVol = v; state.mm.ligaseTouched = true; }, () => { state.mm.ligaseTouched = false; });
    wireMmField(c, 'enzyme', v => { state.mm.enzymeVol = v; state.mm.enzymeTouched = true; }, () => { state.mm.enzymeTouched = false; });
  }

  /* ============ Setup field wiring ============ */
  document.getElementById('in-accFmol').addEventListener('input', e => { state.accFmol = parseFloat(e.target.value)||0; refreshAllResults(); });
  document.getElementById('in-insFmol').addEventListener('input', e => { state.insFmol = parseFloat(e.target.value)||0; refreshAllResults(); });
  document.getElementById('in-totalVol').addEventListener('input', e => {
    state.totalVol = parseFloat(e.target.value)||0;
    renderMasterMixSettings();
    refreshAllResults();
  });
  document.getElementById('in-mmEnabled').addEventListener('change', e => {
    state.mm.enabled = e.target.checked;
    renderAll();
  });
  document.querySelectorAll('input[name="calcMode"]').forEach(r => {
    r.addEventListener('change', e => { state.calcMode = e.target.value; renderAll(); });
  });
  document.getElementById('add-assembly-btn').addEventListener('click', () => {
    state.assemblies.push(mkAssembly());
    renderAll();
  });

  document.getElementById('in-rampRate').addEventListener('input', e => { state.cycling.rampRate = parseFloat(e.target.value)||0; updateCyclingTotal(); });
  document.getElementById('in-heatup').addEventListener('input', e => { state.cycling.heatupMin = parseFloat(e.target.value)||0; updateCyclingTotal(); });

  document.getElementById('clear-page-btn').addEventListener('click', () => {
    if (!confirm('Clear everything and start a new page? This cannot be undone.')) return;
    state.accFmol = 25;
    state.insFmol = 50;
    state.totalVol = 10;
    state.calcMode = 'sequence';
    state.mm = {
      enabled: false, bufferType: 'homemade',
      bufferVol: null, bufferTouched: false,
      bufferAVol: 1, bufferBVol: 1,
      ligaseVol: null, ligaseTouched: false,
      enzymeVol: null, enzymeTouched: false,
    };
    state.mmOverage = 1;
    state.mmChecks = {};
    state.assemblies = [ mkAssembly() ];
    state.cycling = defaultCycling();

    document.getElementById('in-accFmol').value = state.accFmol;
    document.getElementById('in-insFmol').value = state.insFmol;
    document.getElementById('in-totalVol').value = state.totalVol;
    document.getElementById('in-mmEnabled').checked = state.mm.enabled;
    document.querySelector('input[name="calcMode"][value="sequence"]').checked = true;
    document.getElementById('in-rampRate').value = state.cycling.rampRate;
    document.getElementById('in-heatup').value = state.cycling.heatupMin;

    renderMasterMixSettings();
    renderCyclingTable();
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
    details: [a.enz, a.lvl!=null&&a.lvl!=='' ? 'Level '+a.lvl : null, a.pos&&a.pos!=='-' ? 'Position '+a.pos : null, a.mark].filter(Boolean).join(' · ')
  }));
  PARTS.forEach(p => ALL_DB.push({
    kind:'part', raw:p, name:p.n, cat:(p.t||'').trim(), len:p.s.length,
    details: `5′ ${p.s5||''}  ·  3′ ${p.s3||''}`
  }));
  LINKERS.forEach(l => ALL_DB.push({
    kind:'linker', raw:l, name:l.n, cat:(l.t||'').trim(), len:l.s.length,
    details: [l.pos?('Positions '+l.pos):null, l.lvl!=null?('Level '+l.lvl):null].filter(Boolean).join(' · ')
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
      if (!state.assemblies.length) state.assemblies.push(mkAssembly());
      const asm = state.assemblies[state.assemblies.length-1];
      if (d.kind === 'acceptor'){
        asm.acceptor = { mode:'library', part:{kind:'acceptor', raw:d.raw, cat:'Acceptor Vector'}, conc:'', customName:'', customSeq:'', customLen:'' };
      } else {
        const ins = mkInsert();
        ins.part = { kind:d.kind, raw:d.raw, cat:d.cat };
        asm.inserts.push(ins);
      }
      asm.collapsed = false;
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
  renderMasterMixSettings();
  renderCyclingTable();
  renderAll();
}

loadData();
