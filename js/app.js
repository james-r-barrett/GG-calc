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

     Any entry in these three lists may omit `s` (sequence) and give `len` (integer, bp)
     instead, e.g. {n:"pXYZ-001", len:2734, t:"Promoter"}. Such parts are fully selectable
     and usable in volume calculations (length-only mode), but their sequence is never
     fetched, displayed, or exposed via "Copy seq" — only the name, category and length.
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
    'My Parts': '#2B6CA3',
    // Chloroplast (CHLOROMODAS) kit -- categories with no nuclear-kit equivalent
    "5' Homology": '#3E7CB1',
    "3' Homology": '#3E7CB1',
    "5' Connector": '#6B8E4E',
    "3' Connector": '#6B8E4E',
    "Operon Connector (5')": '#8FA85E',
    "Operon Connector (3')": '#8FA85E',
    'IEE': '#B0793E',
    'N-tag': '#9B6EC4',
    'C-tag': '#9B6EC4',
    'Selection Marker': '#C08A2E',
    'E. coli ORI + Resistance': '#7C5C3E',
    'Placeholder': '#9CA3AF',
    'CDS': '#B85C38',
  };
  function catColor(cat){ return CAT_COLORS[cat] || '#8A8F98'; }
  function catDot(cat){ return `<span class="cat-dot" style="background:${catColor(cat)}"></span>`; }
  function catPill(cat){ return `<span class="cat-pill">${catDot(cat)}${escapeHtml(cat||'Custom')}</span>`; }

  /* ============ Kit labels (Part Library kit toggle) ============ */
  const KIT_LABELS = { nuclear: 'Nuclear kit', chloroplast: 'Chloroplast kit' };
  function kitLabel(kit){ return KIT_LABELS[kit] || 'Other'; }

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

  /* ============ Combined searchable acceptor list ============ */
  const ACCEPTABLE = [];
  ACCEPT.forEach(a => ACCEPTABLE.push({ kind:'acceptor', raw:a, cat:'Acceptor Vector' }));

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
  /* ============ User-saved parts (name + optional description + length only, no sequence) ============ */
  // Persisted in this browser only (localStorage) so lab members can remember frequently-reused
  // fragments (e.g. a backbone) by name/length without needing to paste or store the actual sequence.
  const USER_PARTS_KEY = 'gg-calc-user-parts-v1';
  function loadUserParts(){
    try {
      const arr = JSON.parse(localStorage.getItem(USER_PARTS_KEY) || '[]');
      return Array.isArray(arr) ? arr : [];
    } catch(e){ return []; }
  }
  function saveUserParts(){
    try { localStorage.setItem(USER_PARTS_KEY, JSON.stringify(state.userParts)); } catch(e){}
  }

  const state = {
    accFmol: 25, accFmolTouched: false,   // default scales with insert count -- see fmolCurrent()
    insFmol: 50, insFmolTouched: false,
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
    userParts: loadUserParts(),  // [{id, n:name, note:description, len:bp}], no sequence
  };
  // uidCounter starts fresh at 1 each page load, but userParts persist across sessions with ids
  // already minted from a previous counter -- without this, a newly-added part can be assigned an
  // id ("up1", "up2", ...) that collides with an existing saved part, so id-keyed lookups
  // (Use/Delete/role change) silently act on the wrong (first-matching) entry.
  state.userParts.forEach(up => {
    const m = /^up(\d+)$/.exec(up.id || '');
    if (m) uidCounter = Math.max(uidCounter, parseInt(m[1], 10) + 1);
  });
  // Repair any parts already saved with colliding ids from that bug (keep the first, renumber the rest).
  {
    const seenIds = new Set();
    let repaired = false;
    state.userParts.forEach(up => {
      if (seenIds.has(up.id)){ up.id = 'up'+(uidCounter++); repaired = true; }
      seenIds.add(up.id);
    });
    if (repaired) saveUserParts();
  }
  // Recommended digestion temperatures per Type IIS enzyme (NEB Golden Gate guidance).
  // Ligation always runs at 16°C (T4 ligase) regardless of the restriction enzyme, so it isn't listed here.
  // The final digestion step (60°C, 5 min) and hold (4°C) are the same for every enzyme.
  const CYCLING_ENZYMES = [
    { key:'bsai',   label:'BsaI-HF / BsaI-HFv2',         digest:37 },
    { key:'bbsi',   label:'BbsI-HF (BpiI)',              digest:37 },
    { key:'sapi',   label:'SapI / BspQI-HF / PaqCI',     digest:37 },
    { key:'bspqi',  label:'BspQI (non-HF)',              digest:42 },
    { key:'bsmbi',  label:'BsmBI-v2 (Esp3I)',            digest:42 },
    { key:'custom', label:'Custom / other',              digest:null },
  ];
  function enzymeKeyFromName(name){
    const n = (name||'').toLowerCase();
    if (n.includes('bsai')) return 'bsai';
    if (n.includes('bbsi') || n.includes('bpii')) return 'bbsi';
    if (n.includes('sapi') || n.includes('paqci')) return 'sapi';
    if (n.includes('bspqi')) return 'bspqi';
    if (n.includes('bsmbi') || n.includes('esp3i')) return 'bsmbi';
    return null;
  }
  function applyCyclingEnzyme(key){
    const e = CYCLING_ENZYMES.find(x => x.key === key);
    if (!e || e.digest == null) return;
    state.cycling.blocks[0].steps[0].temp = e.digest;      // Digestion
  }
  function defaultCycling(){
    return {
      rampRate: 3, heatupMin: 5,   // fixed rough estimates, not user-editable -- stated in the blurb instead
      enzyme: 'bsai', enzymeTouched: false,   // enzymeTouched: user picked manually, stop auto-syncing to the selected acceptor
      timeBudgetHours: '1', timeBudgetMinutes: '25',   // pre-filled example time -- auto-optimised against on first open, remembered thereafter if the user edits it
      appliedSummary: '',
      timeWarning: null,   // set when a time-optimise fit falls below the standard 30x1min protocol; shown in its own box, not inline in appliedSummary
      summaryKind: 'recommended',   // 'recommended'/'overnight'/'quick': appliedSummary text is regenerated from this + the current enzyme whenever the enzyme changes; 'time' summaries are custom-worded and left alone
      activeMode: 'recommended',   // which of the 3 always-present buttons (recommended/overnight/time) is highlighted + governs whether the time-fill dialog shows
      protocolMode: 'auto', autoFragCount: null,   // 'auto': keep re-applying the recommended tier as fragment count changes; 'manual': user picked time-fill/overnight/single-stage/hand-edited a step, stop tracking
      blocks: [
        { cycles: 30, steps: [ { name:'Digestion', temp:37, time:1 }, { name:'Ligation', temp:16, time:1 } ] },
        { cycles: 1,  steps: [ { name:'Final digestion', temp:60, time:5 } ] },
        { cycles: 1,  steps: [ { name:'Hold', temp:4, time:0 } ] },
      ],
    };
  }

  // Fragment count for an assembly = 1 acceptor/backbone + N inserts.
  function assemblyFragmentCount(asm){ return 1 + asm.inserts.length; }
  function maxFragmentCount(){
    return state.assemblies.reduce((m, a) => Math.max(m, assemblyFragmentCount(a)), 1);
  }
  // Complexity tier by insert count, used to scale default fmol targets and enzyme/ligase volumes:
  // a single insert needs less than a multi-part assembly, and a large (7+) assembly needs more still.
  function insertTier(n){
    if (n <= 1) return 0;
    if (n <= 6) return 1;
    return 2;
  }
  function maxInsertCount(){
    return state.assemblies.reduce((m, a) => Math.max(m, a.inserts.length), 1);
  }
  const FMOL_TIER_DEFAULTS = [ { acc:10, ins:20 }, { acc:20, ins:40 }, { acc:40, ins:80 } ];
  const MM_TIER_DEFAULTS = [ { enzyme:0.3, ligase:0.5 }, { enzyme:0.6, ligase:1 }, { enzyme:1.2, ligase:2 } ];
  function fmolDefaults(){ return FMOL_TIER_DEFAULTS[insertTier(maxInsertCount())]; }
  function fmolCurrent(){
    const d = fmolDefaults();
    return {
      acc: state.accFmolTouched ? state.accFmol : d.acc,
      ins: state.insFmolTouched ? state.insFmol : d.ins,
    };
  }
  // Default cycling scheme by assembly complexity, on the same insert-count tiers as the fmol and
  // enzyme/ligase defaults (see insertTier()): a single insert needs less time than a multi-part
  // assembly, and a large (7+) assembly needs longer steps still.
  function recommendedScheme(insertCount){
    const t = insertTier(insertCount);
    if (t === 0) return { cycles:30, stepMin:1, tier:'1 insert' };
    if (t === 1) return { cycles:30, stepMin:1, tier:'2–6 inserts' };
    return { cycles:30, stepMin:5, tier:'7+ inserts' };
  }
  // The minimal viable protocol per NEB's Golden Gate Assembly Kit manual, offered only as a
  // time-crunch fallback (not the default) when even 30x1min doesn't fit the available time. Fast
  // (37°C) enzymes can skip cycling entirely with a single incubation; slow (42°C) enzymes still
  // need at least 15 cycles.
  function quickScheme(digestTemp){
    return digestTemp === 42
      ? { mode:'cycle', cycles:15, stepMin:1 }
      : { mode:'single', time:15 };
  }
  function applyQuickScheme(scheme, digestTemp){
    if (scheme.mode === 'single') applySingleStageBlocks(digestTemp, scheme.time);
    else applySchemeToBlocks(scheme.cycles, scheme.stepMin);
  }
  // Sets the main cycling block to (cycles, stepMin) and pins the final digestion (60°C, 5 min)
  // and hold (4°C) steps, regardless of enzyme or fragment count.
  function applySchemeToBlocks(cycles, stepMin){
    const b0 = state.cycling.blocks[0];
    const digestTemp = b0.steps[0] ? b0.steps[0].temp : 37;
    b0.cycles = cycles;
    // Rebuilt from scratch (rather than assuming b0.steps[0]/[1] exist) so this also recovers
    // cleanly from the single-stage layout, which collapses this block to one step.
    b0.steps = [ { name:'Digestion', temp: digestTemp, time: stepMin }, { name:'Ligation', temp:16, time: stepMin } ];
    const b1 = state.cycling.blocks[1];
    b1.cycles = 1; b1.steps[0].temp = 60; b1.steps[0].time = 5;
    const b2 = state.cycling.blocks[2];
    b2.cycles = 1; b2.steps[0].temp = 4; b2.steps[0].time = 0;
  }
  // Alternative to cycling for simple, time-pressed assemblies: one combined digestion/ligation
  // incubation instead of repeated cycles, followed by the same fixed final digestion and hold.
  function applySingleStageBlocks(digestTemp, timeMin){
    const b0 = state.cycling.blocks[0];
    b0.cycles = 1;
    b0.steps = [ { name:'Digestion + Ligation', temp: digestTemp, time: timeMin || 15 } ];
    const b1 = state.cycling.blocks[1];
    b1.cycles = 1; b1.steps[0].temp = 60; b1.steps[0].time = 5;
    const b2 = state.cycling.blocks[2];
    b2.cycles = 1; b2.steps[0].temp = 4; b2.steps[0].time = 0;
  }
  // Read-only total-time estimate for a candidate (cycles, stepMin) main-cycling block, followed by
  // the fixed final digestion (60°C/5min) and hold (4°C) tail. Mirrors cyclingTotalMinutes()'s ramp math.
  function estimateTotalMinutes(cycles, stepMin, digestTemp){
    const rampPerMin = (state.cycling.rampRate || 0) * 60;
    let total = state.cycling.heatupMin || 0;
    let currentTemp = null;
    const advance = (temp, time) => {
      if (currentTemp != null && rampPerMin > 0) total += Math.abs(temp - currentTemp) / rampPerMin;
      total += time;
      currentTemp = temp;
    };
    for (let c = 0; c < cycles; c++){
      advance(digestTemp, stepMin);
      advance(16, stepMin);
    }
    advance(60, 5);
    advance(4, 0);
    return total;
  }
  // Fill-the-time search. Anchored on 30 cycles x 5 min -- NEB's own baseline for a standard
  // assembly -- rather than treating cycles and step time as equally free. The cycle cap depends on
  // fragment count: 30x5min is the optimal ceiling for the <20-fragment tiers (see recommendedScheme),
  // so extra time beyond that isn't spent on more cycles unless the assembly is complex enough
  // (>20 fragments) to benefit, per the tested pilot data (see the efficiency/accuracy chart) up to 60.
  //
  // Per cycle, only the step time varies with the fit -- ramp time between steps depends only on the
  // fixed temperatures, not on how long each step runs -- so total time is exactly linear in step
  // time for a given cycle count: estimateTotalMinutes(C, s) = overhead(C) + 2*C*s. That lets step
  // time be solved for exactly (down to the second) instead of only tried at whole-minute increments.
  const FILL_TIME_ANCHOR_CYCLES = 30, FILL_TIME_ANCHOR_STEPMIN = 5;
  function fillTimeCycleCap(fragCount){ return fragCount > 20 ? 60 : FILL_TIME_ANCHOR_CYCLES; }
  function suggestFillTime(availableMinutes, digestTemp, fragCount){
    const cap = fillTimeCycleCap(fragCount);
    const overhead = cycles => estimateTotalMinutes(cycles, 0, digestTemp);
    const timeAt = (cycles, stepMin) => overhead(cycles) + 2 * cycles * stepMin;
    const anchorMinutes = timeAt(FILL_TIME_ANCHOR_CYCLES, FILL_TIME_ANCHOR_STEPMIN);
    const minMinutes = timeAt(FILL_TIME_ANCHOR_CYCLES, 1);

    if (availableMinutes >= anchorMinutes){
      // At or above the 30x5min anchor: grow cycles (5 min each) up to `cap` -- a no-op for the
      // <20-fragment tiers, where cap is already 30, so the result just stays at the 30x5min ceiling.
      let cycles = FILL_TIME_ANCHOR_CYCLES;
      while (cycles < cap && timeAt(cycles + 1, 5) <= availableMinutes) cycles++;
      return { cycles, stepMin: 5, minutes: timeAt(cycles, 5), shortfall: false, belowMinimum: false };
    }

    if (availableMinutes >= minMinutes){
      // Below the anchor but at/above the 30x1min floor: cycles stay fixed at 30 and step time is
      // solved exactly to use up whatever time is left.
      const stepMin = Math.min(5, Math.max(1, (availableMinutes - overhead(FILL_TIME_ANCHOR_CYCLES)) / (2 * FILL_TIME_ANCHOR_CYCLES)));
      return { cycles: FILL_TIME_ANCHOR_CYCLES, stepMin, minutes: timeAt(FILL_TIME_ANCHOR_CYCLES, stepMin), shortfall: false, belowMinimum: false };
    }

    // Even 30 cycles at 1 min doesn't fit -- below the minimum sane protocol. Find the most cycles
    // that fit at the 1 min floor, then spend any leftover slack on step time (still exact), rather
    // than presenting this as a normal fit.
    let cycles = 1;
    while (cycles < FILL_TIME_ANCHOR_CYCLES && timeAt(cycles + 1, 1) <= availableMinutes) cycles++;
    if (timeAt(cycles, 1) > availableMinutes){
      return { cycles, stepMin: 1, minutes: timeAt(cycles, 1), shortfall: true, belowMinimum: true };
    }
    const stepMin = Math.min(5, (availableMinutes - overhead(cycles)) / (2 * cycles));
    return { cycles, stepMin, minutes: timeAt(cycles, stepMin), shortfall: false, belowMinimum: true };
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
  // Length of a library/user part regardless of whether the full sequence is known.
  function partLen(raw){
    if (raw.s) return raw.s.length;
    return Number.isFinite(raw.len) ? raw.len : null;
  }
  // Saved parts (My Parts) only ever record a length, never a sequence, so they can't use the
  // exact sequence-composition method -- force length-only for them regardless of any setting.
  function rowHasSeq(row){
    return row.mode !== 'library' || !!(row.part && row.part.raw.s);
  }
  // A row can override the global sequence/length calculation mode; falls back to the global setting.
  function rowCalcMode(row){
    if (!rowHasSeq(row)) return 'length';
    return row.calcModeOverride || state.calcMode;
  }
  // Fragment length to *display* next to a row, independent of whether it can be used for a volume calc.
  function displayLen(row){
    if (row.mode === 'library') return row.part ? partLen(row.part.raw) : null;
    if (rowCalcMode(row) === 'length') return rowLenBp(row);
    const seq = rowSeq(row);
    return seq ? seq.length : null;
  }

  /* ============ Master mix defaults ============ */
  function mmDefaults(){
    const t = MM_TIER_DEFAULTS[insertTier(maxInsertCount())];
    return { buffer: state.totalVol/10, ligase: t.ligase, enzyme: t.enzyme };
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
    const hasSeq = rowHasSeq(row);
    const lengthMode = rowCalcMode(row) === 'length';

    let html = `<div class="part-row-head">
        <span class="part-row-label">${opts.label}</span>
        <div class="row-actions">
          ${hasSeq ? `<select class="calc-mode-select" title="Volume calculation method for this fragment">
            <option value=""${!row.calcModeOverride?' selected':''}>${state.calcMode==='length'?'Length':'Sequence'} (default)</option>
            <option value="sequence"${row.calcModeOverride==='sequence'?' selected':''}>Sequence composition</option>
            <option value="length"${row.calcModeOverride==='length'?' selected':''}>Length only</option>
          </select>` : `<span class="calc-mode-forced" title="Saved parts only record a length, so the exact sequence-composition method isn't available">Length only</span>`}
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
          <input type="text" class="combo-input" placeholder="Search parts by name or type&hellip;" value="${escapeHtml(inputVal)}" autocomplete="off" role="combobox" aria-expanded="false" aria-autocomplete="list" aria-label="Search parts by name or type">
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
    const calcModeSel = container.querySelector('.calc-mode-select');
    if (calcModeSel){
      calcModeSel.addEventListener('change', e => {
        row.calcModeOverride = e.target.value || null;
        renderAll();
      });
    }
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
      filtered = !q ? items.slice(0, 60) : items.filter(it =>
        it.raw.n.toLowerCase().includes(q) ||
        (it.cat && it.cat.toLowerCase().includes(q)) ||
        (it.raw.note && it.raw.note.toLowerCase().includes(q)) ||
        (it.raw.desc && it.raw.desc.toLowerCase().includes(q))
      ).slice(0, 60);
      if (filtered.length === 0){
        panel.innerHTML = `<div class="combo-empty">No matches</div>`;
      } else {
        panel.innerHTML = filtered.map((it,i) => `
          <div class="combo-item${i===activeIdx?' active':''}" role="option" data-idx="${i}">
            ${catDot(it.cat)}
            <span class="ci-name">${escapeHtml(it.raw.n)}</span>
            <span class="ci-len mono">${partLen(it.raw)} bp</span>
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
      len = row.mode === 'custom' ? rowLenBp(row) : (row.part ? partLen(row.part.raw) : null);
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
    const fm = fmolCurrent();
    const accRow = buildRow(asm.acceptor, fm.acc, 'Acceptor');
    accRow.cat = 'Acceptor Vector';
    accRow.key = 'acceptor';
    rows.push(accRow);
    asm.inserts.forEach((ins, i) => {
      const r = buildRow(ins, fm.ins, 'Insert ' + (i+1));
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
      renderPartRow(accWrap, asm.acceptor, { label:'Acceptor', items: ACCEPTABLE, onRemove:null });

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
    const first = state.assemblies[0];
    const enz = (first && first.acceptor.mode === 'library' && first.acceptor.part) ? first.acceptor.part.raw.enz : null;

    if (enz && !state.cycling.enzymeTouched){
      const key = enzymeKeyFromName(enz);
      if (key && key !== state.cycling.enzyme){
        state.cycling.enzyme = key;
        applyCyclingEnzyme(key);
        renderCyclingTable();
      }
    }

    renderCyclingSuggest();
  }

  /* ============ Cycling protocol suggestion (fragment-count and time-budget based) ============ */
  // Formats a decimal-minute duration to the nearest second, e.g. 1.2833 -> "1m 17s" -- used wherever
  // a time-optimised step time isn't a round number of minutes.
  function fmtMinSec(mins){
    const totalSec = Math.round((mins || 0) * 60);
    const m = Math.floor(totalSec / 60), s = totalSec % 60;
    return s > 0 ? `${m}m ${s}s` : `${m} min`;
  }
  // Compact form for the narrow table cell, e.g. 1.2833 -> "1m17s", 5 -> "5m".
  function fmtMinSecShort(mins){
    const totalSec = Math.round((mins || 0) * 60);
    const m = Math.floor(totalSec / 60), s = totalSec % 60;
    return s > 0 ? `${m}m${s}s` : `${m}m`;
  }
  // Parses "1m 17s", "1m17s", "1:17", "17s", or a plain number (taken as minutes) back to decimal
  // minutes -- the inverse of fmtMinSec/fmtMinSecShort, so the table's time cells can round-trip
  // without forcing users to type decimals.
  function parseMinSec(str){
    const s = String(str || '').trim().toLowerCase();
    if (!s) return 0;
    let m = s.match(/^(\d+(?:\.\d+)?)\s*mi?n?s?\s*(\d+(?:\.\d+)?)?\s*s(?:ec)?s?$/);
    if (m) return parseFloat(m[1]) + (m[2] ? parseFloat(m[2]) / 60 : 0);
    m = s.match(/^(\d+(?:\.\d+)?):(\d+(?:\.\d+)?)$/);
    if (m) return parseFloat(m[1]) + parseFloat(m[2]) / 60;
    m = s.match(/^(\d+(?:\.\d+)?)\s*s(?:ec)?s?$/);
    if (m) return parseFloat(m[1]) / 60;
    m = s.match(/^(\d+(?:\.\d+)?)\s*mi?n?s?$/);
    if (m) return parseFloat(m[1]);
    const n = parseFloat(s);
    return isNaN(n) ? 0 : n;
  }
  function cyclingSchemeText(cycles, stepMin, digestTemp){
    return `${cycles} cycles of ${fmtSmart(digestTemp,1)}&deg;C for ${fmtMinSec(stepMin)} and 16&deg;C for ${fmtMinSec(stepMin)}, then 60&deg;C for 5 min, then hold at 4&deg;C indefinitely.`;
  }
  function singleStageSchemeText(digestTemp, timeMin){
    return `${fmtSmart(digestTemp,1)}&deg;C for ${fmtMinSec(timeMin)}, then 60&deg;C for 5 min, then hold at 4&deg;C indefinitely.`;
  }
  // Text for a quickScheme() result -- the only place that still deals in the single-stage/cycle
  // "mode" distinction, since recommendedScheme() itself only ever returns a cycling scheme now.
  function quickSchemeText(scheme, digestTemp){
    return scheme.mode === 'single' ? singleStageSchemeText(digestTemp, scheme.time) : cyclingSchemeText(scheme.cycles, scheme.stepMin, digestTemp);
  }
  function cyclingEnzymeLabel(){
    const e = CYCLING_ENZYMES.find(x => x.key === state.cycling.enzyme);
    return e ? e.label : 'the selected enzyme';
  }
  // Regenerates the applied-summary text for the "static" kinds (recommended/overnight/quick), all of
  // which name the current enzyme -- called whenever the enzyme selection changes so the text at the
  // top of the box stays in sync. Time-fill summaries are custom-worded per result and left alone.
  function buildAppliedSummaryText(kind){
    const fragCount = maxFragmentCount();
    if (kind === 'recommended') return `Recommended for ${fragCount} fragment${fragCount===1?'':'s'} with ${cyclingEnzymeLabel()}.`;
    if (kind === 'overnight') return `Overnight for ${fragCount} fragment${fragCount===1?'':'s'} with ${cyclingEnzymeLabel()}.`;
    if (kind === 'quick') return `Quick protocol for ${cyclingEnzymeLabel()}.`;
    return null;
  }
  function applyAutoRecommended(){
    const fragCount = maxFragmentCount();
    const scheme = recommendedScheme(maxInsertCount());
    applySchemeToBlocks(scheme.cycles, scheme.stepMin);
    state.cycling.autoFragCount = fragCount;
    state.cycling.timeWarning = null;
    state.cycling.summaryKind = 'recommended';
    state.cycling.appliedSummary = buildAppliedSummaryText('recommended');
  }
  function renderCyclingSuggest(){
    const el = document.getElementById('cycling-suggest');
    if (!el) return;

    // Auto mode keeps the recommended protocol applied and in sync as the fragment count changes,
    // without the user needing to click anything -- mirrors the enzyme auto-sync elsewhere.
    if (state.cycling.protocolMode === 'auto' && state.cycling.autoFragCount !== maxFragmentCount()){
      applyAutoRecommended();
      renderCyclingTable();
    }
    const digestTemp = state.cycling.blocks[0].steps[0].temp;
    const fragCount = maxFragmentCount();
    const activeMode = state.cycling.activeMode;

    // The applied protocol and its action buttons stay visible at all times now -- pressing Time-
    // optimise/Overnight/single-stage just updates the applied text in place rather than swapping to
    // a separate panel that has to be reopened. All 3 buttons are always present, with the active one
    // highlighted so it's clear which protocol is currently applied.
    let html = `<div class="cyc-suggest-box cyc-applied-box">
      <p class="cyc-suggest-line">${state.cycling.appliedSummary}</p>
      <div class="cyc-suggest-actions">
        <button type="button" class="btn-suggest btn-suggest-alt${activeMode==='recommended'?' btn-suggest-active':''}" id="cyc-apply-recommended">Recommended</button>
        <button type="button" class="btn-suggest btn-suggest-alt${activeMode==='overnight'?' btn-suggest-active':''}" id="cyc-suggest-overnight">Overnight</button>
        <button type="button" class="btn-suggest btn-suggest-alt${activeMode==='time'?' btn-suggest-active':''}" id="cyc-toggle-time">Time-optimise</button>
      </div>`;

    if (activeMode === 'time'){
      html += `<div class="cyc-time-box">
        <p class="field-label">Fill the time available before transformation</p>
        <div class="cyc-time-row">
          <label class="mm-field-inline">Hours
            <input type="text" inputmode="decimal" id="cyc-time-hours" value="${escapeHtml(state.cycling.timeBudgetHours||'')}" placeholder="0">
          </label>
          <label class="mm-field-inline">Minutes
            <input type="text" inputmode="decimal" id="cyc-time-minutes" value="${escapeHtml(state.cycling.timeBudgetMinutes||'')}" placeholder="0">
          </label>
          <button type="button" class="btn-suggest" id="cyc-suggest-time">Optimise</button>
        </div>
        <p id="cyc-time-result" class="note" hidden></p>
        ${state.cycling.timeWarning ? `<div class="warning cyc-time-warning">
          <p>${state.cycling.timeWarning}</p>
          <button type="button" class="btn-suggest btn-suggest-alt cyc-time-warning-btn${state.cycling.summaryKind==='quick'?' btn-suggest-active':''}" id="cyc-apply-quick">${state.cycling.summaryKind==='quick'?'Quick protocol applied':'Use quick protocol'}</button>
        </div>` : ''}
      </div>`;
    }

    html += `</div>`;

    el.innerHTML = html;

    // Reads the current (pre-filled or user-edited) time budget, computes the fill-time suggestion,
    // and applies it -- shared by the Time-optimise toggle (so opening/returning to the panel
    // auto-calculates against whatever time is already entered) and the Optimise button itself.
    function runTimeOptimise(){
      const hours = parseFloat(state.cycling.timeBudgetHours) || 0;
      const minutes = parseFloat(state.cycling.timeBudgetMinutes) || 0;
      const totalMinutes = hours * 60 + minutes;
      if (totalMinutes <= 0) return false;
      const suggestion = suggestFillTime(totalMinutes, digestTemp, fragCount);
      applySchemeToBlocks(suggestion.cycles, suggestion.stepMin);
      state.cycling.protocolMode = 'manual';
      state.cycling.activeMode = 'time';
      state.cycling.summaryKind = 'time';
      const quick = quickScheme(digestTemp);
      state.cycling.appliedSummary = `Time-optimised for ${fmtDuration(totalMinutes)} available (~${fmtDuration(suggestion.minutes)} total).`;
      state.cycling.timeWarning = suggestion.belowMinimum
        ? `Not enough time for the standard 30-cycle protocol (30 &times; 1 min needs ~${fmtDuration(estimateTotalMinutes(30,1,digestTemp))}) &mdash; scaled down below (~${fmtDuration(suggestion.minutes)}, doesn't fit ${fmtDuration(totalMinutes)}). Efficiency and accuracy will likely be reduced &mdash; see the efficiency/accuracy graph below. Consider the quicker protocol instead: ${quickSchemeText(quick, digestTemp)}`
        : null;
      return true;
    }

    el.querySelector('#cyc-toggle-time').addEventListener('click', () => {
      const wasAlreadyTime = state.cycling.activeMode === 'time';
      state.cycling.activeMode = 'time';
      if (!wasAlreadyTime) runTimeOptimise();
      renderCyclingSuggest();
      renderCyclingTable();
    });

    el.querySelector('#cyc-apply-recommended').addEventListener('click', () => {
      state.cycling.protocolMode = 'auto';
      state.cycling.activeMode = 'recommended';
      applyAutoRecommended();
      renderCyclingSuggest();
      renderCyclingTable();
    });

    const quickBtn = el.querySelector('#cyc-apply-quick');
    if (quickBtn){
      quickBtn.addEventListener('click', () => {
        const quick = quickScheme(digestTemp);
        applyQuickScheme(quick, digestTemp);
        state.cycling.protocolMode = 'manual';
        state.cycling.summaryKind = 'quick';
        state.cycling.appliedSummary = buildAppliedSummaryText('quick');
        renderCyclingSuggest();
        renderCyclingTable();
      });
    }

    el.querySelector('#cyc-suggest-overnight').addEventListener('click', () => {
      applySchemeToBlocks(40, 5);
      state.cycling.protocolMode = 'manual';
      state.cycling.activeMode = 'overnight';
      state.cycling.summaryKind = 'overnight';
      state.cycling.timeWarning = null;
      state.cycling.appliedSummary = buildAppliedSummaryText('overnight');
      renderCyclingSuggest();
      renderCyclingTable();
    });

    if (activeMode === 'time'){
      el.querySelector('#cyc-time-hours').addEventListener('input', e => { state.cycling.timeBudgetHours = e.target.value; });
      el.querySelector('#cyc-time-minutes').addEventListener('input', e => { state.cycling.timeBudgetMinutes = e.target.value; });
      el.querySelector('#cyc-suggest-time').addEventListener('click', () => {
        const resultEl = el.querySelector('#cyc-time-result');
        if (!runTimeOptimise()){
          resultEl.hidden = false;
          resultEl.textContent = 'Enter the time you have available.';
          return;
        }
        resultEl.hidden = true;
        renderCyclingSuggest();
        renderCyclingTable();
      });
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
  function renderCyclingEnzymeSelect(){
    const el = document.getElementById('cycling-enzyme-select');
    if (!el) return;
    el.innerHTML = `<label class="cycling-enzyme-field">
      <span>Restriction enzyme</span>
      <select id="in-cyclingEnzyme">
        ${CYCLING_ENZYMES.map(e => `<option value="${e.key}"${state.cycling.enzyme===e.key?' selected':''}>${escapeHtml(e.label)}</option>`).join('')}
      </select>
    </label>`;
    el.querySelector('#in-cyclingEnzyme').addEventListener('change', e => {
      state.cycling.enzyme = e.target.value;
      state.cycling.enzymeTouched = true;
      applyCyclingEnzyme(state.cycling.enzyme);
      // The applied-summary text names the enzyme, so keep it in sync for the "static" kinds.
      const regenerated = buildAppliedSummaryText(state.cycling.summaryKind);
      if (regenerated) state.cycling.appliedSummary = regenerated;
      renderCyclingTable();
      renderCyclingSuggest();
    });
  }
  function renderCyclingTable(){
    const wrap = document.getElementById('cycling-table-wrap');
    if (!wrap) return;
    renderCyclingEnzymeSelect();
    let html = `<table class="cycling-table">
      <colgroup><col><col class="col-temp"><col class="col-time"><col class="col-cyc"></colgroup>
      <thead><tr><th>Step</th><th>Temp</th><th>Time</th><th>&times;</th></tr></thead><tbody>`;
    // The hold block (always the last one, 4°C) isn't user-editable -- it's hardcoded to run
    // indefinitely rather than for a set number of minutes, so it's shown as a fixed row below
    // instead of with editable inputs.
    const cyclableBlocks = state.cycling.blocks.filter(block => block.steps[0].name !== 'Hold');
    cyclableBlocks.forEach((block) => {
      const bi = state.cycling.blocks.indexOf(block);
      block.steps.forEach((step, si) => {
        // Flag a hot (42°C) digestion temp in red -- it's the one enzyme-dependent difference in the
        // schedule worth calling out at a glance (e.g. BspQI/BsmBI-v2 vs the more common 37°C enzymes).
        const isHotDigestion = step.name === 'Digestion' && step.temp >= 40;
        html += `<tr>
          <td>${escapeHtml(step.name)}</td>
          <td class="cyc-temp-static${isHotDigestion ? ' cyc-temp-hot' : ''}">${fmtSmart(step.temp,1)}&deg;</td>
          <td><input type="text" class="cyc-input cyc-time" data-block="${bi}" data-step="${si}" value="${fmtMinSecShort(step.time)}"></td>
          ${si===0 ? (block.steps.length>1
            ? `<td rowspan="${block.steps.length}"><input type="text" inputmode="numeric" class="cyc-input cyc-cycles" data-block="${bi}" value="${block.cycles}"></td>`
            : `<td class="cyc-temp-static">${block.cycles}</td>`) : ''}
        </tr>`;
      });
    });
    html += `<tr class="cyc-hold-row"><td>Hold <span class="cyc-hold-flag">see warning</span></td><td class="cyc-temp-static">4&deg;</td><td colspan="2" class="cyc-hold-indefinite">Indefinite</td></tr>`;
    html += `</tbody></table>`;
    wrap.innerHTML = html;

    const holdWarningEl = document.getElementById('cycling-hold-warning');
    if (holdWarningEl){
      holdWarningEl.innerHTML = `<p class="warning cyc-hold-note">Held at 4&deg;C until you're ready to transform. If held for a long period (e.g. overnight), <strong>repeat the final 60&deg;C for 5 min digestion step</strong> immediately before transformation to fully digest any fragments re-ligated during the hold.</p>`;
    }

    wrap.querySelectorAll('.cyc-time').forEach(inp => {
      inp.addEventListener('input', e => {
        state.cycling.blocks[+e.target.dataset.block].steps[+e.target.dataset.step].time = parseMinSec(e.target.value);
        state.cycling.protocolMode = 'manual';   // stop auto-syncing to the tier once the user hand-edits a step
        updateCyclingTotal();
      });
      // Reformat to the compact "1m17s" form once the user's done typing, rather than on every keystroke.
      inp.addEventListener('blur', e => {
        e.target.value = fmtMinSecShort(state.cycling.blocks[+e.target.dataset.block].steps[+e.target.dataset.step].time);
      });
    });
    wrap.querySelectorAll('.cyc-cycles').forEach(inp => inp.addEventListener('input', e => {
      state.cycling.blocks[+e.target.dataset.block].cycles = Math.max(1, parseInt(e.target.value, 10) || 1);
      state.cycling.protocolMode = 'manual';
      updateCyclingTotal();
    }));

    updateCyclingTotal();
  }

  function renderAll(){
    syncSetupFields();
    renderMasterMixSettings();
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
  function mmFieldRow(key, label, value, showReset, defaultHint, unitNote){
    return `<label class="mm-field">
      <span>${label} <span class="unit">&micro;L${unitNote ? `, ${unitNote}` : ''}</span></span>
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
    const d = mmDefaults();
    html += mmFieldRow('ligase', 'T4 DNA ligase', cur.ligase, state.mm.ligaseTouched, roundStr(d.ligase,3), 'assuming 400,000 U/mL');
    html += mmFieldRow('enzyme', 'Restriction enzyme', cur.enzyme, state.mm.enzymeTouched, roundStr(d.enzyme,3));
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
  // Acceptor/insert fmol targets default from the insert-count tier (see fmolCurrent()) until the
  // user edits them directly; syncSetupFields() keeps the displayed value/reset-button in step with
  // the tier as assemblies are added, removed, or resized.
  function syncSetupFields(){
    const d = fmolDefaults();
    const accInput = document.getElementById('in-accFmol');
    const insInput = document.getElementById('in-insFmol');
    const accReset = document.getElementById('reset-accFmol');
    const insReset = document.getElementById('reset-insFmol');
    if (!state.accFmolTouched) accInput.value = d.acc;
    if (!state.insFmolTouched) insInput.value = d.ins;
    accReset.hidden = !state.accFmolTouched;
    insReset.hidden = !state.insFmolTouched;
    accReset.textContent = `reset (${d.acc})`;
    insReset.textContent = `reset (${d.ins})`;
  }
  document.getElementById('in-accFmol').addEventListener('input', e => { state.accFmol = parseFloat(e.target.value)||0; state.accFmolTouched = true; syncSetupFields(); refreshAllResults(); });
  document.getElementById('in-insFmol').addEventListener('input', e => { state.insFmol = parseFloat(e.target.value)||0; state.insFmolTouched = true; syncSetupFields(); refreshAllResults(); });
  document.getElementById('reset-accFmol').addEventListener('click', () => { state.accFmolTouched = false; syncSetupFields(); refreshAllResults(); });
  document.getElementById('reset-insFmol').addEventListener('click', () => { state.insFmolTouched = false; syncSetupFields(); refreshAllResults(); });
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

  document.getElementById('clear-page-btn').addEventListener('click', () => {
    if (!confirm('Clear everything and start a new page? This cannot be undone.')) return;
    state.accFmol = 25; state.accFmolTouched = false;
    state.insFmol = 50; state.insFmolTouched = false;
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

    document.getElementById('in-totalVol').value = state.totalVol;
    document.getElementById('in-mmEnabled').checked = state.mm.enabled;
    document.querySelector('input[name="calcMode"][value="sequence"]').checked = true;

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
    kind:'acceptor', raw:a, name:a.n, cat:'Acceptor Vector', len:partLen(a), kit:a.kit||'nuclear',
    details: [a.desc, a.enz, a.lvl!=null&&a.lvl!=='' ? 'Level '+a.lvl : null, a.pos&&a.pos!=='-' ? 'Position '+a.pos : null,
      (a.s5||a.s3) ? `5′ ${a.s5||''} · 3′ ${a.s3||''}` : null, a.mark, a.well?('Well '+a.well):null].filter(Boolean).join(' · ')
  }));
  PARTS.forEach(p => ALL_DB.push({
    kind:'part', raw:p, name:p.n, cat:(p.t||'').trim(), len:partLen(p), kit:p.kit||'nuclear',
    details: [(p.s5||p.s3) ? `5′ ${p.s5||''} · 3′ ${p.s3||''}` : null, p.pos?('Position '+p.pos):null, p.well?('Well '+p.well):null]
      .filter(Boolean).join(' · ')
  }));
  LINKERS.forEach(l => ALL_DB.push({
    kind:'linker', raw:l, name:l.n, cat:(l.t||'').trim(), len:partLen(l), kit:l.kit||'nuclear',
    details: [l.pos?('Positions '+l.pos):null, l.lvl!=null?('Level '+l.lvl):null].filter(Boolean).join(' · ')
  }));

  // Keeps INSERTABLE (assembly-row combobox) and ALL_DB (Part Library table) in sync with the
  // user's saved parts; rebuilt wholesale on every add/delete rather than tracked incrementally.
  function syncUserPartsIntoLists(){
    for (let i = INSERTABLE.length - 1; i >= 0; i--) if (INSERTABLE[i].kind === 'user') INSERTABLE.splice(i, 1);
    for (let i = ACCEPTABLE.length - 1; i >= 0; i--) if (ACCEPTABLE[i].kind === 'user') ACCEPTABLE.splice(i, 1);
    for (let i = ALL_DB.length - 1; i >= 0; i--) if (ALL_DB[i].kind === 'user') ALL_DB.splice(i, 1);
    state.userParts.forEach(up => {
      const raw = { id: up.id, n: up.n, note: up.note, len: up.len, role: up.role === 'acceptor' ? 'acceptor' : 'insert' };
      INSERTABLE.push({ kind:'user', raw, cat:'My Parts' });
      ACCEPTABLE.push({ kind:'user', raw, cat:'My Parts' });
      ALL_DB.push({ kind:'user', raw, name: up.n, cat:'My Parts', len: up.len, kit:'user', details: up.note || '' });
    });
  }

  let dbActiveFilter = 'All';
  let dbActiveKit = 'All';

  // Items always shown regardless of the kit toggle (the user's own saved parts
  // aren't part of either built-in kit) live under kit:'user'.
  function matchesKit(d){ return dbActiveKit==='All' || d.kit===dbActiveKit || d.kind==='user'; }

  function renderDbKitFilters(){
    const wrap = document.getElementById('db-kit-filters');
    if (!wrap) return;
    const kits = Array.from(new Set(ALL_DB.filter(d => d.kind!=='user').map(d => d.kit || 'other')));
    const order = ['nuclear', 'chloroplast'];
    kits.sort((a,b) => {
      const ia = order.indexOf(a), ib = order.indexOf(b);
      return (ia===-1?99:ia) - (ib===-1?99:ib);
    });
    const chips = ['All', ...kits];
    wrap.innerHTML = chips.map(k => `<button type="button" class="filter-chip kit-chip${k===dbActiveKit?' active':''}" data-kit="${escapeHtml(k)}">${k==='All'?'All kits':escapeHtml(kitLabel(k))}</button>`).join('');
    wrap.querySelectorAll('.filter-chip').forEach(chip => {
      chip.addEventListener('click', () => {
        dbActiveKit = chip.dataset.kit;
        dbActiveFilter = 'All';
        renderDbKitFilters();
        renderDbFilters();
        renderDbTable();
      });
    });
  }

  function renderDbFilters(){
    const wrap = document.getElementById('db-filters');
    const cats = ['All', ...Array.from(new Set(ALL_DB.filter(matchesKit).map(d => d.cat))).sort()];
    if (!cats.includes(dbActiveFilter)) dbActiveFilter = 'All';
    wrap.innerHTML = cats.map(c => `<button type="button" class="filter-chip${c===dbActiveFilter?' active':''}" data-cat="${escapeHtml(c)}">${c==='All'?'':catDot(c)}${escapeHtml(c)}</button>`).join('');
    wrap.querySelectorAll('.filter-chip').forEach(chip => {
      chip.addEventListener('click', () => { dbActiveFilter = chip.dataset.cat; renderDbTable(); });
    });
  }

  // Sends a library/saved-part entry (d: {kind, raw, cat}) to the last assembly, as the acceptor
  // or as a new insert depending on kind/role, then jumps to the Calculator tab.
  function useLibraryItem(d){
    if (!state.assemblies.length) state.assemblies.push(mkAssembly());
    const asm = state.assemblies[state.assemblies.length-1];
    const useAsAcceptor = d.kind === 'acceptor' || (d.kind === 'user' && d.raw.role === 'acceptor');
    if (useAsAcceptor){
      asm.acceptor = { mode:'library', part:{kind:d.kind, raw:d.raw, cat: d.kind==='user' ? 'My Parts' : 'Acceptor Vector'}, conc:'', customName:'', customSeq:'', customLen:'' };
    } else {
      const ins = mkInsert();
      ins.mode = 'library';
      ins.part = { kind:d.kind, raw:d.raw, cat:d.cat };
      asm.inserts.push(ins);
    }
    asm.collapsed = false;
    renderAll();
    document.querySelector('.tab-btn[data-tab="calc"]').click();
  }

  function renderDbTable(){
    const q = document.getElementById('db-search').value.trim().toLowerCase();
    const rows = ALL_DB.filter(d => matchesKit(d) && (dbActiveFilter==='All' || d.cat===dbActiveFilter) && (!q || d.name.toLowerCase().includes(q) || (d.details && d.details.toLowerCase().includes(q))));
    const tbody = document.getElementById('db-tbody');
    tbody.innerHTML = rows.map((d, i) => `
      <tr>
        <td class="dname">${escapeHtml(d.name)}</td>
        <td>${catPill(d.cat)}</td>
        <td class="dlen mono">${d.len} bp</td>
        <td class="ddetails">${escapeHtml(d.details)}</td>
        <td class="dactions">
          ${d.raw.s ? `<button type="button" class="icon-btn" data-copy="${i}">Copy seq</button>` : ''}
          <button type="button" class="icon-btn" data-use="${i}">Use&nbsp;&rarr;</button>
        </td>
      </tr>`).join('');

    tbody.querySelectorAll('[data-copy]').forEach(b => b.addEventListener('click', () => {
      const d = rows[+b.dataset.copy];
      navigator.clipboard.writeText(d.raw.s).then(() => { b.textContent='Copied'; setTimeout(()=>b.textContent='Copy seq',1200); });
    }));
    tbody.querySelectorAll('[data-use]').forEach(b => b.addEventListener('click', () => {
      useLibraryItem(rows[+b.dataset.use]);
    }));

    document.getElementById('db-count').textContent = `${rows.length} of ${ALL_DB.length} parts`;
  }
  document.getElementById('db-search').addEventListener('input', renderDbTable);

  /* ============ My saved parts (name + optional description + length, no sequence) ============ */
  function renderUserPartsList(){
    const countEl = document.getElementById('user-parts-count');
    if (countEl) countEl.textContent = state.userParts.length ? `(${state.userParts.length})` : '';
    const wrap = document.getElementById('my-parts-list');
    if (!wrap) return;
    if (!state.userParts.length){
      wrap.innerHTML = `<p class="note">No saved parts yet &mdash; add one above.</p>`;
      return;
    }
    wrap.innerHTML = `<table class="db-table my-parts-table">
      <thead><tr><th>Name</th><th>Length</th><th>Description</th><th>Use as</th><th></th></tr></thead>
      <tbody>
        ${state.userParts.map(up => `
          <tr>
            <td class="dname">${escapeHtml(up.n)}</td>
            <td class="dlen mono">${up.len} bp</td>
            <td class="ddetails">${escapeHtml(up.note||'')}</td>
            <td class="drole">
              <select class="role-select" data-role="${escapeHtml(up.id)}">
                <option value="insert"${up.role==='acceptor'?'':' selected'}>Insert</option>
                <option value="acceptor"${up.role==='acceptor'?' selected':''}>Acceptor vector</option>
              </select>
            </td>
            <td class="dactions">
              <button type="button" class="icon-btn" data-use="${escapeHtml(up.id)}">Use&nbsp;&rarr;</button>
              <button type="button" class="icon-btn" data-del="${escapeHtml(up.id)}">Delete</button>
            </td>
          </tr>`).join('')}
      </tbody>
    </table>`;
    wrap.querySelectorAll('[data-del]').forEach(b => b.addEventListener('click', () => {
      state.userParts = state.userParts.filter(up => up.id !== b.dataset.del);
      saveUserParts();
      syncUserPartsIntoLists();
      renderUserPartsList();
      renderDbFilters();
      renderDbTable();
    }));
    wrap.querySelectorAll('[data-role]').forEach(sel => sel.addEventListener('change', () => {
      const up = state.userParts.find(u => u.id === sel.dataset.role);
      if (!up) return;
      up.role = sel.value === 'acceptor' ? 'acceptor' : 'insert';
      saveUserParts();
      syncUserPartsIntoLists();
      renderDbTable();
    }));
    wrap.querySelectorAll('[data-use]').forEach(b => b.addEventListener('click', () => {
      const up = state.userParts.find(u => u.id === b.dataset.use);
      if (!up) return;
      useLibraryItem({ kind:'user', raw:{ id: up.id, n: up.n, note: up.note, len: up.len, role: up.role === 'acceptor' ? 'acceptor' : 'insert' }, cat:'My Parts' });
    }));
  }
  document.getElementById('mp-download-btn').addEventListener('click', () => {
    const payload = JSON.stringify(state.userParts.map(up => ({ n: up.n, note: up.note || '', len: up.len, role: up.role === 'acceptor' ? 'acceptor' : 'insert' })), null, 2);
    const blob = new Blob([payload], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'my-parts.json';
    a.click();
    URL.revokeObjectURL(url);
  });
  document.getElementById('mp-restore-btn').addEventListener('click', () => {
    document.getElementById('mp-restore-input').click();
  });
  document.getElementById('mp-restore-input').addEventListener('change', (e) => {
    const file = e.target.files[0];
    e.target.value = '';
    if (!file) return;
    file.text().then(text => {
      let arr;
      try { arr = JSON.parse(text); } catch(err){ alert('That file is not valid JSON.'); return; }
      if (!Array.isArray(arr)){ alert('Expected a JSON array of parts.'); return; }
      const { added, skipped } = importUserParts(arr.map(entry => ({
        n: entry && entry.n, note: entry && entry.note, len: entry && entry.len, role: entry && entry.role
      })));
      if (skipped) alert(`Imported ${added} part(s). Skipped ${skipped} (invalid or already saved).`);
    });
  });

  // Bulk import from a spreadsheet export: name, description, length (bp) columns, header row optional.
  function parseCsv(text){
    const rows = [];
    let row = [], field = '', inQuotes = false;
    for (let i = 0; i < text.length; i++){
      const c = text[i];
      if (inQuotes){
        if (c === '"'){
          if (text[i+1] === '"'){ field += '"'; i++; } else { inQuotes = false; }
        } else field += c;
      } else if (c === '"'){
        inQuotes = true;
      } else if (c === ','){
        row.push(field); field = '';
      } else if (c === '\n' || c === '\r'){
        if (c === '\r' && text[i+1] === '\n') i++;
        row.push(field); field = '';
        rows.push(row); row = [];
      } else field += c;
    }
    if (field !== '' || row.length){ row.push(field); rows.push(row); }
    return rows.filter(r => r.some(f => f.trim() !== ''));
  }
  function importUserParts(entries){
    const existingNames = new Set(state.userParts.map(up => up.n));
    let added = 0, skipped = 0;
    entries.forEach(entry => {
      const name = (entry && entry.n || '').toString().trim();
      const len = parseInt(entry && entry.len, 10);
      if (!name || !Number.isFinite(len) || len <= 0){ skipped++; return; }
      if (existingNames.has(name)){ skipped++; return; }
      const role = (entry && entry.role || '').toString().trim().toLowerCase() === 'acceptor' ? 'acceptor' : 'insert';
      state.userParts.push({ id:'up'+(uidCounter++), n:name, note:((entry && entry.note)||'').toString().trim(), len, role });
      existingNames.add(name);
      added++;
    });
    saveUserParts();
    syncUserPartsIntoLists();
    renderUserPartsList();
    renderDbFilters();
    renderDbTable();
    return { added, skipped };
  }
  document.getElementById('mp-csv-btn').addEventListener('click', () => {
    document.getElementById('mp-csv-input').click();
  });
  document.getElementById('mp-csv-input').addEventListener('change', (e) => {
    const file = e.target.files[0];
    e.target.value = '';
    if (!file) return;
    file.text().then(text => {
      let rows = parseCsv(text);
      if (!rows.length){ alert('That CSV file has no rows.'); return; }
      // Skip an optional header row (its length column won't parse as a positive integer).
      if (!Number.isFinite(parseInt(rows[0][2], 10)) || parseInt(rows[0][2], 10) <= 0) rows = rows.slice(1);
      const entries = rows.map(r => ({ n: r[0], note: r[1], len: r[2], role: r[3] }));
      const { added, skipped } = importUserParts(entries);
      alert(`Imported ${added} part(s) from CSV.${skipped ? ` Skipped ${skipped} (invalid or already saved).` : ''}`);
    });
  });

  // Parses one or more GenBank flat-file records (split on the "//" record terminator), pulling
  // just LOCUS name/length, DEFINITION, and ORIGIN sequence -- the sequence itself is only used to
  // compute a length and is discarded (saved parts never store sequence, see USER_PARTS_KEY above).
  function parseGenbank(text){
    const records = [];
    text.split(/\r?\n\/\/\s*\r?\n?/).map(s => s.trim()).filter(Boolean).forEach(chunk => {
      const locusMatch = /^LOCUS\s+(\S+)\s+(\d+)\s*bp/im.exec(chunk);
      if (!locusMatch) return;
      const defMatch = /^DEFINITION\s+([^\n]*)/im.exec(chunk);
      let desc = defMatch ? defMatch[1].trim() : '';
      if (desc === '.') desc = '';
      const originIdx = chunk.search(/^ORIGIN/im);
      let seq = '';
      if (originIdx !== -1){
        seq = chunk.slice(originIdx).replace(/^ORIGIN.*$/im, '').replace(/[^a-zA-Z]/g, '').toUpperCase();
      }
      records.push({ name: locusMatch[1], desc, seq, len: seq.length || parseInt(locusMatch[2], 10) });
    });
    return records;
  }
  document.getElementById('mp-gb-btn').addEventListener('click', () => {
    document.getElementById('mp-gb-input').click();
  });
  document.getElementById('mp-gb-input').addEventListener('change', (e) => {
    const files = Array.from(e.target.files);
    e.target.value = '';
    if (!files.length) return;
    Promise.all(files.map(f => f.text())).then(texts => {
      const records = texts.flatMap(parseGenbank);
      if (!records.length){ alert('No LOCUS/ORIGIN records found in that file.'); return; }
      // A single file holding exactly one record: fill the form for review before saving.
      // Multiple files and/or multi-record files: nothing to review against, so import directly.
      if (files.length === 1 && records.length === 1){
        const r = records[0];
        document.getElementById('mp-name').value = r.name || '';
        document.getElementById('mp-desc').value = r.desc || '';
        const seqEl = document.getElementById('mp-seq');
        if (r.seq){
          seqEl.value = r.seq;
          seqEl.dispatchEvent(new Event('input'));
        } else {
          document.getElementById('mp-len').value = r.len || '';
        }
        document.getElementById('mp-name').focus();
      } else {
        const entries = records.map(r => ({ n: r.name, note: r.desc, len: r.len }));
        const { added, skipped } = importUserParts(entries);
        alert(`Imported ${added} part(s) from GenBank.${skipped ? ` Skipped ${skipped} (invalid or already saved).` : ''}`);
      }
    });
  });

  // Paste a raw sequence to live-calculate its length (ignores FASTA header lines & whitespace) into the length field.
  document.getElementById('mp-seq').addEventListener('input', (e) => {
    const resultEl = document.getElementById('mp-seq-result');
    const len = e.target.value.split(/\r?\n/).filter(line => !line.startsWith('>')).join('').replace(/\s/g, '').length;
    if (!len){ resultEl.textContent = ''; return; }
    document.getElementById('mp-len').value = len;
    resultEl.textContent = `= ${len} bp`;
  });
  document.getElementById('mp-add-btn').addEventListener('click', () => {
    const nameEl = document.getElementById('mp-name');
    const descEl = document.getElementById('mp-desc');
    const lenEl = document.getElementById('mp-len');
    const roleEl = document.getElementById('mp-role');
    const name = nameEl.value.trim();
    const len = parseInt(lenEl.value, 10);
    if (!name){ nameEl.focus(); return; }
    if (!Number.isFinite(len) || len <= 0){ lenEl.focus(); return; }
    state.userParts.push({ id:'up'+(uidCounter++), n:name, note:descEl.value.trim(), len, role: roleEl.value==='acceptor'?'acceptor':'insert' });
    saveUserParts();
    syncUserPartsIntoLists();
    nameEl.value = ''; descEl.value = ''; lenEl.value = ''; roleEl.value = 'insert';
    document.getElementById('mp-seq').value = '';
    document.getElementById('mp-seq-result').textContent = '';
    renderUserPartsList();
    renderDbFilters();
    renderDbTable();
    nameEl.focus();
  });

  /* ============ Init ============ */
  document.addEventListener('click', (e) => {
    if (openCombo && !openCombo.panel.contains(e.target) && e.target !== openCombo.input) closeOpenCombo();
  });
  syncUserPartsIntoLists();
  renderUserPartsList();
  renderDbKitFilters();
  renderDbFilters();
  renderMasterMixSettings();
  renderCyclingTable();
  renderAll();
}

loadData();
