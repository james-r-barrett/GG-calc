# Golden Gate / MoClo Assembly Volume Calculator

A single-page web tool for planning Golden Gate and MoClo assembly reactions. It turns stock DNA concentrations into pipetting volumes, suggests a thermocycler program for the digestion–ligation cycle, and keeps a searchable library of acceptor vectors, level-0 parts, and linkers on hand while you set up a reaction.

Everything runs client-side in the browser. There is no backend, build step, or server-side storage — the page is opened directly (or served as static files), and nothing entered into it leaves the machine.

## Contents

- [Overview](#overview)
- [Calculator](#calculator)
- [Part library](#part-library)
- [Data files](#data-files)
- [Project structure](#project-structure)
- [Running locally](#running-locally)
- [Notes and assumptions](#notes-and-assumptions)

## Overview

The page has two tabs:

- **Calculator** — build one or more assemblies (an acceptor vector plus its insert fragments), enter stock concentrations, and get the volume of each fragment, buffer, ligase, and restriction enzyme needed for the reaction. Also recommends a thermocycler program for the digestion–ligation step.
- **Part Library** — browse and search the built-in collection of acceptor vectors, level-0 parts, and end linkers, and maintain your own list of frequently-used parts (name and length only) that persists in the browser.

## Calculator

### Assemblies

Each assembly consists of one acceptor vector and one or more insert fragments. Fragments can be:

- picked from the built-in part library or your saved parts (by name), or
- entered as a custom fragment, by pasting a sequence or by giving a length in bp directly.

For each fragment you enter a stock concentration (ng/µL), and the calculator returns the volume to add so that the fragment reaches its target amount (in fmol) in the final reaction. Assemblies can be added, duplicated, or removed independently; multiple assemblies are calculated side by side.

### Volume calculation

Two calculation modes are available:

- **Sequence composition (exact)** — uses the actual A/C/G/T composition of the fragment to compute its molecular weight, giving an exact conversion between fmol and ng. Used automatically whenever a sequence is available.
- **Length only (bp average)** — falls back to an average per-bp molecular weight when only a length in bp is known (no sequence). Less accurate, since it ignores base composition, but works for parts recorded as "name + length" only.

Targets (acceptor fmol, insert fmol per fragment, total reaction volume) are global settings with tier-based defaults that scale automatically with the number of inserts in the largest assembly (a single-insert assembly needs less than a 7+ fragment assembly); any of these can be overridden manually.

### Master mix

Digestion–ligation reagents (buffer, restriction enzyme, T4 ligase) default to per-reaction volumes that scale with insert count, and can be edited individually. A "shared master mix" mode can be switched on to batch these reagents across all assemblies on the page at once, including an adjustable overage for pipetting error, rather than calculating each assembly independently.

### Recommended cycling

Given the number of fragments in each assembly and a choice of Type IIS restriction enzyme (BsaI, BbsI/BpiI, SapI/BspQI/PaqCI, or BsmBI-v2, or a custom enzyme with an editable digestion temperature), the calculator proposes a thermocycler program: a digestion–ligation cycling block, a final digestion step, and an indefinite 4 °C hold. Cycle count and step times are editable, and a "fit to a time budget" option re-optimises the program against a target run length. An expandable panel shows the reasoning (and underlying pilot data) behind the ~30-cycle default: efficiency keeps rising well beyond 30 cycles, but per-colony accuracy plateaus by roughly that point.

### Export

Pipetting volumes for all assemblies can be downloaded as a CSV file.

## Part library

The "Part Library" tab lists all built-in acceptor vectors, level-0 parts, and linkers, with search-by-name and category filters. The built-in library currently covers two kits, switchable with a kit toggle above the category filters (an "Other" bucket appears automatically if a future kit doesn't fit either):

- **Nuclear kit** — the Chlamydomonas MoClo Kit (Crozet et al.), for nuclear engineering.
- **Chloroplast kit** — the CHLOROMODAS chloroplast MoClo kit (Inckemann, Chotel, Burgis et al., *Nat. Plants* 2025; [ChlamyMarburg/ChloroplastTools](https://github.com/ChlamyMarburg/ChloroplastTools)). Category names follow that kit's own naming/position convention (5'/3' Homology, 5'/3' Connector, Operon Connectors, IEE, N-/C-tag, etc.), and part entries carry the 96-well plate position where available.

Alongside the built-in library, you can maintain your own list of saved parts (shown regardless of the kit toggle):

- **Add manually** — record a part by name, optional description, and length in bp, without needing the full sequence. Useful for backbones and fragments the lab reuses often. A helper field will calculate the length of a pasted sequence for you.
- **Add from CSV** — bulk-import a spreadsheet with columns: name, description, length in bp, and an optional fourth column ("acceptor" or "insert", defaulting to insert).
- **Add from GenBank** — read the LOCUS/DEFINITION/ORIGIN fields from one or more `.gb`/`.gbk` files (for example, downloaded from Addgene or a sequencing provider). Selecting a single-record file opens it in the form for review before saving; selecting multiple files, or a multi-record file, imports them all directly. Only the name, description, and calculated length are kept — the sequence itself is never stored.
- **Backup / restore** — saved parts can be downloaded as a JSON backup and restored later.

Saved parts are held in the browser's local storage only. They are not stored on a server or synced between devices or browsers, so clearing browser data will delete them — use the backup feature regularly.

## Data files

The built-in libraries live in `data/` as JSON:

| File | Contents |
|---|---|
| `data/acceptors.json` | Acceptor/resistance vectors — name, sequence, cloning enzyme, assembly level, position, and selection marker. |
| `data/parts.json` | Level-0 parts — name, sequence, category (promoter, 5'/3' UTR, CDS-adjacent elements, reporters, resistance markers, etc.), and 5'/3' overhang syntax. |
| `data/linkers.json` | End linkers and dummy fragments used to fill unused positions in an assembly, with the level and position(s) they occupy. |

Any entry in these files may omit its sequence (`s`) and give a length in bp (`len`) instead — such entries remain fully selectable and usable in volume calculations under length-only mode, but no sequence is fetched, displayed, or exposed for them.

Every entry also carries a `kit` field (`"nuclear"` or `"chloroplast"`) that drives the Part Library kit toggle. Chloroplast entries additionally carry `pos` (the kit's own MoClo position code, e.g. `"3b"` for an N-tag) and, where a confident match to the kit's plate map was possible, `well` (the 96-well plate position).

## Project structure

```
index.html        Page markup and static content (both tabs)
css/styles.css     Styling
js/app.js          Application logic: state, calculations, rendering, import/export
data/*.json        Built-in acceptor, part, and linker libraries
images/logo.svg    Site logo / favicon
```

There is no build process — `index.html` loads `js/app.js` and `css/styles.css` directly, and `js/app.js` fetches the three JSON data files on load.

## Running locally

Because the page fetches the JSON data files with `fetch()`, it needs to be served over HTTP rather than opened as a bare `file://` URL in most browsers. From the project root:

```bash
python3 -m http.server 8000
```

then open `http://localhost:8000` in a browser.

## Notes and assumptions

- Volumes default to using actual base composition (A/C/G/T molecular weights); length-only mode uses an average bp weight instead, and is less accurate.
- Master-mix buffer defaults to 1/10 of the reaction volume; ligase and restriction enzyme default to 0.5/0.3 µL for a single insert, 1/0.6 µL for 2–6 inserts, or 2/1.2 µL for 7+ inserts — all editable, with water making up the remainder.
- Estimated cycling run time assumes a 3 °C/sec block ramp rate and a 5-minute initial heat-up; these are rough estimates, not measured values.
- Nothing entered into the calculator or the part library leaves the browser — there is no server component.
