<p align="center">
  <img src=".github/logo_final.png" alt="StrucTTY Logo" width="320"/>
</p>

<h1 align="center">StrucTTY</h1>

<p align="center">
  <b>Interactive, Terminal-Native Protein Structure Viewer</b>
</p>

<p align="center">
  <a href="#-installation"><img src="https://img.shields.io/badge/C%2B%2B-17-blue.svg" alt="C++17"/></a>
  <a href="#-installation"><img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg" alt="Platform"/></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-GPLv3-green.svg" alt="License"/></a>
  <a href="https://github.com/steineggerlab/StrucTTY"><img src="https://img.shields.io/badge/build-CMake-orange.svg" alt="Build"/></a>
  <a href="https://doi.org/10.64898/2026.03.17.712308"><img src="https://img.shields.io/badge/bioRxiv-10.64898%2F2026.03.17.712308-b31b1b.svg?logo=biorxiv" alt="bioRxiv"/></a>
</p>

---

**StrucTTY** is a lightweight, terminal-based protein structure visualizer built in C++17. It renders 3D protein structures directly in the terminal using **Unicode Braille sub-pixel rendering**, providing 8x resolution compared to standard character-based rendering.

StrucTTY supports simultaneous visualization of up to 9 proteins, 9 color modes with 3-band depth fog, and integrates with **Foldseek** and **FoldMason** for structural search and multiple structure alignment.

## Features

- **Braille sub-pixel rendering** — each terminal cell maps to a 2×4 logical pixel grid
- **Up to 9 proteins** rendered simultaneously with independent controls
- **9 color modes** — `protein`, `chain`, `rainbow`, `plddt`, `interface`, `conservation`, `align`, `align-fs`, `align-near`
- **3-band depth fog** — near (bright), mid (normal), far (dark with hue retention) for depth perception
- **Secondary structure visualization** — helix cylinders and sheet ribbons
- **Foldseek integration** — load `.m8`/`_report` results (`-fsr`), navigate hits, and take targets from a Foldseek DB, a local directory, or automatic download (`-fst`)
- **FoldMason integration** — MSA superposition with conservation coloring
- **MSA conservation scoring** — Shannon entropy from FASTA/A3M alignments
- **Interface detection** — inter-chain contact residue highlighting
- **Alignment visualization** — structural alignment region highlighting
- **Screenshot export** — PNG output via stb_image_write
- **Chain selection** — filter specific chains per protein

## Installation

### Requirements

| Dependency | Version |
|-----------|---------|
| C++ compiler | GCC ≥ 7.1 or Clang ≥ 5.0 (C++17) |
| CMake | ≥ 3.15 |

Supported platforms: **Linux**, **macOS**

### Build (Linux / macOS)

```bash
git clone --recurse-submodules https://github.com/steineggerlab/StrucTTY.git
cd StrucTTY
mkdir build && cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make -j $(nproc)
```

> The output binary will be generated at `build/StrucTTY`.

### Example data

The structures used below are distributed separately. Unpack them in the
repository root, next to `build/`:

```bash
cd StrucTTY
curl -L -o example.zip https://github.com/user-attachments/files/31012856/structty_example.zip
unzip example.zip && rm example.zip
```

## Quick Start

### Single structure

```bash
./StrucTTY ../example/1NPL-assembly1.cif --mode chain
```
https://github.com/user-attachments/assets/64c37f49-7ae8-4dea-bec3-37575583a767

### Secondary structure visualization

```bash
./StrucTTY ../example/3HGM-assembly1.cif --mode chain
./StrucTTY ../example/3HGM-assembly1.cif --mode chain -s
```
https://github.com/user-attachments/assets/01d91bc2-4c49-46f8-80c5-d736ccdeea07

### Color modes

```bash
./StrucTTY ../example/3A0C-assembly1.cif                  # protein (default)
./StrucTTY ../example/3A0C-assembly1.cif --mode chain
./StrucTTY ../example/3A0C-assembly1.cif --mode rainbow
./StrucTTY  ../example/3A0C-assembly1_colabfold.pdb --mode plddt
./StrucTTY ../example/3A0C-assembly1.cif --mode interface
./StrucTTY ../example/3A0C-assembly1.cif --mode conservation \
  --msa ../example/msa_result/query.a3m
```
https://github.com/user-attachments/assets/d4fea46a-2dc5-4e92-85cc-aa47ebcdb2d1

### Multiple structures

```bash
./StrucTTY ../example/1CJK-assembly1.cif \
  ../example/1NPL-assembly1.cif \
  ../example/3A0C-assembly1.cif \
  ../example/3HGM-assembly1.cif \
  ../example/3OAG-assembly1.cif \
  ../example/9FL9-assembly1.cif \
  ../example/AF-A0A233SAX3-F1-model_v6.cif \
  ../example/9N47-assembly1.cif \
  ../example/8KGM-assembly1.cif
```
https://github.com/user-attachments/assets/eb2fae4c-4f64-489e-b891-26505d55179c

### Chain selection

```bash
./StrucTTY ../example/9N47-assembly1.cif -m chain
./StrucTTY ../example/9N47-assembly1.cif -c ../example/chainfile_9N47.tsv -m chain
```
https://github.com/user-attachments/assets/f9cfca51-bba3-4090-b021-a93ad1e671bc

### Foldseek hit navigation

```bash
./StrucTTY ../example/foldseek_result/DB1/ \
  -fst ../example/foldseek_result/DB2/ \
  -fsr ../example/foldseek_result/result \
  -m align-fs -s
```
```bash
./StrucTTY ../example/foldseek_result/DB1/3cna-assembly1.cif \
  -fst ../example/foldseek_result/DB2/ \
  -fsr ../example/foldseek_result/result \
  -m align-fs -s
```

https://github.com/user-attachments/assets/9bf428eb-4506-40f2-9c86-b92fab0d47b7

### FoldMason MSA superposition

```bash
./StrucTTY ../example/3A0C-assembly1.cif  ../example/L7RCY6.pdb 
./StrucTTY ../example/3A0C-assembly1.cif  ../example/L7RCY6.pdb \
  -fm ../example/foldmason_result/foldmason.json -m align
```
![FoldMason_alignment](https://github.com/user-attachments/assets/3237329b-b0bf-4fa1-b580-c12df4df203e)

## Usage

```
./StrucTTY <query...> [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-m, --mode <MODE>` | Color mode: `protein` (default), `chain`, `rainbow`, `plddt`, `interface`, `conservation`, `align`, `align-fs`, `align-near` |
| `-c, --chains <FILE>` | Chain selection file (TSV: index + chain IDs) |
| `-s, --structure` | Show secondary structure (helix/sheet) |
| `--msa <FILE>` | MSA file for conservation scoring (FASTA/A3M) |
| `-fst, --foldseek-target <PATH>` | Target source for Foldseek hits: Foldseek DB, structure directory, structure file, or `auto` (download from public DBs) |
| `-fsr, --foldseek-result <FILE>` | Foldseek result: `.m8` (12/17/21/29 columns) or multimer `_report` (14 columns) |
| `-fm, --foldmason <FILE>` | FoldMason result (JSON or FASTA MSA) |
| `-n, --nopanel` | Hide info panel |

`-fst` and `-fsr` must be given together — one without the other is an error.

### Supported inputs

The kind of every input is detected automatically (no format flags):

| Kind | What it is | Accepted as |
|---|---|---|
| Structure file | `.pdb` / `.cif` / `.ent` (+ `.gz`) | query, `-fst` |
| Structure directory | a directory of those files | query, `-fst` |
| Foldseek DB | base path of a DB built **from structures** (needs `<db>_ca`) | query, `-fst` |
| Foldseek result | `.m8` (12/17/21/29 columns) or multimer `_report` (14 columns) | `-fsr` |

A `_report` (14 columns) enters the multimer path and requires a Foldseek **query DB** as the query,
because the per-complex chains are read from that DB.

With `-fst <directory>`, hits are looked up by accession. Foldseek splits multimers per chain
(`1dci-assembly1_B-2`), so the trailing `_<chain>` is stripped until a file matches
(`1dci-assembly1.cif`) and that chain alone is drawn — a plain directory of structures works as a
target, no `createdb` needed.

**Sequence FASTA is not supported.** It carries no 3D coordinates, and `foldseek createdb
--prostt5-model` predicts 3Di (`_ss`) without writing any `_ca`, so a sequence-derived DB cannot be
rendered either. Both cases fail before rendering starts, with the reason printed.

## Keyboard Controls

| Key | Action |
|-----|--------|
| `0` | Control all proteins |
| `1`–`9` | Control individual protein |
| `W` / `A` / `S` / `D` | Move up / left / down / right |
| `X` / `Y` / `Z` | Rotate around X / Y / Z axis |
| `R` / `F` | Zoom in / out |
| `N` / `P` | Next hit / Previous hit |
| `[` / `]` | Next query / Previous query |
| `Q` | Quit |

> Mouse hover displays residue information in the info panel.

## Color Modes

| Mode | Description |
|------|-------------|
| `protein` | One color per protein (9 distinct colors, cycling) |
| `chain` | One color per chain (15 colors) |
| `rainbow` | N→C gradient (20-step hue spectrum) |
| `plddt` | AlphaFold confidence: blue (≥90), cyan (70–90), yellow (50–70), orange (<50) |
| `interface` | Inter-chain contacts (CA–CA < 8 Å): magenta vs. dim |
| `conservation` | MSA Shannon entropy: blue (variable) → red (conserved) |
| `align` | Structurally aligned regions: bright vs. dim gray. Uses the alignment when the result has one, otherwise falls back to distance |
| `align-fs` | Only what Foldseek aligned (`qaln`/`taln` columns). Never falls back — errors out if the result carries no alignment |
| `align-near` | Distance only: residues with a counterpart within the cutoff, whatever the result says |

The panel's `Align:` line names the source actually used: `aln-string` (Foldseek
alignment), `msa-col` (FoldMason MSA columns) or `nearest-nbr` (distance).

> The alignment-based modes need `qaln`/`taln`, i.e. a 17/21/29-column `.m8`.
> Foldseek's default 12-column output has none.
> Generate a usable result with `--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln`
> (the search itself must run with `-a`). `foldseek ... --view-structty` does this for you.
> Multimer `_report` files carry no alignment strings at all, so they only work with `align-near`.

All modes support **3-band depth fog**: near (vivid), mid (normal), far (dark, hue-retaining).

## Integrations

### Foldseek

StrucTTY reads Foldseek `easy-search` output (`.m8` format) with support for 12, 17, 21, and 29 column formats. Features include:

- Interactive hit navigation with automatic structure downloading
- Direct Foldseek DB reading (`-fst <DB>`) read Cα coordinates directly from Foldseek `_ca` database, eliminating network dependency. Uses hit-based selective scanning for minimal memory usage (~152KB for 1000 hits, even on AFDB50)
- Structural superposition using U/T rotation-translation matrices
- Alignment string visualization (`qaln`/`taln`)
- Multi-database support: PDB, AlphaFold DB, ESMAtlas, CATH, BFVD, and more
- Multi-query navigation (`]`/`[`) across queries in a single `.m8`. The query can be a Foldseek DB
  **or a plain directory** — each accession in the result is looked up inside it, so a standalone
  run reproduces what `foldseek --view-structty` shows
- With a structure file as the query, hits are filtered to that file — one `.m8` covering a whole
  query directory no longer walks you through other queries' hits. Foldseek indexes a multimer per
  chain, so `]`/`[` steps through the query's chains and `N`/`P` through that chain's hits; the
  panel shows `Q[3/6][2 / 14]`
- Multimer (complex-level) report viewing with per-complex superposition

#### Launch from Foldseek

StrucTTY is **embedded directly into Foldseek as a static library** (`add_subdirectory(lib/structty)`), so no external binary or `PATH` lookup is required — Foldseek calls `structty::run()` in-process. The viewer opens automatically once the search finishes, reading the query and target structures directly from the search's temporary databases (folder/tar/gz inputs supported). Temporary DBs are kept alive for the viewer and cleaned up after it closes. Supported workflows: `easy-search` and `easy-multimersearch`.

**Automatic launch after a search** — add the `--view-structty` flag (it takes no value):

```bash
foldseek easy-search query.cif targetDir result.m8 tmp --view-structty \
  -a --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln"
```

The viewer builds that 17-column layout internally either way; passing
`-a --format-output` keeps a copy in `result.m8`, so the same hits can be
reopened later with `-fsr` and still show lDDT, TM scores and the aligned
regions. Without it the file keeps Foldseek's 12-column default, which carries
no alignment — `align-fs` then refuses to run and `align` falls back to
distance.

**Multimer (complex-level) search** — the viewer works out of the box, since the per-complex report (`--multimer-report-mode 1`) is the default; setting `--multimer-report-mode 0` skips the launch. The viewer reads the `_report` file the workflow already writes, so it stays after the viewer closes:

```bash
foldseek easy-multimersearch queryDir targetDir result tmp --view-structty
```

This writes `result_report`, the 14-column file the viewer reads. Both sides must
be complexes — a single-chain query or target yields an empty report, and the
viewer then refuses it.

#### Standalone usage

First produce a result that carries the alignment. The two columns that matter
are `qaln` and `taln`; `-a` makes the search keep the backtraces they come from:

```bash
foldseek easy-search query.pdb targetDir result.m8 tmp \
  -a --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln"
```

Add `--alignment-type 1` to align with TM-align instead of 3Di+AA — slower, but
it usually reports a tighter set of aligned residues.

For a result that already exists as a Foldseek database, convert it instead of
searching again:

```bash
foldseek convertalis queryDB targetDB resultDB result.m8 \
  --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln"
```

```bash
# Read target Cα coordinates straight from a Foldseek DB (offline)
./StrucTTY query.pdb -fst /path/to/targetDB -fsr result.m8

# Look targets up by hit accession in a local structure directory
./StrucTTY query.pdb -fst /path/to/pdbs/ -fsr result.m8

# Download hit structures from public DBs (PDB, AFDB, ESMAtlas, CATH, BFVD, ...)
./StrucTTY query.pdb -fst auto -fsr result.m8

# Query from a Foldseek DB — multi-query navigation with ]/[
./StrucTTY /path/to/queryDB -fst /path/to/targetDB -fsr result.m8 -m align-fs

# Colour by distance instead, whatever the result says (works on 12-column files)
./StrucTTY query.pdb -fst /path/to/targetDB -fsr result.m8 -m align-near

# Multimer: 14-column _report, chains read per complex from the query DB.
# A _report carries no alignment strings, so align-near is the only align mode.
./StrucTTY /path/to/queryDB -fst /path/to/targetDB -fsr result_report -m align-near
```

### FoldMason

StrucTTY loads FoldMason MSA results (JSON with Cα coordinates or FASTA) for:

- Kabsch-based structural superposition
- Column-wise conservation scoring
- Gap-aware alignment visualization

### MSA Conservation

Load FASTA or A3M multiple sequence alignments to compute per-residue conservation scores via Shannon entropy, visualized with the `conservation` color mode.

## Performance

StrucTTY renders interactively even for large complexes. The table below measures load time, time-to-first-frame (TTFF), per-frame render time, and input-to-frame latency across structures of increasing size:

| Structure | Cα | Load (ms) | TTFF (ms) | Frame mean (ms) | Frame p95 (ms) | Latency mean (ms) | Latency p95 (ms) |
|-----------|----:|----------:|----------:|----------------:|---------------:|------------------:|-----------------:|
| 1CRN | 46 | 5.4 | 44.6 | <0.5 | <0.5 | <0.5 | <0.5 |
| 1STP | 121 | 12.4 | 106.0 | <0.5 | <0.5 | <0.5 | <0.5 |
| 3BIK | 446 | 47.0 | 205.8 | 0.01 | <0.5 | 0.01 | <0.5 |
| 6VXX | 2916 | 255.0 | 894.4 | 2.97 | 3.00 | 2.99 | 3.00 |
| 4V4Q | 11463 | 3336.8 | 6199.2 | 13.36 | 14.00 | 13.68 | 14.80 |

- **Frame time stays under the 16 ms (60 fps) budget** even at 11,463 Cα — interaction remains smooth for structures spanning three orders of magnitude in size.
- **Input latency tracks frame time closely**, so rotation/zoom feels responsive with no perceptible input lag.
- Load and TTFF scale with atom count; the immediate z-test rasterizer keeps per-frame cost roughly linear in Cα count.

Benchmarks are reproducible via the built-in benchmark mode, which replays a fixed key script and logs per-frame timings to a `structty_bench_*.csv` file.

## Third-Party Libraries

| Library | License | Purpose |
|---------|---------|---------|
| [Gemmi](https://gemmi.readthedocs.io/) | MPL-2.0 | mmCIF/PDB file parsing |
| [stb_image_write](https://github.com/nothings/stb) | MIT / public domain | PNG screenshot encoding |

See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for detailed license information.

## License

This project is licensed under the [MIT License](LICENSE).

---

<p align="center">
  Developed by Luna Sung-eun Jang, Soo Young Cha — <a href="https://github.com/steineggerlab">Steinegger Lab</a>
</p>
