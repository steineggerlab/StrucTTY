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
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="License"/></a>
  <a href="https://github.com/steineggerlab/StrucTTY/actions/workflows/build.yml"><img src="https://github.com/steineggerlab/StrucTTY/actions/workflows/build.yml/badge.svg" alt="Build"/></a>
  <a href="https://doi.org/10.64898/2026.03.17.712308"><img src="https://img.shields.io/badge/bioRxiv-10.64898%2F2026.03.17.712308-b31b1b.svg?logo=biorxiv" alt="bioRxiv"/></a>
</p>

---

**StrucTTY** is a lightweight, terminal-based protein structure visualizer built in C++17. It renders 3D protein structures directly in the terminal using **Unicode Braille sub-pixel rendering**, providing 8x resolution compared to standard character-based rendering.

StrucTTY supports simultaneous visualization of up to 9 proteins, 11 color modes with 3-band depth fog, and integrates with **Foldseek** and **FoldMason** for structural search and multiple structure alignment.

## Features

- **Braille sub-pixel rendering** — each terminal cell maps to a 2×4 logical pixel grid
- **Up to 9 proteins** rendered simultaneously with independent controls
- **11 color modes** — `protein`, `chain`, `rainbow`, `bfactor`, `plddt`, `interface`, `conservation`, `scale`, `align`, `align-fs`, `align-near`
- **3-band depth fog** — near (bright), mid (normal), far (dark with hue retention) for depth perception
- **Secondary structure visualization** — helices as helical ribbons, strands as flat zig-zag ribbons. HELIX/SHEET records are used when the file carries them, otherwise [P-SEA](https://doi.org/10.1093/bioinformatics/13.3.291) is applied to the Cα trace — the only class of method that works on the Cα-only coordinates of a Foldseek database. Cross-checked against the reference implementation in Biotite (92.4% of residues, exact on continuous chains)
- **Foldseek integration** — load `.m8`/`_report` results (`-fsr`), navigate hits, and take targets from a Foldseek DB, a local directory, or automatic download (`-fst`)
- **FoldMason integration** — MSA superposition with conservation coloring
- **MSA conservation scoring** — Shannon entropy from FASTA/A3M alignments
- **Interface detection** — inter-chain contact residue highlighting
- **Alignment visualization** — structural alignment region highlighting
- **Biological assemblies** — `--assembly <ID>` expands the requested assembly instead of the deposited asymmetric unit
- **Light and dark themes** — `--theme light` remaps the palette for light terminal backgrounds
- **Aspect correction** — `--aspect <FLOAT>` for terminals whose character cells are not in the usual 1:2 ratio
- **Custom residue scales** — `-m scale --scale <FILE>` colors by any per-residue property (hydrophobicity, charge, ...)
- **Screenshot export** — PNG output via stb_image_write (`C` key)
- **Chain selection** — filter specific chains per protein

## Installation

### Requirements

| Dependency | Version |
|-----------|---------|
| C++ compiler | GCC ≥ 7.1 or Clang ≥ 5.0 (C++17) |
| CMake | ≥ 3.15 |

The terminal needs UTF-8 and 256-colour support. Mouse hover additionally needs SGR
mouse reporting (`1006`) and is simply inactive without it. If the terminal size cannot
be queried, an 80×24 layout is assumed. No ncurses and no graphical subsystem is used.

#### Verified configurations

Every push builds from a clean checkout and runs a rendering smoke test on:

| OS | Architecture | Compiler |
|---|---|---|
| macOS (GitHub `macos-latest` image) | Apple Silicon (arm64) | Apple Clang (system `c++`) |
| Ubuntu (GitHub `ubuntu-latest` image) | x86-64 | GCC (system `c++`) |

The exact image and compiler versions of each run are printed in the
[build log](https://github.com/steineggerlab/StrucTTY/actions/workflows/build.yml).
The smoke test covers single-structure rendering, Foldseek hit rendering with
superposition, and conservation colouring from an A3M alignment, using the published
example archive — so a broken example link fails the build as well.

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

The structures used below are distributed separately, from the
[`data` release](https://github.com/steineggerlab/StrucTTY/releases/tag/data).
Unpack them in the repository root, next to `build/`:

```bash
cd StrucTTY
curl -L -o example.zip https://github.com/user-attachments/files/31627970/structty_example.zip
unzip -d example example.zip && rm example.zip
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
./StrucTTY ../example/3A0C-assembly1_colabfold.pdb --mode plddt    # predicted model
./StrucTTY ../example/3A0C-assembly1.cif            --mode bfactor  # experimental
./StrucTTY ../example/3A0C-assembly1.cif --mode interface
./StrucTTY ../example/3A0C-assembly1.cif --mode conservation \
  --msa ../example/msa_result/query.a3m
```

Any per-residue scale can be used as a colour source. The file is one
`<one-letter amino acid> <value>` per line (`#` comments allowed); values are
rescaled to the same ten-step ramp, and residues missing from the file stay uncoloured:

```bash
# Kyte-Doolittle hydrophobicity
printf 'A 1.8\nR -4.5\nN -3.5\nD -3.5\nC 2.5\nQ -3.5\nE -3.5\nG -0.4\nH -3.2\nI 4.5\nL 3.8\nK -3.9\nM 1.9\nF 2.8\nP -1.6\nS -0.8\nT -0.7\nW -0.9\nY -1.3\nV 4.2\n' > kd.txt
./StrucTTY ../example/3A0C-assembly1.cif -m scale --scale kd.txt
```

### Biological assemblies

Crystal structures are deposited as an asymmetric unit, which may hold more or fewer
chains than the biological unit. `--assembly` expands the one you ask for:

```bash
curl -fsSL -o 1STP.cif https://files.rcsb.org/download/1STP.cif

./StrucTTY 1STP.cif -m chain                 # asymmetric unit — chain A, 121 residues
./StrucTTY 1STP.cif -m chain --assembly 1    # biological assembly — A1..A4, 4 x 121
```

Duplicated chains get a numeric suffix (`A1`, `A2`, ...) so they stay distinguishable in
the panel and in `-m chain`. An unknown ID lists the ones the file actually defines, and a
structure carrying no assembly records says so; in both cases StrucTTY exits with an error
rather than opening a viewer that does not show what was asked for.
The option applies to structure files; coordinates taken from a Foldseek database are
already per chain.

### Light terminal backgrounds

```bash
./StrucTTY ../example/3A0C-assembly1.cif -m chain --theme light
```

The far depth band fades toward white instead of black, so depth stays readable on a
light background. Use `--aspect` if your terminal's character cells are not close to
1:2 and the structure looks horizontally squashed or stretched.
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

The same hit under both alignment modes — `align-fs` colours only what Foldseek
aligned, `align-near` colours whatever ends up close in space:

```bash
./StrucTTY ../example/foldseek_result/DB1/3cna-assembly1.cif \
  -fst ../example/foldseek_result/DB2/ \
  -fsr ../example/foldseek_result/result \
  -m align-fs -s
```
```bash
./StrucTTY ../example/foldseek_result/DB1/3cna-assembly1.cif \
  -fst ../example/foldseek_result/DB2/ \
  -fsr ../example/foldseek_result/result \
  -m align-near -s
```

In the recording, the right pane is `align-fs` and the left one is `align-near`.

https://github.com/user-attachments/assets/0118ba8f-e77e-455c-87c8-197f008ffe7a

Pass the query directory instead of one structure to walk every query in the
result with `]`/`[`:

```bash
./StrucTTY ../example/foldseek_result/DB1/ \
  -fst ../example/foldseek_result/DB2/ \
  -fsr ../example/foldseek_result/result \
  -m align-fs -s
```

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
| `-m, --mode <MODE>` | Color mode: `protein` (default), `chain`, `rainbow`, `bfactor`, `plddt`, `interface`, `conservation`, `scale`, `align`, `align-fs`, `align-near` |
| `-c, --chains <FILE>` | Chain selection file (TSV: index + chain IDs) |
| `-s, --structure` | Show secondary structure (helix/sheet) |
| `--msa <FILE>` | MSA file for conservation scoring (FASTA/A3M) |
| `-fst, --foldseek-target <PATH>` | Target source for Foldseek hits: Foldseek DB, structure directory, structure file, or `auto` (download from public DBs) |
| `-fsr, --foldseek-result <FILE>` | Foldseek result: `.m8` (12/17/21/29 columns) or multimer `_report` (14 columns) |
| `-fm, --foldmason <FILE>` | FoldMason result (JSON or FASTA MSA) |
| `-n, --nopanel` | Hide info panel |
| `--aspect <FLOAT>` | Horizontal stretch of the scene (default `1.0`), for terminals whose cells are not ~1:2 |
| `--theme <dark\|light>` | Colour scheme for a dark (default) or light terminal background |
| `--assembly <ID>` | Load the given biological assembly (e.g. `1`) instead of the asymmetric unit |
| `--scale <FILE>` | Per-residue values used by `-m scale`; each line is `<one-letter AA> <value>` |
| `-h, --help` | Show the option list, the accepted input kinds, and example invocations |

`-fst` and `-fsr` must be given together — one without the other is an error.

### Supported inputs

The kind of every input is detected automatically (no format flags):

| Kind | What it is | Accepted as |
|---|---|---|
| Structure file | `.pdb` / `.cif` / `.ent` (+ `.gz`) | query, `-fst` |
| Structure directory | a directory of those files | query, `-fst` |
| Foldseek DB | base path of a DB built **from structures** (needs `<db>_ca`) | query, `-fst` |
| Foldseek result | `.m8` (12/17/21/29 columns) or multimer `_report` (14 columns) | `-fsr` |

Only residues carrying a Cα atom are read. Waters, ligands and ions therefore never
appear, with no flag needed; modified residues that do carry a Cα (MSE, for instance)
are kept. A content/extension mismatch — mmCIF text in a `.pdb`, say — is reported as a
diagnostic instead of aborting.

A `_report` (14 columns) enters the multimer path. The query is where the per-complex chains are
read from, so it must be a Foldseek **query DB** or a **directory** of structures.

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
| `]` / `[` | Next query / Previous query |
| `C` | Save a PNG screenshot of the current view |
| `Q` | Quit |

> Mouse hover displays residue information in the info panel.

## Color Modes

| Mode | Description |
|------|-------------|
| `protein` | One color per protein (9 distinct colors, cycling) |
| `chain` | One color per chain (15 colors) |
| `rainbow` | N→C gradient (20-step hue spectrum) |
| `bfactor` | B-factor column on a blue–white–red ramp spanning the observed range of the loaded structures, as in ChimeraX. The panel shows that range next to the hovered value |
| `plddt` | Same column, but binned on the AlphaFold confidence scale that predicted models store there: blue (≥90), cyan (70–90), yellow (50–70), orange (<50) |
| `interface` | Inter-chain contacts (CA–CA < 8 Å): magenta vs. dim |
| `conservation` | MSA Shannon entropy: blue (variable) → red (conserved) |
| `scale` | Any per-residue scale given with `--scale` (e.g. hydrophobicity), rescaled to the same ten-step ramp |
| `align` | Structurally aligned regions: bright vs. dim gray. Uses the alignment when the result has one, otherwise falls back to distance |
| `align-fs` | Only what Foldseek aligned (`qaln`/`taln` columns). Never falls back — errors out if the result carries no alignment |
| `align-near` | Distance only: residues with a counterpart within the cutoff, whatever the result says |

The panel's `Align:` line names the colour source actually used: `aln-string`
(Foldseek alignment), `msa-col` (FoldMason MSA columns) or `nearest-nbr`
(distance). The `Superpose:` line in the hit block is a different thing — how the
structures were laid on top of each other — and does not change with the colour
mode. It is chosen automatically in this order, and names the step that ran:

| `Superpose:` | Source |
|---|---|
| `ut-matrix` | The U/T rotation-translation matrices carried by a 29-column result, applied directly |
| `kabsch-alns` | Kabsch fit on the aligned coordinates of a 21-column `alis` result |
| `kabsch-qaln` | Kabsch fit on the residue pairs implied by `qaln`/`taln` (17 columns) |
| `kabsch-seq` | Kabsch fit on residues in file order, when the result carries none of the above |

`RMSD` is reported over the residue pairs the fit used; with `ut-matrix` there is
no explicit pair list, so the alignment length is shown instead. Multimer
`_report` files carry their own U/T matrices and always use them.

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
- Structural superposition, chosen automatically in four steps and named in the panel's
  `Superpose:` line: the U/T rotation-translation matrices if the result carries them
  (29 columns), else a Kabsch fit on the aligned coordinates of an `alis` result
  (21 columns), else on the residue pairs implied by `qaln`/`taln` (17 columns), else on
  residues in file order. A 12-column default result therefore still opens, but only
  reaches that last step
- Per-hit `RMSD: <value> A over <N> res` in the panel, over the residue pairs the fit
  actually used. When the superposition comes from U/T matrices there is no explicit pair
  list, so the alignment length is shown instead
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

StrucTTY is **embedded directly into Foldseek as a static library** (`add_subdirectory(lib/structty)`), so no external binary and no `PATH` lookup is involved — the viewer lives inside the `foldseek` binary itself and is reached through the `foldseek structty` subcommand, which the `easy-search` and `easy-multimersearch` workflows invoke for you. The viewer opens automatically once the search finishes, reading the query and target structures directly from the search's temporary databases (folder/tar/gz inputs supported). Temporary DBs are kept alive for the viewer and cleaned up after it closes. Supported workflows: `easy-search` and `easy-multimersearch`.

**Automatic launch after a search** — add the `--view-structty` flag (it takes no value):

```bash
foldseek easy-search query.cif targetDir result.m8 tmp --view-structty \
  -a --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln,u,t"
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
  -a --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln,u,t"
```

Add `--alignment-type 1` to align with TM-align instead of 3Di+AA — slower, but
it usually reports a tighter set of aligned residues.

For a result that already exists as a Foldseek database, convert it instead of
searching again:

```bash
foldseek convertalis queryDB targetDB resultDB result.m8 \
  --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,qaln,taln,u,t"
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

# Multimer: 14-column _report, chains read per complex from the query DB or directory.
# A _report carries no alignment strings, so align-near is the only align mode.
./StrucTTY ../example/foldseek_result/DB1/ \
  -fst ../example/foldseek_result/DB2/ \
  -fsr ../example/foldseek_result/result_multimer_report -m align-near
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

| Structure | Cα | Load (ms) | TTFF (ms) | Frame mean (ms) | Frame p95 (ms) | Latency mean (ms) | Latency p95 (ms) | Peak RSS (MB) |
|-----------|----:|----------:|----------:|----------------:|---------------:|------------------:|-----------------:|--------------:|
| 1CRN | 46 | 4.8 | 39.6 | 0.16 | 0.18 | <1 | <1 | 6.0 |
| 1STP | 121 | 11.0 | 74.4 | 0.28 | 0.32 | <1 | <1 | 7.3 |
| 3BIK | 446 | 45.2 | 175.6 | 0.62 | 0.67 | <1 | <1 | 13.4 |
| 6VXX | 2916 | 201.4 | 856.8 | 3.25 | 3.37 | 3.00 | 3.00 | 44.7 |
| 4V4Q | 11463 | 2575.6 | 4744.8 | 10.58 | 10.78 | 10.34 | 11.00 | 379.0 |

Measured on a single core (pinned with `taskset`) of an Arm Cortex-X925, memory capped at
4 GB, Ubuntu 24.04 LTS, g++ 13.3, CMake Release. Mean of 5 independent runs; each run is
200 warm-up iterations followed by 2000 timed ones. Frame time is the sum of the six
microsecond-resolution stage timers; the other metrics are logged at millisecond
resolution, so latencies below that limit are shown as `<1`. Peak resident set size comes
from `/usr/bin/time -v` over a load-and-render run.

- **Frame time stays under the 16 ms (60 fps) budget** even at 11,463 Cα — interaction remains smooth for structures spanning three orders of magnitude in size.
- **Input latency tracks frame time closely**, so rotation/zoom feels responsive with no perceptible input lag.
- Load and TTFF scale with atom count; the immediate z-test rasterizer keeps per-frame cost roughly linear in Cα count.
- **Peak memory is set by the largest single file being parsed, not by how much is on screen** — opening the same 446-residue structure nine times at once moves the peak from 13.4 to 14.4 MB. Loading a Foldseek result and stepping through its hits stays around 22 MB whether the targets come from a directory or straight from a Foldseek DB, and conservation colouring from an A3M alignment is indistinguishable from the same structure without it (9.2 MB either way).

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
