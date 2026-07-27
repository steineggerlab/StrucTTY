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

StrucTTY supports simultaneous visualization of up to 9 proteins, 7 color modes with 3-band depth fog, and integrates with **Foldseek** and **FoldMason** for structural search and multiple structure alignment.

## Features

- **Braille sub-pixel rendering** — each terminal cell maps to a 2×4 logical pixel grid
- **Up to 9 proteins** rendered simultaneously with independent controls
- **7 color modes** — `protein`, `chain`, `rainbow`, `plddt`, `interface`, `conservation`, `aligned`
- **3-band depth fog** — near (bright), mid (normal), far (dark with hue retention) for depth perception
- **Secondary structure visualization** — helix cylinders and sheet ribbons
- **Foldseek integration** — load `.m8` results, navigate hits, auto-download structures or read directly from Foldseek DB
- **FoldMason integration** — MSA superposition with conservation coloring
- **MSA conservation scoring** — Shannon entropy from FASTA/A3M alignments
- **Interface detection** — inter-chain contact residue highlighting
- **Alignment visualization** — structural alignment region highlighting
- **Screenshot export** — PNG output via lodepng
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

## Quick Start

### Single structure

```bash
./StrucTTY ../example/1NPL-assembly1.cif --mode chain
```
![1NPL chain mode](.github/Single_structure.gif)

### Secondary structure visualization

```bash
./StrucTTY ../example/3HGM-assembly1.cif --mode chain
./StrucTTY ../example/3HGM-assembly1.cif --mode chain -s
```
![Secondary structure](.github/Secondary_structure_visualization.gif)

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
![Color modes](.github/Color_modes.gif)

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
![Multi-structure comparison](.github/Multiple_structures.gif)

### Chain selection

```bash
./StrucTTY ../example/9N47-assembly1.cif -m chain
./StrucTTY ../example/9N47-assembly1.cif -c ../example/chainfile_9N47.tsv -m chain
```
![Chain_selection](.github/Chain_selection.gif)

### Foldseek hit navigation

```bash
./StrucTTY ../example/3A0C-assembly1.cif \
  -fs ../example/foldseek_result/alis_afdb50.m8 \
  --db-path /path/to/pdb/
```
![Foldseek_navigation](.github/Foldseek_hit_navigation.gif)

### FoldMason MSA superposition

```bash
./StrucTTY ../example/3A0C-assembly1.cif  ../example/L7RCY6.pdb 
./StrucTTY ../example/3A0C-assembly1.cif  ../example/L7RCY6.pdb \
  -fm ../example/foldmason_result/foldmason.json -m aligned
```
![FoldMason_alignment](.github/FoldMason_MSA_superposition.png)

## Usage

```
./StrucTTY <input_files...> [OPTIONS]
```

| Option | Description |
|--------|-------------|
| `-m, --mode <MODE>` | Color mode: `protein` (default), `chain`, `rainbow`, `plddt`, `interface`, `conservation`, `aligned` |
| `-c, --chains <FILE>` | Chain selection file (TSV: index + chain IDs) |
| `-s, --structure` | Show secondary structure (helix/sheet) |
| `-ut, --utmatrix <FILE>` | Apply rotation/translation matrix for alignment |
| `--msa <FILE>` | MSA file for conservation scoring (FASTA/A3M) |
| `-fs, --foldseek <FILE>` | Foldseek `.m8` result for hit navigation |
| `--db <PATH>` | Foldseek structure database path for offline Cα coordinate reading |
| `--db-path <DIR>` | PDB directory for Foldseek hit loading |
| `-fm, --foldmason <FILE>` | FoldMason result (JSON or FASTA MSA) |
| `-n, --nopanel` | Hide info panel |

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
| `aligned` | Structurally aligned regions: bright vs. dim gray |

All modes support **3-band depth fog**: near (vivid), mid (normal), far (dark, hue-retaining).

## Integrations

### Foldseek

StrucTTY reads Foldseek `easy-search` output (`.m8` format) with support for 12, 17, 21, and 29 column formats. Features include:

- Interactive hit navigation with automatic structure downloading
- Direct Foldseek DB reading (`--db`) read Cα coordinates directly from Foldseek `_ca` database, eliminating network dependency. Uses hit-based selective scanning for minimal memory usage (~152KB for 1000 hits, even on AFDB50)
- Structural superposition using U/T rotation-translation matrices
- Alignment string visualization (`qaln`/`taln`)
- Multi-database support: PDB, AlphaFold DB, ESMAtlas, CATH, BFVD, and more
- Multi-query navigation (`]`/`[`) across queries in a single `.m8`
- Multimer (complex-level) report viewing with per-complex superposition

#### Launch from Foldseek

StrucTTY is **embedded directly into Foldseek as a static library** (`add_subdirectory(lib/structty)`), so no external binary or `PATH` lookup is required — Foldseek calls `structty::run()` in-process.

**Automatic launch after a search** — add the `--view-structty` flag (it takes no value) to any structure search workflow:

```bash
foldseek easy-search query.cif targetDir result.m8 tmp --view-structty
foldseek search queryDB targetDB result tmp --view-structty
foldseek easy-rbh query.cif targetDir result.m8 tmp --view-structty
```

The viewer opens automatically once the search finishes, reading the query and target structures directly from the search's temporary databases (folder/tar/gz inputs supported). Temporary DBs are kept alive for the viewer and cleaned up after it closes. Supported workflows: `easy-search`, `easy-rbh`, `search`, and `easy-multimersearch`.

**Multimer (complex-level) search** — the viewer works out of the box, since the per-complex report (`--multimer-report-mode 1`) is the default; setting `--multimer-report-mode 0` skips the launch:

```bash
foldseek easy-multimersearch queryDir targetDir result tmp --view-structty
```

> If a search runs without `--view-structty`, results are written normally and no viewer is launched.

#### Standalone usage

```bash
# Online mode (download structures)
./StrucTTY query.pdb --foldseek result.m8

# Offline mode (read from Foldseek DB)
./StrucTTY query.pdb --foldseek result.m8 --db /path/to/targetDB

# Hybrid (DB first, fallback to download)
./StrucTTY query.pdb --foldseek result.m8 --db /path/to/targetDB --db-path /path/to/pdbs/
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
| [LodePNG](https://lodev.org/lodepng/) | zlib | PNG screenshot encoding |

See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for detailed license information.

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).

---

<p align="center">
  Developed by Luna Sung-eun Jang, Soo Young Cha — <a href="https://github.com/steineggerlab">Steinegger Lab</a>
</p>
