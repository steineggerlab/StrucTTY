# StrucTTY

An Interactive, Terminal-Native Protein Structure Viewer

**Structty** is a lightweight, terminal-based protein structure visualizer designed for fast and interactive comparison of protein structures. It supports simultaneous visualization of two proteins, customizable chain selection, and even external transformation matrices (u, t) for alignment. It is built for extensibility and will support Foldseek integration in the near future.

## ✨ Features

* Render up to six proteins in the terminal
* Move and rotate proteins independently or simultaneously
* Support for chain-specific views
* Helix/sheet secondary structure visualization
* Adjustable screen width and height

## 📦 Installation

### Requirements

* **C++17 compiler** 
  * requires **GCC ≥ 7.1** or **Clang ≥ 5.0**
  * Tested with **GCC 12.2.0 (Debian 12.2.0-14+deb12u1)**
* **CMake ≥ 3.15**
* **Linux or macOS** (Terminal-compatible)

### Build

```bash
git clone --recurse-submodules https://github.com/steineggerlab/StrucTTY.git
cd StrucTTY
mkdir build && cd build
cmake ..
make -j 10
```

```bash
git clone --recurse-submodules https://github.com/steineggerlab/StrucTTY.git
cd StrucTTY
mkdir build && cd build
brew install ncurses
cmake .. \
  -DCURSES_INCLUDE_PATH=/opt/homebrew/opt/ncurses/include \
  -DCURSES_LIBRARY=/opt/homebrew/opt/ncurses/lib/libncursesw.dylib \
  -DCMAKE_EXE_LINKER_FLAGS="-L/opt/homebrew/opt/ncurses/lib -lncursesw -Wl,-rpath,/opt/homebrew/opt/ncurses/lib" \
  -DCMAKE_CXX_FLAGS="-I/opt/homebrew/opt/ncurses/include"
make
```

> Output binary will be generated at `build/StrucTTY`.

## 🚀 Quick Start

### Run a single PDB file (show secondary structure & rainbow color):

```bash
./StrucTTY 1mh1.pdb -s --mode rainbow
```
![example](.github/Structty_example_single.gif)

### Compare two PDB files:

```bash
 ./StrucTTY 7DN3.pdb 2Y0S.pdb 7DU2.pdb
```

### Compare two PDB files using rotation and transform matrix:

```bash
 ./StrucTTY 7DN3.pdb 2Y0S.pdb 7DU2.pdb -ut utfile_ex 
```
![example](.github/Structty_example_multiple.gif)

### With options:

```bash
./StrucTTY 1A52.pdb 1ERR.pdb \
  -c A,B \                     # select chains A and B
  -m chain \                   # color mode: chain / rainbow / default
  -w 3 -h 2 \                  # terminal screen size (width x height units, 1~5)
  -ut utmatrix.txt \            # optional u, t matrix from Foldseek
  -s                          # show secondary structure (helix/sheet)
```

> You can move/rotate proteins individually or together during runtime with keyboard controls.

## 🎮 Keyboard Controls

* `0` — Control all proteins
* `1` — Control only the first protein
* `2` — Control only the second protein
* ...
* `6` — Control only the sixth protein

### Movement

* `W` / `w`: move up (+Y)
* `A` / `a`: move left (-X)
* `S` / `s`: move down (-Y)
* `D` / `d`: move right (+X)

### Rotation

* `X` / `x`: rotate around X-axis
* `Y` / `y`: rotate around Y-axis
* `Z` / `z`: rotate around Z-axis

### Zoom

* `R` / `r`: zoom in
* `F` / `f`: zoom out

### Exit

* `Q` / `q`: quit program

## 🔗 Integration with Foldseek

Structty can accept external transformation matrices (`u`, `t`) output by **Foldseek** for protein alignment visualization.

Future releases will include:

* Direct Foldseek output parser
* Visual alignment validation
* Foldseek-based pairwise comparison interface


## License

This project is licensed under the GNU General Public License v3.0.

This software uses the following third-party libraries:
- Gemmi (MPL-2.0)
- LodePNG (zlib)

See `THIRD_PARTY_NOTICES.md` for detailed license information.


> Developed by Luna Sung-eun Jang and Soo young Cha and C++
