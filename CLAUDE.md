# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

UMD (Universal Molecular Dynamics) is a Python package for post-processing ab initio and classical molecular dynamics simulations. It works around a central ASCII file format (`.umd.dat`) that stores atomic trajectories and thermodynamic data extracted from simulation codes (VASP, QBox, LAMMPS).

## Building C Extensions

Performance-critical routines are implemented as C shared libraries loaded at runtime via `ctypes`. They **must be compiled before running any analysis script**. All commands run from `src/`:

**macOS:**
```bash
gcc -dynamiclib -o c_UMDprocess.dylib c_UMDprocess.c
gcc -dynamiclib -o c_autocorrelation_vib.dylib c_autocorrelation_vib.c
gcc -dynamiclib -o c_bonds_full.dylib c_bonds_full.c
gcc -dynamiclib -o c_clusters.dylib c_clusters.c
gcc -dynamiclib -o c_msd.dylib c_msd.c
gcc -dynamiclib -o c_gofr_new_alt.dylib c_gofr_new_alt.c
```

**Linux / Windows (`.so` / `.dll`):** See `compilC-mac.txt` and `compilC-windows.txt` for the analogous commands.

`umd_processes_fast.py` auto-detects the OS and loads the correct library extension (`.dylib` / `.so` / `.dll`). The library path is resolved relative to the script's own `__file__` location, so scripts must be run from the `src/` directory or the libraries must be on the system path.

## Dependencies

```
python3, numpy, scipy, matplotlib
```

No package installation is required — scripts are run directly.

## Running Scripts

All scripts are standalone executables. Typical invocation pattern:
```bash
python3 src/VaspParser.py -f <OUTCAR_file> [-i InitialStep] [-s SkipStep]
python3 src/gofr_umd.py -f <file.umd.dat> [options]
python3 src/averages.py -f <file.umd.dat> [-s SkipSteps]
```

Use `-h` or read the argument-parsing block at the bottom of each script to discover its options.

## Architecture

### Core Libraries
- **`crystallography.py`** — Defines the two central data structures used everywhere:
  - `Atom`: position (`xred`, `xcart`, `absxcart`), velocities, forces, charge, magnetization
  - `Lattice`: simulation cell (vectors, volume, angles), atom list, and all thermodynamic scalars (energy, pressure, temperature, stress tensor, etc.). Also contains `Elements2rest()` (element → atomic number / mass lookup) and coordinate-transformation methods.
  - `Lattice.__init__`'s `acell`/`angles`/`rprim`/`rprimd`/`gprimd`/`stress` used to be mutable default arguments, so every `Lattice()` created without explicit values silently shared (and could corrupt) the same list — fixed to build a fresh list per instance. No script constructs `Lattice()`/`Atom()` with explicit keyword arguments, so this was a pure latent-bug fix.
- **`umd_processes_fast.py`** — Python wrapper around `c_UMDprocess.dylib` (ctypes). Provides fast reading of UMD snapshots and bond-definition routines used by most analysis scripts. This is the main I/O layer for `.umd.dat` files.

### UMD File Format (`.umd.dat`)
Human-readable ASCII: one header section per snapshot followed by per-atom lines. Fields include `timestep`, `time`, `InternalEnergy`, `Enthalpy`, `Temperature`, `Pressure`, `StressTensor`, `acell`, `rprim`, and per-atom `xred`/`xcart`/`vels`/`forces`. Parsers write this format; analysis scripts read it.

### Parsers (simulation code → `.umd.dat`)
| Script | Source format |
|---|---|
| `VaspParser.py` | VASP OUTCAR |
| `QBoxParser.py` | QBox XML output |
| `LAMMPSParser2umd.py` | LAMMPS dump files |

### Analysis Scripts
Each script imports `crystallography` and `umd_processes_fast`; most also load a dedicated C library for the heavy computation.

**Structural:**
- `gofr_umd.py` — pair distribution functions (uses `c_gofr_new_alt`)
- `analyze_gofr.py` / `analyze_gofr_semi_automatic.py` — post-process gofr output, find bond lengths
- `bonding_umd.py` — bond analysis (uses `c_bonds_full`)
- `speciation_and_angles.py` — chemical speciation and interatomic angles (uses `c_clusters`)
- `lifetimes.py` — bond/cluster lifetimes from population files

**Transport:**
- `msd_umd.py` / `msd_all_umd.py` / `msd_clusters_umd.py` — mean-square displacements and diffusion coefficients (uses `c_msd`)
- `viscosity_umd.py` — shear viscosity via stress autocorrelation
- `analyze_msd.py` — post-process MSD output

**Thermodynamic:**
- `averages.py` / `fullaverages.py` — thermodynamic averages with plots
- `vibr_spectrum_umd_fast.py` — vibrational spectra via velocity autocorrelation (uses `c_autocorrelation_vib`)

**Utilities:**
- `umd2out.py` — convert UMD snapshots to XYZ or POSCAR
- `build_supercell.py` — build an `nx × ny × nz` supercell from a unit cell read from a VASP5-style POSCAR/CONTCAR file (e.g. one exported by VESTA), in reduced or cartesian coordinates. Replication is done in reduced coordinates (exact for any cell shape), writing both `.vasp` (POSCAR) and `.xyz` output. Kept as a standalone script rather than folded into `insert_umd.py`, since it reads a different input format (POSCAR vs. the `molecules.dat` format) and serves a different purpose (replicating a periodic cell vs. inserting molecules into a UMD trajectory); a shared "cell manipulation" library backing both may be considered later.
- `insert_umd.py` — insert molecules into a UMD trajectory (or an empty box); output format is selectable with `-t` (1 = xyz, 2 = vasp poscar, 3 = umd file). `-s -1` (default) builds an empty cubic box of side `-a`; `-s 0` inserts into the last snapshot of the UMD file given by `-f`; `-s Ksteps>0` inserts into every Ksteps-th snapshot. `insert_umd_xred.py` was a broken, unused duplicate and has been removed.
- `stat-concentrate.py` — aggregate statistics from multiple speciation runs
- `check_overlap.py` — detect atomic overlaps
- `crystallography.py` also contains `Elements2rest()`, used widely for element → mass/atomic-number lookup

### Data Flow
```
OUTCAR / QBox / LAMMPS dump
        ↓  Parser scripts
    <name>.umd.dat
        ↓  Analysis scripts (read via umd_processes_fast / C libs)
    <name>.gofr.dat / .bonding.dat / .msd.dat / .vibr.dat / …
        ↓  Analyze/plot scripts
    summary files + matplotlib figures
```

## Naming Conventions
Output files inherit the UMD filename as prefix: `<base>.umd.dat` → `<base>.gofr.dat`, `<base>.bonding.dat`, `<base>.r0.popul.dat`, etc. Test data in `data/` follows this convention with `glycine.*` files.
