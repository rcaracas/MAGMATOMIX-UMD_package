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
python3 src/averages.py -f <file.umd.dat> [-p Pattern] [-s SkipSteps] [-l] [-v]
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
- `msd_umd.py` — mean-square displacements per element (`-m elements`) or per atom and element (`-m atoms`), written to `<base>.msd.dat` (`#` comment lines, then a `time_(fs)` column header). Computed in numpy with the FFT algorithm (Calandrini/nMoldyn), averaged over all time origins, for lags up to half the trajectory; this replaced the former `c_msd` C library (removed), which took ~100x more CPU time and only used origins from the first half. `-a cart` / `-a cryst` (`Auto` is a synonym) / `-a [x1,y1,z1,...]` add the MSD projected on unit vectors along x,y,z / the a,b,c vectors of the first snapshot's cell / arbitrary axes, as extra columns `El_x`, `El_a`, `El_1`… For a non-orthogonal cell the a,b,c components do not sum to the total. `-b` drops lags shorter than the ballistic regime, `-v` writes every v-th lag, `-z` is accepted but ignored (kept because `UMD_GUI.py` and `MSDdistrib.py` pass it). `main()` returns `[msdfile, MSD, Instants, elements, Axes]` with `MSD[species][0]` the total and `MSD[species][1+k]` axis k, normalized per atom; `UMD_GUI.py` and `MSDdistrib.py` rely on this. `msd_all_umd.py` was an older pure-Python duplicate and has been removed.
- `msd_clusters_umd.py` — MSD of the mass-weighted centre of mass of each cluster listed in a `.popul.dat` file from `speciation_and_angles.py`, one `<popul base>__<cluster>.msd.dat` file per cluster type (one column per cluster instance, plus their mean). Reuses `msd_umd.msd_fft`. The popul Begin/End steps are converted to UMD snapshots using the file's lifetime column, so a speciation run with sampling > 1 is handled.
- `viscosity_umd.py` — shear viscosity via stress autocorrelation
- `analyze_msd.py` — fits every column of all `*.msd.dat` files below the current folder and writes `diffusivities.txt` per folder; D = slope/6 for total-MSD columns and slope/2 for single-axis columns (any name containing `_`).

**Thermodynamic:**
- `averages.py` — mean and stdev of any numerical field of one UMD file. Default: Pressure, Temperature and InternalEnergy on one tab-separated line (`file #PstdTstdEstd nsteps P stdev_P T stdev_T E stdev_E`, a format that may be parsed by downstream scripts — keep it stable). `-p <Pattern>` averages a single field instead; multi-valued fields (`StressTensor`, `acell`, …) are averaged component by component. `-l` gives the long per-component output (mean, variance, stdev, unit), `-v` opens a plot with Previous/Next buttons stepping through the quantities/components one at a time, each on its own y-scale (default: no display). `-i` and `-s` are synonyms for the number of initial snapshots to skip (default 0 = all). `averages_detailed.py` and `averages_visual.py` were merged into it and removed; `averages_forGUI.py` is a separate copy kept only because `UMD_GUI.py` imports it.
- `fullaverages.py` — batch thermodynamics over every `*.umd.dat` in the current folder (blocking-method errors, bootstrap Cv), written to `fullthermo.dat`
- `vibr_spectrum_umd_fast.py` — vibrational spectra via velocity autocorrelation (uses `c_autocorrelation_vib`). Writes `<base>.vels.scf.dat` (VACF), `<base>.vibr.dat` (per-species + total vDOS and diffusion coefficients), and `<base>.vibr_win.dat` (a raised-cosine-tapered version of the total vDOS, on the normalized VACF scale rather than the physical units of `.vibr.dat`, meant to reduce spectral ringing from truncating the correlation function). `vibr_hann.py` was an in-development fork with the same windowing feature (plus a broken `-tau` flag and a bug that always wrote the windowed output to the current directory as a fixed `vib_win.dat` instead of alongside the input file); its working parts were merged in here and it has been removed.

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
