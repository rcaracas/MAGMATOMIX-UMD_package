# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

UMD (Universal Molecular Dynamics) is a Python package for post-processing ab initio and classical molecular dynamics simulations. It works around a central ASCII file format (`.umd.dat`) that stores atomic trajectories and thermodynamic data extracted from simulation codes (VASP, QBox, LAMMPS).

## Building C Extensions

Performance-critical routines are implemented as C shared libraries loaded at runtime via `ctypes`. They **must be compiled before running any analysis script**. All commands run from `src/`:

**macOS:**
```bash
gcc -dynamiclib -o c_UMDprocess.dylib c_UMDprocess.c
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

`VaspParser.py` unwraps the positions (`absxcart`) at every MD step from `-i` on, even those not written with `-s`, so no crossing of a cell face is missed. It finds each crossing as a whole-number change of the displacement in reduced coordinates, which works for any cell shape. The `timestep` it writes is the time between written snapshots (`-s` × POTIM), as the analysis scripts expect. On-the-fly ML learning steps print the positions twice (ML and DFT blocks); the previous positions are taken only from the first block of each ionic step.

### Analysis Scripts
Each script imports `crystallography` and `umd_processes_fast`; most also load a dedicated C library for the heavy computation.

**Structural:**
- `gofr_umd.py` — pair distribution functions (uses `c_gofr_new_alt`)
- `analyze_gofr.py` / `analyze_gofr_semi_automatic.py` — post-process gofr output, find bond lengths
- `bonding_umd.py` — bond analysis (uses `c_bonds_full`)
- `speciation_and_angles.py` — chemical speciation and interatomic angles (uses `c_clusters`)
- `lifetimes.py` — bond/cluster lifetimes from population files

**Transport:**
- `msd_umd.py` — mean-square displacements per element (`-m elements`) or per atom and element (`-m atoms`), written to `<base>.msd.dat` (`#` comment lines, then a `time_(fs)` column header). Computed in numpy with the FFT algorithm (Calandrini/nMoldyn), averaged over all time origins, for lags up to half the trajectory; this replaced the former `c_msd` C library (removed), which took ~100x more CPU time and only used origins from the first half. `-a cart` / `-a cryst` (`Auto` is a synonym) / `-a [x1,y1,z1,...]` add the MSD projected on unit vectors along x,y,z / the a,b,c vectors of the first snapshot's cell / arbitrary axes, as extra columns `El_x`, `El_a`, `El_1`… For a non-orthogonal cell the a,b,c components do not sum to the total. `-b` drops lags shorter than the ballistic regime, `-v` writes every v-th lag, `-z` is accepted but ignored (kept because `UMD_GUI.py` and `MSDdistrib.py` pass it). It also fits D (slope/6 total, slope/2 per axis; `analyze_msd.calculation_diffusivities` is imported, so both tools share one fit) between `-b` and `-e` (default 0.75 of the last lag of a block, i.e. 7.5 % of the trajectory for 5 blocks; with `-n` 0/1 the 5-block length is used; fitting to long lags made D ~3× noisier in random-walk tests) and writes `<base>.diffusion.dat` (D, error, R², alpha, with a warning when alpha is outside 0.9–1.1). Errors come from `-n` consecutive blocks (default 5; MSDs written to `<base>.msd.blocks.dat`, which `analyze_msd.py` also uses): error = std(D_blocks)·√(t_end,full / (n·t_end,block)), because the variance of a fitted D grows with the fit-window end and falls with trajectory length. The plain std/√n would underestimate it 2–3× when the full fit is longer than the blocks. Calibrated on synthetic random walks (reported/actual scatter 0.83–1.08); the linear-regression slope error is ~200× too small. `main()` returns `[msdfile, MSD, Instants, elements, Axes]` with `MSD[species][0]` the total and `MSD[species][1+k]` axis k, normalized per atom; `UMD_GUI.py` and `MSDdistrib.py` rely on this. `msd_all_umd.py` was an older pure-Python duplicate and has been removed.
- `msd_clusters_umd.py` — MSD of the mass-weighted centre of mass of each cluster listed in a `.popul.dat` file from `speciation_and_angles.py`, one `<popul base>__<cluster>.msd.dat` file per cluster type (one column per cluster instance, plus their mean). Reuses `msd_umd.msd_fft`. The popul Begin/End steps are converted to UMD snapshots using the file's lifetime column, so a speciation run with sampling > 1 is handled.
- `viscosity_umd.py` — shear viscosity by Green–Kubo from the UMD `StressTensor`; the only viscosity tool (the viscosity modes of `spectra_umd.py` and `spectra_umd_corrected.py` were removed). By default uses the 5 independent components of the traceless stress tensor N = σ − pI (Daivis–Evans weights: 1/5 per off-diagonal, 2/15 per diagonal, the diagonal running integrals ×3/4); `-c offdiag` uses only xy, yz, zx. The mean of the whole run is removed once (blocks included), and η(t) is corrected for the bias this causes: η = η_raw/(1 − 2 t_end/T_run). FFT autocorrelation, cumulative trapezoid integral. η is read at the end of the integral, `-e` (steps; default and max = half a block, 1/8 of the run with the default `-n 4` blocks); the error is the standard error over the blocks. A warning flags η(t) still rising by > 2 errors over the second half (not converged → lower bound). Automatic windows (Sokal) and the double-exponential time-decomposition fit (Zhang–Otani–Maginn) were tried and rejected: melts have a small but long stress-autocorrelation tail that Sokal cuts (hydrous andesite: 2.0 vs ~5.5 mPa·s), and the fit was unstable with a few blocks. Synthetic tests (known η): unbiased within noise, error ~0.7–0.85 of the real scatter. `-i` skips initial snapshots, `-v` plots autocorrelations and running viscosities. Output `<base>.visc.dat`.
- `analyze_msd.py` — fits all element/atom `*.msd.dat` files (from `msd_umd.py`) below the current folder and writes one `diffusivities.txt` per folder. Cluster MSD files are recognized by their `cl0 … mean` header and skipped. Columns: `file`, `fit_start(fs)`, `fit_end(fs)`, `error` (blocks/fit), then 6 per species (D, D_stdev, R², slope, intercept, alpha). The species columns are the union over all files of the folder, `nan` where a file lacks one. D = slope/6 for total-MSD columns, slope/2 for single-axis columns (any name containing `_`). Fit window `-s` to `-e` in fs (default end 0.75 of the last block time if a blocks file exists, else 0.15 of the last MSD time — the same point as for 5 blocks). D_stdev is the block error when `<base>.msd.blocks.dat` exists, otherwise the (underestimated) regression error. alpha is the log-log slope of the MSD over the window (MSD ~ t^alpha): a warning is printed outside 0.9–1.1 (sub-diffusive / ballistic). `-v` shows each MSD with its fit and shaded window, with Previous/Next buttons (`plot_fits`, matplotlib imported lazily so `msd_umd.py`, which imports this module, does not need it).
- `analyze_msd_clusters.py` — the same for the cluster MSD files of `msd_clusters_umd.py`, one `diffusivities_clusters.txt` per folder (one row per cluster file: instances, instances used, window, D, error, R², slope, intercept, alpha). Instances live for different times, so only those living through the whole fit window are averaged and fitted; default end = 0.75 of the median instance lifetime; error = std of the per-instance D / √n. Short-lived species often stay ballistic (alpha ≈ 2) — check alpha. Reuses the fit, alpha and plot functions of `analyze_msd.py`; `-v` as there.

**Thermodynamic:**
- `averages.py` — mean and stdev of any numerical field of one UMD file. Default: Pressure, Temperature and InternalEnergy on one tab-separated line (`file #PstdTstdEstd nsteps P stdev_P T stdev_T E stdev_E`, a format that may be parsed by downstream scripts — keep it stable). `-p <Pattern>` averages a single field instead; multi-valued fields (`StressTensor`, `acell`, …) are averaged component by component. `-l` gives the long per-component output (mean, variance, stdev, unit), `-v` opens a plot with Previous/Next buttons stepping through the quantities/components one at a time, each on its own y-scale (default: no display). `-i` and `-s` are synonyms for the number of initial snapshots to skip (default 0 = all). `averages_detailed.py` and `averages_visual.py` were merged into it and removed; `averages_forGUI.py` is a separate copy kept only because `UMD_GUI.py` imports it.
- `fullaverages.py` — batch thermodynamics over every `*.umd.dat` in the current folder (blocking-method errors, bootstrap Cv), written to `fullthermo.dat`
- `vibrations_umd.py` — velocity autocorrelation functions (VACF), vibrational densities of states (vDOS) and diffusion coefficients from the UMD velocities, per element and, with `-a cart` / `-a cryst` / `-a [x1,y1,z1,...]` (as in `msd_umd.py`), per element along each axis (not per atom). The only vibration tool: it replaced `vibr_spectrum_umd_fast.py` (its C library `c_autocorrelation_vib`, O(N²), took ~10x more CPU time), `spectra_umd.py` (broken with current scipy) and `spectra_umd_corrected.py`, all removed. VACFs by FFT over all time origins. vDOS = 2 m/(kB T) × the cosine transform of the VACF, mass-weighted, in states per cm⁻¹ (the old script wrote per fs⁻¹ against a cm⁻¹ axis), so the vDOS of an element integrates to 3 N_El · T_kin/T (N_El per axis). T is the mean `Temperature` of the UMD file (`-t` overrides it; the old default was 5000 K); the kinetic temperature of the velocities is also printed. With `-i 0` the first snapshot is dropped when its ⟨v²⟩ is 0 or more than 10× from the median of the next 100: the parsers have no previous step there (current `VaspParser.py` writes 0; older versions wrote the positions, as in the pyrolite and glycine test files, which gave a spike at t=0 in the VACF and a flat vDOS up to Nyquist). `-w` is the VACF length used for the spectra (steps, default half a block, max half the run), `-p hann` (default) / `none` a one-sided Hann taper applied to all spectra (suppresses truncation ringing), `-g` a Gaussian smoothing in cm⁻¹ (reflection at 0, exact for an even spectrum), zero-padding ×4 for a finer frequency grid. D (Green–Kubo) is the running integral 1/3∫⟨v(0)·v(t)⟩ (∫⟨v_k v_k⟩ per axis), read at `-e` (steps; default and max half a block), with the standard error over `-n` blocks (default 4, min 3), and a warning when D changes by > 2 errors between half the integral and its end (solids: D oscillates around 0, which the block error does not include). Checks: the MSD of `msd_umd.py` equals 2∫∫VACF within 1 % up to 1 ps on the hydrous andesite; synthetic OU velocities give unbiased D with error/scatter 0.98–1.03; synthetic oscillators give the right peaks and exact normalization. D from the VACF and from the MSD fit differ when the system is not yet diffusive over the windows (they use different times). Outputs: `<base>.vels.scf.dat` (normalized VACFs, mass-weighted `Average`, running `D_El`), `<base>.vibr.dat` (one header line `Frequency(cm^-1)`, `VDos_El`…, `Total_DOS`, the format read by `UMD_GUI.py`; keep it), `<base>.vibr.axes.dat` (with `-a`), `<base>.vacf.diffusion.dat`. `-v` plots the VACFs and vDOS with their error band (the running D is in `<base>.vels.scf.dat`). vDOS errors (`<base>.vibr.err.dat`: each vDOS, also along axes, followed by `_err`; `main` returns them last): spread of the vDOS of min(`-n`, run/window) blocks at least one window long, scaled to the whole run by the variance of a Blackman–Tukey spectrum on the unbiased VACF, 1/L + Σ_{t≥1} w(t)²/(L−t) (lag t has variance ~1/(L−t)); an earlier L·Σw²/(L−t)² overweighted the long lags of short blocks and gave errors ~0.8× too small, and 20× too small without taper. With Hann the relative variance var/g² is averaged over the FWHM (~24000/T cm⁻¹, printed), without taper the variance itself. Checks: synthetic Gaussian velocities give error/scatter median 0.92–1.03, 0.55–0.75 at the top of peaks whose lifetime (~1 ps) is comparable to the blocks; real-data halves test std(z) 0.9–1.2 for pyrolite and hydrous andesite, ~1.8 for ZnS at 1000 K (long-lived phonons), 1.7 for glycine (not stationary: ⟨v²⟩ rises 26 % over the run while its `Temperature` field stays 300.0). Three checks are printed in the terminal (`check_zero_velocities`, `check_temperatures`, `check_halves`): (0) snapshots whose velocities are all exactly 0, the signature of `VaspParser.py` before 48975f1 on on-the-fly ML learning steps (`data/al2o3.umd.dat`, hexagonal corundum + 1 He: 151 such snapshots, 4.5 % of block 1, which made block 1 ~60 K colder and triggered the drift warning; its absxcart is also wrong, from the missed crossings and the old acell/2 test in a hexagonal cell — re-parse the OUTCAR or use `repair_umd.py`; `data/al2o3.new.umd.dat` is the re-parsed version); (1) the kinetic temperature of each element per `-n` block and over the run, with an error from the integrated autocorrelation time of the linearly detrended series (`mean_error`; without detrending a drift inflates τ and hides itself: glycine H τ 676 → 6 steps); warnings for a block deviating by > 4 block errors (drift) and for an element differing from all atoms by > 4 errors (no equipartition). Glycine: O 408→312 K and H 227→297 K over the run while the total stays ~301 K (the earlier "⟨v²⟩ rises 26 %" was the unweighted mean, dominated by H); ZnS: S 968 K vs Zn 1024 K (±4), steady. Melts pass (max 2.2 block, 3.2 equipartition). (2) the vDOS of the two halves, same window, each with its block errors: rms(diff)/rms(√(e₁²+e₂²)) per element and total, warning above 1.5 for the total (melts 0.86–0.96; glycine 1.6, ZnS 2.3 with `-n 4` and only 2.1 with `-n 3`, so fewer blocks do not fix crystals: the `-n` help says to scale the band by this ratio; pyrolite with `-p none` 1.7); skipped when the window exceeds 1/4 of the run. Non-orthogonal cell check on al2o3 with `-a cryst`: the equivalent a and b axes give the same vDOS (Al integrals 92.8/94.5, mean frequency 550/546 cm⁻¹), c differs (518 cm⁻¹); halves ratio 1.19. `UMD_GUI.py` calls `main(["-f", file])`, plus `-t` only when its temperature field is filled.

**Utilities:**
- `umd2out.py` — convert UMD snapshots to XYZ or POSCAR
- `build_supercell.py` — build an `nx × ny × nz` supercell from a unit cell read from a VASP5-style POSCAR/CONTCAR file (e.g. one exported by VESTA), in reduced or cartesian coordinates. Replication is done in reduced coordinates (exact for any cell shape), writing both `.vasp` (POSCAR) and `.xyz` output. Kept as a standalone script rather than folded into `insert_umd.py`, since it reads a different input format (POSCAR vs. the `molecules.dat` format) and serves a different purpose (replicating a periodic cell vs. inserting molecules into a UMD trajectory); a shared "cell manipulation" library backing both may be considered later.
- `repair_umd.py` — repairs a UMD file written by `VaspParser.py` before 48975f1 (23 Oct 2026) when the OUTCAR is gone: recomputes the velocities, the unwrapped positions (abs.diff.) and xred from the wrapped xcart (the OUTCAR values), with the same arithmetic and summation order as the current parser, and copies everything else; writes `<base>.repaired.umd.dat` (`-o`). Checks: `data/al2o3.umd.dat` repaired is byte-identical to the re-parsed `data/al2o3.new.umd.dat` (152 zero-velocity snapshots, 14561 crossings); a file from the current parser (ZnS) comes out byte-identical; old pyrolite: only the first-snapshot velocities change, its abs.diff. differ by a constant cell vector (old (0.5,0.5,0.5) start), so its MSD was right. The first snapshot gets v = 0. With snapshots s > 1 steps apart (old files: Δtime/timestep = s) the velocities can only be means over s steps (warning); it warns too when most nonzero velocities change, e.g. a current-parser file written with `-s` > 1, whose velocities were right.
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
