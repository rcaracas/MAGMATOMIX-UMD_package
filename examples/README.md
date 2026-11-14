# Example runs

The example runs are UMD files, attached to the GitHub release **examples** (they are too large for the repository; their links do not change from one version to the next):

https://github.com/rcaracas/MAGMATOMIX-UMD_package/releases/tag/examples

| File | Run | Size |
|---|---|---|
| `pyrolite.a13.T3000K.umd.zip` | pyrolite melt, 153 atoms (Na Ca₂ Fe₄ Mg₃₀ Al₃ Si₂₄ O₈₉), cubic cell of 13 Å (2.445 g/cm³), 3000 K, 3000 snapshots of 1 fs (3 ps), ab initio molecular dynamics (VASP, NVT), with velocities, forces and stresses | 24 MB (64 MB unzipped) |

sha256 of `pyrolite.a13.T3000K.umd.zip`: `42110d4f597106fcd58ba9d4f3598c7559442d49432ef94c41fee81c6ac63bfc`

## Tutorial

Download and unzip the example in a folder of your own (here `work`, next to the package), compile the C libraries (see `README.txt`), then run the scripts from that folder. They write their outputs next to the UMD file.

```bash
cd work
python3 ../MAGMATOMIX-UMD_package/src/check_umd.py -f pyrolite.a13.T3000K.umd.dat         # quality checks of the run
python3 ../MAGMATOMIX-UMD_package/src/averages.py -a -i 500 -f pyrolite.a13.T3000K.umd.dat # density, P, T, E, Cv with errors, 500 first snapshots skipped
python3 ../MAGMATOMIX-UMD_package/src/gofr_umd.py -f pyrolite.a13.T3000K.umd.dat          # pair distribution functions
python3 ../MAGMATOMIX-UMD_package/src/analyze_gofr.py -f pyrolite.a13.T3000K.gofr.dat     # bond lengths (cutoffs) and coordination numbers
python3 ../MAGMATOMIX-UMD_package/src/bonding_umd.py -f pyrolite.a13.T3000K.umd.dat -c pyrolite.a13.T3000K.maximum_bonds.inp
python3 ../MAGMATOMIX-UMD_package/src/speciation_and_angles.py -f pyrolite.a13.T3000K.umd.dat -b pyrolite.a13.T3000K.bonding.dat -c Si,Al -a O -r 1 -t 1
python3 ../MAGMATOMIX-UMD_package/src/msd_umd.py -f pyrolite.a13.T3000K.umd.dat           # mean square displacements and diffusion coefficients
python3 ../MAGMATOMIX-UMD_package/src/umd2out.py -f pyrolite.a13.T3000K.umd.dat -r xyz -s 10   # trajectory for OVITO or VMD
python3 ../MAGMATOMIX-UMD_package/src/umd_workflow.py -f pyrolite.a13.T3000K.umd.dat      # or all the analyses at once
```

Each script runs in 1–2 s on a laptop, the whole workflow (14 steps) in about 15 s. What you should find:

- `averages.py -a -i 500`: density 2.445 g/cm³, P 1.29 ± 0.16 GPa, T 2999.5 ± 16 K, Cv 3.46 ± 0.22 N kB
- `analyze_gofr.py` (`.coordinations.inp`): Si-O 4.00, Al-O 4.09, Fe-O 4.14, Mg-O 4.93, Na-O 5.90, Ca-O 7.72; Si-O cutoff 2.42 Å
- `speciation_and_angles.py` (`.r1.stat.dat`): SiO₄ 69 % of the time, then free O, AlO₄, SiO₅, SiO₃; mean O-Si-O angle in SiO₄ 108.5°
- `msd_umd.py` (`.diffusion.dat`): D(Si) 6.3 ± 0.6 × 10⁻⁹ m²/s, D(O) 10.6 ± 0.9 × 10⁻⁹ m²/s
- `check_umd.py`: one warning, the energy drift (22 jumps of the conserved energy: electronic steps not converged in this run)

Every script prints its options with `-h`; `documentation/umd-quick-help.pdf` has one card per script.
