# UMD_package
Post-processing package that performs analysis of structural, transport, and thermodynamic properties from ab initio and classical molecular dynamics simulations.

UMD stores the atomic trajectories and the thermodynamic information of a run in one ASCII format, the `.umd.dat` file. Parsers write it from VASP, QBox and LAMMPS outputs; the analysis scripts read it.

## Quick start
Python 3 with numpy, scipy and matplotlib, and a C compiler. Compile the C libraries in `src/` (macOS; Linux and Windows in `README.txt`):
```bash
cd src
gcc -O2 -dynamiclib -o c_UMDprocess.dylib c_UMDprocess.c
gcc -dynamiclib -o c_bonds_full.dylib c_bonds_full.c
gcc -dynamiclib -o c_clusters.dylib c_clusters.c
gcc -O2 -dynamiclib -o c_gofr.dylib c_gofr.c
```
Then, for example:
```bash
python3 VaspParser.py -f OUTCAR            # writes OUTCAR.umd.dat
python3 check_umd.py -f OUTCAR.umd.dat     # quality checks of the run
python3 gofr_umd.py -f OUTCAR.umd.dat      # pair distribution functions
python3 msd_umd.py -f OUTCAR.umd.dat       # mean square displacements and diffusion coefficients
python3 averages.py -a -f OUTCAR.umd.dat   # thermodynamics with errors
```
Every script prints its options with `-h`. The list of the scripts and the changelog are in `README`.

## Citation
Caracas, R. et al. Analyzing Melts and Fluids from Ab Initio Molecular Dynamics Simulations with the UMD Package. J. Vis. Exp. (175), e61534, doi:10.3791/61534 (2021)

## License
GNU GPL v3 and CeCILL v2 (mutually compatible). Authors: see `AUTHORS`.
