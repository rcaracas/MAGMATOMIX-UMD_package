# UMD_package
Post-processing package that performs analysis of structural, transport, and thermodynamic properties from ab initio and classical molecular dynamics simulations.

UMD stores the atomic trajectories and the thermodynamic information of a run in one ASCII format, the `.umd.dat` file. Parsers write it from VASP, ABINIT, CP2K, QBox and LAMMPS outputs; the analysis scripts read it.

## Quick start
Python 3 with numpy, scipy and matplotlib, and a C compiler (gcc). Compile the four C libraries once, on macOS, Linux or Windows (the commands themselves are in `README.txt`):
```bash
cd src
python3 compile_libraries.py
```
**Windows**: see [`WINDOWS.md`](WINDOWS.md); every release has the libraries compiled for Windows (`umd-c-libraries-windows.zip`), so no compiler is needed. The package is tested on Linux, Windows and macOS at every change (GitHub Actions).
Then, for example:
```bash
python3 VaspParser.py -f OUTCAR            # writes OUTCAR.umd.dat
python3 check_umd.py -f OUTCAR.umd.dat     # quality checks of the run
python3 gofr_umd.py -f OUTCAR.umd.dat      # pair distribution functions
python3 msd_umd.py -f OUTCAR.umd.dat       # mean square displacements and diffusion coefficients
python3 averages.py -a -f OUTCAR.umd.dat   # thermodynamics with errors
python3 umd_workflow.py -f OUTCAR.umd.dat  # or all the analyses at once (-g structural, -t transport, -d thermodynamics), summary in umd_workflow.summary.dat
```
Every script prints its options with `-h`. The list of the scripts and the changelog are in `README`; `documentation/umd-quick-help.pdf` has one card per script.

An example run (a pyrolite melt at 3000 K) with a tutorial and the results to expect: `examples/README.md`, the file in the release [examples](https://github.com/rcaracas/MAGMATOMIX-UMD_package/releases/tag/examples).

`validation/run_validation.py` runs every script on reference runs and on synthetic files with exact answers, and compares the results with committed references (see `validation/README`; the reference runs themselves are not distributed).

## Citation
Caracas, R. et al. Analyzing Melts and Fluids from Ab Initio Molecular Dynamics Simulations with the UMD Package. J. Vis. Exp. (175), e61534, doi:10.3791/61534 (2021)

## License
GNU GPL v3 and CeCILL v2 (mutually compatible). Authors: see `AUTHORS`.
