# MAGMATOMIX/UMD on Windows

The package runs natively on Windows 10 and 11: Python, and four small C libraries that come compiled with every release (no compiler needed). The commands below are typed in a **PowerShell** or **Command Prompt** window.

## 1. Python

Install Python 3.12 or 3.13 from [python.org](https://www.python.org/downloads/windows/). In the first window of the installer, tick **"Add python.exe to PATH"**. Then, in a new PowerShell window:

```
python --version
python -m pip install numpy scipy matplotlib
```

## 2. The package

From the [releases page](https://github.com/rcaracas/MAGMATOMIX-UMD_package/releases), take the release you want and download two files:

- **Source code (zip)**: the package. Unzip it, for example in your `Documents` folder, and rename the folder `MAGMATOMIX-UMD_package` (GitHub adds the version to its name), so that the commands of the tutorial work as written.
- **umd-c-libraries-windows.zip**: the compiled C libraries. Unzip its four files (`c_UMDprocess.dll`, `c_bonds_full.dll`, `c_clusters.dll`, `c_gofr.dll`) into the `src` folder of the package, next to the `.py` scripts.

If you have the gcc compiler of MinGW-w64 (e.g. from [MSYS2](https://www.msys2.org) or [WinLibs](https://winlibs.com)), you can compile the libraries yourself instead: `python compile_libraries.py` in the `src` folder.

Check that everything is in place:

```
cd Documents\MAGMATOMIX-UMD_package\src
python check_umd.py -h
```

It prints the banner of the package and the options of the script.

## 3. The example run and the tutorial

Download [pyrolite.a13.T3000K.umd.zip](https://github.com/rcaracas/MAGMATOMIX-UMD_package/releases/download/examples/pyrolite.a13.T3000K.umd.zip) and unzip it in a folder `work` next to the package (`Documents\work`). Then follow the tutorial of `examples\README.md`, with two changes: `python` instead of `python3`, and `\` instead of `/` in the paths:

```
cd Documents\work
python ..\MAGMATOMIX-UMD_package\src\check_umd.py -f pyrolite.a13.T3000K.umd.dat
python ..\MAGMATOMIX-UMD_package\src\gofr_umd.py -f pyrolite.a13.T3000K.umd.dat
python ..\MAGMATOMIX-UMD_package\src\umd_workflow.py -f pyrolite.a13.T3000K.umd.dat
```

The last command runs all the analyses (about 20 s) and writes `umd_workflow.summary.dat`; compare your numbers with those of the tutorial. Scripts with `-v` open their plots in a window.

## Problems

- `Could not find module 'c_UMDprocess.dll'` (or another `.dll`): the four `.dll` files are not in the `src` folder of the package.
- `python` opens the Microsoft Store, or is not found: Python was installed without "Add python.exe to PATH"; run the installer again and tick it (or use `py` instead of `python`).
- A file you downloaded is blocked: right click, Properties, tick **Unblock**.

Please report any other problem as an issue on [GitHub](https://github.com/rcaracas/MAGMATOMIX-UMD_package/issues), with the command you typed, what it printed, your Windows and Python versions (`python --version`).

The Windows Subsystem for Linux (WSL) also works: in a PowerShell opened as administrator, `wsl --install`, restart, open *Ubuntu*, and follow the Linux instructions of `README.txt`.
