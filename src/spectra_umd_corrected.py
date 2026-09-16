#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spectra_umd_fixed.py
Rewritten and corrected version of your UMD spectra tool.

- Robust handling of snapshots assumed in flat arrays: [vx0, vy0, vz0, vx1, vy1, vz1, ...]
- VACF computed with multiple time origins using FFT convolution (fast) and unbiased normalization.
- Hann taper applied before rFFT to remove dependence on arbitrary truncation/window_size.
- Produces VDOS (power spectrum) in frequency units (Hz, THz) and wavenumbers (cm^-1).
- Compatible with UMD read_values/read_stresses_4visc routines (umd_processes_fast).
"""

import sys
import getopt
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.ndimage import gaussian_filter1d
# Import UMD helper routines (assumed available)
import umd_processes_fast as umdpf
import crystallography as cr

# ---------------------------
# Utility & physics constants
# ---------------------------
HBAR = 1.054571817e-34  # J*s
KB = 1.3806503e-23      # J/K
EV_TO_J = 1.602176634e-19
FS_TO_S = 1e-15
CM1_TO_HZ = 3.3333333333333335e-11  # 1 cm^-1 in Hz



# ---------------------------
# Low-level helpers
# ---------------------------
def _safe_snapshot_to_array(snapshot):
    """
    Convert snapshot into a 1D numpy array. We support:
      - list/tuple/np.ndarray -> returned as np.array
      - objects with attributes: values, velocity, velocities, vel, data
      - dict-like with keys 'velocities', 'velocity', 'values', 'vel'
    If none match, raises ValueError.
    """
    if isinstance(snapshot, (list, tuple, np.ndarray)):
        return np.asarray(snapshot, dtype=float)
    for attr in ("values", "velocity", "velocities", "vel", "data"):
        if hasattr(snapshot, attr):
            return np.asarray(getattr(snapshot, attr), dtype=float)
    if hasattr(snapshot, "get") and callable(snapshot.get):
        for key in ("velocities", "velocity", "values", "vel"):
            if key in snapshot:
                return np.asarray(snapshot.get(key), dtype=float)
    raise ValueError("Unrecognized snapshot format; cannot extract numeric array.")


def _next_pow2(x):
    return 1 << (int(x) - 1).bit_length()


# ---------------------------
# VACF via FFT convolution
# ---------------------------
def autocorr_fft_unbiased(x):
    """
    Compute the unbiased autocorrelation for 1D real array x using FFT convolution.
    Returns autocorrelation for lags 0..(n-1) and the array of denominators (n - tau).
    """
    x = np.asarray(x, dtype=float)
    n = x.shape[0]
    # subtract mean (stationarity assumption helps)
    x = x - x.mean()
    # pad to 2*n for circular convolution
    m = 1 << ((2 * n - 1).bit_length())
    fx = np.fft.rfft(x, n=m)
    acf_full = np.fft.irfft(fx * np.conjugate(fx), n=m)
    acf = acf_full[:n]
    # unbiased normalization: divide by (n - tau)
    denom = np.arange(n, 0, -1)
    return acf, denom


def compute_vacf_from_velocities(velocities):
    """
    Compute normalized VACF from velocities array with shape (nsteps, natom, 3).
    Returns vacf (length nsteps) normalized by vacf[0] so vacf[0] == 1.
    Uses unbiased averaging over time origins and atoms & directions.
    """
    velocities = np.asarray(velocities, dtype=float)
    nsteps, natom, _ = velocities.shape
    if nsteps < 2:
        raise ValueError("Need at least 2 time frames to compute VACF.")

    # Convert to shape (nsteps, natom*3) for easier per-component handling
    flat = velocities.reshape(nsteps, natom * 3)

    # We'll compute autocorr for each of the M = natom*3 series using FFT,
    # then average them with unbiased normalization.
    M = flat.shape[1]
    acfs = np.zeros((M, nsteps))
    denoms = np.zeros((M, nsteps))
    for j in range(M):
        acf_j, denom_j = autocorr_fft_unbiased(flat[:, j])
        acfs[j, :] = acf_j
        denoms[j, :] = denom_j

    # average across components with correct weighting: sum(acf)/sum(denom)
    # but we want average of expectation <v(0)v(t)> so perform:
    # vacf(t) = (1/M) * sum_j ( acf_j(t) / (n - t) ) * (n - t)
    # Simpler: average the unbiased estimates: acf_j / denom_j * denom_common and then divide by M.
    unbiased = acfs / denoms  # shape M x nsteps
    vacf = unbiased.mean(axis=0)

    # scale back by vacf[0] (should be variance); normalize to 1
    if vacf[0] == 0:
        raise ValueError("VACF[0] evaluated to zero (zero velocities?), cannot normalize.")
    vacf /= vacf[0]
    return vacf


# ---------------------------
# VDOS (power spectrum) from VACF
# ---------------------------
def compute_vdos_from_vacf(vacf, dt_fs, apply_taper=True):
    """
    Compute VDOS (power spectrum) from normalized VACF.
    - vacf: 1D array length n
    - dt_fs: timestep in femtoseconds
    Returns frequencies (Hz) and spectral power (real, >=0) for rfft bins.
    """
    n = len(vacf)
    dt = dt_fs * FS_TO_S

    # taper to reduce ringing (Hann)
    if apply_taper:
        if n > 1:
            hann = 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(n) / (n - 1)))
        else:
            hann = np.ones(1)
        vacf_tap = vacf * hann
    else:
        vacf_tap = vacf

    # rFFT (real-input FFT) -> gives n//2 + 1 frequency bins
    spec = np.fft.rfft(vacf_tap)
    # Because VACF is real and even, the rfft should be real-valued; but numerical noise may give small imag parts
    spec_real = np.real(spec)

    freqs = np.fft.rfftfreq(n, d=dt)  # in Hz
    # Make sure spec is non-negative (numerical rounding)
    spec_real[spec_real < 0] = 0.0
    return freqs, spec_real

def smooth_vdos(freq, vdos, sigma, sigma_unit='Hz'):
    """
    Smooth a VDOS curve with a Gaussian kernel.

    Parameters
    ----------
    freq : 1D array_like
        Frequency grid corresponding to `vdos`. **Units must match `sigma_unit`**:
        - if sigma_unit == 'Hz' : freq in Hz
        - if sigma_unit == 'THz': freq in THz
        - if sigma_unit == 'cm-1': freq in cm^-1
    vdos : 1D array_like
        Spectral power (same length as freq).
    sigma : float
        Gaussian kernel width expressed in the units indicated by `sigma_unit`.
        Interpreted as the Gaussian standard deviation (not FWHM).
    sigma_unit : {'Hz','THz','cm-1'}, optional
        Unit for `sigma`. Default 'Hz'.

    Returns
    -------
    vdos_smooth : ndarray
        Smoothed VDOS (same shape as `vdos`).
    """
    freq = np.asarray(freq, dtype=float)
    vdos = np.asarray(vdos, dtype=float)

    if freq.ndim != 1 or vdos.ndim != 1:
        raise ValueError("freq and vdos must be 1D arrays.")
    if freq.size != vdos.size:
        raise ValueError("freq and vdos must have the same length.")
    n = freq.size
    if n < 2 or np.allclose(vdos, 0.0):
        # Nothing sensible to do; return a copy
        return vdos.copy()

    # convert sigma to same units as freq (Hz)
    if sigma_unit == 'Hz':
        sigma_hz = float(sigma)
    elif sigma_unit == 'THz':
        sigma_hz = float(sigma) * 1e12
    elif sigma_unit == 'cm-1':
        sigma_hz = float(sigma) / CM1_TO_HZ  # CM1_TO_HZ is 1 cm^-1 in Hz (your constant)
    else:
        raise ValueError("Unsupported sigma_unit. Use 'Hz', 'THz' or 'cm-1'.")

    # get typical spacing (use median to be robust to small nonuniformities)
    df = np.median(np.diff(freq))
    if df <= 0:
        raise ValueError("Frequency grid must be strictly increasing.")

    # convert sigma (Hz) to index units for gaussian_filter1d
    sigma_idx = sigma_hz / df

    # guard against extremely small sigma_idx
    if sigma_idx <= 1.0e-6:
        # effectively no smoothing
        return vdos.copy()

    # perform Gaussian smoothing in index space
    # use reflect mode to avoid large edge artifacts
    vdos_smooth = gaussian_filter1d(vdos, sigma=sigma_idx, mode='reflect')

    # ensure non-negative (safety)
    vdos_smooth[vdos_smooth < 0] = 0.0
    return vdos_smooth


# ---------------------------
# Plotting and writing
# ---------------------------
def plot_vdos(freqs_hz, spec, title=None, show_cm1=True):
    """
    Plot VDOS on THz axis and optionally cm^-1 axis on top.
    """
    freqs_thz = freqs_hz * 1e-12
    plt.figure(figsize=(8, 5))
    plt.plot(freqs_thz, spec)
    plt.xlabel("Frequency (THz)")
    plt.ylabel("Spectral Power (arb. units)")
    if title:
        plt.title(title)
    plt.show()


def write_spectrum_to_file(filename, freqs_hz, spec, units="THz"):
    """
    Write spectrum file with frequency and power columns.
    By default frequencies are written in THz.
    """
    freqs_thz = freqs_hz * 1e-12
    with open(filename, "w") as f:
        f.write("# freq(THz)\tpower\n")
        for fthz, p in zip(freqs_thz, spec):
            f.write(f"{fthz:.8e}\t{p:.12e}\n")


# ---------------------------
# Process velocities -> per-atom autocorrs (for older-style outputs)
# ---------------------------
def build_velocity_array_from_AllSnapshots(AllSnapshots, MyCrystal):
    """
    Given AllSnapshots (list-like of flat arrays/objects) and a MyCrystal with natom,
    return velocities in shape (nsteps, natom, 3).
    Assumes snapshot ordering [vx0, vy0, vz0, vx1, vy1, vz1, ...].
    """
    nsteps = len(AllSnapshots)
    natom = MyCrystal.natom
    velocities = np.zeros((nsteps, natom, 3), dtype=float)
    flat_snapshots = []
    for isnap, snap in enumerate(AllSnapshots):
        arr = _safe_snapshot_to_array(snap)
        flat_snapshots.append(arr)
    # validate length
    if len(flat_snapshots) == 0:
        raise ValueError("No snapshots found in AllSnapshots.")
    example_len = flat_snapshots[0].size
    expected_len = natom * 3
    if example_len < expected_len:
        raise ValueError(f"Snapshot length ({example_len}) < expected 3*natom ({expected_len}). Check AllSnapshots format.")

    for isnap in range(nsteps):
        arr = flat_snapshots[isnap]
        # parse into velocities array
        for iatom in range(natom):
            base = iatom * 3
            velocities[isnap, iatom, 0] = arr[base]
            velocities[isnap, iatom, 1] = arr[base + 1]
            velocities[isnap, iatom, 2] = arr[base + 2]
    return velocities


# ---------------------------
# Main CLI / control flow
# ---------------------------
def prefactor(option, volume, temperature, timestep_fs):
    """
    compute:
        
        the prefactor using the formula B= V* dt / (kb * T) of the GK formula for viscosity 
    
    OR 
    
        the prefactor for scaling of the VDOS
    
    """
    kb = 1.3806503  # this appears to be in different units in your original; preserve their formula
    if option == 0:
        # original: 1.0E-04 * volume * timestep / (kb * temperature)
        print('Prefactor for viscosity is ', volume * timestep_fs / (kb * temperature) * 1.0E-04)
        return 1.0E-04 * volume * timestep_fs / (kb * temperature)
    elif option in [1, 11]:
        Boltzmann = 8.6173303e-5  # eV / K
        Avogadro = 6.022140857e23
        au_angstrom_square = 1.0 / Avogadro * 1.0e-3 * 1.0e10 / (1.602176634e-19)
        print(' au_angstrom_square is ', au_angstrom_square)
        return (2.0 * au_angstrom_square / (temperature * Boltzmann))
    return 1.0


def main():
    # defaults
    UMDfile = 'OUTCAR.umd.dat'
    firststep = 0
    window_size = None 
    property = 0  
    flaggraph = 0
    mean_volume = 1.0
    mean_temperature = 300.0

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:i:t:w:r:g:")
    except getopt.GetoptError:
        print('spectra_umd_fixed.py -f <umdfile> -i <InitialStep> -w <WindowSize> -r <Property> -g <Graphical>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('spectra_umd_fixed.py program to compute velocity auto-correlation and the related viscosity or power spectra')
            print('spectra_umd_fixed.py -f <umdfile> -i <InitialStep> -w <WindowSize> -r <Property> -g <Graphical>')
            print('default values: -f OUTCAR.umd.dat -i 0 -w None -t 300 -r 0 -g 0')
            print('-r options:')
            print('0 : viscosity  (Default)')
            print('1: combined VDOS analysis')
            print('11 : separated VDOS by atomic type')
            sys.exit()
        elif opt == '-f':
            UMDfile = str(arg)
        elif opt == '-i':
            firststep = int(arg)
        elif opt == '-t':
            mean_temperature = float(arg)
        elif opt == '-w':
            window_size = int(arg)
        elif opt == '-r':
            property = int(arg)
        elif opt == '-g':
            flaggraph = int(arg)

    if not os.path.isfile(UMDfile):
        print('the umdfile ', UMDfile, ' does not exist')
        sys.exit(1)

    # Branch for viscosity (kept largely as before)
    if property == 0:
        MyCrystal, AllSnapshots, TimeStep = umdpf.read_stresses_4visc(UMDfile)
        TimeStep_fs = TimeStep  # assume returned in fs
        print('after reading the stresses, total snapshots =', len(AllSnapshots))
        temps = [s.temperature for s in AllSnapshots]
        mean_temperature = np.mean(temps)
        volumes = [s.cellvolume for s in AllSnapshots]
        mean_volume = np.mean(volumes)
        pref = prefactor(property, mean_volume, mean_temperature, TimeStep_fs)

        stresses_xy = np.array([snapshot.stress[3] for snapshot in AllSnapshots])
        stresses_yz = np.array([snapshot.stress[4] for snapshot in AllSnapshots])
        stresses_zx = np.array([snapshot.stress[5] for snapshot in AllSnapshots])

        autocorr_xy, _ = autocorr_fft_unbiased(stresses_xy - np.mean(stresses_xy))
        autocorr_yz, _ = autocorr_fft_unbiased(stresses_yz - np.mean(stresses_yz))
        autocorr_zx, _ = autocorr_fft_unbiased(stresses_zx - np.mean(stresses_zx))

        intvaluesxy = np.array([simpson(autocorr_xy[:i], dx=1) for i in range(1, len(autocorr_xy) + 1)]) * pref
        intvaluesyz = np.array([simpson(autocorr_yz[:i], dx=1) for i in range(1, len(autocorr_yz) + 1)]) * pref
        intvalueszx = np.array([simpson(autocorr_zx[:i], dx=1) for i in range(1, len(autocorr_zx) + 1)]) * pref
        visc_mean = (intvaluesxy + intvaluesyz + intvalueszx) / 3.0

        # plot and write
        plot_vdos(np.fft.rfftfreq(len(autocorr_xy), TimeStep_fs * FS_TO_S), intvaluesxy, title="Viscosity integral (xy)")
        viscfile = UMDfile.replace(".umd.dat", "") + ".visc.dat"
        write_spectrum_to_file(viscfile, np.fft.rfftfreq(len(autocorr_xy), TimeStep_fs * FS_TO_S), visc_mean)
        print("Wrote viscosity integral to", viscfile)
        return

    # Velocity-based spectra (10..13)
    if property in [1, 11]:
        # read velocities
        # umdpf.read_values returns (MyCrystal, AllSnapshots, TimeStep, length)
        MyCrystal, AllSnapshots, TimeStep, length = umdpf.read_values(UMDfile, "velocity", "line", 1, firststep, laststep=None, cutoff="all", nCores=None)
        TimeStep_fs = TimeStep  # assumed in fs
        natom = MyCrystal.natom
        print("Spectra from velocities: size of AllSnapshots", len(AllSnapshots), " by ", len(AllSnapshots[0]), " where natom =", natom)

        # Build velocities array shape (nsteps, natom, 3)
        velocities = build_velocity_array_from_AllSnapshots(AllSnapshots, MyCrystal)

        # Optionally truncate VACF length by window_size (but compute vacf from full traj then trim for FFT if desired)
        vacf = compute_vacf_from_velocities(velocities)
        if window_size is not None and window_size > 0 and window_size < len(vacf):
            vacf = vacf[:window_size]

        # Compute VDOS
        freqs_hz, spec = compute_vdos_from_vacf(vacf, TimeStep_fs, apply_taper=True)
        # Apply prefactor if desired (e.g., for units conversion)
        pref = prefactor(property, 1.0, mean_temperature, TimeStep_fs)
        spec_scaled = spec * pref

        base_out = UMDfile.replace(".umd.dat", "") + f".r{property}"
        specfile = base_out + ".vdos.dat"
        write_spectrum_to_file(specfile, freqs_hz, spec_scaled)
        print("Wrote spectrum to", specfile)

        # Plot combined VDOS (sum over all atoms/directions already averaged inside compute_vacf_from_velocities)
        plot_vdos(freqs_hz, spec_scaled, title=f"VDOS (property={property})")

        # Additional outputs / branches to mimic original behaviour:

        if property == 11:
            # separated by atomic types and directions
            # We'll compute per-type averaged VACF then VDOS
            types = MyCrystal.typat  # assuming typat is list of length natom with integer type indices
            ntypes = MyCrystal.ntypat
            print(ntypes)
            vacf_bytype = []
            for itype in range(ntypes):
                indices = [i for i in range(natom) if types[i] == itype]
                if len(indices) == 0:
                    continue
                # average velocities across atoms of this type
                vel_sub = velocities[:, indices, :].copy()
                vacf_t = compute_vacf_from_velocities(vel_sub)
                freqs_t, spec_t = compute_vdos_from_vacf(vacf_t, TimeStep_fs, apply_taper=True)
                write_spectrum_to_file(base_out + f".type{itype}.vdos.dat", freqs_t, spec_t * pref)
                vacf_bytype.append((itype, freqs_t, spec_t))
            # plot first few
            for (itype, f_hz, s) in vacf_bytype[:8]:
                plt.plot(f_hz * CM1_TO_HZ , s * pref, label=f"type {itype}")
            plt.xlabel("Frequency (cm-1)")
            plt.ylabel("Spectral power (arb.)")
            plt.legend()
            plt.title("VDOS by atomic type")
            plt.tight_layout()
            plt.savefig('power_spectra.pdf')
            plt.show()


if __name__ == "__main__":
    main()
