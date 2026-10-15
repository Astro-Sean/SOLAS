#!/usr/bin/env python3
"""Convert CMFGEN SN_HYDRO_DATA to SEDONA .mod format."""

import sys
import numpy as np

# Atomic number mapping for CMFGEN element names
cmfgen_to_z = {
    'HYD': 1, 'HE': 2, 'LI': 3, 'BE': 4, 'B': 5, 'CARB': 6, 'NIT': 7, 'OXY': 8,
    'F': 9, 'NEON': 10, 'NA': 11, 'MG': 12, 'AL': 13, 'SIL': 14, 'P': 15,
    'SUL': 16, 'CL': 17, 'AR': 18, 'K': 19, 'CA': 20, 'CAL': 20, 'SC': 21,
    'TI': 22, 'V': 23, 'CR': 24, 'MN': 25, 'FE': 26, 'IRON': 26, 'CO': 27,
    'COB': 27, 'NI': 28, 'NICK': 28
}

# Dominant isotope mass number per atomic number, written as the "A" in
# the .mod "Z.A" species tokens.  SEDONA's ascii parser requires Z.A; a
# bare "Z" is silently parsed as A=Z, halving the particle mass for all
# metals and He (grid_1D_sphere.cpp:read_ascii_file).
z_to_a = {
    1: 1, 2: 4, 3: 7, 4: 9, 5: 11, 6: 12, 7: 14, 8: 16, 9: 19, 10: 20,
    11: 23, 12: 24, 13: 27, 14: 28, 15: 31, 16: 32, 17: 35, 18: 40,
    19: 39, 20: 40, 21: 45, 22: 48, 23: 51, 24: 52, 25: 55, 26: 56,
    27: 59, 28: 58
}

def read_sn_hydro_data(filepath):
    """Read CMFGEN SN_HYDRO_DATA file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    nd = 156  # default
    time_days = 0.0
    lum_lsun = 1e9
    n_mass_frac = 0
    
    for line in lines[:10]:
        stripped = line.strip()
        if 'Number of data points' in stripped:
            parts = stripped.split()
            for i, p in enumerate(parts):
                if p == 'points:':
                    nd = int(parts[i+1])
                    break
        elif 'Time(days) since explosion' in stripped:
            parts = stripped.split()
            for i, p in enumerate(parts):
                if p == 'explosion:':
                    time_days = float(parts[i+1])
                    break
        elif 'Energy Lum [lsun]' in stripped:
            parts = stripped.split()
            for i, p in enumerate(parts):
                if p == '[lsun]:':
                    lum_lsun = float(parts[i+1])
                    break
        elif 'Number of mass fractions' in stripped:
            parts = stripped.split()
            for i, p in enumerate(parts):
                if p == 'fractions:':
                    n_mass_frac = int(parts[i+1])
                    break
    
    # Parse data sections
    sections = {}
    current_section = None
    current_data = []
    
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        
        # Check for section headers
        if 'Radius grid' in stripped:
            current_section = 'radius'
            continue
        elif 'Velocity' in stripped and 'km/s' in stripped:
            current_section = 'velocity'
            continue
        elif 'Density' in stripped and 'gm/cm^3' in stripped:
            current_section = 'density'
            continue
        elif 'Temperature' in stripped:
            current_section = 'temperature'
            continue
        elif 'Sigma' in stripped:
            current_section = 'sigma'
            continue
        elif 'Electron density' in stripped:
            current_section = 'electron_density'
            continue
        elif 'Atom density' in stripped:
            current_section = 'atom_density'
            continue
        elif 'Kappa' in stripped:
            current_section = 'kappa'
            continue
        elif 'Rosseland' in stripped:
            current_section = 'rosseland'
            continue
        elif 'mass fraction' in stripped:
            parts = stripped.split()
            elem_name = parts[0]
            # Isotope sections look like "NICK 56 mass fraction": keep the
            # mass number in the key so isotopes do not overwrite the
            # elemental column of the same name (mass_frac_NICK vs
            # mass_frac_NICK56).
            if len(parts) > 1 and parts[1].isdigit():
                elem_name += parts[1]
            current_section = f'mass_frac_{elem_name}'
            continue
        elif 'Isotope' in stripped or 'isotope' in stripped:
            current_section = 'isotope'
            continue
        
        # Check if this is a header line (starts with text)
        if stripped and not stripped[0].isdigit() and not stripped[0] == '-':
            if 'mass fraction' in stripped:
                parts = stripped.split()
                elem_name = parts[0]
                if len(parts) > 1 and parts[1].isdigit():
                    elem_name += parts[1]
                current_section = f'mass_frac_{elem_name}'
                continue
            else:
                current_section = None
                continue
        
        # Parse numeric data
        if current_section and stripped:
            try:
                values = [float(x) for x in stripped.split()]
                current_data.extend(values)
                
                # If we have enough values for this section, save it
                if len(current_data) >= nd:
                    sections[current_section] = np.array(current_data[:nd])
                    current_data = []
                    current_section = None
            except ValueError:
                pass
    
    # Handle any remaining data
    if current_section and current_data:
        sections[current_section] = np.array(current_data[:nd])
    
    # Extract mass fraction elements
    mass_fracs = {}
    for key in sections:
        if key.startswith('mass_frac_'):
            elem = key.replace('mass_frac_', '')
            mass_fracs[elem] = sections[key]
    
    return {
        'nd': nd,
        'time_days': time_days,
        'lum_lsun': lum_lsun,
        'n_mass_frac': n_mass_frac,
        'radius': sections.get('radius', np.zeros(nd)),  # in 10^10 cm
        'velocity': sections.get('velocity', np.zeros(nd)),  # in km/s
        'density': sections.get('density', np.zeros(nd)),  # in g/cm^3
        'temperature': sections.get('temperature', np.zeros(nd)),  # in K
        'mass_fracs': mass_fracs,
        'rosseland': sections.get('rosseland', np.zeros(nd)),  # chi_R in 10^-10 cm^-1
    }

def convert_to_sedona(sn_data, output_file, species_list=None, time_days=None,
                        homologous=False, tau_max=None):
    """Convert SN_HYDRO_DATA to SEDONA .mod format.

    Parameters
    ----------
    sn_data : dict
        Parsed SN_HYDRO_DATA dictionary.
    output_file : str
        Path to output .mod file.
    species_list : list[int] or None
        Atomic numbers to include. Default: H, He, C, O, Si, Ca, Fe, Co, Ni.
    time_days : float or None
        Time since explosion in days. Overrides sn_data['time_days'] if given.
    homologous : bool
        If True, writes the correct time and outer radius fields for use
        with hydro_module="homologous" in SEDONA (Sobolev line expansion).
        If False, writes 0.0 for time and v_max for outer radius for backward
        compatibility with hydro_module="none".
    tau_max : float or None
        If set, drop innermost zones beyond this cumulative Rosseland tau
        (measured inward from the outer edge). Uses the 'Rosseland mean
        opacity' column of SN_HYDRO_DATA.
    """
    nd = sn_data['nd']

    # Default species list (atomic numbers); zero-abundance entries are dropped
    if species_list is None:
        species_list = [1, 2, 6, 8, 14, 20, 26, 27, 28]  # H, He, C, O, Si, Ca, Fe, Co, Ni

    # Extract data
    r = sn_data['radius']  # in 10^10 cm
    v = sn_data['velocity']  # in km/s
    rho = sn_data['density']  # in g/cm^3

    # Convert units
    r_cm = r * 1e10  # convert to cm
    v_cms = v * 1e5  # convert to cm/s
    
    # Calculate temperature if not present
    if 'temperature' in sn_data and len(sn_data['temperature']) == nd and np.any(sn_data['temperature'] > 0):
        # SN_HYDRO_DATA temperature is in 10^4 K units
        T = sn_data['temperature'] * 1e4  # convert to K
        print(f"  Using temperature from SN_HYDRO_DATA: {np.min(T):.1f} - {np.max(T):.1f} K")
    else:
        # Estimate temperature from luminosity and radius
        # For a blackbody: L = 4*pi*R^2*sigma*T^4
        L = sn_data['lum_lsun'] * 3.828e33  # erg/s
        sigma = 5.6704e-5  # Stefan-Boltzmann constant
        
        # Calculate temperature at inner boundary (photosphere)
        r_inner = np.min(r_cm) if np.any(r_cm > 0) else 1e15
        T_phot = (L / (4 * np.pi * r_inner**2 * sigma))**0.25
        print(f"  L = {L:.3e} erg/s, r_inner = {r_inner:.3e} cm, T_phot = {T_phot:.1f} K")
        
        # Scale temperature with radius (approximate)
        # T ~ R^(-1/2) for radiation-dominated atmospheres
        T = np.zeros(nd)
        for i in range(nd):
            if r_cm[i] > 0:
                T[i] = T_phot * (r_inner / r_cm[i])**0.5
            else:
                T[i] = T_phot
        
        # Cap temperature to reasonable values
        T = np.clip(T, 3000, 100000)
        print(f"  Calculated temperature range: {np.min(T):.1f} - {np.max(T):.1f} K")
    
    # Get mass fractions
    mass_fracs = sn_data['mass_fracs']
    
    # Map CMFGEN element names to atomic numbers.  Isotope keys
    # (e.g. "NICK56" from "NICK 56 mass fraction") do not match
    # cmfgen_to_z and are handled separately below.
    elem_to_z = {}
    for elem in mass_fracs:
        z = cmfgen_to_z.get(elem.upper(), 0)
        if z > 0:
            elem_to_z[elem] = z

    # Build mass fraction array for each species.  species_tokens holds
    # the "Z.A" strings written to the .mod header, kept in parallel with
    # species_list so the zero-abundance filter stays aligned.
    n_species = len(species_list)
    species_mass_fracs = np.zeros((nd, n_species))
    species_tokens = [f"{z}.{z_to_a.get(z, 2 * z)}" for z in species_list]

    for i, z in enumerate(species_list):
        # Find the element with this atomic number
        found = False
        for elem, elem_z in elem_to_z.items():
            if elem_z == z and elem in mass_fracs:
                species_mass_fracs[:, i] = mass_fracs[elem]
                found = True
                break
        if not found:
            # Set to zero if not found
            species_mass_fracs[:, i] = 0.0

    # Append explicit isotope columns (e.g. "NICK56" -> species "28.56").
    # These are separate .mod species so SEDONA's radioactive.decay can
    # track chains like Ni56 -> Co56 -> Fe56.
    import re
    for key, values in mass_fracs.items():
        m = re.match(r'^([A-Z]+)(\d+)$', key.upper())
        if not m:
            continue
        z = cmfgen_to_z.get(m.group(1), 0)
        if z == 0:
            continue
        token = f"{z}.{int(m.group(2))}"
        col = np.asarray(values, dtype=float).reshape(-1, 1)
        if token in species_tokens:
            # Same Z.A as an elemental column (e.g. decay-chain IRON 56
            # duplicating elemental IRON): merge into the existing column
            # so the .mod has no duplicate species.
            j = species_tokens.index(token)
            species_mass_fracs[:, j] += col[:, 0]
            print(f"  Merged isotope column {key} into species {token}")
        else:
            species_mass_fracs = np.hstack([species_mass_fracs, col])
            species_tokens.append(token)
            print(f"  Added isotope column {key} as species {token}")

    # Remove species columns that are all zeros (no atoms of this element)
    non_zero_mask = np.any(species_mass_fracs > 0, axis=0)
    if not np.all(non_zero_mask):
        removed = [t for i, t in enumerate(species_tokens) if not non_zero_mask[i]]
        print(f"  Removing zero-abundance species from .mod: {removed}")
        species_tokens = [t for i, t in enumerate(species_tokens) if non_zero_mask[i]]
        species_mass_fracs = species_mass_fracs[:, non_zero_mask]
        n_species = len(species_tokens)
        print(f"  Remaining species: {species_tokens}")
    else:
        n_species = len(species_tokens)

    # Normalize mass fractions to sum to 1
    for i in range(nd):
        total = np.sum(species_mass_fracs[i, :])
        if total > 0:
            species_mass_fracs[i, :] /= total

    # Optional tau-trim: SN_HYDRO_DATA rows run outer -> inner (r decreasing),
    # so the cumulative Rosseland tau grows monotonically inward.
    rosseland = sn_data.get('rosseland')
    if tau_max is not None:
        if rosseland is None or len(rosseland) != nd or not np.any(rosseland > 0):
            print("  WARNING: tau_max requested but no Rosseland opacity column;"
                  " skipping trim")
        else:
            chi = rosseland * 1e-10  # 10^-10 cm^-1 -> cm^-1
            tau = np.zeros(nd)
            for i in range(1, nd):
                dr = abs(r_cm[i - 1] - r_cm[i])
                tau[i] = tau[i - 1] + 0.5 * (chi[i - 1] + chi[i]) * dr
            print(f"  Rosseland tau at inner edge: {tau[-1]:.1f}")
            keep = tau <= tau_max
            if not np.any(keep):
                keep[0] = True  # always keep the outermost zone
            n_drop = nd - int(keep.sum())
            r_cm = r_cm[keep]
            v_cms = v_cms[keep]
            rho = rho[keep]
            T = T[keep]
            species_mass_fracs = species_mass_fracs[keep, :]
            nd = int(keep.sum())
            print(f"  Tau-trimmed to tau_max={tau_max:g}: dropped {n_drop} inner "
                  f"zones, nd={nd}, new inner edge r={r_cm[-1]:.3e} cm, "
                  f"tau={tau[keep][-1]:.1f}")

    # Sort by radius in ascending order (Sedona requires ascending r_out)
    sort_idx = np.argsort(r_cm)
    r_cm = r_cm[sort_idx]
    v_cms = v_cms[sort_idx]
    rho = rho[sort_idx]
    T = T[sort_idx]
    species_mass_fracs = species_mass_fracs[sort_idx, :]

    # Determine header values
    t_model = time_days if time_days is not None else sn_data['time_days']
    t_sec = t_model * 24.0 * 3600.0
    r_max = r_cm[-1]
    v_max = np.max(v_cms) if len(v_cms) > 0 else 1e9

    # Verify data integrity
    print(f"  Data verification before writing:")
    print(f"    r_cm: {r_cm.min():.3e} - {r_cm.max():.3e} cm, first: {r_cm[0]:.3e}")
    print(f"    v_cms: {v_cms.min():.3e} - {v_cms.max():.3e} cm/s, first: {v_cms[0]:.3e}")
    print(f"    rho: {rho.min():.3e} - {rho.max():.3e} g/cm^3, first: {rho[0]:.3e}")
    print(f"    T: {T.min():.1f} - {T.max():.1f} K, first: {T[0]:.1f}")

    # Write output file
    with open(output_file, 'w') as f:
        # SEDONA 1D_sphere standard header: nzon  rmin  texp  n_elems
        # (see grid_1D_sphere.cpp:read_ascii_file, lines 288-296)
        #   rmin: inner boundary radius. SEDONA creates zone 0 from rmin to r_out[0].
        #         Must NOT be 0 for a trimmed CSM shell, because rmin=0 creates a huge
        #         fake zone 0 (volume 4/3*pi*r[0]^3) adding ~16 Msun of fake mass and
        #         tau~1200, trapping all photons.
        #   texp: time since explosion in seconds (used for t_now in SEDONA).
        f.write("1D_sphere standard\n")
        if homologous:
            # Homologous: rmin is not used (SEDONA computes from t and v), write t for both
            rmin = r_cm[0] * 0.99  # just inside first zone
            f.write(f"{nd}\t{rmin:.8e}\t{t_sec:.1f}\t{n_species}\n")
            print(f"  Header (homologous): {nd} zones, rmin={rmin:.3e} cm, t_now={t_sec:.1f}s ({t_model:.2f}d)")
        else:
            # Static snapshot: set rmin just inside r[0] (half a zone width inward)
            if nd > 1:
                dr0 = abs(r_cm[1] - r_cm[0])
            else:
                dr0 = r_cm[0] * 0.01
            rmin = max(r_cm[0] - dr0, 0.0)
            f.write(f"{nd}\t{rmin:.8e}\t{t_sec:.1f}\t{n_species}\n")
            print(f"  Header (static): {nd} zones, rmin={rmin:.3e} cm, t_now={t_sec:.1f}s ({t_model:.2f}d)")
        # Species as "Z.A" isotope tokens (required by
        # grid_1D_sphere::read_ascii_file; bare "Z" parses as A=Z).
        f.write(" ".join(species_tokens) + "\n")

        # Data lines: r, v, rho, T, mass_fracs...
        # Use .8e precision for r, v, rho, T to avoid zero-thickness zones
        # when adjacent zones differ by less than 1e-4 * r (e.g. inner CSM zones
        # with dr ~ 2e9 cm at r ~ 4e14 cm would round to the same value with .4e).
        for i in range(nd):
            line = f"{r_cm[i]:.8e} {v_cms[i]:.8e} {rho[i]:.8e} {T[i]:.8e}"
            for j in range(n_species):
                line += f" {species_mass_fracs[i, j]:.4f}"
            line += "\n"
            f.write(line)

    print(f"Written {output_file}: {nd} zones, {n_species} species")
    print(f"  Time: {t_model:.3f} days")
    print(f"  Luminosity: {sn_data['lum_lsun']:.3e} Lsun")
    print(f"  Radius range: [{r_cm[0]:.3e}, {r_cm[-1]:.3e}] cm")
    print(f"  Velocity range: [{v_cms[0]:.3e}, {v_cms[-1]:.3e}] cm/s")
    print(f"  Density range: [{rho[0]:.3e}, {rho[-1]:.3e}] g/cm^3")
    print(f"  Temperature range: [{T[0]:.1f}, {T[-1]:.1f}] K")

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Convert CMFGEN SN_HYDRO_DATA to SEDONA .mod format."
    )
    parser.add_argument('input_file', help='Path to SN_HYDRO_DATA file')
    parser.add_argument('output_file', help='Path to output .mod file')
    parser.add_argument('--species', type=str, default=None,
                        help='Comma-separated atomic numbers (e.g. 1,2,6,8,14,26)')
    parser.add_argument('--time-days', type=float, default=None,
                        help='Time since explosion in days (overrides SN_HYDRO_DATA header)')
    parser.add_argument('--homologous', action='store_true',
                        help='Write correct time and outer radius for hydro_module="homologous" '
                             '(Sobolev line expansion). Otherwise writes legacy header for '
                             'hydro_module="none".')
    parser.add_argument('--tau-max', type=float, default=None,
                        help='Drop innermost zones beyond this cumulative Rosseland tau')
    args = parser.parse_args()

    species_list = None
    if args.species:
        species_list = [int(x.strip()) for x in args.species.split(',')]

    print(f"Reading {args.input_file}...")
    sn_data = read_sn_hydro_data(args.input_file)

    print(f"Converting to SEDONA format...")
    convert_to_sedona(sn_data, args.output_file, species_list,
                        time_days=args.time_days, homologous=args.homologous,
                        tau_max=args.tau_max)

if __name__ == '__main__':
    main()
