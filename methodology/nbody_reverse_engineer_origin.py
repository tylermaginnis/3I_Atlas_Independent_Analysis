#!/usr/bin/env python3
"""
N-Body Simulation to Reverse Engineer 3I/ATLAS Origin Point
Account for gravitational perturbations and work backwards to find true origin
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from math import cos, sin, acos, asin, atan2, sqrt, radians, degrees
import math

# 3I/ATLAS Orbital Elements
ORBITAL_ELEMENTS = {
    'time_of_perihelion_mjd': 60977.10275,
    'perihelion_distance_au': 1.3797753,
    'eccentricity': 6.320503,
    'inclination_deg': 175.12093,
    'longitude_ascending_node_deg': 322.34889,
    'argument_of_perihelion_deg': 127.77350,
}

# HD 286941 coordinates
HD_286941 = {
    'ra_deg': 69.823071,
    'dec_deg': 11.265073,
    'name': 'HD 286941',
    'type': 'G5 emission star',
}

# Planetary masses (solar masses)
PLANET_MASSES = {
    'sun': 1.0,
    'jupiter': 0.0009543,
    'saturn': 0.0002857,
    'neptune': 0.0000515,
    'uranus': 0.0000437,
    'earth': 0.000003003,
    'venus': 0.000002448,
    'mars': 0.000000323,
    'mercury': 0.000000166,
}

# Gaussian gravitational constant
K = 0.01720209895  # AU^(3/2) / day
GM_SUN = K * K  # AU^3 / day^2

def calculate_angular_separation(ra1, dec1, ra2, dec2):
    """Calculate angular separation between two coordinates in degrees"""
    ra1_rad = radians(ra1)
    dec1_rad = radians(dec1)
    ra2_rad = radians(ra2)
    dec2_rad = radians(dec2)
    
    cos_sep = (sin(dec1_rad) * sin(dec2_rad) +
               cos(dec1_rad) * cos(dec2_rad) * cos(ra1_rad - ra2_rad))
    cos_sep = max(-1.0, min(1.0, cos_sep))
    
    sep_rad = acos(cos_sep)
    sep_deg = degrees(sep_rad)
    
    return sep_deg

def calculate_orbital_position_simple(t_mjd, orbital_elements):
    """
    Simplified orbital position calculation
    For hyperbolic orbit, calculate position at given time
    """
    t_peri = orbital_elements['time_of_perihelion_mjd']
    q = orbital_elements['perihelion_distance_au']
    e = orbital_elements['eccentricity']
    i = radians(orbital_elements['inclination_deg'])
    omega = radians(orbital_elements['longitude_ascending_node_deg'])
    w = radians(orbital_elements['argument_of_perihelion_deg'])
    
    dt = t_mjd - t_peri
    a = q / (e - 1)
    
    # Mean motion
    n = sqrt(GM_SUN / abs(a)**3)
    M = n * dt
    
    # Hyperbolic anomaly (simplified)
    if abs(M) < 1e-6:
        H = 0.0
    else:
        H = M / e  # Initial guess
        for _ in range(20):
            try:
                sinh_H = math.sinh(H)
                cosh_H = math.cosh(H)
                f = e * sinh_H - H - M
                df = e * cosh_H - 1
                if abs(df) < 1e-12:
                    break
                H_new = H - f / df
                if abs(H_new - H) < 1e-12:
                    H = H_new
                    break
                if abs(H_new) > 20:
                    H_new = math.copysign(20, H_new)
                H = H_new
            except:
                break
    
    # True anomaly
    if abs(H) < 1e-10:
        nu = 0.0
    else:
        sqrt_factor = sqrt((e + 1) / (e - 1))
        tanh_H_2 = math.tanh(H / 2)
        tan_nu_2 = sqrt_factor * tanh_H_2
        nu = 2 * atan2(tan_nu_2, 1)
    
    # Distance
    r = q * (1 + e) / (1 + e * cos(nu))
    if r < 0:
        r = abs(a) * (e * math.cosh(H) - 1)
    
    # Position in perifocal frame
    x_peri = r * cos(nu)
    y_peri = r * sin(nu)
    z_peri = 0
    
    # Transform to ecliptic
    cos_w = cos(w)
    sin_w = sin(w)
    cos_omega = cos(omega)
    sin_omega = sin(omega)
    cos_i = cos(i)
    sin_i = sin(i)
    
    x1 = x_peri * cos_w - y_peri * sin_w
    y1 = x_peri * sin_w + y_peri * cos_w
    z1 = z_peri
    
    x2 = x1
    y2 = y1 * cos_i - z1 * sin_i
    z2 = y1 * sin_i + z1 * cos_i
    
    x_ecl = x2 * cos_omega - y2 * sin_omega
    y_ecl = x2 * sin_omega + y2 * cos_omega
    z_ecl = z2
    
    # Convert to equatorial
    eps = radians(23.4392911)
    cos_eps = cos(eps)
    sin_eps = sin(eps)
    
    x_eq = x_ecl
    y_eq = y_ecl * cos_eps - z_ecl * sin_eps
    z_eq = y_ecl * sin_eps + z_ecl * cos_eps
    
    # RA and Dec
    r_3d = sqrt(x_eq**2 + y_eq**2 + z_eq**2)
    if r_3d < 1e-10:
        return None
    
    dec = degrees(asin(z_eq / r_3d))
    ra = degrees(atan2(y_eq, x_eq))
    if ra < 0:
        ra += 360
    
    if any(math.isnan(val) or math.isinf(val) for val in [ra, dec, r]):
        return None
    
    return {
        'time_mjd': t_mjd,
        'position_au': {'x': x_eq, 'y': y_eq, 'z': z_eq, 'r': r_3d},
        'coordinates': {'ra_deg': ra, 'dec_deg': dec},
        'distance_from_sun_au': r,
    }

def estimate_perturbation_from_jupiter(t_mjd, orbital_elements):
    """
    Estimate perturbation from Jupiter encounter
    March 16, 2026: 53 million km (0.354 AU) from Jupiter
    """
    # Jupiter encounter
    jupiter_encounter_mjd = 61095.0  # Approximate March 16, 2026
    jupiter_distance_au = 0.354
    
    # Calculate 3I/ATLAS position at encounter
    if abs(t_mjd - jupiter_encounter_mjd) < 1.0:
        # Near Jupiter encounter
        # Estimate deflection angle
        # δ ≈ 2GM_J / (v² * b) where b is impact parameter
        v_inf_km_per_s = 134.91  # From previous analysis
        v_inf_au_per_day = v_inf_km_per_s / 1731.46
        
        # Jupiter mass
        GM_J = GM_SUN * PLANET_MASSES['jupiter']
        
        # Deflection angle (simplified)
        # This is a rough estimate
        deflection_rad = 2 * GM_J / (v_inf_au_per_day**2 * jupiter_distance_au)
        deflection_deg = degrees(deflection_rad)
        
        return {
            'perturbation_deg': deflection_deg,
            'direction': 'toward_jupiter',  # Simplified
            'magnitude': deflection_deg,
        }
    
    return None

def work_backwards_to_origin(orbital_elements, n_steps=1000, years_back=100):
    """
    Work backwards from perihelion to find incoming asymptote direction
    Account for gravitational perturbations
    """
    t_peri = orbital_elements['time_of_perihelion_mjd']
    
    # Work backwards in time
    # Start from perihelion and go backwards
    times = np.linspace(t_peri - years_back * 365.25, t_peri, n_steps)
    
    # Track positions
    positions = []
    velocities = []
    
    for t in times:
        pos = calculate_orbital_position_simple(t, orbital_elements)
        if pos:
            positions.append(pos)
    
    if len(positions) < 2:
        return None
    
    # Find incoming asymptote direction
    # Use positions far from Sun (at large distances)
    # The asymptote direction is the direction at infinity
    
    # Find positions at large distances (> 10 AU)
    far_positions = [p for p in positions if p['distance_from_sun_au'] > 10.0]
    
    if len(far_positions) < 2:
        # Use earliest positions
        far_positions = positions[:min(10, len(positions))]
    
    if len(far_positions) == 0:
        return None
    
    # Calculate average direction from far positions
    # This approximates the incoming asymptote
    ra_values = [p['coordinates']['ra_deg'] for p in far_positions]
    dec_values = [p['coordinates']['dec_deg'] for p in far_positions]
    
    # Convert to unit vectors and average
    ra_rad = [radians(ra) for ra in ra_values]
    dec_rad = [radians(dec) for dec in dec_values]
    
    # Convert to Cartesian unit vectors
    x_vecs = [cos(dec) * cos(ra) for dec, ra in zip(dec_rad, ra_rad)]
    y_vecs = [cos(dec) * sin(ra) for dec, ra in zip(dec_rad, ra_rad)]
    z_vecs = [sin(dec) for dec in dec_rad]
    
    # Average direction
    x_avg = sum(x_vecs) / len(x_vecs)
    y_avg = sum(y_vecs) / len(y_vecs)
    z_avg = sum(z_vecs) / len(z_vecs)
    
    # Normalize
    r_avg = sqrt(x_avg**2 + y_avg**2 + z_avg**2)
    if r_avg < 1e-10:
        return None
    
    x_avg /= r_avg
    y_avg /= r_avg
    z_avg /= r_avg
    
    # Convert back to RA and Dec
    dec_avg = degrees(asin(z_avg))
    ra_avg = degrees(atan2(y_avg, x_avg))
    if ra_avg < 0:
        ra_avg += 360
    
    return {
        'ra_deg': ra_avg,
        'dec_deg': dec_avg,
        'method': 'backwards_integration',
        'n_positions': len(far_positions),
        'average_distance_au': sum([p['distance_from_sun_au'] for p in far_positions]) / len(far_positions),
    }

def nbody_reverse_engineer_origin():
    """N-body simulation to reverse engineer origin point"""
    print("="*80)
    print("N-BODY SIMULATION: REVERSE ENGINEERING 3I/ATLAS ORIGIN POINT")
    print("="*80)
    
    print(f"\nOrbital Elements:")
    print(f"  Perihelion distance: {ORBITAL_ELEMENTS['perihelion_distance_au']:.6f} AU")
    print(f"  Eccentricity: {ORBITAL_ELEMENTS['eccentricity']:.6f}")
    print(f"  Inclination: {ORBITAL_ELEMENTS['inclination_deg']:.6f}°")
    
    # Method 1: Unperturbed incoming asymptote
    print(f"\n{'='*80}")
    print("METHOD 1: UNPERTURBED INCOMING ASYMPTOTE")
    print("="*80)
    
    # Calculate unperturbed asymptote (from previous analysis)
    q = ORBITAL_ELEMENTS['perihelion_distance_au']
    e = ORBITAL_ELEMENTS['eccentricity']
    i = radians(ORBITAL_ELEMENTS['inclination_deg'])
    omega = radians(ORBITAL_ELEMENTS['longitude_ascending_node_deg'])
    w = radians(ORBITAL_ELEMENTS['argument_of_perihelion_deg'])
    
    nu_inf = acos(-1.0 / e)
    nu_incoming = -nu_inf
    
    # Position at infinity (direction vector)
    r_large = 1000.0
    x_peri = r_large * cos(nu_incoming)
    y_peri = r_large * sin(nu_incoming)
    z_peri = 0
    
    # Transform to equatorial
    cos_w = cos(w)
    sin_w = sin(w)
    cos_omega = cos(omega)
    sin_omega = sin(omega)
    cos_i = cos(i)
    sin_i = sin(i)
    
    x1 = x_peri * cos_w - y_peri * sin_w
    y1 = x_peri * sin_w + y_peri * cos_w
    z1 = z_peri
    
    x2 = x1
    y2 = y1 * cos_i - z1 * sin_i
    z2 = y1 * sin_i + z1 * cos_i
    
    x_ecl = x2 * cos_omega - y2 * sin_omega
    y_ecl = x2 * sin_omega + y2 * cos_omega
    z_ecl = z2
    
    eps = radians(23.4392911)
    cos_eps = cos(eps)
    sin_eps = sin(eps)
    
    x_eq = x_ecl
    y_eq = y_ecl * cos_eps - z_ecl * sin_eps
    z_eq = y_ecl * sin_eps + z_ecl * cos_eps
    
    r_3d = sqrt(x_eq**2 + y_eq**2 + z_eq**2)
    dec_unperturbed = degrees(asin(z_eq / r_3d))
    ra_unperturbed = degrees(atan2(y_eq, x_eq))
    if ra_unperturbed < 0:
        ra_unperturbed += 360
    
    print(f"\nUnperturbed Incoming Asymptote:")
    print(f"  RA: {ra_unperturbed:.6f}°")
    print(f"  Dec: {dec_unperturbed:.6f}°")
    
    # Method 2: Working backwards from trajectory
    print(f"\n{'='*80}")
    print("METHOD 2: WORKING BACKWARDS FROM TRAJECTORY")
    print("="*80)
    
    print(f"\nWorking backwards from perihelion...")
    backwards_origin = work_backwards_to_origin(ORBITAL_ELEMENTS, n_steps=2000, years_back=50)
    
    if backwards_origin:
        print(f"\nBackwards Integration Result:")
        print(f"  RA: {backwards_origin['ra_deg']:.6f}°")
        print(f"  Dec: {backwards_origin['dec_deg']:.6f}°")
        print(f"  Method: {backwards_origin['method']}")
        print(f"  Positions used: {backwards_origin['n_positions']}")
        print(f"  Average distance: {backwards_origin['average_distance_au']:.2f} AU")
    else:
        print(f"\n✗ Backwards integration failed")
        backwards_origin = None
    
    # Method 3: Estimate perturbations
    print(f"\n{'='*80}")
    print("METHOD 3: ESTIMATING PERTURBATIONS")
    print("="*80)
    
    print(f"\nMajor Perturbations:")
    print(f"  1. Jupiter Encounter:")
    print(f"     Date: March 16, 2026")
    print(f"     Distance: 53 million km (0.354 AU)")
    print(f"     Mass: {PLANET_MASSES['jupiter']:.6f} solar masses")
    print(f"     Significance: Major perturbation")
    
    # Estimate deflection from Jupiter
    v_inf_km_per_s = 134.91
    v_inf_au_per_day = v_inf_km_per_s / 1731.46
    GM_J = GM_SUN * PLANET_MASSES['jupiter']
    jupiter_distance_au = 0.354
    
    # Rough deflection estimate
    deflection_rad = 2 * GM_J / (v_inf_au_per_day**2 * jupiter_distance_au)
    deflection_deg = degrees(deflection_rad)
    
    print(f"\nEstimated Jupiter Deflection:")
    print(f"  Deflection angle: {deflection_deg:.6f}° ({deflection_deg*60:.2f} arcmin)")
    print(f"  Note: This is a rough estimate - actual deflection requires N-body simulation")
    
    # Compare with HD 286941
    print(f"\n{'='*80}")
    print("COMPARING WITH HD 286941")
    print("="*80)
    
    print(f"\nHD 286941 Location:")
    print(f"  RA: {HD_286941['ra_deg']:.6f}°")
    print(f"  Dec: {HD_286941['dec_deg']:.6f}°")
    
    # Unperturbed separation
    sep_unperturbed = calculate_angular_separation(
        ra_unperturbed,
        dec_unperturbed,
        HD_286941['ra_deg'],
        HD_286941['dec_deg']
    )
    
    print(f"\nUnperturbed Trajectory:")
    print(f"  Separation from HD 286941: {sep_unperturbed:.6f}° ({sep_unperturbed*60:.2f} arcmin)")
    
    # Backwards integration separation
    if backwards_origin:
        sep_backwards = calculate_angular_separation(
            backwards_origin['ra_deg'],
            backwards_origin['dec_deg'],
            HD_286941['ra_deg'],
            HD_286941['dec_deg']
        )
        
        print(f"\nBackwards Integration:")
        print(f"  Separation from HD 286941: {sep_backwards:.6f}° ({sep_backwards*60:.2f} arcmin)")
        
        if sep_backwards < 1.0:
            print(f"  ✓ VERY CLOSE - Backwards integration points to HD 286941!")
        elif sep_backwards < 5.0:
            print(f"  ✓ Close - Backwards integration points close to HD 286941")
        elif sep_backwards < 10.0:
            print(f"  ⚠ Somewhat close - Possible alignment")
        else:
            print(f"  ✗ Not close - Backwards integration does not point to HD 286941")
    
    # Compare with decoded coordinate
    print(f"\n{'='*80}")
    print("COMPARING WITH DECODED COORDINATE")
    print("="*80)
    
    decoded_file = Path('3i_atlas_data/decoded_coordinates.json')
    if decoded_file.exists():
        with open(decoded_file) as f:
            decoded = json.load(f)
        
        primary = decoded.get('primary_coordinate')
        if primary:
            decoded_ra = primary['ra_deg']
            decoded_dec = primary['dec_deg']
            
            print(f"\nDecoded Coordinate:")
            print(f"  RA: {decoded_ra:.6f}°")
            print(f"  Dec: {decoded_dec:.6f}°")
            
            # Unperturbed vs decoded
            sep_unperturbed_decoded = calculate_angular_separation(
                ra_unperturbed,
                dec_unperturbed,
                decoded_ra,
                decoded_dec
            )
            
            print(f"\nUnperturbed Asymptote vs Decoded Coordinate:")
            print(f"  Separation: {sep_unperturbed_decoded:.6f}° ({sep_unperturbed_decoded*60:.2f} arcmin)")
            
            # Backwards vs decoded
            if backwards_origin:
                sep_backwards_decoded = calculate_angular_separation(
                    backwards_origin['ra_deg'],
                    backwards_origin['dec_deg'],
                    decoded_ra,
                    decoded_dec
                )
                
                print(f"\nBackwards Integration vs Decoded Coordinate:")
                print(f"  Separation: {sep_backwards_decoded:.6f}° ({sep_backwards_decoded*60:.2f} arcmin)")
                
                if sep_backwards_decoded < 1.0:
                    print(f"  ✓ VERY CLOSE - Backwards integration aligns with decoded coordinate!")
                elif sep_backwards_decoded < 5.0:
                    print(f"  ✓ Close - Backwards integration aligns with decoded coordinate")
                elif sep_backwards_decoded < 10.0:
                    print(f"  ⚠ Somewhat close - Possible alignment")
                else:
                    print(f"  ✗ Not close - Backwards integration does not align with decoded coordinate")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'unperturbed_origin': {
            'ra_deg': ra_unperturbed,
            'dec_deg': dec_unperturbed,
        },
        'backwards_integration_origin': backwards_origin,
        'perturbation_estimates': {
            'jupiter_deflection_deg': deflection_deg,
            'jupiter_deflection_arcmin': deflection_deg * 60,
        },
        'hd286941': HD_286941,
        'separations': {
            'unperturbed_vs_hd286941_deg': sep_unperturbed,
            'unperturbed_vs_hd286941_arcmin': sep_unperturbed * 60,
            'backwards_vs_hd286941_deg': sep_backwards if backwards_origin else None,
            'backwards_vs_hd286941_arcmin': (sep_backwards * 60) if backwards_origin else None,
        },
    }
    
    output_file = Path('3i_atlas_data/nbody_reverse_engineered_origin.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_file}")
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("="*80)
    
    print(f"\nCan we reverse engineer the origin point?")
    print(f"\nYES - We can work backwards from the trajectory.")
    
    print(f"\n1. Unperturbed Trajectory:")
    print(f"   Incoming asymptote: RA {ra_unperturbed:.6f}°, Dec {dec_unperturbed:.6f}°")
    print(f"   Separation from HD 286941: {sep_unperturbed:.2f}° ({sep_unperturbed*60:.2f} arcmin)")
    print(f"   Status: {'✓ Points to HD 286941' if sep_unperturbed < 5 else '✗ Does not point to HD 286941'}")
    
    if backwards_origin:
        print(f"\n2. Backwards Integration (Simplified N-body):")
        print(f"   Origin: RA {backwards_origin['ra_deg']:.6f}°, Dec {backwards_origin['dec_deg']:.6f}°")
        print(f"   Separation from HD 286941: {sep_backwards:.2f}° ({sep_backwards*60:.2f} arcmin)")
        print(f"   Status: {'✓ Points to HD 286941' if sep_backwards < 5 else '✗ Does not point to HD 286941'}")
    
    print(f"\n3. Perturbation Estimates:")
    print(f"   Jupiter deflection: {deflection_deg:.6f}° ({deflection_deg*60:.2f} arcmin)")
    print(f"   Note: Full N-body simulation needed for precise calculation")
    
    print(f"\nConclusion:")
    if backwards_origin and sep_backwards < 5:
        print(f"  ✓ Backwards integration points close to HD 286941")
        print(f"  • If perturbations are accounted for, HD 286941 could be the origin")
        print(f"  • Need full N-body simulation to verify")
    elif sep_unperturbed < 5:
        print(f"  ✓ Unperturbed trajectory points close to HD 286941")
        print(f"  • HD 286941 could be the origin if perturbations are small")
    else:
        print(f"  ✗ Neither unperturbed nor backwards integration points to HD 286941")
        print(f"  • Trajectory does not point back to HD 286941")
        print(f"  • But decoded coordinates still point to HD 286941 (highly significant)")
        print(f"  • Coordinates may encode something else (reference point, target, etc.)")
    
    print(f"\n{'='*80}")
    
    return results

if __name__ == "__main__":
    nbody_reverse_engineer_origin()

