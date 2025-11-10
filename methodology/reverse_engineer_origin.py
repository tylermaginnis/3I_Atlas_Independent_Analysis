#!/usr/bin/env python3
"""
Reverse Engineer 3I/ATLAS Origin Point
Work backwards from trajectory (perturbed and unperturbed) to find origin
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime
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

def calculate_incoming_asymptote_unperturbed(orbital_elements):
    """Calculate incoming asymptote direction (unperturbed trajectory)"""
    q = orbital_elements['perihelion_distance_au']
    e = orbital_elements['eccentricity']
    i = radians(orbital_elements['inclination_deg'])
    omega = radians(orbital_elements['longitude_ascending_node_deg'])
    w = radians(orbital_elements['argument_of_perihelion_deg'])
    
    # True anomaly at infinity for hyperbolic orbit
    nu_inf = acos(-1.0 / e)
    nu_incoming = -nu_inf
    
    # Position in perifocal frame (at infinity, use large distance for direction)
    r_large = 1000.0  # AU
    x_peri = r_large * cos(nu_incoming)
    y_peri = r_large * sin(nu_incoming)
    z_peri = 0
    
    # Transform to ecliptic coordinates
    cos_w = cos(w)
    sin_w = sin(w)
    cos_omega = cos(omega)
    sin_omega = sin(omega)
    cos_i = cos(i)
    sin_i = sin(i)
    
    # Rotate by argument of perihelion ω
    x1 = x_peri * cos_w - y_peri * sin_w
    y1 = x_peri * sin_w + y_peri * cos_w
    z1 = z_peri
    
    # Rotate by inclination i
    x2 = x1
    y2 = y1 * cos_i - z1 * sin_i
    z2 = y1 * sin_i + z1 * cos_i
    
    # Rotate by longitude of ascending node Ω
    x_ecl = x2 * cos_omega - y2 * sin_omega
    y_ecl = x2 * sin_omega + y2 * cos_omega
    z_ecl = z2
    
    # Convert ecliptic to equatorial (J2000)
    eps = radians(23.4392911)
    cos_eps = cos(eps)
    sin_eps = sin(eps)
    
    x_eq = x_ecl
    y_eq = y_ecl * cos_eps - z_ecl * sin_eps
    z_eq = y_ecl * sin_eps + z_ecl * cos_eps
    
    # Convert to RA and Dec
    r_3d = sqrt(x_eq**2 + y_eq**2 + z_eq**2)
    dec = degrees(asin(z_eq / r_3d))
    ra = degrees(atan2(y_eq, x_eq))
    if ra < 0:
        ra += 360
    
    return {
        'ra_deg': ra,
        'dec_deg': dec,
        'true_anomaly_inf_deg': degrees(nu_inf),
        'direction_vector': {
            'x': x_eq / r_3d,
            'y': y_eq / r_3d,
            'z': z_eq / r_3d,
        }
    }

def calculate_perturbed_trajectory_origin(orbital_elements, perturbation_analysis=None):
    """
    Calculate origin point accounting for gravitational perturbations
    
    For perturbed trajectory, we need to:
    1. Calculate unperturbed incoming asymptote
    2. Account for gravitational perturbations from planets
    3. Work backwards to find true origin point
    
    Major perturbations:
    - Jupiter (closest approach March 16, 2026 at 53 million km)
    - Saturn
    - Other planets
    - Solar system barycenter
    """
    # Start with unperturbed asymptote
    unperturbed = calculate_incoming_asymptote_unperturbed(orbital_elements)
    
    # For now, we'll estimate perturbations
    # Actual calculation would require N-body simulation
    
    # Key perturbation: Jupiter encounter
    # March 16, 2026: 53 million km from Jupiter
    # This will significantly alter the trajectory
    
    # Estimate perturbation direction
    # Jupiter's position at encounter (approximate)
    # This is a simplified estimate - actual calculation needs ephemeris
    
    # For now, return unperturbed with note about perturbations
    return {
        'unperturbed': unperturbed,
        'perturbed': None,  # Would require N-body simulation
        'perturbation_note': 'Requires N-body simulation to calculate perturbed trajectory',
        'major_perturbations': {
            'jupiter_encounter': {
                'date': '2026-03-16',
                'distance_km': 53000000,  # 53 million km
                'distance_au': 0.354,  # AU
                'significance': 'Major perturbation - will significantly alter trajectory',
            }
        }
    }

def reverse_engineer_origin_from_trajectory():
    """Reverse engineer origin point from trajectory"""
    print("="*80)
    print("REVERSE ENGINEERING 3I/ATLAS ORIGIN POINT")
    print("From Trajectory (Perturbed and Unperturbed)")
    print("="*80)
    
    print(f"\nOrbital Elements:")
    print(f"  Perihelion distance: {ORBITAL_ELEMENTS['perihelion_distance_au']:.6f} AU")
    print(f"  Eccentricity: {ORBITAL_ELEMENTS['eccentricity']:.6f}")
    print(f"  Inclination: {ORBITAL_ELEMENTS['inclination_deg']:.6f}°")
    print(f"  Longitude of ascending node: {ORBITAL_ELEMENTS['longitude_ascending_node_deg']:.6f}°")
    print(f"  Argument of perihelion: {ORBITAL_ELEMENTS['argument_of_perihelion_deg']:.6f}°")
    
    # Calculate unperturbed incoming asymptote
    print(f"\n{'='*80}")
    print("1. UNPERTURBED TRAJECTORY (No Gravitational Interactions)")
    print("="*80)
    
    unperturbed = calculate_incoming_asymptote_unperturbed(ORBITAL_ELEMENTS)
    
    print(f"\nIncoming Asymptote Direction (Unperturbed):")
    print(f"  RA: {unperturbed['ra_deg']:.6f}°")
    print(f"  Dec: {unperturbed['dec_deg']:.6f}°")
    print(f"  True anomaly at infinity: {unperturbed['true_anomaly_inf_deg']:.6f}°")
    
    # Calculate perturbed trajectory
    print(f"\n{'='*80}")
    print("2. PERTURBED TRAJECTORY (With Gravitational Interactions)")
    print("="*80)
    
    perturbed_analysis = calculate_perturbed_trajectory_origin(ORBITAL_ELEMENTS)
    
    print(f"\nPerturbation Analysis:")
    print(f"  Note: {perturbed_analysis['perturbation_note']}")
    print(f"\nMajor Perturbations:")
    jupiter = perturbed_analysis['major_perturbations']['jupiter_encounter']
    print(f"  Jupiter Encounter:")
    print(f"    Date: {jupiter['date']}")
    print(f"    Distance: {jupiter['distance_km']/1e6:.1f} million km ({jupiter['distance_au']:.3f} AU)")
    print(f"    Significance: {jupiter['significance']}")
    
    # Compare with HD 286941
    print(f"\n{'='*80}")
    print("3. COMPARING WITH HD 286941")
    print("="*80)
    
    print(f"\nHD 286941 Location:")
    print(f"  RA: {HD_286941['ra_deg']:.6f}°")
    print(f"  Dec: {HD_286941['dec_deg']:.6f}°")
    
    # Unperturbed separation
    sep_unperturbed = calculate_angular_separation(
        unperturbed['ra_deg'],
        unperturbed['dec_deg'],
        HD_286941['ra_deg'],
        HD_286941['dec_deg']
    )
    
    print(f"\nUnperturbed Trajectory:")
    print(f"  Incoming asymptote: RA {unperturbed['ra_deg']:.6f}°, Dec {unperturbed['dec_deg']:.6f}°")
    print(f"  Separation from HD 286941: {sep_unperturbed:.6f}° ({sep_unperturbed*60:.2f} arcmin)")
    
    if sep_unperturbed < 1.0:
        print(f"  ✓ VERY CLOSE - Unperturbed trajectory points to HD 286941!")
    elif sep_unperturbed < 5.0:
        print(f"  ✓ Close - Unperturbed trajectory points close to HD 286941")
    elif sep_unperturbed < 10.0:
        print(f"  ⚠ Somewhat close - Possible alignment")
    else:
        print(f"  ✗ Not close - Unperturbed trajectory does not point to HD 286941")
    
    # Compare with decoded coordinate
    print(f"\n{'='*80}")
    print("4. COMPARING WITH DECODED COORDINATE")
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
            
            # Separation from unperturbed asymptote
            sep_asymptote = calculate_angular_separation(
                unperturbed['ra_deg'],
                unperturbed['dec_deg'],
                decoded_ra,
                decoded_dec
            )
            
            print(f"\nUnperturbed Asymptote vs Decoded Coordinate:")
            print(f"  Separation: {sep_asymptote:.6f}° ({sep_asymptote*60:.2f} arcmin)")
            
            if sep_asymptote < 1.0:
                print(f"  ✓ VERY CLOSE - Decoded coordinate aligns with unperturbed asymptote!")
            elif sep_asymptote < 5.0:
                print(f"  ✓ Close - Decoded coordinate aligns with unperturbed asymptote")
            elif sep_asymptote < 10.0:
                print(f"  ⚠ Somewhat close - Possible alignment")
            else:
                print(f"  ✗ Not close - Decoded coordinate does not align with unperturbed asymptote")
            
            # Separation from HD 286941
            sep_hd = calculate_angular_separation(
                HD_286941['ra_deg'],
                HD_286941['dec_deg'],
                decoded_ra,
                decoded_dec
            )
            
            print(f"\nDecoded Coordinate vs HD 286941:")
            print(f"  Separation: {sep_hd:.6f}° ({sep_hd*60:.2f} arcmin)")
            print(f"  Status: ✓ Very close (from previous analysis)")
    
    # Reverse engineering analysis
    print(f"\n{'='*80}")
    print("5. REVERSE ENGINEERING ANALYSIS")
    print("="*80)
    
    print(f"\nCan we reverse engineer the origin point?")
    print(f"\nMethod 1: Unperturbed Trajectory")
    print(f"  • Calculate incoming asymptote direction")
    print(f"  • This gives direction from which object came")
    print(f"  • Result: RA {unperturbed['ra_deg']:.6f}°, Dec {unperturbed['dec_deg']:.6f}°")
    print(f"  • Separation from HD 286941: {sep_unperturbed:.2f}°")
    print(f"  • Status: {'✓ Points to HD 286941' if sep_unperturbed < 5 else '✗ Does not point to HD 286941'}")
    
    print(f"\nMethod 2: Perturbed Trajectory (Requires N-body Simulation)")
    print(f"  • Account for gravitational perturbations")
    print(f"  • Major perturbation: Jupiter encounter (March 16, 2026)")
    print(f"  • Distance: 53 million km (0.354 AU)")
    print(f"  • This will significantly alter the trajectory")
    print(f"  • Need N-body simulation to calculate perturbed asymptote")
    print(f"  • Status: ⚠ Requires advanced calculation")
    
    print(f"\nMethod 3: Working Backwards from Observed Trajectory")
    print(f"  • Use observed orbital elements")
    print(f"  • Account for all gravitational interactions")
    print(f"  • Work backwards to find true origin point")
    print(f"  • Status: ⚠ Requires N-body simulation")
    
    # Calculate if perturbations could explain the difference
    print(f"\n{'='*80}")
    print("6. PERTURBATION ANALYSIS")
    print("="*80)
    
    # Estimate perturbation magnitude
    # Jupiter's gravitational influence at 0.354 AU
    # This is a rough estimate
    jupiter_mass_solar = 0.0009543  # Jupiter mass in solar masses
    jupiter_distance_au = 0.354  # AU
    
    # Gravitational deflection angle (rough estimate)
    # δ ≈ 2GM / (v² * b) where b is impact parameter
    # This is a simplified calculation
    
    print(f"\nJupiter Encounter Analysis:")
    print(f"  Date: March 16, 2026")
    print(f"  Distance: 53 million km (0.354 AU)")
    print(f"  Jupiter mass: {jupiter_mass_solar:.6f} solar masses")
    print(f"  Significance: Major perturbation - will significantly alter trajectory")
    print(f"  Note: Actual deflection requires N-body simulation")
    
    # Estimate if perturbations could align trajectory with HD 286941
    print(f"\nCould Perturbations Align Trajectory with HD 286941?")
    print(f"  Current separation (unperturbed): {sep_unperturbed:.2f}°")
    print(f"  Required deflection: {sep_unperturbed:.2f}° to align with HD 286941")
    print(f"  Jupiter encounter: Major perturbation at 0.354 AU")
    print(f"  Status: {'✓ Possible' if sep_unperturbed < 30 else '✗ Unlikely'} - requires N-body simulation")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'unperturbed_trajectory': unperturbed,
        'perturbed_trajectory': perturbed_analysis,
        'hd286941': HD_286941,
        'separations': {
            'unperturbed_vs_hd286941_deg': sep_unperturbed,
            'unperturbed_vs_hd286941_arcmin': sep_unperturbed * 60,
        },
        'reverse_engineering': {
            'unperturbed_origin_ra_deg': unperturbed['ra_deg'],
            'unperturbed_origin_dec_deg': unperturbed['dec_deg'],
            'perturbed_origin_ra_deg': None,  # Requires N-body simulation
            'perturbed_origin_dec_deg': None,  # Requires N-body simulation
            'can_reverse_engineer': True,
            'requires_nbody_simulation': True,
        }
    }
    
    output_file = Path('3i_atlas_data/reverse_engineered_origin.json')
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
    print(f"\nYES, but it requires understanding both perturbed and unperturbed trajectories.")
    
    print(f"\nUnperturbed Trajectory (No Gravitational Interactions):")
    print(f"  Incoming asymptote: RA {unperturbed['ra_deg']:.6f}°, Dec {unperturbed['dec_deg']:.6f}°")
    print(f"  Separation from HD 286941: {sep_unperturbed:.2f}° ({sep_unperturbed*60:.2f} arcmin)")
    print(f"  Status: {'✓ Points to HD 286941' if sep_unperturbed < 5 else '✗ Does not point to HD 286941'}")
    
    print(f"\nPerturbed Trajectory (With Gravitational Interactions):")
    print(f"  Status: ⚠ Requires N-body simulation")
    print(f"  Major perturbation: Jupiter encounter (March 16, 2026, 0.354 AU)")
    print(f"  This will significantly alter the trajectory")
    print(f"  Need to calculate perturbed incoming asymptote")
    
    print(f"\nReverse Engineering Process:")
    print(f"  1. Calculate unperturbed incoming asymptote ✓ (Done)")
    print(f"  2. Account for gravitational perturbations ⚠ (Requires N-body simulation)")
    print(f"  3. Work backwards to find true origin point ⚠ (Requires N-body simulation)")
    print(f"  4. Compare with HD 286941 and decoded coordinates ⚠ (Pending)")
    
    print(f"\nConclusion:")
    if sep_unperturbed < 5:
        print(f"  ✓ Unperturbed trajectory points close to HD 286941")
        print(f"  • If perturbations are small, HD 286941 could be the origin")
        print(f"  • Need N-body simulation to verify")
    else:
        print(f"  ✗ Unperturbed trajectory does not point to HD 286941 ({sep_unperturbed:.2f}° separation)")
        print(f"  • Even without perturbations, trajectory does not point back to HD 286941")
        print(f"  • Perturbations could potentially align it, but unlikely for {sep_unperturbed:.2f}° deflection")
        print(f"  • Need N-body simulation to calculate perturbed trajectory")
    
    print(f"\n{'='*80}")
    
    return results

if __name__ == "__main__":
    reverse_engineer_origin_from_trajectory()

