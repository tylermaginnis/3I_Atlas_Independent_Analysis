#!/usr/bin/env python3
"""
Analyze Unperturbed 3I/ATLAS Trajectory
Calculate incoming asymptote direction (no gravitational interactions)
Check if it points to HD 286941
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

def calculate_incoming_asymptote_direction(orbital_elements):
    """
    Calculate incoming asymptote direction for hyperbolic orbit
    (direction from which object came, assuming no gravitational interactions)
    
    For hyperbolic orbit:
    - True anomaly at infinity: ν_inf = arccos(-1/e)
    - Incoming asymptote is at ν = -ν_inf (before perihelion)
    - Outgoing asymptote is at ν = +ν_inf (after perihelion)
    """
    q = orbital_elements['perihelion_distance_au']
    e = orbital_elements['eccentricity']
    i = radians(orbital_elements['inclination_deg'])
    omega = radians(orbital_elements['longitude_ascending_node_deg'])
    w = radians(orbital_elements['argument_of_perihelion_deg'])
    
    # True anomaly at infinity for hyperbolic orbit
    # cos(ν_inf) = -1/e
    nu_inf = acos(-1.0 / e)
    
    # Incoming asymptote: ν = -ν_inf (before perihelion)
    nu_incoming = -nu_inf
    
    # Distance at infinity (very large)
    # For hyperbolic: r = a(e*cosh(H) - 1) → ∞ as H → ∞
    # At asymptote, we use the direction, not distance
    
    # Position in orbital plane (perifocal coordinates) at incoming asymptote
    # At infinity, the direction is given by the asymptote angle
    # We use a large distance to get the direction
    r_large = 1000.0  # AU (large distance for direction)
    
    # Position in perifocal frame
    x_peri = r_large * cos(nu_incoming)
    y_peri = r_large * sin(nu_incoming)
    z_peri = 0
    
    # Transform to ecliptic coordinates
    # Rotation order: ω (argument of perihelion), i (inclination), Ω (ascending node)
    cos_w = cos(w)
    sin_w = sin(w)
    cos_omega = cos(omega)
    sin_omega = sin(omega)
    cos_i = cos(i)
    sin_i = sin(i)
    
    # Step 1: Rotate by argument of perihelion ω (around z-axis in perifocal frame)
    x1 = x_peri * cos_w - y_peri * sin_w
    y1 = x_peri * sin_w + y_peri * cos_w
    z1 = z_peri
    
    # Step 2: Rotate by inclination i (around x-axis)
    x2 = x1
    y2 = y1 * cos_i - z1 * sin_i
    z2 = y1 * sin_i + z1 * cos_i
    
    # Step 3: Rotate by longitude of ascending node Ω (around z-axis)
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
    
    # Convert to RA and Dec (direction vector)
    r_3d = sqrt(x_eq**2 + y_eq**2 + z_eq**2)
    dec = degrees(asin(z_eq / r_3d))
    ra = degrees(atan2(y_eq, x_eq))
    if ra < 0:
        ra += 360
    
    return {
        'ra_deg': ra,
        'dec_deg': dec,
        'true_anomaly_inf_deg': degrees(nu_inf),
        'true_anomaly_incoming_deg': degrees(nu_incoming),
        'direction_vector': {
            'x': x_eq / r_3d,
            'y': y_eq / r_3d,
            'z': z_eq / r_3d,
        }
    }

def analyze_unperturbed_trajectory():
    """Analyze unperturbed trajectory (no gravitational interactions)"""
    print("="*80)
    print("ANALYZING UNPERTURBED 3I/ATLAS TRAJECTORY")
    print("(Assuming NO gravitational interactions)")
    print("="*80)
    
    print(f"\nOrbital Elements:")
    print(f"  Perihelion distance: {ORBITAL_ELEMENTS['perihelion_distance_au']:.6f} AU")
    print(f"  Eccentricity: {ORBITAL_ELEMENTS['eccentricity']:.6f}")
    print(f"  Inclination: {ORBITAL_ELEMENTS['inclination_deg']:.6f}°")
    print(f"  Longitude of ascending node: {ORBITAL_ELEMENTS['longitude_ascending_node_deg']:.6f}°")
    print(f"  Argument of perihelion: {ORBITAL_ELEMENTS['argument_of_perihelion_deg']:.6f}°")
    
    # Calculate incoming asymptote direction
    print(f"\n{'='*80}")
    print("CALCULATING INCOMING ASYMPTOTE DIRECTION")
    print("="*80)
    
    incoming = calculate_incoming_asymptote_direction(ORBITAL_ELEMENTS)
    
    print(f"\nIncoming Asymptote Direction (unperturbed):")
    print(f"  RA: {incoming['ra_deg']:.6f}°")
    print(f"  Dec: {incoming['dec_deg']:.6f}°")
    print(f"  True anomaly at infinity: {incoming['true_anomaly_inf_deg']:.6f}°")
    print(f"  True anomaly (incoming): {incoming['true_anomaly_incoming_deg']:.6f}°")
    
    # Compare with HD 286941
    print(f"\n{'='*80}")
    print("COMPARING WITH HD 286941")
    print("="*80)
    
    print(f"\nHD 286941 Location:")
    print(f"  RA: {HD_286941['ra_deg']:.6f}°")
    print(f"  Dec: {HD_286941['dec_deg']:.6f}°")
    
    # Calculate angular separation
    sep = calculate_angular_separation(
        incoming['ra_deg'],
        incoming['dec_deg'],
        HD_286941['ra_deg'],
        HD_286941['dec_deg']
    )
    
    print(f"\nAngular Separation:")
    print(f"  Separation: {sep:.6f}° ({sep*60:.2f} arcmin)")
    
    if sep < 1.0:
        print(f"  ✓ VERY CLOSE - Strong evidence for origin from HD 286941!")
    elif sep < 5.0:
        print(f"  ✓ Close - Good evidence for origin from HD 286941")
    elif sep < 10.0:
        print(f"  ⚠ Somewhat close - Possible origin from HD 286941")
    elif sep < 30.0:
        print(f"  ⚠ Not very close - Weak evidence for origin")
    else:
        print(f"  ✗ Not close - Unlikely origin from HD 286941")
    
    # Compare with decoded coordinate
    print(f"\n{'='*80}")
    print("COMPARING WITH DECODED COORDINATE")
    print("="*80)
    
    # Load decoded coordinates
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
            
            # Separation from incoming asymptote
            sep_asymptote = calculate_angular_separation(
                incoming['ra_deg'],
                incoming['dec_deg'],
                decoded_ra,
                decoded_dec
            )
            
            print(f"\nSeparation from Incoming Asymptote:")
            print(f"  Separation: {sep_asymptote:.6f}° ({sep_asymptote*60:.2f} arcmin)")
            
            if sep_asymptote < 1.0:
                print(f"  ✓ VERY CLOSE - Decoded coordinate aligns with incoming asymptote!")
            elif sep_asymptote < 5.0:
                print(f"  ✓ Close - Decoded coordinate aligns with incoming asymptote")
            elif sep_asymptote < 10.0:
                print(f"  ⚠ Somewhat close - Possible alignment")
            else:
                print(f"  ✗ Not close - Decoded coordinate does not align with incoming asymptote")
            
            # Separation from HD 286941
            sep_hd = calculate_angular_separation(
                HD_286941['ra_deg'],
                HD_286941['dec_deg'],
                decoded_ra,
                decoded_dec
            )
            
            print(f"\nSeparation from HD 286941:")
            print(f"  Separation: {sep_hd:.6f}° ({sep_hd*60:.2f} arcmin)")
            print(f"  Status: ✓ Very close (from previous analysis)")
    
    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print("="*80)
    
    print(f"\nIncoming Asymptote Direction (unperturbed trajectory):")
    print(f"  RA: {incoming['ra_deg']:.6f}°")
    print(f"  Dec: {incoming['dec_deg']:.6f}°")
    
    print(f"\nHD 286941 Location:")
    print(f"  RA: {HD_286941['ra_deg']:.6f}°")
    print(f"  Dec: {HD_286941['dec_deg']:.6f}°")
    
    print(f"\nAngular Separation: {sep:.6f}° ({sep*60:.2f} arcmin)")
    
    if sep < 1.0:
        print(f"\n✓ STRONG EVIDENCE: Incoming asymptote points directly to HD 286941!")
        print(f"  If gravitational interactions did not affect the trajectory,")
        print(f"  3I/ATLAS likely originated from HD 286941.")
    elif sep < 5.0:
        print(f"\n✓ GOOD EVIDENCE: Incoming asymptote points close to HD 286941")
        print(f"  If gravitational interactions did not affect the trajectory,")
        print(f"  3I/ATLAS may have originated from HD 286941.")
    elif sep < 10.0:
        print(f"\n⚠ POSSIBLE: Incoming asymptote points somewhat close to HD 286941")
        print(f"  If gravitational interactions did not affect the trajectory,")
        print(f"  3I/ATLAS might have originated from HD 286941.")
    else:
        print(f"\n✗ WEAK EVIDENCE: Incoming asymptote does not point to HD 286941")
        print(f"  Even without gravitational interactions,")
        print(f"  the trajectory does not point back to HD 286941.")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'assumption': 'No gravitational interactions (unperturbed trajectory)',
        'incoming_asymptote': incoming,
        'hd286941': HD_286941,
        'angular_separation_deg': sep,
        'angular_separation_arcmin': sep * 60,
        'interpretation': {
            'very_close': sep < 1.0,
            'close': sep < 5.0,
            'somewhat_close': sep < 10.0,
            'not_close': sep >= 10.0,
        }
    }
    
    output_file = Path('3i_atlas_data/unperturbed_trajectory_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_file}")
    
    return results

if __name__ == "__main__":
    analyze_unperturbed_trajectory()

