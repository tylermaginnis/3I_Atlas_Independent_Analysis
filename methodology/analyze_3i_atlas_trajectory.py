#!/usr/bin/env python3
"""
Analyze 3I/ATLAS Trajectory from Orbital Elements
Cross-reference with decoded coordinates and HD 286941
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from math import cos, sin, acos, asin, atan2, sqrt, radians, degrees
import math

# For hyperbolic orbits
try:
    from math import asinh
except ImportError:
    # Fallback for older Python versions
    def asinh(x):
        return math.log(x + math.sqrt(x*x + 1))

# 3I/ATLAS Orbital Elements (from Seligman et al. 2025)
ORBITAL_ELEMENTS = {
    'time_of_perihelion_mjd': 60977.10275,  # MJD
    'perihelion_distance_au': 1.3797753,  # AU
    'eccentricity': 6.320503,
    'inclination_deg': 175.12093,  # degrees
    'longitude_ascending_node_deg': 322.34889,  # Ω
    'argument_of_perihelion_deg': 127.77350,  # ω
    'last_updated': 'July 6, 2025',
    'source': 'Seligman et al. (2025)',
}

def mjd_to_datetime(mjd):
    """Convert Modified Julian Date to datetime"""
    # MJD = JD - 2400000.5
    # JD epoch is January 1, 4713 BC
    # MJD epoch is November 17, 1858
    mjd_epoch = datetime(1858, 11, 17)
    return mjd_epoch + timedelta(days=mjd)

def calculate_orbital_position(t_mjd, orbital_elements):
    """
    Calculate 3D position of 3I/ATLAS at given time
    
    Args:
        t_mjd: Time in Modified Julian Date
        orbital_elements: Dictionary with orbital elements
    
    Returns:
        Dictionary with position (x, y, z) in AU and (ra, dec) in degrees
    """
    # Extract orbital elements
    t_peri = orbital_elements['time_of_perihelion_mjd']
    q = orbital_elements['perihelion_distance_au']  # perihelion distance
    e = orbital_elements['eccentricity']  # eccentricity
    i = radians(orbital_elements['inclination_deg'])  # inclination
    omega = radians(orbital_elements['longitude_ascending_node_deg'])  # Ω
    w = radians(orbital_elements['argument_of_perihelion_deg'])  # ω
    
    # Calculate time since perihelion (in days)
    dt = t_mjd - t_peri
    
    # For hyperbolic orbit (e > 1), we need to solve for true anomaly
    # Semi-major axis for hyperbolic orbit: a = q / (e - 1) (negative)
    a = q / (e - 1)
    
    # Mean motion for hyperbolic orbit
    # n = sqrt(GM / |a|^3) where GM = 0.01720209895^2 AU^3 / day^2 (for Sun)
    # Actually, using k = 0.01720209895 (Gaussian gravitational constant)
    k = 0.01720209895  # AU^(3/2) / day
    GM = k * k  # AU^3 / day^2
    n = sqrt(GM / abs(a)**3)
    
    # Mean anomaly (preserve sign for hyperbolic)
    M = n * dt
    
    # For hyperbolic orbit, solve M = e * sinh(H) - H for H (hyperbolic anomaly)
    # Use Newton-Raphson method with better initial guess
    if abs(M) < 1e-6:
        H = 0.0
    else:
        # Initial guess: preserve sign of M
        # For small M: H ≈ M/e
        # For larger M: H ≈ asinh(M/e) (preserves sign)
        if abs(M) < 1.0:
            H = M / e
        else:
            # asinh preserves sign
            H = asinh(M / e) if M >= 0 else -asinh(-M / e)
        
        # Newton-Raphson iteration
        for iteration in range(50):
            try:
                sinh_H = math.sinh(H)
                cosh_H = math.cosh(H)
                f = e * sinh_H - H - M
                df = e * cosh_H - 1
                
                if abs(df) < 1e-12:
                    break
                
                H_new = H - f / df
                
                # Check for convergence
                if abs(H_new - H) < 1e-12:
                    H = H_new
                    break
                
                # Limit H to prevent overflow
                if abs(H_new) > 20:
                    H_new = math.copysign(20, H_new)
                
                H = H_new
            except (OverflowError, ZeroDivisionError, ValueError):
                # Fallback: use approximation
                H = asinh(M / e)
                break
    
    # True anomaly from hyperbolic anomaly
    # tan(ν/2) = sqrt((e+1)/(e-1)) * tanh(H/2)
    if abs(H) < 1e-10:
        nu = 0.0
    else:
        sqrt_factor = sqrt((e + 1) / (e - 1))
        tanh_H_2 = math.tanh(H / 2)
        tan_nu_2 = sqrt_factor * tanh_H_2
        nu = 2 * atan2(tan_nu_2, 1)  # true anomaly
    
    # Distance from Sun (hyperbolic orbit)
    # r = a(1 - e*cosh(H)) for hyperbolic, or r = q(1+e)/(1+e*cos(nu))
    r = q * (1 + e) / (1 + e * cos(nu))
    
    # Check for valid distance
    if r < 0 or not (r > 0 and r < 1e6):
        # Fallback calculation
        r = abs(a) * (e * math.cosh(H) - 1)
        if r < 0:
            r = q  # Use perihelion distance as fallback
    
    # Position in orbital plane (perifocal coordinates)
    x_peri = r * cos(nu)
    y_peri = r * sin(nu)
    z_peri = 0
    
    # Transform to ecliptic coordinates
    # Standard rotation order: Ω (ascending node), i (inclination), ω (argument of perihelion)
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
    # Ecliptic obliquity at J2000: 23.4392911 degrees
    eps = radians(23.4392911)
    cos_eps = cos(eps)
    sin_eps = sin(eps)
    
    x_eq = x_ecl
    y_eq = y_ecl * cos_eps - z_ecl * sin_eps
    z_eq = y_ecl * sin_eps + z_ecl * cos_eps
    
    # Convert to RA and Dec
    r_3d = sqrt(x_eq**2 + y_eq**2 + z_eq**2)
    
    # Check for valid distance
    if r_3d < 1e-10 or math.isnan(r_3d) or math.isinf(r_3d):
        # Return None to indicate calculation failure
        return None
    
    dec = degrees(asin(z_eq / r_3d))
    ra = degrees(atan2(y_eq, x_eq))
    if ra < 0:
        ra += 360
    
    # Final NaN checks
    if any(math.isnan(val) or math.isinf(val) for val in [ra, dec, r, r_3d]):
        # Return None to indicate calculation failure
        return None
    
    return {
        'time_mjd': t_mjd,
        'time_datetime': mjd_to_datetime(t_mjd),
        'position_au': {
            'x': x_eq,
            'y': y_eq,
            'z': z_eq,
            'r': r_3d,
        },
        'coordinates': {
            'ra_deg': ra,
            'dec_deg': dec,
        },
        'true_anomaly_deg': degrees(nu),
        'distance_from_sun_au': r,
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

def find_closest_approach_to_target(orbital_elements, target_ra, target_dec, 
                                     t_start_mjd, t_end_mjd, n_points=1000):
    """
    Find closest approach of 3I/ATLAS to a target coordinate
    
    Args:
        orbital_elements: Dictionary with orbital elements
        target_ra: Target RA in degrees
        target_dec: Target Dec in degrees
        t_start_mjd: Start time in MJD
        t_end_mjd: End time in MJD
        n_points: Number of points to sample
    
    Returns:
        Dictionary with closest approach information
    """
    times = np.linspace(t_start_mjd, t_end_mjd, n_points)
    min_separation = float('inf')
    closest_time = None
    closest_position = None
    
    for t in times:
        pos = calculate_orbital_position(t, orbital_elements)
        
        # Skip if calculation failed
        if pos is None:
            continue
        
        ra = pos['coordinates']['ra_deg']
        dec = pos['coordinates']['dec_deg']
        
        # Skip if coordinates are invalid
        if math.isnan(ra) or math.isnan(dec):
            continue
        
        separation = calculate_angular_separation(ra, dec, target_ra, target_dec)
        
        if separation < min_separation:
            min_separation = separation
            closest_time = t
            closest_position = pos
    
    return {
        'closest_time_mjd': closest_time,
        'closest_time_datetime': mjd_to_datetime(closest_time) if closest_time else None,
        'closest_separation_deg': min_separation,
        'closest_separation_arcmin': min_separation * 60,
        'closest_position': closest_position,
        'target_ra': target_ra,
        'target_dec': target_dec,
    }

def analyze_3i_atlas_trajectory():
    """Analyze 3I/ATLAS trajectory and cross-reference with decoded coordinates"""
    print("="*80)
    print("3I/ATLAS TRAJECTORY ANALYSIS")
    print("="*80)
    
    print(f"\nOrbital Elements (from Seligman et al. 2025):")
    print(f"  Time of perihelion: MJD {ORBITAL_ELEMENTS['time_of_perihelion_mjd']:.5f}")
    print(f"  Perihelion distance: {ORBITAL_ELEMENTS['perihelion_distance_au']:.6f} AU")
    print(f"  Eccentricity: {ORBITAL_ELEMENTS['eccentricity']:.6f}")
    print(f"  Inclination: {ORBITAL_ELEMENTS['inclination_deg']:.6f}°")
    print(f"  Longitude of ascending node: {ORBITAL_ELEMENTS['longitude_ascending_node_deg']:.6f}°")
    print(f"  Argument of perihelion: {ORBITAL_ELEMENTS['argument_of_perihelion_deg']:.6f}°")
    
    perihelion_time = mjd_to_datetime(ORBITAL_ELEMENTS['time_of_perihelion_mjd'])
    print(f"\nPerihelion time: {perihelion_time}")
    
    # Calculate position at perihelion
    perihelion_pos = calculate_orbital_position(
        ORBITAL_ELEMENTS['time_of_perihelion_mjd'],
        ORBITAL_ELEMENTS
    )
    
    print(f"\nPosition at perihelion:")
    print(f"  RA: {perihelion_pos['coordinates']['ra_deg']:.6f}°")
    print(f"  Dec: {perihelion_pos['coordinates']['dec_deg']:.6f}°")
    print(f"  Distance from Sun: {perihelion_pos['distance_from_sun_au']:.6f} AU")
    
    # Load decoded coordinates
    decoded_file = Path('3i_atlas_data/decoded_coordinates.json')
    if not decoded_file.exists():
        print("\nError: Decoded coordinates not found")
        return
    
    with open(decoded_file) as f:
        decoded = json.load(f)
    
    # Targets to analyze
    targets = [
        {
            'name': 'HD 286941',
            'ra_deg': 69.823071,
            'dec_deg': 11.265073,
            'type': 'Primary target (G5 emission star)',
        },
        {
            'name': 'Primary Decoded Coordinate',
            'ra_deg': decoded['primary_coordinate']['ra_deg'],
            'dec_deg': decoded['primary_coordinate']['dec_deg'],
            'type': 'Decoded coordinate',
        },
        {
            'name': 'HD 80838',
            'ra_deg': 139.650260,
            'dec_deg': -67.536354,
            'type': 'Structure 2 (K0IV star)',
        },
        {
            'name': 'HD 171547',
            'ra_deg': 279.421567,
            'dec_deg': -45.047173,
            'type': 'Structure 3 (A8/9IV star)',
        },
    ]
    
    # Time range: 6 months before to 6 months after perihelion
    t_peri = ORBITAL_ELEMENTS['time_of_perihelion_mjd']
    t_start = t_peri - 180  # 6 months before
    t_end = t_peri + 180    # 6 months after
    
    print(f"\n{'='*80}")
    print("FINDING CLOSEST APPROACH TO TARGETS")
    print("="*80)
    
    results = {
        'analysis_date': datetime.now().isoformat(),
        'orbital_elements': ORBITAL_ELEMENTS,
        'perihelion_position': perihelion_pos,
        'closest_approaches': [],
    }
    
    for target in targets:
        print(f"\n{target['name']} ({target['type']}):")
        print(f"  Target: RA {target['ra_deg']:.6f}°, Dec {target['dec_deg']:.6f}°")
        
        closest = find_closest_approach_to_target(
            ORBITAL_ELEMENTS,
            target['ra_deg'],
            target['dec_deg'],
            t_start,
            t_end,
            n_points=2000
        )
        
        print(f"  Closest approach:")
        print(f"    Time: {closest['closest_time_datetime']}")
        print(f"    Separation: {closest['closest_separation_deg']:.6f}° ({closest['closest_separation_arcmin']:.2f} arcmin)")
        
        if closest['closest_position']:
            pos = closest['closest_position']
            print(f"    3I/ATLAS position:")
            print(f"      RA: {pos['coordinates']['ra_deg']:.6f}°")
            print(f"      Dec: {pos['coordinates']['dec_deg']:.6f}°")
            print(f"      Distance from Sun: {pos['distance_from_sun_au']:.6f} AU")
        
        # Check if very close
        if closest['closest_separation_arcmin'] < 10:
            print(f"    ⚠ VERY CLOSE APPROACH!")
        elif closest['closest_separation_arcmin'] < 60:
            print(f"    ⚠ Close approach")
        
        results['closest_approaches'].append({
            'target': target,
            'closest_approach': closest,
        })
    
    # Save results
    output_file = Path('3i_atlas_data/trajectory_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n{'='*80}")
    print("TRAJECTORY ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_file}")
    
    # Summary
    print(f"\n{'='*80}")
    print("TRAJECTORY ANALYSIS SUMMARY")
    print("="*80)
    
    hd286941_approach = next(
        (ca for ca in results['closest_approaches'] 
         if ca['target']['name'] == 'HD 286941'),
        None
    )
    
    if hd286941_approach:
        sep_arcmin = hd286941_approach['closest_approach']['closest_separation_arcmin']
        time = hd286941_approach['closest_approach']['closest_time_datetime']
        
        print(f"\nHD 286941 Closest Approach:")
        print(f"  Separation: {sep_arcmin:.2f} arcmin")
        print(f"  Time: {time}")
        
        if sep_arcmin < 10:
            print(f"  ✓ VERY CLOSE - Highly significant!")
        elif sep_arcmin < 60:
            print(f"  ✓ Close approach - Significant")
        else:
            print(f"  ✗ Not particularly close")
    
    print(f"\n{'='*80}")
    
    return results

if __name__ == "__main__":
    analyze_3i_atlas_trajectory()

