#!/usr/bin/env python3
"""
Comprehensive Extragalactic Analysis for 3I/ATLAS
1. Query LEDA 1363602 distance (redshift, distance measurement)
2. Calculate trajectory from LEDA 1363602 direction (if extragalactic)
3. Analyze velocity requirements for intergalactic travel
4. Investigate other objects at Diophantine v3 coordinates
"""

import json
import requests
import math
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt, asinh

# Constants
C_LIGHT = 299792.458  # km/s (speed of light)
H0 = 70.0  # km/s/Mpc (Hubble constant, approximate)
AU_TO_KM = 149597870.7  # km
LY_TO_KM = 9.461e12  # km
PC_TO_KM = 3.086e13  # km
MPC_TO_KM = 3.086e19  # km

# Diophantine v3 coordinates (pointing to LEDA 1363602 galaxy)
DIOPHANTINE_V3_COORDINATES = {
    'ra_deg': 72.888794,
    'dec_deg': 9.364471,
    'ra_hms': '04h51m33.311s',
    'dec_dms': '+09°21\'52.096"',
}

# LEDA 1363602 (galaxy at Diophantine v3 coordinates)
LEDA_1363602 = {
    'identifier': 'LEDA 1363602',
    'type': 'Galaxy (G)',
    'ra_deg': 72.939583,  # 04h51m45.5s
    'dec_deg': 9.334444,  # +09°20'04"
    'ra_hms': '04h51m45.5s',
    'dec_dms': '+09°20\'04"',
    'distance_arcsec': 210.32,
    'distance_arcmin': 3.51,
}

# 3I/ATLAS orbital elements
ORBITAL_ELEMENTS = {
    'time_of_perihelion_mjd': 60977.10275,
    'perihelion_distance_au': 1.3797753,
    'eccentricity': 6.320503,
    'inclination_deg': 175.12093,
    'longitude_ascending_node_deg': 322.34889,
    'argument_perihelion_deg': 127.77350,
}

# HIP 21684 (HD 286941) - star in our galaxy
HIP_21684 = {
    'ra_deg': 69.822875,
    'dec_deg': 11.265222,
    'distance_ly': 491.95,
    'distance_pc': 150.83,
}

def query_leda_1363602_detailed():
    """Query SIMBAD for LEDA 1363602 detailed information including redshift"""
    url = "https://simbad.cds.unistra.fr/simbad/sim-id"
    
    params = {
        'Ident': 'LEDA 1363602',
        'output.format': 'ASCII',
    }
    
    try:
        print(f"  Querying SIMBAD for LEDA 1363602 detailed information...")
        response = requests.get(url, params=params, timeout=15)
        
        if response.status_code == 200:
            result = response.text
            return {'success': True, 'data': result}
        else:
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def parse_redshift(data):
    """Parse redshift from SIMBAD data"""
    if not data:
        return None
    
    lines = data.split('\n')
    redshift = None
    redshift_error = None
    
    for line in lines:
        if 'z' in line.lower() and 'redshift' in line.lower():
            # Try to extract redshift value
            parts = line.split()
            for i, part in enumerate(parts):
                try:
                    z = float(part)
                    if 0 < z < 10:  # Reasonable redshift range
                        redshift = z
                        # Try to get error if available
                        if i + 1 < len(parts):
                            try:
                                redshift_error = float(parts[i + 1])
                            except:
                                pass
                        break
                except:
                    continue
    
    return {'redshift': redshift, 'redshift_error': redshift_error}

def calculate_distance_from_redshift(z, h0=H0):
    """Calculate distance from redshift using Hubble's law"""
    if z is None or z <= 0:
        return None
    
    # Hubble's law: v = H0 * d
    # For small redshifts: z ≈ v/c
    # For larger redshifts, need relativistic correction
    
    if z < 0.1:
        # Non-relativistic approximation
        v = z * C_LIGHT  # km/s
        d_mpc = v / h0  # Mpc
        d_km = d_mpc * MPC_TO_KM  # km
        d_ly = d_km / LY_TO_KM  # light-years
    else:
        # Relativistic correction
        # z = sqrt((1 + v/c) / (1 - v/c)) - 1
        # For large z, use cosmological distance
        # Simplified: d ≈ (c/H0) * z * (1 + z/2) for z < 1
        d_mpc = (C_LIGHT / h0) * z * (1 + z / 2)  # Mpc (approximate)
        d_km = d_mpc * MPC_TO_KM  # km
        d_ly = d_km / LY_TO_KM  # light-years
    
    return {
        'redshift': z,
        'distance_mpc': d_mpc,
        'distance_km': d_km,
        'distance_ly': d_ly,
        'distance_pc': d_mpc * 1e6,  # parsecs
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

def calculate_trajectory_from_leda():
    """Calculate trajectory from LEDA 1363602 direction"""
    # Calculate direction vector from LEDA 1363602 to solar system
    # (assuming 3I/ATLAS came from LEDA 1363602)
    
    leda_ra = radians(LEDA_1363602['ra_deg'])
    leda_dec = radians(LEDA_1363602['dec_deg'])
    
    # Convert to unit vector
    leda_x = cos(leda_dec) * cos(leda_ra)
    leda_y = cos(leda_dec) * sin(leda_ra)
    leda_z = sin(leda_dec)
    
    # Calculate incoming asymptote direction from orbital elements
    # For hyperbolic orbit, incoming asymptote is in direction of -v_infinity
    
    # Simplified: use orbital elements to estimate incoming direction
    # This is a rough approximation - full calculation would need N-body simulation
    
    return {
        'leda_direction': {
            'ra_deg': LEDA_1363602['ra_deg'],
            'dec_deg': LEDA_1363602['dec_deg'],
            'unit_vector': [leda_x, leda_y, leda_z],
        },
        'note': 'Full trajectory calculation requires N-body simulation',
    }

def analyze_velocity_requirements(distance_ly, travel_time_years=None):
    """Analyze velocity requirements for intergalactic travel"""
    if distance_ly is None:
        return None
    
    distance_km = distance_ly * LY_TO_KM
    
    # Calculate required velocities for different travel times
    travel_times = [
        1e6,  # 1 million years
        10e6,  # 10 million years
        100e6,  # 100 million years
        1e9,  # 1 billion years
        10e9,  # 10 billion years
    ]
    
    velocities = []
    
    for time_years in travel_times:
        time_seconds = time_years * 365.25 * 24 * 3600
        velocity_km_s = distance_km / time_seconds
        velocity_c = velocity_km_s / C_LIGHT
        
        velocities.append({
            'travel_time_years': time_years,
            'velocity_km_s': velocity_km_s,
            'velocity_c': velocity_c,
            'velocity_fraction_of_c': velocity_c,
        })
    
    # Calculate minimum velocity (if travel time is specified)
    if travel_time_years:
        time_seconds = travel_time_years * 365.25 * 24 * 3600
        min_velocity_km_s = distance_km / time_seconds
        min_velocity_c = min_velocity_km_s / C_LIGHT
    else:
        min_velocity_km_s = None
        min_velocity_c = None
    
    # Galactic escape velocity
    galactic_escape_velocity_km_s = 550  # km/s (approximate)
    galactic_escape_velocity_c = galactic_escape_velocity_km_s / C_LIGHT
    
    return {
        'distance_ly': distance_ly,
        'distance_km': distance_km,
        'galactic_escape_velocity_km_s': galactic_escape_velocity_km_s,
        'galactic_escape_velocity_c': galactic_escape_velocity_c,
        'velocities_for_travel_times': velocities,
        'min_velocity_km_s': min_velocity_km_s,
        'min_velocity_c': min_velocity_c,
    }

def investigate_other_objects():
    """Investigate other objects at Diophantine v3 coordinates"""
    # Objects found at Diophantine v3 coordinates
    objects = [
        {
            'identifier': 'LEDA 1363602',
            'type': 'Galaxy (G)',
            'distance_arcsec': 210.32,
            'distance_arcmin': 3.51,
        },
        {
            'identifier': 'UCAC4 497-008460',
            'type': 'Emission star (Em*)',
            'distance_arcsec': 276.27,
            'distance_arcmin': 4.60,
            'v_magnitude': 13.739,
            'b_magnitude': 15.133,
            'spectral_type': 'K7',
        },
        {
            'identifier': 'NVSS J045150+092332',
            'type': 'Radio source (Rad)',
            'distance_arcsec': 277.71,
            'distance_arcmin': 4.63,
        },
        {
            'identifier': 'Gaia DR3 3292871264573792640',
            'type': 'Star (*)',
            'distance_arcsec': 296.58,
            'distance_arcmin': 4.94,
        },
    ]
    
    return objects

def comprehensive_analysis():
    """Perform comprehensive extragalactic analysis"""
    print("="*80)
    print("COMPREHENSIVE EXTRAGALACTIC ANALYSIS FOR 3I/ATLAS")
    print("="*80)
    
    # 1. Query LEDA 1363602 distance
    print("\n1. QUERYING LEDA 1363602 DISTANCE (REDSHIFT, DISTANCE MEASUREMENT)")
    print("="*80)
    
    simbad_result = query_leda_1363602_detailed()
    
    if simbad_result.get('success'):
        print(f"  ✓ SIMBAD query successful")
        
        # Parse redshift
        redshift_data = parse_redshift(simbad_result['data'])
        
        if redshift_data and redshift_data['redshift']:
            print(f"\n  Redshift Information:")
            print(f"    Redshift (z): {redshift_data['redshift']:.6f}")
            if redshift_data['redshift_error']:
                print(f"    Redshift error: {redshift_data['redshift_error']:.6f}")
            
            # Calculate distance from redshift
            distance_data = calculate_distance_from_redshift(redshift_data['redshift'])
            
            if distance_data:
                print(f"\n  Distance from Redshift:")
                print(f"    Distance: {distance_data['distance_ly']:.2e} light-years")
                print(f"    Distance: {distance_data['distance_mpc']:.2f} Mpc")
                print(f"    Distance: {distance_data['distance_pc']:.2e} parsecs")
                
                leda_distance = distance_data['distance_ly']
            else:
                print(f"    ✗ Could not calculate distance from redshift")
                leda_distance = None
        else:
            print(f"  ✗ Redshift not found in SIMBAD data")
            print(f"  Data preview:")
            print(f"  {simbad_result['data'][:1000]}...")
            leda_distance = None
    else:
        print(f"  ✗ SIMBAD query failed: {simbad_result.get('error', 'Unknown')}")
        leda_distance = None
    
    # 2. Calculate trajectory from LEDA 1363602 direction
    print(f"\n2. CALCULATING TRAJECTORY FROM LEDA 1363602 DIRECTION")
    print("="*80)
    
    trajectory_data = calculate_trajectory_from_leda()
    
    print(f"\n  LEDA 1363602 Direction:")
    print(f"    RA: {trajectory_data['leda_direction']['ra_deg']:.6f}° ({LEDA_1363602['ra_hms']})")
    print(f"    Dec: {trajectory_data['leda_direction']['dec_deg']:.6f}° ({LEDA_1363602['dec_dms']})")
    print(f"    Unit vector: {trajectory_data['leda_direction']['unit_vector']}")
    
    print(f"\n  Note: {trajectory_data['note']}")
    print(f"    Full trajectory calculation requires N-body simulation")
    print(f"    Would need to account for galactic rotation, dark matter, etc.")
    
    # 3. Analyze velocity requirements
    print(f"\n3. ANALYZING VELOCITY REQUIREMENTS FOR INTERGALACTIC TRAVEL")
    print("="*80)
    
    if leda_distance:
        velocity_analysis = analyze_velocity_requirements(leda_distance)
        
        print(f"\n  Distance: {velocity_analysis['distance_ly']:.2e} light-years")
        print(f"  Distance: {velocity_analysis['distance_km']:.2e} km")
        
        print(f"\n  Galactic Escape Velocity:")
        print(f"    Velocity: {velocity_analysis['galactic_escape_velocity_km_s']:.2f} km/s")
        print(f"    Velocity: {velocity_analysis['galactic_escape_velocity_c']:.6f} c")
        
        print(f"\n  Required Velocities for Different Travel Times:")
        for v in velocity_analysis['velocities_for_travel_times']:
            print(f"\n    Travel time: {v['travel_time_years']:.2e} years")
            print(f"      Required velocity: {v['velocity_km_s']:.2e} km/s")
            print(f"      Required velocity: {v['velocity_c']:.6f} c")
            print(f"      Required velocity: {v['velocity_fraction_of_c']*100:.4f}% of speed of light")
            
            if v['velocity_c'] > 1.0:
                print(f"      ⚠ EXCEEDS SPEED OF LIGHT (impossible)")
            elif v['velocity_c'] > 0.1:
                print(f"      ⚠ Requires relativistic speeds (very difficult)")
            elif v['velocity_c'] > 0.01:
                print(f"      ⚠ Requires very high speeds (difficult)")
            else:
                print(f"      ✓ Achievable with advanced propulsion")
    else:
        print(f"  ✗ Cannot analyze velocity requirements without distance")
        print(f"    LEDA 1363602 distance not available")
        velocity_analysis = None
    
    # 4. Investigate other objects
    print(f"\n4. INVESTIGATING OTHER OBJECTS AT DIOPHANTINE V3 COORDINATES")
    print("="*80)
    
    other_objects = investigate_other_objects()
    
    print(f"\n  Objects Found at Diophantine v3 Coordinates:")
    print(f"    Total: {len(other_objects)} objects within 5 arcmin")
    
    for i, obj in enumerate(other_objects, 1):
        print(f"\n  Object {i}: {obj['identifier']}")
        print(f"    Type: {obj['type']}")
        print(f"    Distance: {obj['distance_arcmin']:.2f} arcmin ({obj['distance_arcsec']:.2f} arcsec)")
        
        if 'v_magnitude' in obj:
            print(f"    V Magnitude: {obj['v_magnitude']}")
        if 'b_magnitude' in obj:
            print(f"    B Magnitude: {obj['b_magnitude']}")
        if 'spectral_type' in obj:
            print(f"    Spectral Type: {obj['spectral_type']}")
        
        if obj['type'] == 'Galaxy (G)':
            print(f"    ⚠ EXTRAGALACTIC OBJECT")
        elif obj['type'] == 'Radio source (Rad)':
            print(f"    ⚠ RADIO SOURCE (could be extragalactic)")
        elif obj['type'] == 'Emission star (Em*)':
            print(f"    ✓ Emission star (in our galaxy)")
        elif obj['type'] == 'Star (*)':
            print(f"    ✓ Star (in our galaxy)")
    
    # Summary
    print(f"\n5. SUMMARY")
    print("="*80)
    
    print(f"\nLEDA 1363602 Analysis:")
    if leda_distance:
        print(f"  Distance: {leda_distance:.2e} light-years")
        print(f"  Status: ✓ Extragalactic (confirmed by redshift)")
    else:
        print(f"  Distance: Unknown (redshift not found)")
        print(f"  Status: ⚠ Extragalactic (assumed, but not confirmed)")
    
    print(f"\nTrajectory Analysis:")
    print(f"  LEDA 1363602 Direction: RA {LEDA_1363602['ra_deg']:.6f}°, Dec {LEDA_1363602['dec_deg']:.6f}°")
    print(f"  Note: Full trajectory calculation requires N-body simulation")
    
    if velocity_analysis:
        print(f"\nVelocity Requirements:")
        print(f"  Distance: {velocity_analysis['distance_ly']:.2e} light-years")
        print(f"  Galactic escape velocity: {velocity_analysis['galactic_escape_velocity_km_s']:.2f} km/s")
        
        # Find minimum reasonable travel time
        min_reasonable_velocity = None
        for v in velocity_analysis['velocities_for_travel_times']:
            if v['velocity_c'] < 1.0 and v['velocity_c'] > 0.01:
                min_reasonable_velocity = v
                break
        
        if min_reasonable_velocity:
            print(f"  Minimum reasonable travel time: {min_reasonable_velocity['travel_time_years']:.2e} years")
            print(f"    Required velocity: {min_reasonable_velocity['velocity_km_s']:.2e} km/s ({min_reasonable_velocity['velocity_c']:.6f} c)")
    
    print(f"\nOther Objects:")
    print(f"  Total: {len(other_objects)} objects within 5 arcmin")
    print(f"  Extragalactic: {sum(1 for obj in other_objects if 'Galaxy' in obj['type'] or 'Radio' in obj['type'])}")
    print(f"  Galactic: {sum(1 for obj in other_objects if 'Star' in obj['type'] or 'Emission' in obj['type'])}")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'leda_1363602': {
            'identifier': LEDA_1363602['identifier'],
            'type': LEDA_1363602['type'],
            'coordinates': {
                'ra_deg': LEDA_1363602['ra_deg'],
                'dec_deg': LEDA_1363602['dec_deg'],
            },
            'redshift': redshift_data if redshift_data else None,
            'distance': distance_data if leda_distance else None,
        },
        'trajectory_analysis': trajectory_data,
        'velocity_analysis': velocity_analysis,
        'other_objects': other_objects,
        'simbad_query': simbad_result,
    }
    
    output_file = Path('3i_atlas_data/comprehensive_extragalactic_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    return results

if __name__ == "__main__":
    comprehensive_analysis()

