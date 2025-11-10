#!/usr/bin/env python3
"""
Query NED (NASA/IPAC Extragalactic Database) for LEDA 1363602
NED is specifically designed for extragalactic objects and should have redshift information
"""

import json
import requests
from pathlib import Path
from datetime import datetime

# LEDA 1363602 coordinates
LEDA_1363602 = {
    'identifier': 'LEDA 1363602',
    'ra_deg': 72.939583,  # 04h51m45.5s
    'dec_deg': 9.334444,  # +09°20'04"
    'ra_hms': '04h51m45.5s',
    'dec_dms': '+09°20\'04"',
}

# Constants
C_LIGHT = 299792.458  # km/s (speed of light)
H0 = 70.0  # km/s/Mpc (Hubble constant, approximate)
LY_TO_KM = 9.461e12  # km
MPC_TO_KM = 3.086e19  # km

def query_ned_by_name(name):
    """Query NED by object name"""
    url = "https://ned.ipac.caltech.edu/cgi-bin/objsearch"
    
    params = {
        'objname': name,
        'extend': 'no',
        'hconst': '73',
        'omegam': '0.27',
        'omegav': '0.73',
        'corr_z': '1',
        'out_csys': 'Equatorial',
        'out_equinox': 'J2000.0',
        'obj_sort': 'RA or Longitude',
        'of': 'pre_text',
        'zv_breaker': '30000.0',
        'list_limit': '5',
        'img_stamp': 'YES',
    }
    
    try:
        print(f"  Querying NED for {name}...")
        response = requests.get(url, params=params, timeout=15)
        
        if response.status_code == 200:
            result = response.text
            return {'success': True, 'data': result}
        else:
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def query_ned_by_coordinates(ra_deg, dec_deg, radius_arcmin=5.0):
    """Query NED by coordinates"""
    url = "https://ned.ipac.caltech.edu/cgi-bin/objsearch"
    
    params = {
        'search_type': 'Coordinates',
        'objname': f"{ra_deg} {dec_deg}",
        'extend': 'no',
        'hconst': '73',
        'omegam': '0.27',
        'omegav': '0.73',
        'corr_z': '1',
        'out_csys': 'Equatorial',
        'out_equinox': 'J2000.0',
        'obj_sort': 'Distance to search center',
        'of': 'pre_text',
        'zv_breaker': '30000.0',
        'list_limit': '5',
        'img_stamp': 'YES',
        'radius': f"{radius_arcmin}",
    }
    
    try:
        print(f"  Querying NED for coordinates RA {ra_deg:.6f}°, Dec {dec_deg:.6f}°...")
        response = requests.get(url, params=params, timeout=15)
        
        if response.status_code == 200:
            result = response.text
            return {'success': True, 'data': result}
        else:
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def parse_ned_redshift(data):
    """Parse redshift from NED data"""
    if not data:
        return None
    
    lines = data.split('\n')
    redshift = None
    redshift_error = None
    velocity = None
    distance = None
    
    for i, line in enumerate(lines):
        # Look for redshift information
        if 'Redshift' in line or 'z =' in line or 'z=' in line:
            # Try to extract redshift value
            parts = line.split()
            for j, part in enumerate(parts):
                try:
                    # Look for z value
                    if 'z' in part.lower() or 'redshift' in part.lower():
                        # Try next part as value
                        if j + 1 < len(parts):
                            try:
                                z = float(parts[j + 1])
                                if 0 < z < 10:  # Reasonable redshift range
                                    redshift = z
                                    # Try to get error if available
                                    if j + 2 < len(parts):
                                        try:
                                            redshift_error = float(parts[j + 2])
                                        except:
                                            pass
                                    break
                            except:
                                pass
                except:
                    continue
        
        # Look for velocity information
        if 'Velocity' in line or 'km/s' in line:
            parts = line.split()
            for j, part in enumerate(parts):
                try:
                    v = float(part)
                    if 0 < v < 300000:  # Reasonable velocity range (km/s)
                        velocity = v
                        break
                except:
                    continue
        
        # Look for distance information
        if 'Distance' in line or 'Mpc' in line or 'Mly' in line:
            parts = line.split()
            for j, part in enumerate(parts):
                try:
                    d = float(part)
                    if 0 < d < 10000:  # Reasonable distance range (Mpc)
                        distance = d
                        break
                except:
                    continue
    
    return {
        'redshift': redshift,
        'redshift_error': redshift_error,
        'velocity_km_s': velocity,
        'distance_mpc': distance,
    }

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

def query_leda_1363602_ned():
    """Query NED for LEDA 1363602"""
    print("="*80)
    print("QUERYING NED FOR LEDA 1363602")
    print("="*80)
    
    # Try querying by name first
    print("\n1. Querying NED by Object Name:")
    print("="*80)
    
    result_by_name = query_ned_by_name('LEDA 1363602')
    
    if result_by_name.get('success'):
        print(f"  ✓ NED query by name successful")
        print(f"  Data preview:")
        print(f"  {result_by_name['data'][:1000]}...")
        
        # Parse redshift
        redshift_data = parse_ned_redshift(result_by_name['data'])
        
        if redshift_data and redshift_data['redshift']:
            print(f"\n  Redshift Information:")
            print(f"    Redshift (z): {redshift_data['redshift']:.6f}")
            if redshift_data['redshift_error']:
                print(f"    Redshift error: {redshift_data['redshift_error']:.6f}")
            if redshift_data['velocity_km_s']:
                print(f"    Velocity: {redshift_data['velocity_km_s']:.2f} km/s")
            if redshift_data['distance_mpc']:
                print(f"    Distance: {redshift_data['distance_mpc']:.2f} Mpc")
            
            # Calculate distance from redshift
            distance_data = calculate_distance_from_redshift(redshift_data['redshift'])
            
            if distance_data:
                print(f"\n  Distance from Redshift:")
                print(f"    Distance: {distance_data['distance_ly']:.2e} light-years")
                print(f"    Distance: {distance_data['distance_mpc']:.2f} Mpc")
                print(f"    Distance: {distance_data['distance_pc']:.2e} parsecs")
                
                return {
                    'success': True,
                    'redshift_data': redshift_data,
                    'distance_data': distance_data,
                    'raw_data': result_by_name['data'],
                }
        else:
            print(f"  ✗ Redshift not found in NED data")
    else:
        print(f"  ✗ NED query by name failed: {result_by_name.get('error', 'Unknown')}")
    
    # Try querying by coordinates
    print("\n2. Querying NED by Coordinates:")
    print("="*80)
    
    result_by_coords = query_ned_by_coordinates(
        LEDA_1363602['ra_deg'],
        LEDA_1363602['dec_deg'],
        radius_arcmin=5.0
    )
    
    if result_by_coords.get('success'):
        print(f"  ✓ NED query by coordinates successful")
        print(f"  Data preview:")
        print(f"  {result_by_coords['data'][:1000]}...")
        
        # Parse redshift
        redshift_data = parse_ned_redshift(result_by_coords['data'])
        
        if redshift_data and redshift_data['redshift']:
            print(f"\n  Redshift Information:")
            print(f"    Redshift (z): {redshift_data['redshift']:.6f}")
            if redshift_data['redshift_error']:
                print(f"    Redshift error: {redshift_data['redshift_error']:.6f}")
            if redshift_data['velocity_km_s']:
                print(f"    Velocity: {redshift_data['velocity_km_s']:.2f} km/s")
            if redshift_data['distance_mpc']:
                print(f"    Distance: {redshift_data['distance_mpc']:.2f} Mpc")
            
            # Calculate distance from redshift
            distance_data = calculate_distance_from_redshift(redshift_data['redshift'])
            
            if distance_data:
                print(f"\n  Distance from Redshift:")
                print(f"    Distance: {distance_data['distance_ly']:.2e} light-years")
                print(f"    Distance: {distance_data['distance_mpc']:.2f} Mpc")
                print(f"    Distance: {distance_data['distance_pc']:.2e} parsecs")
                
                return {
                    'success': True,
                    'redshift_data': redshift_data,
                    'distance_data': distance_data,
                    'raw_data': result_by_coords['data'],
                }
        else:
            print(f"  ✗ Redshift not found in NED data")
    else:
        print(f"  ✗ NED query by coordinates failed: {result_by_coords.get('error', 'Unknown')}")
    
    return {
        'success': False,
        'error': 'Redshift not found in NED data',
    }

if __name__ == "__main__":
    result = query_leda_1363602_ned()
    
    if result.get('success'):
        print(f"\n{'='*80}")
        print("QUERY COMPLETE")
        print(f"{'='*80}")
        print(f"\nRedshift: {result['redshift_data']['redshift']:.6f}")
        print(f"Distance: {result['distance_data']['distance_ly']:.2e} light-years")
        
        # Save results
        output_file = Path('3i_atlas_data/ned_leda_1363602_query.json')
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=2)
        
        print(f"\nResults saved to: {output_file}")
    else:
        print(f"\n{'='*80}")
        print("QUERY FAILED")
        print(f"{'='*80}")
        print(f"\nError: {result.get('error', 'Unknown')}")

