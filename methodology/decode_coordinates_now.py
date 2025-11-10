#!/usr/bin/env python3
"""
Decode 3I/ATLAS Coordinates NOW
Convert 37-bit structures to celestial coordinates and cross-reference
"""

import json
import numpy as np
from pathlib import Path
from datetime import datetime
import math

def load_decoding_results():
    """Load message decoding results"""
    results_file = Path('3i_atlas_data/message_decoding_analysis.json')
    if not results_file.exists():
        return None
    
    with open(results_file) as f:
        return json.load(f)

def convert_ra_dec_to_standard(ra_deg, dec_deg):
    """Convert RA/Dec from degrees to standard format"""
    # RA: degrees to hours, minutes, seconds
    ra_hours = ra_deg / 15.0
    ra_h = int(ra_hours)
    ra_m = int((ra_hours - ra_h) * 60)
    ra_s = ((ra_hours - ra_h) * 60 - ra_m) * 60
    
    # Dec: degrees to degrees, arcminutes, arcseconds
    dec_sign = '+' if dec_deg >= 0 else '-'
    dec_abs = abs(dec_deg)
    dec_d = int(dec_abs)
    dec_m = int((dec_abs - dec_d) * 60)
    dec_s = ((dec_abs - dec_d) * 60 - dec_m) * 60
    
    return {
        'ra_hms': f"{ra_h:02d}h {ra_m:02d}m {ra_s:06.3f}s",
        'ra_deg': ra_deg,
        'ra_hours': ra_hours,
        'dec_dms': f"{dec_sign}{dec_d:02d}° {dec_m:02d}' {dec_s:06.3f}\"",
        'dec_deg': dec_deg,
    }

def get_distance_info(ra_deg, dec_deg):
    """Get distance information for coordinate (requires astrometric data)"""
    # Distance calculation requires proper astrometric data from catalogs
    # No distance assumption is made
    return None

def decode_all_coordinates():
    """Decode all 37-bit structures to celestial coordinates"""
    print("="*80)
    print("DECODING 3I/ATLAS COORDINATES NOW")
    print("="*80)
    
    results = load_decoding_results()
    if not results:
        print("Error: Decoding results not found")
        return
    
    # Get all 37-bit structures
    structures = results.get('37_bit_structures', {}).get('atlas_spectrum', [])
    
    # Get decoded messages
    decoded_messages = results.get('decoded_messages', {})
    
    all_coordinates = []
    
    print(f"\n{'='*80}")
    print("37-BIT STRUCTURE COORDINATES")
    print(f"{'='*80}\n")
    
    for i, struct in enumerate(structures, 1):
        ra_deg = struct['ra_degrees']
        dec_deg = struct['dec_degrees']
        
        # Convert to standard format
        std_format = convert_ra_dec_to_standard(ra_deg, dec_deg)
        
        # Distance info requires astrometric data from catalogs
        distance_info = get_distance_info(ra_deg, dec_deg)
        
        coord_info = {
            'structure_number': i,
            'bit_string': struct['bit_string'],
            'coordinate_value': struct['coordinate_value'],
            'ra_deg': ra_deg,
            'dec_deg': dec_deg,
            'ra_hms': std_format['ra_hms'],
            'dec_dms': std_format['dec_dms'],
        }
        if distance_info:
            coord_info['distance_info'] = distance_info
        
        all_coordinates.append(coord_info)
        
        print(f"Structure {i}:")
        print(f"  Bit String: {struct['bit_string'][:50]}...")
        print(f"  RA: {std_format['ra_hms']} ({ra_deg:.6f}°)")
        print(f"  Dec: {std_format['dec_dms']} ({dec_deg:.6f}°)")
        print(f"  Coordinate Value: {struct['coordinate_value']:,}")
        print()
    
    print(f"\n{'='*80}")
    print("ENCODING SCHEME COORDINATES")
    print(f"{'='*80}\n")
    
    for key, message in decoded_messages.items():
        ra_deg = message['ra_degrees']
        dec_deg = message['dec_degrees']
        
        # Convert to standard format
        std_format = convert_ra_dec_to_standard(ra_deg, dec_deg)
        
        # Distance info requires astrometric data from catalogs
        distance_info = get_distance_info(ra_deg, dec_deg)
        
        coord_info = {
            'encoding_scheme': key,
            'ra_deg': ra_deg,
            'dec_deg': dec_deg,
            'ra_hms': std_format['ra_hms'],
            'dec_dms': std_format['dec_dms'],
        }
        if distance_info:
            coord_info['distance_info'] = distance_info
        
        all_coordinates.append(coord_info)
        
        print(f"{key.replace('_', ' ').title()}:")
        print(f"  RA: {std_format['ra_hms']} ({ra_deg:.6f}°)")
        print(f"  Dec: {std_format['dec_dms']} ({dec_deg:.6f}°)")
        print()
    
    # Find convergence region
    print(f"\n{'='*80}")
    print("COORDINATE CONVERGENCE ANALYSIS")
    print(f"{'='*80}\n")
    
    # Group coordinates by proximity
    convergence_groups = []
    
    for coord in all_coordinates:
        ra = coord['ra_deg']
        dec = coord['dec_deg']
        
        # Find if it belongs to an existing group (within 10 degrees)
        found_group = False
        for group in convergence_groups:
            group_ra = group['center_ra']
            group_dec = group['center_dec']
            
            # Calculate angular separation
            ra_diff = abs(ra - group_ra)
            if ra_diff > 180:
                ra_diff = 360 - ra_diff
            
            dec_diff = abs(dec - group_dec)
            
            # Simple distance check (within 10 degrees)
            if ra_diff < 10 and dec_diff < 10:
                group['coordinates'].append(coord)
                # Update center (average)
                group['center_ra'] = (group['center_ra'] * len(group['coordinates']) + ra) / (len(group['coordinates']) + 1)
                group['center_dec'] = (group['center_dec'] * len(group['coordinates']) + dec) / (len(group['coordinates']) + 1)
                found_group = True
                break
        
        if not found_group:
            convergence_groups.append({
                'center_ra': ra,
                'center_dec': dec,
                'coordinates': [coord]
            })
    
    # Sort by number of coordinates
    convergence_groups.sort(key=lambda x: len(x['coordinates']), reverse=True)
    
    print(f"Found {len(convergence_groups)} convergence regions:\n")
    
    for i, group in enumerate(convergence_groups, 1):
        center_std = convert_ra_dec_to_standard(group['center_ra'], group['center_dec'])
        
        print(f"Convergence Region {i}:")
        print(f"  Center: RA {center_std['ra_hms']}, Dec {center_std['dec_dms']}")
        print(f"  Number of coordinates: {len(group['coordinates'])}")
        print(f"  Coordinates:")
        for coord in group['coordinates']:
            if 'structure_number' in coord:
                print(f"    - Structure {coord['structure_number']}: RA {coord['ra_hms']}, Dec {coord['dec_dms']}")
            else:
                print(f"    - {coord['encoding_scheme']}: RA {coord['ra_hms']}, Dec {coord['dec_dms']}")
        print()
    
    # Primary coordinate (most converged)
    if convergence_groups:
        primary = convergence_groups[0]
        primary_std = convert_ra_dec_to_standard(primary['center_ra'], primary['center_dec'])
        
        print(f"\n{'='*80}")
        print("PRIMARY DECODED COORDINATE")
        print(f"{'='*80}\n")
        print(f"Right Ascension: {primary_std['ra_hms']} ({primary['center_ra']:.6f}°)")
        print(f"Declination: {primary_std['dec_dms']} ({primary['center_dec']:.6f}°)")
        print(f"\nConvergence: {len(primary['coordinates'])} coordinates converge on this region")
        print(f"\nThis coordinate is the most consistent across multiple encoding schemes.")
        print(f"It may represent:")
        print(f"  • The source location of 3I/ATLAS")
        print(f"  • A target object to observe")
        print(f"  • A reference point in the coordinate system")
    
    # Save decoded coordinates
    output_file = Path('3i_atlas_data/decoded_coordinates.json')
    with open(output_file, 'w') as f:
        json.dump({
            'decoding_date': datetime.now().isoformat(),
            'all_coordinates': all_coordinates,
            'convergence_groups': [
                {
                    'center_ra': g['center_ra'],
                    'center_dec': g['center_dec'],
                    'center_ra_hms': convert_ra_dec_to_standard(g['center_ra'], g['center_dec'])['ra_hms'],
                    'center_dec_dms': convert_ra_dec_to_standard(g['center_ra'], g['center_dec'])['dec_dms'],
                    'num_coordinates': len(g['coordinates']),
                    'coordinates': g['coordinates']
                }
                for g in convergence_groups
            ],
            'primary_coordinate': {
                'ra_deg': primary['center_ra'],
                'dec_deg': primary['center_dec'],
                'ra_hms': primary_std['ra_hms'],
                'dec_dms': primary_std['dec_dms'],
            } if convergence_groups else None
        }, f, indent=2)
    
    print(f"\n{'='*80}")
    print("DECODING COMPLETE")
    print(f"{'='*80}")
    print(f"\nAll decoded coordinates saved to: {output_file}")
    print(f"\nNext steps:")
    print(f"  1. Cross-reference with astronomical catalogs (SIMBAD, Gaia, etc.)")
    print(f"  2. Determine distance from astrometric data (parallax measurements)")
    print(f"  3. Check for known objects at these coordinates")
    print(f"  4. Analyze trajectory from multiple coordinate points")
    print(f"\n{'='*80}")
    
    return all_coordinates, convergence_groups

if __name__ == "__main__":
    decode_all_coordinates()

