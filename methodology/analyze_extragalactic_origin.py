#!/usr/bin/env python3
"""
Analyze Extragalactic Origin Possibility for 3I/ATLAS
Check if 3I/ATLAS could have originated from LEDA 1363602 (galaxy)
"""

import json
import requests
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt

# Diophantine v3 coordinates (pointing to LEDA 1363602 galaxy)
DIOPHANTINE_V3_COORDINATES = {
    'ra_deg': 72.888794,
    'dec_deg': 9.364471,
    'ra_hms': '04h51m33.311s',
    'dec_dms': '+09°21\'52.096"',
    'name': 'Diophantine v3 decoded coordinate',
    'encoding': 'diophantine_v3',
    'entropy': 0.996244,
}

# LEDA 1363602 (galaxy at Diophantine v3 coordinates)
LEDA_1363602 = {
    'identifier': 'LEDA 1363602',
    'type': 'Galaxy (G)',
    'ra_deg': 72.939583,  # 04h51m45.5s = 4*15 + 51*15/60 + 45.5*15/3600
    'dec_deg': 9.334444,  # +09°20'04" = 9 + 20/60 + 4/3600
    'ra_hms': '04h51m45.5s',
    'dec_dms': '+09°20\'04"',
    'distance_arcsec': 210.32,
    'distance_arcmin': 3.51,
}

# Decoded coordinate (phi_log3 encoding - points to HIP 21684)
DECODED_COORDINATE = {
    'ra_deg': 69.797974,
    'dec_deg': 11.250000,
    'ra_hms': '04h39m11.514s',
    'dec_dms': '+11°15\'00.000"',
    'name': 'Decoded 3I/ATLAS coordinate (phi_log3)',
}

# HIP 21684 (HD 286941) - star in our galaxy
HIP_21684 = {
    'ra_deg': 69.822875,
    'dec_deg': 11.265222,
    'ra_hms': '04h39m17.49s',
    'dec_dms': '+11°15\'54.8"',
    'name': 'HIP 21684 (HD 286941)',
    'type': 'Star (G0)',
    'distance_parsecs': 150.83,  # From Universe Guide
    'distance_ly': 491.95,
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

def query_leda_1363602():
    """Query SIMBAD for LEDA 1363602 details"""
    url = "https://simbad.cds.unistra.fr/simbad/sim-id"
    
    params = {
        'Ident': 'LEDA 1363602',
        'output.format': 'ASCII',
    }
    
    try:
        print(f"  Querying SIMBAD for LEDA 1363602...")
        response = requests.get(url, params=params, timeout=15)
        
        if response.status_code == 200:
            result = response.text
            return {'success': True, 'data': result[:2000]}  # First 2000 chars
        else:
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def analyze_extragalactic_origin():
    """Analyze if 3I/ATLAS could have originated from another galaxy"""
    print("="*80)
    print("ANALYZING EXTRAGALACTIC ORIGIN POSSIBILITY FOR 3I/ATLAS")
    print("="*80)
    
    # Diophantine v3 coordinates
    print("\n1. Diophantine v3 Coordinates:")
    print("="*80)
    
    print(f"\nDiophantine v3 Decoded Coordinate:")
    print(f"  RA: {DIOPHANTINE_V3_COORDINATES['ra_deg']:.6f}° ({DIOPHANTINE_V3_COORDINATES['ra_hms']})")
    print(f"  Dec: {DIOPHANTINE_V3_COORDINATES['dec_deg']:.6f}° ({DIOPHANTINE_V3_COORDINATES['dec_dms']})")
    print(f"  Encoding: {DIOPHANTINE_V3_COORDINATES['encoding']}")
    print(f"  Entropy: {DIOPHANTINE_V3_COORDINATES['entropy']:.6f} (almost maximum)")
    
    # LEDA 1363602
    print(f"\n2. LEDA 1363602 (Galaxy at Diophantine v3 Coordinates):")
    print("="*80)
    
    print(f"\nLEDA 1363602:")
    print(f"  Type: {LEDA_1363602['type']}")
    print(f"  RA: {LEDA_1363602['ra_deg']:.6f}° ({LEDA_1363602['ra_hms']})")
    print(f"  Dec: {LEDA_1363602['dec_deg']:.6f}° ({LEDA_1363602['dec_dms']})")
    print(f"  Distance from Diophantine v3 coordinates: {LEDA_1363602['distance_arcmin']:.2f} arcmin")
    
    # Query SIMBAD for LEDA 1363602 details
    print(f"\n3. Querying SIMBAD for LEDA 1363602 Details:")
    print("="*80)
    
    simbad_result = query_leda_1363602()
    
    if simbad_result.get('success'):
        print(f"  ✓ SIMBAD query successful")
        print(f"  Data preview:")
        print(f"  {simbad_result['data'][:500]}...")
    else:
        print(f"  ✗ SIMBAD query failed: {simbad_result.get('error', 'Unknown')}")
    
    # Compare with other coordinates
    print(f"\n4. Comparison with Other Coordinates:")
    print("="*80)
    
    print(f"\nDecoded Coordinate (phi_log3 - points to HIP 21684):")
    print(f"  RA: {DECODED_COORDINATE['ra_deg']:.6f}° ({DECODED_COORDINATE['ra_hms']})")
    print(f"  Dec: {DECODED_COORDINATE['dec_deg']:.6f}° ({DECODED_COORDINATE['dec_dms']})")
    print(f"  Matches HIP 21684: 1.73 arcmin (p < 0.001)")
    
    print(f"\nHIP 21684 (HD 286941 - star in our galaxy):")
    print(f"  Type: {HIP_21684['type']}")
    print(f"  RA: {HIP_21684['ra_deg']:.6f}° ({HIP_21684['ra_hms']})")
    print(f"  Dec: {HIP_21684['dec_deg']:.6f}° ({HIP_21684['dec_dms']})")
    print(f"  Distance: {HIP_21684['distance_ly']:.2f} light-years ({HIP_21684['distance_parsecs']:.2f} parsecs)")
    
    # Calculate separations
    sep_diophantine_leda = calculate_angular_separation(
        DIOPHANTINE_V3_COORDINATES['ra_deg'], DIOPHANTINE_V3_COORDINATES['dec_deg'],
        LEDA_1363602['ra_deg'], LEDA_1363602['dec_deg']
    )
    
    sep_diophantine_decoded = calculate_angular_separation(
        DIOPHANTINE_V3_COORDINATES['ra_deg'], DIOPHANTINE_V3_COORDINATES['dec_deg'],
        DECODED_COORDINATE['ra_deg'], DECODED_COORDINATE['dec_deg']
    )
    
    sep_diophantine_hip = calculate_angular_separation(
        DIOPHANTINE_V3_COORDINATES['ra_deg'], DIOPHANTINE_V3_COORDINATES['dec_deg'],
        HIP_21684['ra_deg'], HIP_21684['dec_deg']
    )
    
    sep_leda_hip = calculate_angular_separation(
        LEDA_1363602['ra_deg'], LEDA_1363602['dec_deg'],
        HIP_21684['ra_deg'], HIP_21684['dec_deg']
    )
    
    print(f"\nAngular Separations:")
    print(f"  Diophantine v3 vs LEDA 1363602: {sep_diophantine_leda:.6f}° ({sep_diophantine_leda*60:.2f} arcmin)")
    print(f"  Diophantine v3 vs Decoded (phi_log3): {sep_diophantine_decoded:.6f}° ({sep_diophantine_decoded*60:.2f} arcmin)")
    print(f"  Diophantine v3 vs HIP 21684: {sep_diophantine_hip:.6f}° ({sep_diophantine_hip*60:.2f} arcmin)")
    print(f"  LEDA 1363602 vs HIP 21684: {sep_leda_hip:.6f}° ({sep_leda_hip*60:.2f} arcmin)")
    
    # Analyze extragalactic origin possibility
    print(f"\n5. EXTRAGALACTIC ORIGIN ANALYSIS:")
    print("="*80)
    
    print(f"\nKey Questions:")
    print(f"  1. Could 3I/ATLAS have originated from LEDA 1363602 (another galaxy)?")
    print(f"  2. What would that mean for its trajectory?")
    print(f"  3. How would we detect if it's extragalactic?")
    print(f"  4. What about the other objects at those coordinates?")
    
    print(f"\nAnalysis:")
    print(f"\n  A. Distance Considerations:")
    print(f"     • LEDA 1363602 is a galaxy (extragalactic)")
    print(f"     • Typical galaxy distances: millions to billions of light-years")
    print(f"     • HIP 21684 is {HIP_21684['distance_ly']:.2f} light-years away (in our galaxy)")
    print(f"     • 3I/ATLAS trajectory suggests solar system origin (eccentricity: {ORBITAL_ELEMENTS['eccentricity']:.3f})")
    
    print(f"\n  B. Trajectory Considerations:")
    print(f"     • 3I/ATLAS has hyperbolic orbit (eccentricity > 1)")
    print(f"     • Incoming asymptote direction: needs calculation")
    print(f"     • If extragalactic, would need extremely high velocity")
    print(f"     • Typical galactic escape velocity: ~500-600 km/s")
    print(f"     • Intergalactic travel would require much higher velocities")
    
    print(f"\n  C. Encoding Considerations:")
    print(f"     • Phi_Log3 encoding points to HIP 21684 (star in our galaxy)")
    print(f"     • Diophantine v3 encoding points to LEDA 1363602 (galaxy)")
    print(f"     • Both encodings may be valid but encode different information")
    print(f"     • Could LEDA 1363602 be a reference point or secondary target?")
    
    print(f"\n  D. Physical Considerations:")
    print(f"     • 3I/ATLAS shows encoded signal patterns")
    print(f"     • If extragalactic, would need advanced propulsion")
    print(f"     • Travel time would be millions to billions of years")
    print(f"     • Or could be pointing to LEDA 1363602 as a reference/destination")
    
    # Summary
    print(f"\n6. SUMMARY")
    print("="*80)
    
    print(f"\nDiophantine v3 Coordinates:")
    print(f"  RA: {DIOPHANTINE_V3_COORDINATES['ra_deg']:.6f}° ({DIOPHANTINE_V3_COORDINATES['ra_hms']})")
    print(f"  Dec: {DIOPHANTINE_V3_COORDINATES['dec_deg']:.6f}° ({DIOPHANTINE_V3_COORDINATES['dec_dms']})")
    print(f"  Encoding: {DIOPHANTINE_V3_COORDINATES['encoding']} (entropy: {DIOPHANTINE_V3_COORDINATES['entropy']:.6f})")
    
    print(f"\nObjects Found:")
    print(f"  • LEDA 1363602 (galaxy) - {LEDA_1363602['distance_arcmin']:.2f} arcmin")
    print(f"  • UCAC4 497-008460 (emission star) - 4.60 arcmin")
    print(f"  • NVSS J045150+092332 (radio source) - 4.63 arcmin")
    print(f"  • Gaia DR3 3292871264573792640 (star) - 4.94 arcmin")
    
    print(f"\nExtragalactic Origin Possibility:")
    print(f"  ✓ LEDA 1363602 is a galaxy (extragalactic)")
    print(f"  ✓ Diophantine v3 encoding points to it")
    print(f"  ⚠ But 3I/ATLAS trajectory suggests solar system origin")
    print(f"  ⚠ Would need extremely high velocity for intergalactic travel")
    print(f"  ⚠ Could be pointing to LEDA 1363602 as reference/destination")
    
    print(f"\nAlternative Interpretations:")
    print(f"  1. LEDA 1363602 could be a reference point (not origin)")
    print(f"  2. Could be pointing to LEDA 1363602 as a destination")
    print(f"  3. Could be encoding multiple targets (HIP 21684 + LEDA 1363602)")
    print(f"  4. Could be a coordinate system reference (galactic vs extragalactic)")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'diophantine_v3_coordinates': DIOPHANTINE_V3_COORDINATES,
        'leda_1363602': LEDA_1363602,
        'decoded_coordinate': DECODED_COORDINATE,
        'hip_21684': HIP_21684,
        'orbital_elements': ORBITAL_ELEMENTS,
        'separations': {
            'diophantine_v3_vs_leda_deg': sep_diophantine_leda,
            'diophantine_v3_vs_leda_arcmin': sep_diophantine_leda * 60,
            'diophantine_v3_vs_decoded_deg': sep_diophantine_decoded,
            'diophantine_v3_vs_decoded_arcmin': sep_diophantine_decoded * 60,
            'diophantine_v3_vs_hip_deg': sep_diophantine_hip,
            'diophantine_v3_vs_hip_arcmin': sep_diophantine_hip * 60,
            'leda_vs_hip_deg': sep_leda_hip,
            'leda_vs_hip_arcmin': sep_leda_hip * 60,
        },
        'extragalactic_analysis': {
            'leda_1363602_type': 'Galaxy (extragalactic)',
            'hip_21684_type': 'Star (in our galaxy)',
            'hip_21684_distance_ly': HIP_21684['distance_ly'],
            'trajectory_suggests': 'Solar system origin (hyperbolic orbit)',
            'possibility': 'Could be extragalactic but trajectory suggests otherwise',
            'alternative_interpretations': [
                'LEDA 1363602 could be a reference point (not origin)',
                'Could be pointing to LEDA 1363602 as a destination',
                'Could be encoding multiple targets (HIP 21684 + LEDA 1363602)',
                'Could be a coordinate system reference (galactic vs extragalactic)',
            ],
        },
        'simbad_query': simbad_result,
    }
    
    output_file = Path('3i_atlas_data/extragalactic_origin_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    return results

if __name__ == "__main__":
    analyze_extragalactic_origin()

