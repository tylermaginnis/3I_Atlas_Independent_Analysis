#!/usr/bin/env python3
"""
Analyze NVSS Radio Source at Diophantine v3 Coordinates
Investigate NVSS J045150+092332 and its relationship to LEDA 1363602
"""

import json
import requests
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt

# NVSS Radio Source at Diophantine v3 coordinates
NVSS_J045150_092332 = {
    'identifier': 'NVSS J045150+092332',
    'type': 'Radio source (Rad)',
    'ra_deg': 72.961750,  # 04h51m50.82s = 4*15 + 51*15/60 + 50.82*15/3600
    'dec_deg': 9.392222,  # +09°23'32.0" = 9 + 23/60 + 32/3600
    'ra_hms': '04h51m50.82s',
    'dec_dms': '+09°23\'32.0"',
    'distance_arcsec': 277.71,
    'distance_arcmin': 4.63,
}

# LEDA 1363602 (galaxy at Diophantine v3 coordinates)
LEDA_1363602 = {
    'identifier': 'LEDA 1363602',
    'type': 'Galaxy (G)',
    'ra_deg': 72.939583,  # 04h51m45.5s
    'dec_deg': 9.334444,  # +09°20'04"
    'ra_hms': '04h51m45.5s',
    'dec_dms': '+09°20\'04"',
    'distance_ly': 720e6,  # 720 million light-years
    'distance_mpc': 220.8,
    'redshift': 0.053915,
}

# Diophantine v3 coordinates
DIOPHANTINE_V3_COORDINATES = {
    'ra_deg': 72.888794,
    'dec_deg': 9.364471,
    'ra_hms': '04h51m33.311s',
    'dec_dms': '+09°21\'52.096"',
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

def query_nvss_source(name):
    """Query NVSS for radio source information"""
    # NVSS catalog browser URL
    url = "http://www.cv.nrao.edu/nvss/cgi-bin/NVSSlist.cgi"
    
    params = {
        'RA': name.split('J')[1][:6] if 'J' in name else '',
        'DEC': name.split('+')[1] if '+' in name else name.split('-')[1] if '-' in name else '',
        'Radius': '0.1',  # 0.1 degrees = 6 arcmin
        'RadiusUnits': 'deg',
    }
    
    try:
        print(f"  Querying NVSS for {name}...")
        response = requests.get(url, params=params, timeout=15)
        
        if response.status_code == 200:
            result = response.text
            return {'success': True, 'data': result[:2000]}  # First 2000 chars
        else:
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def query_ned_for_nvss(name):
    """Query NED for NVSS source information"""
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
            return {'success': True, 'data': result[:2000]}  # First 2000 chars
        else:
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def analyze_nvss_radio_source():
    """Analyze NVSS radio source at Diophantine v3 coordinates"""
    print("="*80)
    print("ANALYZING NVSS RADIO SOURCE AT DIOPHANTINE V3 COORDINATES")
    print("="*80)
    
    # NVSS Radio Source
    print("\n1. NVSS Radio Source Information:")
    print("="*80)
    
    print(f"\nNVSS J045150+092332:")
    print(f"  Type: {NVSS_J045150_092332['type']}")
    print(f"  RA: {NVSS_J045150_092332['ra_deg']:.6f}° ({NVSS_J045150_092332['ra_hms']})")
    print(f"  Dec: {NVSS_J045150_092332['dec_deg']:.6f}° ({NVSS_J045150_092332['dec_dms']})")
    print(f"  Distance from Diophantine v3 coordinates: {NVSS_J045150_092332['distance_arcmin']:.2f} arcmin")
    
    # LEDA 1363602
    print(f"\n2. LEDA 1363602 (Galaxy at Diophantine v3 Coordinates):")
    print("="*80)
    
    print(f"\nLEDA 1363602:")
    print(f"  Type: {LEDA_1363602['type']}")
    print(f"  RA: {LEDA_1363602['ra_deg']:.6f}° ({LEDA_1363602['ra_hms']})")
    print(f"  Dec: {LEDA_1363602['dec_deg']:.6f}° ({LEDA_1363602['dec_dms']})")
    print(f"  Distance: {LEDA_1363602['distance_ly']/1e6:.0f} million light-years ({LEDA_1363602['distance_mpc']:.1f} Mpc)")
    print(f"  Redshift: z = {LEDA_1363602['redshift']:.6f}")
    
    # Calculate separations
    print(f"\n3. Angular Separations:")
    print("="*80)
    
    sep_nvss_leda = calculate_angular_separation(
        NVSS_J045150_092332['ra_deg'], NVSS_J045150_092332['dec_deg'],
        LEDA_1363602['ra_deg'], LEDA_1363602['dec_deg']
    )
    
    sep_nvss_diophantine = calculate_angular_separation(
        NVSS_J045150_092332['ra_deg'], NVSS_J045150_092332['dec_deg'],
        DIOPHANTINE_V3_COORDINATES['ra_deg'], DIOPHANTINE_V3_COORDINATES['dec_deg']
    )
    
    print(f"\n  NVSS J045150+092332 vs LEDA 1363602:")
    print(f"    Angular separation: {sep_nvss_leda:.6f}° ({sep_nvss_leda*60:.2f} arcmin)")
    print(f"    RA difference: {abs(NVSS_J045150_092332['ra_deg'] - LEDA_1363602['ra_deg'])*60:.2f} arcmin")
    print(f"    Dec difference: {abs(NVSS_J045150_092332['dec_deg'] - LEDA_1363602['dec_deg'])*60:.2f} arcmin")
    
    print(f"\n  NVSS J045150+092332 vs Diophantine v3 coordinates:")
    print(f"    Angular separation: {sep_nvss_diophantine:.6f}° ({sep_nvss_diophantine*60:.2f} arcmin)")
    
    # Check if NVSS source is associated with LEDA 1363602
    if sep_nvss_leda * 60 < 5.0:  # Within 5 arcmin
        print(f"\n    ✓ NVSS source is VERY CLOSE to LEDA 1363602!")
        print(f"      Could be radio emission from the galaxy")
    elif sep_nvss_leda * 60 < 10.0:  # Within 10 arcmin
        print(f"\n    ✓ NVSS source is close to LEDA 1363602")
        print(f"      Could be related to the galaxy")
    else:
        print(f"\n    ⚠ NVSS source is not close to LEDA 1363602")
        print(f"      May be unrelated")
    
    # Query NED for NVSS source
    print(f"\n4. Querying NED for NVSS Source:")
    print("="*80)
    
    ned_result = query_ned_for_nvss('NVSS J045150+092332')
    
    if ned_result.get('success'):
        print(f"  ✓ NED query successful")
        print(f"  Data preview:")
        print(f"  {ned_result['data'][:500]}...")
    else:
        print(f"  ✗ NED query failed: {ned_result.get('error', 'Unknown')}")
    
    # Analysis
    print(f"\n5. ANALYSIS:")
    print("="*80)
    
    print(f"\n  A. NVSS Radio Source Properties:")
    print(f"     • Type: Radio source (Rad)")
    print(f"     • Position: RA {NVSS_J045150_092332['ra_deg']:.6f}°, Dec {NVSS_J045150_092332['dec_deg']:.6f}°")
    print(f"     • Distance from Diophantine v3: {NVSS_J045150_092332['distance_arcmin']:.2f} arcmin")
    print(f"     • Distance from LEDA 1363602: {sep_nvss_leda*60:.2f} arcmin")
    
    print(f"\n  B. Relationship to LEDA 1363602:")
    if sep_nvss_leda * 60 < 5.0:
        print(f"     • NVSS source is VERY CLOSE to LEDA 1363602 ({sep_nvss_leda*60:.2f} arcmin)")
        print(f"     • Could be radio emission from the galaxy")
        print(f"     • LEDA 1363602 is 720 million light-years away")
        print(f"     • If associated, NVSS source is also extragalactic")
    else:
        print(f"     • NVSS source is not close to LEDA 1363602 ({sep_nvss_leda*60:.2f} arcmin)")
        print(f"     • May be unrelated to the galaxy")
    
    print(f"\n  C. Relationship to Diophantine v3 Coordinates:")
    print(f"     • NVSS source is {sep_nvss_diophantine*60:.2f} arcmin from Diophantine v3 coordinates")
    print(f"     • One of 4 objects found at Diophantine v3 coordinates")
    print(f"     • Could be part of the encoded message")
    
    print(f"\n  D. Extragalactic Possibility:")
    print(f"     • If NVSS source is associated with LEDA 1363602, it is extragalactic")
    print(f"     • Distance: 720 million light-years")
    print(f"     • Could be radio emission from active galactic nucleus (AGN)")
    print(f"     • Or could be unrelated foreground/background source")
    
    # Summary
    print(f"\n6. SUMMARY")
    print("="*80)
    
    print(f"\nNVSS J045150+092332:")
    print(f"  Type: Radio source (Rad)")
    print(f"  Position: RA {NVSS_J045150_092332['ra_deg']:.6f}°, Dec {NVSS_J045150_092332['dec_deg']:.6f}°")
    print(f"  Distance from Diophantine v3: {NVSS_J045150_092332['distance_arcmin']:.2f} arcmin")
    print(f"  Distance from LEDA 1363602: {sep_nvss_leda*60:.2f} arcmin")
    
    if sep_nvss_leda * 60 < 5.0:
        print(f"\n  ✓ NVSS source is VERY CLOSE to LEDA 1363602!")
        print(f"    Could be radio emission from the galaxy (extragalactic)")
    else:
        print(f"\n  ⚠ NVSS source is not close to LEDA 1363602")
        print(f"    May be unrelated to the galaxy")
    
    print(f"\nRelationship to Diophantine v3 Coordinates:")
    print(f"  • One of 4 objects found at Diophantine v3 coordinates")
    print(f"  • Could be part of the encoded message")
    print(f"  • If extragalactic, could be pointing to extragalactic origin/destination")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'nvss_source': NVSS_J045150_092332,
        'leda_1363602': LEDA_1363602,
        'diophantine_v3_coordinates': DIOPHANTINE_V3_COORDINATES,
        'separations': {
            'nvss_vs_leda_deg': sep_nvss_leda,
            'nvss_vs_leda_arcmin': sep_nvss_leda * 60,
            'nvss_vs_diophantine_deg': sep_nvss_diophantine,
            'nvss_vs_diophantine_arcmin': sep_nvss_diophantine * 60,
        },
        'analysis': {
            'nvss_close_to_leda': sep_nvss_leda * 60 < 5.0,
            'nvss_could_be_galaxy_emission': sep_nvss_leda * 60 < 5.0,
            'nvss_extragalactic_possibility': sep_nvss_leda * 60 < 5.0,
        },
        'ned_query': ned_result,
    }
    
    output_file = Path('3i_atlas_data/nvss_radio_source_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    return results

if __name__ == "__main__":
    analyze_nvss_radio_source()

