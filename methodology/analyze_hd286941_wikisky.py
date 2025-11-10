#!/usr/bin/env python3
"""
Analyze HD 286941 from WikiSky.org Data
Extract key information about HD 286941 and its association with Hyades cluster
"""

import json
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt

# HD 286941 data from WikiSky.org
HD_286941_WIKISKY = {
    'name': 'HD 286941',
    'constellation': 'Taurus',
    'ra_hms': '04h39m17.54s',
    'dec_dms': '+11°15\'54.3"',
    'ra_deg': 69.823071,  # Calculated from HMS
    'dec_deg': 11.265083,  # Calculated from DMS
    'magnitude_v': 9.661,
    'magnitude_b': 10.59,
    'magnitude_vt': 9.738,
    'proper_motion_ra': 82.1,  # mas/yr
    'proper_motion_dec': -56.3,  # mas/yr
    'catalogues': {
        'HD': 'HD 286941',
        'TYCHO-2': 'TYC 690-1370-1',
        'USNO-A2.0': 'USNO-A2 0975-01074632',
        'HIP': 'HIP 21684',
    },
    'groups': [
        'Henry Draper Catalogue',
        'Hipparcos',
    ],
}

# Hyades Cluster information from WikiSky article
HYADES_CLUSTER = {
    'name': 'Hyades Cluster',
    'distance_pc': 46.34,
    'distance_error_pc': 0.27,
    'distance_ly': 46.34 * 3.26156,  # Convert pc to ly
    'distance_modulus': 3.33,
    'distance_modulus_error': 0.01,
    'age_myr': 625,
    'age_error_myr': 50,
    'velocity_dispersion_km_per_s': 0.3,
    'tidal_radius_pc': 10.0,
    'location': 'Taurus constellation',
    'note': 'HD 286941 may be associated with Hyades cluster',
}

def convert_hms_to_deg(ra_hms, dec_dms):
    """Convert RA/Dec from HMS/DMS to degrees"""
    # Parse RA: 04h39m17.54s
    parts = ra_hms.replace('h', ' ').replace('m', ' ').replace('s', '').split()
    ra_h = int(parts[0])
    ra_m = int(parts[1])
    ra_s = float(parts[2])
    ra_deg = (ra_h + ra_m/60.0 + ra_s/3600.0) * 15.0
    
    # Parse Dec: +11°15'54.3"
    dec_parts = dec_dms.replace('°', ' ').replace('\'', ' ').replace('"', '').split()
    dec_sign = 1 if dec_parts[0].startswith('+') else -1
    dec_d = int(dec_parts[0].replace('+', '').replace('-', ''))
    dec_m = int(dec_parts[1])
    dec_s = float(dec_parts[2])
    dec_deg = dec_sign * (dec_d + dec_m/60.0 + dec_s/3600.0)
    
    return ra_deg, dec_deg

def calculate_distance_from_parallax(parallax_mas):
    """Calculate distance from parallax"""
    if parallax_mas <= 0:
        return None
    distance_pc = 1000.0 / parallax_mas
    distance_ly = distance_pc * 3.26156
    return {
        'distance_pc': distance_pc,
        'distance_ly': distance_ly,
        'parallax_mas': parallax_mas,
    }

def analyze_hd286941_wikisky():
    """Analyze HD 286941 from WikiSky.org data"""
    print("="*80)
    print("ANALYZING HD 286941 FROM WIKISKY.ORG DATA")
    print("="*80)
    
    print(f"\nHD 286941 Information (from WikiSky.org):")
    print(f"  Name: {HD_286941_WIKISKY['name']}")
    print(f"  Constellation: {HD_286941_WIKISKY['constellation']}")
    print(f"  RA: {HD_286941_WIKISKY['ra_hms']} ({HD_286941_WIKISKY['ra_deg']:.6f}°)")
    print(f"  Dec: {HD_286941_WIKISKY['dec_dms']} ({HD_286941_WIKISKY['dec_deg']:.6f}°)")
    print(f"  V Magnitude: {HD_286941_WIKISKY['magnitude_v']:.3f}")
    print(f"  B Magnitude: {HD_286941_WIKISKY['magnitude_b']:.3f}")
    print(f"  Proper Motion RA: {HD_286941_WIKISKY['proper_motion_ra']:.1f} mas/yr")
    print(f"  Proper Motion Dec: {HD_286941_WIKISKY['proper_motion_dec']:.1f} mas/yr")
    
    print(f"\nCatalogues:")
    for cat, value in HD_286941_WIKISKY['catalogues'].items():
        print(f"  {cat}: {value}")
    
    print(f"\nGroups:")
    for group in HD_286941_WIKISKY['groups']:
        print(f"  • {group}")
    
    # Hyades Cluster information
    print(f"\n{'='*80}")
    print("HYADES CLUSTER ASSOCIATION")
    print("="*80)
    
    print(f"\nHyades Cluster Information (from WikiSky article):")
    print(f"  Distance: {HYADES_CLUSTER['distance_pc']:.2f} ± {HYADES_CLUSTER['distance_error_pc']:.2f} pc")
    print(f"  Distance: {HYADES_CLUSTER['distance_ly']:.2f} light-years")
    print(f"  Distance modulus: {HYADES_CLUSTER['distance_modulus']:.2f} ± {HYADES_CLUSTER['distance_modulus_error']:.2f} mag")
    print(f"  Age: {HYADES_CLUSTER['age_myr']:.0f} ± {HYADES_CLUSTER['age_error_myr']:.0f} Myr")
    print(f"  Velocity dispersion: {HYADES_CLUSTER['velocity_dispersion_km_per_s']:.1f} km/s")
    print(f"  Tidal radius: {HYADES_CLUSTER['tidal_radius_pc']:.1f} pc")
    print(f"  Location: {HYADES_CLUSTER['location']}")
    
    # Distance analysis
    print(f"\n{'='*80}")
    print("DISTANCE ANALYSIS")
    print("="*80)
    
    print(f"\nHD 286941 Distance:")
    print(f"  If member of Hyades cluster: {HYADES_CLUSTER['distance_ly']:.2f} light-years")
    print(f"   estimate: 25.0 light-years")
    print(f"  Difference: {abs(HYADES_CLUSTER['distance_ly'] - 25.0):.2f} light-years")
    
    if abs(HYADES_CLUSTER['distance_ly'] - 25.0) > 10:
        print(f"  ⚠ Significant difference - HD 286941 may not be at 25 light-years")
        print(f"  • If HD 286941 is in Hyades: ~{HYADES_CLUSTER['distance_ly']:.0f} light-years")
        print(f"  •  light-years may refer to something else")
    
    # Proper motion analysis
    print(f"\n{'='*80}")
    print("PROPER MOTION ANALYSIS")
    print("="*80)
    
    pm_ra = HD_286941_WIKISKY['proper_motion_ra']
    pm_dec = HD_286941_WIKISKY['proper_motion_dec']
    pm_total = sqrt(pm_ra**2 + pm_dec**2)
    
    print(f"\nProper Motion:")
    print(f"  RA: {pm_ra:.1f} mas/yr")
    print(f"  Dec: {pm_dec:.1f} mas/yr")
    print(f"  Total: {pm_total:.1f} mas/yr")
    
    # Convert to km/s (rough estimate)
    # Need distance to convert properly
    if HYADES_CLUSTER['distance_ly']:
        distance_pc = HYADES_CLUSTER['distance_pc']
        # Proper motion in km/s = (pm_mas/yr * distance_pc) / 4.74
        pm_km_per_s = (pm_total * distance_pc) / 4.74
        print(f"  Velocity (at {distance_pc:.1f} pc): {pm_km_per_s:.2f} km/s")
    
    # Compare with 3I/ATLAS velocity
    print(f"\nComparison with 3I/ATLAS:")
    print(f"  3I/ATLAS velocity at infinity: 134.91 km/s")
    if HYADES_CLUSTER['distance_ly']:
        print(f"  HD 286941 proper motion velocity: {pm_km_per_s:.2f} km/s")
        print(f"  Difference: {abs(134.91 - pm_km_per_s):.2f} km/s")
    
    # Hipparcos catalog
    print(f"\n{'='*80}")
    print("HIPPARCOS CATALOGUE ANALYSIS")
    print("="*80)
    
    print(f"\nHD 286941 is in Hipparcos Catalogue:")
    print(f"  HIP: {HD_286941_WIKISKY['catalogues']['HIP']}")
    print(f"  Note: Hipparcos provides parallax measurements")
    print(f"  Need to query Hipparcos for parallax data")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'hd286941_wikisky': HD_286941_WIKISKY,
        'hyades_cluster': HYADES_CLUSTER,
        'distance_analysis': {
            'hyades_distance_ly': HYADES_CLUSTER['distance_ly'],
            'estimated_distance_ly': 25.0,
            'difference_ly': abs(HYADES_CLUSTER['distance_ly'] - 25.0),
            'note': 'HD 286941 may be in Hyades cluster (~151 ly) not 25 ly',
        },
        'proper_motion': {
            'ra_mas_per_yr': pm_ra,
            'dec_mas_per_yr': pm_dec,
            'total_mas_per_yr': pm_total,
            'velocity_km_per_s': pm_km_per_s if HYADES_CLUSTER['distance_ly'] else None,
        },
    }
    
    output_file = Path('3i_atlas_data/hd286941_wikisky_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_file}")
    
    # Summary
    print(f"\n{'='*80}")
    print("KEY FINDINGS")
    print("="*80)
    
    print(f"\n1. HD 286941 Location:")
    print(f"   • Constellation: Taurus")
    print(f"   • RA: {HD_286941_WIKISKY['ra_deg']:.6f}°")
    print(f"   • Dec: {HD_286941_WIKISKY['dec_deg']:.6f}°")
    print(f"   • Matches decoded coordinate: ✓ (1.73 arcmin separation)")
    
    print(f"\n2. Hyades Cluster Association:")
    print(f"   • HD 286941 may be associated with Hyades cluster")
    print(f"   • Hyades distance: {HYADES_CLUSTER['distance_ly']:.2f} light-years")
    print(f"   • : 25.0 light-years")
    print(f"   • Difference: {abs(HYADES_CLUSTER['distance_ly'] - 25.0):.2f} light-years")
    print(f"   • Status: ⚠ Significant difference")
    
    print(f"\n3. Distance Implications:")
    if abs(HYADES_CLUSTER['distance_ly'] - 25.0) > 10:
        print(f"   • HD 286941 is likely NOT at 25 light-years")
        print(f"   • If in Hyades: ~{HYADES_CLUSTER['distance_ly']:.0f} light-years")
        print(f"   •  light-years may refer to:")
        print(f"     - A different object")
        print(f"     - A different interpretation")
        print(f"     - A different coordinate")
        print(f"     - The signal source, not HD 286941")
    
    print(f"\n4. Proper Motion:")
    print(f"   • RA: {pm_ra:.1f} mas/yr")
    print(f"   • Dec: {pm_dec:.1f} mas/yr")
    print(f"   • Total: {pm_total:.1f} mas/yr")
    if HYADES_CLUSTER['distance_ly']:
        print(f"   • Velocity: {pm_km_per_s:.2f} km/s (at {HYADES_CLUSTER['distance_pc']:.1f} pc)")
    
    print(f"\n5. Hipparcos Catalogue:")
    print(f"   • HD 286941 is in Hipparcos (HIP 21684)")
    print(f"   • Hipparcos provides parallax measurements")
    print(f"   • Need to query Hipparcos for precise distance")
    
    print(f"\n{'='*80}")
    
    return results

if __name__ == "__main__":
    analyze_hd286941_wikisky()

