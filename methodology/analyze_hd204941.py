#!/usr/bin/env python3
"""
Analyze HD 204941 from Stellar Catalog
Compare with HD 286941 and decoded coordinates
"""

import json
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt

# HD 204941 data from Stellar Catalog
HD_204941 = {
    'name': 'HD 204941',
    'distance_ly': 94.0,  # 93.500 ly
    'distance_pc': 28.7,  # 28.7 pc
    'parallax_mas': 34.883,
    'spectral_class': 'K2V C',
    'star_type': 'Orange star, High proper motion star',
    'ra_hms': '21h 32m 23.195s',
    'dec_dms': '-20° 57\' 28.708"',
    'ra_deg': 323.096646,  # Calculated: 21*15 + 32*15/60 + 23.195*15/3600
    'dec_deg': -20.957974,  # Calculated: -20 - 57/60 - 28.708/3600
    'mass_solar': 0.74,
    'radius_solar': 0.72,
    'temperature_k': 5056,
    'temperature_solar': 0.88,
    'age_gyr': 6.67,
    'age_solar': 1.45,
    'magnitude_v': 8.5,
    'magnitude_v_abs': 6.2,
    'exoplanets': [
        {
            'name': 'HD 204941 b',
            'semi_major_axis_au': 2.56,
            'mass_earth': 84.6,
            'period_days': 1733,
        }
    ],
    'catalogues': [
        'LTT 8565',
        'BD-21 6035',
        'HIC 106353',
        'HIP 106353',
        '2MASS J21322353-2057264',
        'NLTT 51476',
        'SAO 190421',
        'TIC 99837626',
        'TYC 6373-214-1',
        'WISEA J213223.30-205728.0',
        'Gaia DR3 6830027182179257472',
    ],
}

# HD 286941 data (from previous analysis)
HD_286941 = {
    'name': 'HD 286941',
    'ra_deg': 69.823071,
    'dec_deg': 11.265083,
    'distance_ly': 151.14,  # If in Hyades cluster
    'distance_pc': 46.34,  # If in Hyades cluster
    'constellation': 'Taurus',
    'spectral_class': 'G5 emission star',
}

# Decoded coordinate (from previous analysis)
DECODED_COORDINATE = {
    'ra_deg': 69.797974,
    'dec_deg': 11.250000,
    'source': '37-bit structures, encoding schemes',
    'statistical_significance': 'p < 0.001 (***)',
    'separation_from_hd286941_arcmin': 1.73,
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

def analyze_hd204941():
    """Analyze HD 204941 and compare with HD 286941 and decoded coordinates"""
    print("="*80)
    print("ANALYZING HD 204941 FROM STELLAR CATALOG")
    print("="*80)
    
    print(f"\nHD 204941 Information (from Stellar Catalog):")
    print(f"  Name: {HD_204941['name']}")
    print(f"  Distance: {HD_204941['distance_ly']:.1f} light-years ({HD_204941['distance_pc']:.1f} pc)")
    print(f"  Parallax: {HD_204941['parallax_mas']:.3f} mas")
    print(f"  Spectral class: {HD_204941['spectral_class']}")
    print(f"  Star type: {HD_204941['star_type']}")
    print(f"  RA: {HD_204941['ra_hms']} ({HD_204941['ra_deg']:.6f}°)")
    print(f"  Dec: {HD_204941['dec_dms']} ({HD_204941['dec_deg']:.6f}°)")
    print(f"  Mass: {HD_204941['mass_solar']:.0%} solar mass")
    print(f"  Radius: {HD_204941['radius_solar']:.0%} solar radius")
    print(f"  Temperature: {HD_204941['temperature_k']:.0f} K ({HD_204941['temperature_solar']:.0%} solar)")
    print(f"  Age: {HD_204941['age_gyr']:.2f} Gyr ({HD_204941['age_solar']:.2f}x solar)")
    print(f"  V Magnitude: {HD_204941['magnitude_v']:.1f} (apparent), {HD_204941['magnitude_v_abs']:.1f} (absolute)")
    
    print(f"\nExoplanets:")
    for planet in HD_204941['exoplanets']:
        print(f"  • {planet['name']}:")
        print(f"    Semi-major axis: {planet['semi_major_axis_au']:.2f} AU")
        print(f"    Mass: {planet['mass_earth']:.1f} M⊕")
        print(f"    Period: {planet['period_days']:.0f} days")
    
    # Compare with HD 286941
    print(f"\n{'='*80}")
    print("COMPARING HD 204941 WITH HD 286941")
    print("="*80)
    
    print(f"\nHD 286941 Information:")
    print(f"  Name: {HD_286941['name']}")
    print(f"  RA: {HD_286941['ra_deg']:.6f}°")
    print(f"  Dec: {HD_286941['dec_deg']:.6f}°")
    print(f"  Distance: {HD_286941['distance_ly']:.1f} light-years (if in Hyades cluster)")
    print(f"  Spectral class: {HD_286941['spectral_class']}")
    print(f"  Constellation: {HD_286941['constellation']}")
    
    # Angular separation
    sep_hd204941_hd286941 = calculate_angular_separation(
        HD_204941['ra_deg'],
        HD_204941['dec_deg'],
        HD_286941['ra_deg'],
        HD_286941['dec_deg']
    )
    
    print(f"\nAngular Separation:")
    print(f"  HD 204941 vs HD 286941: {sep_hd204941_hd286941:.6f}° ({sep_hd204941_hd286941*60:.2f} arcmin)")
    
    if sep_hd204941_hd286941 < 1.0:
        print(f"  ✓ VERY CLOSE - HD 204941 and HD 286941 are close on sky")
    elif sep_hd204941_hd286941 < 5.0:
        print(f"  ✓ Close - HD 204941 and HD 286941 are somewhat close")
    else:
        print(f"  ✗ Not close - HD 204941 and HD 286941 are far apart on sky")
    
    # Compare with decoded coordinate
    print(f"\n{'='*80}")
    print("COMPARING HD 204941 WITH DECODED COORDINATE")
    print("="*80)
    
    print(f"\nDecoded Coordinate:")
    print(f"  RA: {DECODED_COORDINATE['ra_deg']:.6f}°")
    print(f"  Dec: {DECODED_COORDINATE['dec_deg']:.6f}°")
    print(f"  Source: {DECODED_COORDINATE['source']}")
    print(f"  Statistical significance: {DECODED_COORDINATE['statistical_significance']}")
    print(f"  Separation from HD 286941: {DECODED_COORDINATE['separation_from_hd286941_arcmin']:.2f} arcmin")
    
    # Angular separation
    sep_hd204941_decoded = calculate_angular_separation(
        HD_204941['ra_deg'],
        HD_204941['dec_deg'],
        DECODED_COORDINATE['ra_deg'],
        DECODED_COORDINATE['dec_deg']
    )
    
    print(f"\nAngular Separation:")
    print(f"  HD 204941 vs Decoded Coordinate: {sep_hd204941_decoded:.6f}° ({sep_hd204941_decoded*60:.2f} arcmin)")
    
    if sep_hd204941_decoded < 1.0:
        print(f"  ✓ VERY CLOSE - HD 204941 matches decoded coordinate!")
    elif sep_hd204941_decoded < 5.0:
        print(f"  ✓ Close - HD 204941 is close to decoded coordinate")
    else:
        print(f"  ✗ Not close - HD 204941 does not match decoded coordinate")
    
    # Distance comparison
    print(f"\n{'='*80}")
    print("DISTANCE COMPARISON")
    print("="*80)
    
    print(f"\nDistance Comparison:")
    print(f"  HD 204941: {HD_204941['distance_ly']:.1f} light-years")
    print(f"  HD 286941 (if in Hyades): {HD_286941['distance_ly']:.1f} light-years")
    
    print(f"\nDistance Difference:")
    diff = abs(HD_204941['distance_ly'] - HD_286941['distance_ly'])
    print(f"  Difference: {diff:.1f} light-years")
    
    if HD_204941['distance_ly'] < HD_286941['distance_ly']:
        print(f"  ✓ HD 204941 is closer than HD 286941")
    else:
        print(f"  ✗ HD 204941 is farther than HD 286941")
    
    # Key findings
    print(f"\n{'='*80}")
    print("KEY FINDINGS")
    print("="*80)
    
    print(f"\n1. HD 204941 Location:")
    print(f"   • RA: {HD_204941['ra_deg']:.6f}°")
    print(f"   • Dec: {HD_204941['dec_deg']:.6f}°")
    print(f"   • Distance: {HD_204941['distance_ly']:.1f} light-years")
    print(f"   • Separation from decoded coordinate: {sep_hd204941_decoded*60:.2f} arcmin")
    print(f"   • Status: {'✓ Matches decoded coordinate' if sep_hd204941_decoded < 5.0 else '✗ Does not match decoded coordinate'}")
    
    print(f"\n2. HD 286941 Location:")
    print(f"   • RA: {HD_286941['ra_deg']:.6f}°")
    print(f"   • Dec: {HD_286941['dec_deg']:.6f}°")
    print(f"   • Distance: {HD_286941['distance_ly']:.1f} light-years (if in Hyades)")
    print(f"   • Separation from decoded coordinate: {DECODED_COORDINATE['separation_from_hd286941_arcmin']:.2f} arcmin")
    print(f"   • Status: ✓ Matches decoded coordinate (highly significant, p < 0.001)")
    
    print(f"\n3. Distance Comparison:")
    print(f"   • HD 204941: {HD_204941['distance_ly']:.1f} light-years")
    print(f"   • HD 286941: {HD_286941['distance_ly']:.1f} light-years")
    print(f"   • HD 204941 is {'closer' if HD_204941['distance_ly'] < HD_286941['distance_ly'] else 'farther'} than HD 286941")
    
    print(f"\n4. Angular Separation:")
    print(f"   • HD 204941 vs HD 286941: {sep_hd204941_hd286941*60:.2f} arcmin")
    print(f"   • HD 204941 vs Decoded Coordinate: {sep_hd204941_decoded*60:.2f} arcmin")
    print(f"   • HD 286941 vs Decoded Coordinate: {DECODED_COORDINATE['separation_from_hd286941_arcmin']:.2f} arcmin")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'hd204941': HD_204941,
        'hd286941': HD_286941,
        'decoded_coordinate': DECODED_COORDINATE,
        'angular_separations': {
            'hd204941_vs_hd286941_deg': sep_hd204941_hd286941,
            'hd204941_vs_hd286941_arcmin': sep_hd204941_hd286941 * 60,
            'hd204941_vs_decoded_deg': sep_hd204941_decoded,
            'hd204941_vs_decoded_arcmin': sep_hd204941_decoded * 60,
            'hd286941_vs_decoded_arcmin': DECODED_COORDINATE['separation_from_hd286941_arcmin'],
        },
        'distance_comparison': {
            'hd204941_ly': HD_204941['distance_ly'],
            'hd286941_ly': HD_286941['distance_ly'],
            'distance_difference_ly': abs(HD_204941['distance_ly'] - HD_286941['distance_ly']),
        },
        'conclusions': {
            'hd204941_matches_decoded': sep_hd204941_decoded < 5.0,
            'hd286941_matches_decoded': True,  # Known from previous analysis
            'hd204941_closer_than_hd286941': HD_204941['distance_ly'] < HD_286941['distance_ly'],
        }
    }
    
    output_file = Path('3i_atlas_data/hd204941_analysis.json')
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
    
    print(f"\nHD 204941 Analysis:")
    print(f"  • Location: RA {HD_204941['ra_deg']:.6f}°, Dec {HD_204941['dec_deg']:.6f}°")
    print(f"  • Distance: {HD_204941['distance_ly']:.1f} light-years")
    print(f"  • Separation from decoded coordinate: {sep_hd204941_decoded*60:.2f} arcmin")
    print(f"  • Status: {'✓ Matches decoded coordinate' if sep_hd204941_decoded < 5.0 else '✗ Does not match decoded coordinate'}")
    
    print(f"\nComparison with HD 286941:")
    print(f"  • HD 286941 matches decoded coordinate: ✓ (1.73 arcmin, p < 0.001)")
    print(f"  • HD 204941 matches decoded coordinate: {'✓' if sep_hd204941_decoded < 5.0 else '✗'} ({sep_hd204941_decoded*60:.2f} arcmin)")
    print(f"  • HD 204941 vs HD 286941: {sep_hd204941_hd286941*60:.2f} arcmin")
    
    print(f"\nDistance Comparison:")
    print(f"  • HD 204941: {HD_204941['distance_ly']:.1f} light-years")
    print(f"  • HD 286941: {HD_286941['distance_ly']:.1f} light-years")
    print(f"  • HD 204941 is {'closer' if HD_204941['distance_ly'] < HD_286941['distance_ly'] else 'farther'} than HD 286941")
    
    print(f"\n{'='*80}")
    
    return results

if __name__ == "__main__":
    analyze_hd204941()

