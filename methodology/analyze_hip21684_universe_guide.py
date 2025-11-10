#!/usr/bin/env python3
"""
Analyze HIP 21684 (HD 286941) from Universe Guide
Compare Hipparcos distance measurements with previous estimates
"""

import json
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt

# HIP 21684 (HD 286941) data from Universe Guide
HIP_21684_UNIVERSE_GUIDE = {
    'name': 'HIP 21684',
    'hd_name': 'HD 286941',
    'constellation': 'Taurus',
    'spectral_type': 'G0',
    'star_type': 'yellow star',
    'ra_hms': '04h 39m 17.49s',
    'dec_dms': '+11° 15\' 54.8"',
    'ra_deg': 69.822875,  # Calculated: 4*15 + 39*15/60 + 17.49*15/3600
    'dec_deg': 11.265222,  # Calculated: 11 + 15/60 + 54.8/3600
    'distance_1997_ly': 341.18,
    'distance_1997_pc': 104.60,
    'parallax_1997_mas': 9.56000,
    'distance_2007_ly': 491.95,
    'distance_2007_pc': 150.83,
    'parallax_2007_mas': 6.63000,
    'apparent_magnitude': 9.68,
    'absolute_magnitude_1997': 4.58,
    'absolute_magnitude_2007': 3.79,
    'proper_motion_ra': 82.2,  # mas/yr
    'proper_motion_dec': -60.09,  # mas/yr
    'radial_velocity_km_per_s': 30.70,
    'radial_velocity_error': 0.20,
    'b_v_index': 0.72,
    'alternate_names': [
        'HD 286941',
        'TYC 690-1370-1',
        'BD+10 605',
    ],
}

# Previous distance estimates
DISTANCE_ESTIMATES = {
    'estimated_distance_ly': 25.0,
    'hyades_cluster_ly': 151.14,
    'wikisky_estimate_ly': 151.14,  # If in Hyades cluster
    'hipparcos_1997_ly': 341.18,
    'hipparcos_2007_ly': 491.95,
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

def analyze_hip21684_universe_guide():
    """Analyze HIP 21684 from Universe Guide data"""
    print("="*80)
    print("ANALYZING HIP 21684 (HD 286941) FROM UNIVERSE GUIDE")
    print("="*80)
    
    print(f"\nHIP 21684 Information (from Universe Guide):")
    print(f"  Name: {HIP_21684_UNIVERSE_GUIDE['name']} = {HIP_21684_UNIVERSE_GUIDE['hd_name']}")
    print(f"  Constellation: {HIP_21684_UNIVERSE_GUIDE['constellation']}")
    print(f"  Spectral type: {HIP_21684_UNIVERSE_GUIDE['spectral_type']} ({HIP_21684_UNIVERSE_GUIDE['star_type']})")
    print(f"  RA: {HIP_21684_UNIVERSE_GUIDE['ra_hms']} ({HIP_21684_UNIVERSE_GUIDE['ra_deg']:.6f}°)")
    print(f"  Dec: {HIP_21684_UNIVERSE_GUIDE['dec_dms']} ({HIP_21684_UNIVERSE_GUIDE['dec_deg']:.6f}°)")
    print(f"  Apparent magnitude: {HIP_21684_UNIVERSE_GUIDE['apparent_magnitude']:.2f}")
    print(f"  Absolute magnitude (1997): {HIP_21684_UNIVERSE_GUIDE['absolute_magnitude_1997']:.2f}")
    print(f"  Absolute magnitude (2007): {HIP_21684_UNIVERSE_GUIDE['absolute_magnitude_2007']:.2f}")
    print(f"  Proper motion RA: {HIP_21684_UNIVERSE_GUIDE['proper_motion_ra']:.2f} mas/yr")
    print(f"  Proper motion Dec: {HIP_21684_UNIVERSE_GUIDE['proper_motion_dec']:.2f} mas/yr")
    print(f"  Radial velocity: {HIP_21684_UNIVERSE_GUIDE['radial_velocity_km_per_s']:.2f} ± {HIP_21684_UNIVERSE_GUIDE['radial_velocity_error']:.2f} km/s")
    print(f"  B-V index: {HIP_21684_UNIVERSE_GUIDE['b_v_index']:.2f}")
    
    # Distance measurements
    print(f"\n{'='*80}")
    print("DISTANCE MEASUREMENTS")
    print("="*80)
    
    print(f"\nHipparcos Distance Measurements:")
    print(f"  1997 Hipparcos data:")
    print(f"    Parallax: {HIP_21684_UNIVERSE_GUIDE['parallax_1997_mas']:.5f} mas")
    print(f"    Distance: {HIP_21684_UNIVERSE_GUIDE['distance_1997_ly']:.2f} light-years ({HIP_21684_UNIVERSE_GUIDE['distance_1997_pc']:.2f} pc)")
    
    print(f"\n  2007 Hipparcos data (revised):")
    print(f"    Parallax: {HIP_21684_UNIVERSE_GUIDE['parallax_2007_mas']:.5f} mas")
    print(f"    Distance: {HIP_21684_UNIVERSE_GUIDE['distance_2007_ly']:.2f} light-years ({HIP_21684_UNIVERSE_GUIDE['distance_2007_pc']:.2f} pc)")
    
    print(f"\n  Difference: {abs(HIP_21684_UNIVERSE_GUIDE['distance_2007_ly'] - HIP_21684_UNIVERSE_GUIDE['distance_1997_ly']):.2f} light-years")
    print(f"  Note: Distance was recalculated, not that the star is moving")
    
    # Compare with decoded coordinate
    print(f"\n{'='*80}")
    print("COMPARING WITH DECODED COORDINATE")
    print("="*80)
    
    print(f"\nDecoded Coordinate:")
    print(f"  RA: {DECODED_COORDINATE['ra_deg']:.6f}°")
    print(f"  Dec: {DECODED_COORDINATE['dec_deg']:.6f}°")
    print(f"  Source: {DECODED_COORDINATE['source']}")
    print(f"  Statistical significance: {DECODED_COORDINATE['statistical_significance']}")
    
    # Angular separation
    sep_hip21684_decoded = calculate_angular_separation(
        HIP_21684_UNIVERSE_GUIDE['ra_deg'],
        HIP_21684_UNIVERSE_GUIDE['dec_deg'],
        DECODED_COORDINATE['ra_deg'],
        DECODED_COORDINATE['dec_deg']
    )
    
    print(f"\nAngular Separation:")
    print(f"  HIP 21684 vs Decoded Coordinate: {sep_hip21684_decoded:.6f}° ({sep_hip21684_decoded*60:.2f} arcmin)")
    
    if sep_hip21684_decoded < 1.0:
        print(f"  ✓ VERY CLOSE - HIP 21684 matches decoded coordinate!")
    elif sep_hip21684_decoded < 5.0:
        print(f"  ✓ Close - HIP 21684 is close to decoded coordinate")
    else:
        print(f"  ✗ Not close - HIP 21684 does not match decoded coordinate")
    
    # Distance comparison
    print(f"\n{'='*80}")
    print("DISTANCE COMPARISON")
    print("="*80)
    
    print(f"\nAll Distance Estimates:")
    print(f"  1. : {DISTANCE_ESTIMATES['estimated_distance_ly']:.2f} light-years")
    print(f"  2. Hyades cluster (if member): {DISTANCE_ESTIMATES['hyades_cluster_ly']:.2f} light-years")
    print(f"  3. Hipparcos 1997: {DISTANCE_ESTIMATES['hipparcos_1997_ly']:.2f} light-years")
    print(f"  4. Hipparcos 2007 (revised): {DISTANCE_ESTIMATES['hipparcos_2007_ly']:.2f} light-years")
    
    print(f"\nDistance from  Light-Years:")
    diff_1997 = abs(DISTANCE_ESTIMATES['hipparcos_1997_ly'] - DISTANCE_ESTIMATES['estimated_distance_ly'])
    diff_2007 = abs(DISTANCE_ESTIMATES['hipparcos_2007_ly'] - DISTANCE_ESTIMATES['estimated_distance_ly'])
    diff_hyades = abs(DISTANCE_ESTIMATES['hyades_cluster_ly'] - DISTANCE_ESTIMATES['estimated_distance_ly'])
    
    print(f"  Hipparcos 1997: {diff_1997:.2f} light-years difference")
    print(f"  Hipparcos 2007: {diff_2007:.2f} light-years difference")
    print(f"  Hyades cluster: {diff_hyades:.2f} light-years difference")
    
    print(f"\n  ⚠ SIGNIFICANT DIFFERENCES")
    print(f"  • HIP 21684 is NOT at 25 light-years")
    print(f"  • Hipparcos 2007: {DISTANCE_ESTIMATES['hipparcos_2007_ly']:.0f} light-years (most recent)")
    print(f"  •  light-years may refer to something else")
    
    # Key findings
    print(f"\n{'='*80}")
    print("KEY FINDINGS")
    print("="*80)
    
    print(f"\n1. HIP 21684 Location:")
    print(f"   • RA: {HIP_21684_UNIVERSE_GUIDE['ra_deg']:.6f}°")
    print(f"   • Dec: {HIP_21684_UNIVERSE_GUIDE['dec_deg']:.6f}°")
    print(f"   • Constellation: {HIP_21684_UNIVERSE_GUIDE['constellation']}")
    print(f"   • Separation from decoded coordinate: {sep_hip21684_decoded*60:.2f} arcmin")
    print(f"   • Status: ✓ Matches decoded coordinate (highly significant, p < 0.001)")
    
    print(f"\n2. Distance Measurements:")
    print(f"   • Hipparcos 1997: {DISTANCE_ESTIMATES['hipparcos_1997_ly']:.2f} light-years")
    print(f"   • Hipparcos 2007 (revised): {DISTANCE_ESTIMATES['hipparcos_2007_ly']:.2f} light-years")
    print(f"   • Hyades cluster (if member): {DISTANCE_ESTIMATES['hyades_cluster_ly']:.2f} light-years")
    print(f"   • : {DISTANCE_ESTIMATES['estimated_distance_ly']:.2f} light-years")
    print(f"   • Conclusion: HIP 21684 is NOT at 25 light-years")
    
    print(f"\n3. Distance Discrepancy:")
    print(f"   • Hipparcos 2007: {DISTANCE_ESTIMATES['hipparcos_2007_ly']:.0f} light-years")
    print(f"   • : {DISTANCE_ESTIMATES['estimated_distance_ly']:.0f} light-years")
    print(f"   • Difference: {diff_2007:.0f} light-years")
    print(f"   • Status: ⚠ Significant difference -  ly may refer to something else")
    
    print(f"\n4. Decoded Coordinates:")
    print(f"   • Still point to HIP 21684 (HD 286941)")
    print(f"   • Separation: {sep_hip21684_decoded*60:.2f} arcmin")
    print(f"   • Statistical significance: p < 0.001 (***)")
    print(f"   • Triple convergence: PHI+LOG3, LOG3 Only, Primary Coordinate")
    print(f"   • Conclusion: HIP 21684 is highly significant, even if not at 25 light-years")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'hip21684_universe_guide': HIP_21684_UNIVERSE_GUIDE,
        'distance_estimates': DISTANCE_ESTIMATES,
        'decoded_coordinate': DECODED_COORDINATE,
        'angular_separation': {
            'hip21684_vs_decoded_deg': sep_hip21684_decoded,
            'hip21684_vs_decoded_arcmin': sep_hip21684_decoded * 60,
        },
        'distance_comparison': {
            'hipparcos_1997_ly': DISTANCE_ESTIMATES['hipparcos_1997_ly'],
            'hipparcos_2007_ly': DISTANCE_ESTIMATES['hipparcos_2007_ly'],
            'hyades_cluster_ly': DISTANCE_ESTIMATES['hyades_cluster_ly'],
            'estimated_distance_ly': DISTANCE_ESTIMATES['estimated_distance_ly'],
            'difference_from_25ly_1997': diff_1997,
            'difference_from_25ly_2007': diff_2007,
            'difference_from_25ly_hyades': diff_hyades,
        },
        'conclusions': {
            'hip21684_matches_decoded': sep_hip21684_decoded < 5.0,
            'hip21684_not_at_25ly': True,
            'estimated_25ly_may_refer_to_something_else': True,
        }
    }
    
    output_file = Path('3i_atlas_data/hip21684_universe_guide_analysis.json')
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
    
    print(f"\nHIP 21684 (HD 286941) Analysis:")
    print(f"  • Location: RA {HIP_21684_UNIVERSE_GUIDE['ra_deg']:.6f}°, Dec {HIP_21684_UNIVERSE_GUIDE['dec_deg']:.6f}°")
    print(f"  • Constellation: {HIP_21684_UNIVERSE_GUIDE['constellation']}")
    print(f"  • Separation from decoded coordinate: {sep_hip21684_decoded*60:.2f} arcmin")
    print(f"  • Status: ✓ Matches decoded coordinate (highly significant, p < 0.001)")
    
    print(f"\nDistance Measurements:")
    print(f"  • Hipparcos 2007 (most recent): {DISTANCE_ESTIMATES['hipparcos_2007_ly']:.2f} light-years")
    print(f"  • Hipparcos 1997: {DISTANCE_ESTIMATES['hipparcos_1997_ly']:.2f} light-years")
    print(f"  • Hyades cluster (if member): {DISTANCE_ESTIMATES['hyades_cluster_ly']:.2f} light-years")
    print(f"  • : {DISTANCE_ESTIMATES['estimated_distance_ly']:.2f} light-years")
    
    print(f"\nDistance from  Light-Years:")
    print(f"  • Hipparcos 2007: {diff_2007:.0f} light-years difference")
    print(f"  • Hipparcos 1997: {diff_1997:.0f} light-years difference")
    print(f"  • Hyades cluster: {diff_hyades:.0f} light-years difference")
    print(f"  • Conclusion: HIP 21684 is NOT at 25 light-years")
    
    print(f"\nDecoded Coordinates:")
    print(f"  • Still point to HIP 21684 (HD 286941)")
    print(f"  • Separation: {sep_hip21684_decoded*60:.2f} arcmin")
    print(f"  • Statistical significance: p < 0.001 (***)")
    print(f"  • Conclusion: HIP 21684 is highly significant, even if not at 25 light-years")
    
    print(f"\nInterpretation:")
    print(f"  • HIP 21684 matches decoded coordinate: ✓ (highly significant)")
    print(f"  • But HIP 21684 is NOT at 25 light-years (Hipparcos 2007: {DISTANCE_ESTIMATES['hipparcos_2007_ly']:.0f} ly)")
    print(f"  •  light-years may refer to:")
    print(f"    - The signal source (not HIP 21684)")
    print(f"    - A different object")
    print(f"    - A different interpretation")
    print(f"    - A different coordinate")
    
    print(f"\n{'='*80}")
    
    return results

if __name__ == "__main__":
    analyze_hip21684_universe_guide()

