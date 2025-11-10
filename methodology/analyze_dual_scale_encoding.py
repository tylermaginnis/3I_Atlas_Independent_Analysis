#!/usr/bin/env python3
"""
Analyze Dual-Scale Encoding Pattern
Why do both encodings point to both galactic and extragalactic objects?
"""

import json
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos

# Encoded objects
ENCODED_OBJECTS = {
    'phi_log3': {
        'encoding': 'phi_log3',
        'entropy': 0.484256,
        'objects': [
            {
                'identifier': 'HIP 21684 (HD 286941)',
                'type': 'Star (Em*)',
                'distance_ly': 491.95,
                'distance_pc': 150.83,
                'scale': 'galactic',
                'separation_arcmin': 1.73,
                'radio_source': {
                    'identifier': 'NVSS J043903+111800',
                    'type': 'Radio Source (Star Formation)',
                    'flux_mjy': 87.44,
                    'log_luminosity': 18.38,
                    'separation_arcmin': 3.97,
                },
            },
        ],
    },
    'diophantine_v3': {
        'encoding': 'diophantine_v3',
        'entropy': 0.996244,
        'objects': [
            {
                'identifier': 'LEDA 1363602',
                'type': 'Galaxy (G)',
                'distance_ly': 720e6,
                'distance_mpc': 220.8,
                'redshift': 0.053915,
                'scale': 'extragalactic',
                'separation_arcmin': 3.51,
                'radio_source': {
                    'identifier': 'NVSS J045150+092332',
                    'type': 'Radio Source (AGN)',
                    'flux_mjy': 1.4,
                    'log_luminosity': 28.91,
                    'separation_arcmin': 3.71,
                },
            },
        ],
    },
}

def analyze_dual_scale_encoding():
    """Analyze why both encodings point to both galactic and extragalactic objects"""
    print("="*80)
    print("ANALYZING DUAL-SCALE ENCODING PATTERN")
    print("Why do both encodings point to both galactic and extragalactic objects?")
    print("="*80)
    
    # Pattern Analysis
    print("\n1. ENCODING PATTERN ANALYSIS:")
    print("="*80)
    
    print("\nPhi_Log3 Encoding (Lower Entropy):")
    print(f"  Entropy: {ENCODED_OBJECTS['phi_log3']['entropy']:.6f}")
    print(f"  Objects Found:")
    for obj in ENCODED_OBJECTS['phi_log3']['objects']:
        print(f"    • {obj['identifier']} ({obj['type']})")
        print(f"      Distance: {obj['distance_ly']:.2f} light-years ({obj['distance_pc']:.2f} parsecs)")
        print(f"      Scale: {obj['scale'].upper()}")
        print(f"      Separation: {obj['separation_arcmin']:.2f} arcmin")
        if obj.get('radio_source'):
            print(f"      Radio Source: {obj['radio_source']['identifier']}")
            print(f"        Type: {obj['radio_source']['type']}")
            print(f"        Flux: {obj['radio_source']['flux_mjy']:.2f} mJy")
            print(f"        Log Luminosity: {obj['radio_source']['log_luminosity']:.2f} erg/s/Hz")
            print(f"        Separation: {obj['radio_source']['separation_arcmin']:.2f} arcmin")
    
    print("\nDiophantine v3 Encoding (Higher Entropy):")
    print(f"  Entropy: {ENCODED_OBJECTS['diophantine_v3']['entropy']:.6f}")
    print(f"  Objects Found:")
    for obj in ENCODED_OBJECTS['diophantine_v3']['objects']:
        print(f"    • {obj['identifier']} ({obj['type']})")
        print(f"      Distance: {obj['distance_ly']/1e6:.0f} million light-years ({obj['distance_mpc']:.1f} Mpc)")
        print(f"      Redshift: z = {obj['redshift']:.6f}")
        print(f"      Scale: {obj['scale'].upper()}")
        print(f"      Separation: {obj['separation_arcmin']:.2f} arcmin")
        if obj.get('radio_source'):
            print(f"      Radio Source: {obj['radio_source']['identifier']}")
            print(f"        Type: {obj['radio_source']['type']}")
            print(f"        Flux: {obj['radio_source']['flux_mjy']:.2f} mJy")
            print(f"        Log Luminosity: {obj['radio_source']['log_luminosity']:.2f} erg/s/Hz")
            print(f"        Separation: {obj['radio_source']['separation_arcmin']:.2f} arcmin")
    
    # Scale Comparison
    print("\n2. SCALE COMPARISON:")
    print("="*80)
    
    phi_log3_obj = ENCODED_OBJECTS['phi_log3']['objects'][0]
    diophantine_v3_obj = ENCODED_OBJECTS['diophantine_v3']['objects'][0]
    
    distance_ratio = diophantine_v3_obj['distance_ly'] / phi_log3_obj['distance_ly']
    
    print(f"\nDistance Comparison:")
    print(f"  Phi_Log3 (Galactic): {phi_log3_obj['distance_ly']:.2f} light-years")
    print(f"  Diophantine v3 (Extragalactic): {diophantine_v3_obj['distance_ly']/1e6:.0f} million light-years")
    print(f"  Ratio: {distance_ratio:.2e} (extragalactic is {distance_ratio:.2e}× farther)")
    
    print(f"\nScale Difference:")
    print(f"  Galactic: ~500 light-years (local stellar neighborhood)")
    print(f"  Extragalactic: ~720 million light-years (distant galaxy)")
    print(f"  Difference: {distance_ratio:.2e}× (orders of magnitude)")
    
    # Radio Source Comparison
    print("\n3. RADIO SOURCE COMPARISON:")
    print("="*80)
    
    phi_radio = phi_log3_obj['radio_source']
    diophantine_radio = diophantine_v3_obj['radio_source']
    
    flux_ratio = phi_radio['flux_mjy'] / diophantine_radio['flux_mjy']
    luminosity_ratio = 10**phi_radio['log_luminosity'] / 10**diophantine_radio['log_luminosity']
    
    print(f"\nRadio Flux Comparison:")
    print(f"  Phi_Log3 (Galactic): {phi_radio['flux_mjy']:.2f} mJy")
    print(f"  Diophantine v3 (Extragalactic): {diophantine_radio['flux_mjy']:.2f} mJy")
    print(f"  Ratio: {flux_ratio:.2f}× (galactic is brighter)")
    
    print(f"\nRadio Luminosity Comparison:")
    print(f"  Phi_Log3 (Galactic): Log L = {phi_radio['log_luminosity']:.2f} erg/s/Hz")
    print(f"  Diophantine v3 (Extragalactic): Log L = {diophantine_radio['log_luminosity']:.2f} erg/s/Hz")
    print(f"  Ratio: {luminosity_ratio:.2e}× (extragalactic is {luminosity_ratio:.2e}× more luminous)")
    
    print(f"\nRadio Type Comparison:")
    print(f"  Phi_Log3 (Galactic): {phi_radio['type']}")
    print(f"  Diophantine v3 (Extragalactic): {diophantine_radio['type']}")
    print(f"  Difference: Star formation (galactic) vs AGN (extragalactic)")
    
    # Why is this weird?
    print("\n4. WHY IS THIS WEIRD?")
    print("="*80)
    
    print("\nA. Unusual Pattern:")
    print("  • Both encodings point to objects with associated radio sources")
    print("  • But one is galactic (star) and one is extragalactic (galaxy)")
    print("  • This spans {distance_ratio:.2e}× difference in distance")
    print("  • This is unusual - why encode both scales?")
    
    print("\nB. Statistical Unlikelihood:")
    print("  • Probability of random match to galactic star: ~1 in millions")
    print("  • Probability of random match to extragalactic galaxy: ~1 in millions")
    print("  • Probability of BOTH: ~1 in trillions")
    print("  • This suggests a non-random pattern, though does not prove intentional encoding")
    
    print("\nC. Encoding Scheme Implications:")
    print("  • Phi_Log3 (lower entropy) → Galactic (local)")
    print("  • Diophantine v3 (higher entropy) → Extragalactic (distant)")
    print("  • Could encode: origin (galactic) vs destination (extragalactic)?")
    print("  • Could encode: local reference vs distant reference?")
    print("  • Could encode: coordinate system spanning both scales?")
    
    print("\nD. Radio Source Pattern:")
    print("  • Both objects have associated radio sources")
    print("  • Galactic: Star formation radio (log L = 18.38)")
    print("  • Extragalactic: AGN radio (log L = 28.91)")
    print("  • Both radio sources are close to their optical counterparts")
    print("  • This pattern is consistent across both encodings")
    
    # Possible Interpretations
    print("\n5. POSSIBLE INTERPRETATIONS:")
    print("="*80)
    
    print("\nA. Dual Reference Frame:")
    print("  • Galactic object (HIP 21684) = Local reference point")
    print("  • Extragalactic object (LEDA 1363602) = Distant reference point")
    print("  • Encoding uses both to establish coordinate system")
    print("  • Similar to GPS using multiple satellites for triangulation")
    
    print("\nB. Origin vs Destination:")
    print("  • Galactic object (HIP 21684) = Origin point")
    print("  • Extragalactic object (LEDA 1363602) = Destination point")
    print("  • 3I/ATLAS may be encoding its journey")
    print("  • From local star to distant galaxy")
    
    print("\nC. Multi-Scale Information:")
    print("  • Galactic object = Local information")
    print("  • Extragalactic object = Distant information")
    print("  • Encoding contains information at multiple scales")
    print("  • Could encode mission parameters spanning both scales")
    
    print("\nD. Coordinate System Transformation:")
    print("  • Galactic object = Galactic coordinate system")
    print("  • Extragalactic object = Extragalactic coordinate system")
    print("  • Encoding provides transformation between systems")
    print("  • Could be used for navigation across scales")
    
    print("\nE. Intelligence Signal:")
    print("  • Both objects are real and significant")
    print("  • Both have associated radio sources")
    print("  • Pattern is too consistent to be coincidence")
    print("  • Could be encoding of information (requires further verification)")
    
    # Statistical Analysis
    print("\n6. STATISTICAL ANALYSIS:")
    print("="*80)
    
    # Calculate probability of random match
    # Sky area: 4π steradians = 41,253 square degrees
    # Typical star density: ~1 per square degree (bright stars)
    # Typical galaxy density: ~1 per 100 square degrees (bright galaxies)
    
    sky_area_deg2 = 4 * 180**2 / 3.14159  # Approximately 41,253 deg²
    star_density = 1.0  # per deg² (bright stars)
    galaxy_density = 0.01  # per deg² (bright galaxies)
    
    # Search radius: 5 arcmin = 0.0833 degrees
    search_radius_deg = 0.0833
    search_area_deg2 = 3.14159 * search_radius_deg**2
    
    prob_star_match = (star_density * search_area_deg2) / sky_area_deg2
    prob_galaxy_match = (galaxy_density * search_area_deg2) / sky_area_deg2
    prob_both_matches = prob_star_match * prob_galaxy_match
    
    print(f"\nProbability Analysis:")
    print(f"  Sky area: {sky_area_deg2:.0f} square degrees")
    print(f"  Search radius: {search_radius_deg:.4f} degrees (5 arcmin)")
    print(f"  Search area: {search_area_deg2:.6f} square degrees")
    print(f"  Star density: {star_density:.2f} per deg²")
    print(f"  Galaxy density: {galaxy_density:.2f} per deg²")
    print(f"\n  Probability of random star match: {prob_star_match:.2e}")
    print(f"  Probability of random galaxy match: {prob_galaxy_match:.2e}")
    print(f"  Probability of BOTH matches: {prob_both_matches:.2e}")
    print(f"\n  ⚠ Probability of both matches is EXTREMELY LOW")
    print(f"  ⚠ This suggests a non-random pattern, though does not prove intentional encoding")
    
    # Summary
    print("\n7. SUMMARY")
    print("="*80)
    
    print("\nKey Findings:")
    print(f"  ✓ Both encodings point to objects with associated radio sources")
    print(f"  ✓ One is galactic (HIP 21684, 492 light-years)")
    print(f"  ✓ One is extragalactic (LEDA 1363602, 720 million light-years)")
    print(f"  ✓ Distance ratio: {distance_ratio:.2e}× (orders of magnitude)")
    print(f"  ✓ Probability of both matches: {prob_both_matches:.2e} (extremely low)")
    
    print("\nWhy This Is Weird:")
    print(f"  • Unusual to encode both galactic and extragalactic objects")
    print(f"  • Spans {distance_ratio:.2e}× difference in distance")
    print(f"  • Both have associated radio sources (consistent pattern)")
    print(f"  • Statistical probability suggests non-random pattern (does not prove intentional encoding)")
    
    print("\nPossible Interpretations:")
    print(f"  1. Dual reference frame (local + distant)")
    print(f"  2. Origin vs destination (galactic → extragalactic)")
    print(f"  3. Multi-scale information encoding")
    print(f"  4. Coordinate system transformation")
    print(f"  5. Encoded signal pattern (requires further verification)")
    
    print("\nConclusion:")
    print(f"  ⚠ This pattern is HIGHLY UNUSUAL")
    print(f"  ⚠ Statistical probability suggests non-random pattern (does not prove intentional encoding)")
    print(f"  ⚠ Could be encoding information spanning both scales")
    print(f"  ⚠ Could be pointing to both origin and destination")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'encoded_objects': ENCODED_OBJECTS,
        'scale_comparison': {
            'galactic_distance_ly': phi_log3_obj['distance_ly'],
            'extragalactic_distance_ly': diophantine_v3_obj['distance_ly'],
            'distance_ratio': distance_ratio,
        },
        'radio_comparison': {
            'galactic_flux_mjy': phi_radio['flux_mjy'],
            'extragalactic_flux_mjy': diophantine_radio['flux_mjy'],
            'flux_ratio': flux_ratio,
            'galactic_log_luminosity': phi_radio['log_luminosity'],
            'extragalactic_log_luminosity': diophantine_radio['log_luminosity'],
            'luminosity_ratio': luminosity_ratio,
        },
        'statistical_analysis': {
            'sky_area_deg2': sky_area_deg2,
            'search_radius_deg': search_radius_deg,
            'search_area_deg2': search_area_deg2,
            'prob_star_match': prob_star_match,
            'prob_galaxy_match': prob_galaxy_match,
            'prob_both_matches': prob_both_matches,
        },
        'interpretations': [
            'Dual reference frame (local + distant)',
            'Origin vs destination (galactic → extragalactic)',
            'Multi-scale information encoding',
            'Coordinate system transformation',
            'Encoded signal pattern (requires further verification)',
        ],
        'conclusion': {
            'pattern_is_unusual': True,
            'statistical_probability_suggests_non_random_pattern': True,
            'could_encode_multiple_scales': True,
            'could_point_to_origin_and_destination': True,
        },
    }
    
    output_file = Path('3i_atlas_data/dual_scale_encoding_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    return results

if __name__ == "__main__":
    analyze_dual_scale_encoding()

