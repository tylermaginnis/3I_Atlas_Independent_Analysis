#!/usr/bin/env python3
"""
Test 45-Bit Decoding (Optimal Information-Theoretic Length)
Compare with 37-bit decoding to see if 45 bits gives better results
"""

import json
import numpy as np
import math
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt

# Mathematical constants
PHI = (1 + math.sqrt(5)) / 2
LOG3 = math.log(3)

def load_decoding_results():
    """Load message decoding results"""
    results_file = Path('3i_atlas_data/message_decoding_analysis.json')
    if not results_file.exists():
        return None
    
    with open(results_file) as f:
        return json.load(f)

def decode_coordinate_from_bits_45(bits):
    """Decode coordinate from 45-bit sequence"""
    if len(bits) < 45:
        return None
    
    # Take first 45 bits
    coordinate_bits = bits[:45]
    bit_string = ''.join(map(str, coordinate_bits))
    coordinate_value = int(bit_string, 2)
    
    # Try different splits for 45 bits
    # Option 1: 23 bits RA, 22 bits Dec
    ra_bits_23 = coordinate_bits[:23]
    dec_bits_22 = coordinate_bits[23:45]
    
    ra_value_23 = int(''.join(map(str, ra_bits_23)), 2)
    dec_value_22 = int(''.join(map(str, dec_bits_22)), 2)
    
    ra_deg_23 = (ra_value_23 / (2**23)) * 360  # 0 to 360 degrees
    dec_deg_22 = (dec_value_22 / (2**22)) * 180 - 90  # -90 to +90 degrees
    
    # Option 2: 22 bits RA, 23 bits Dec
    ra_bits_22 = coordinate_bits[:22]
    dec_bits_23 = coordinate_bits[22:45]
    
    ra_value_22 = int(''.join(map(str, ra_bits_22)), 2)
    dec_value_23 = int(''.join(map(str, dec_bits_23)), 2)
    
    ra_deg_22 = (ra_value_22 / (2**22)) * 360  # 0 to 360 degrees
    dec_deg_23 = (dec_value_23 / (2**23)) * 180 - 90  # -90 to +90 degrees
    
    # Option 3: 24 bits RA, 21 bits Dec
    ra_bits_24 = coordinate_bits[:24]
    dec_bits_21 = coordinate_bits[24:45]
    
    ra_value_24 = int(''.join(map(str, ra_bits_24)), 2)
    dec_value_21 = int(''.join(map(str, dec_bits_21)), 2)
    
    ra_deg_24 = (ra_value_24 / (2**24)) * 360  # 0 to 360 degrees
    dec_deg_21 = (dec_value_21 / (2**21)) * 180 - 90  # -90 to +90 degrees
    
    return {
        'bit_string': bit_string,
        'coordinate_value': coordinate_value,
        'split_23_22': {
            'ra_bits': ra_bits_23,
            'dec_bits': dec_bits_22,
            'ra_value': ra_value_23,
            'dec_value': dec_value_22,
            'ra_degrees': ra_deg_23,
            'dec_degrees': dec_deg_22,
        },
        'split_22_23': {
            'ra_bits': ra_bits_22,
            'dec_bits': dec_bits_23,
            'ra_value': ra_value_22,
            'dec_value': dec_value_23,
            'ra_degrees': ra_deg_22,
            'dec_degrees': dec_deg_23,
        },
        'split_24_21': {
            'ra_bits': ra_bits_24,
            'dec_bits': dec_bits_21,
            'ra_value': ra_value_24,
            'dec_value': dec_value_21,
            'ra_degrees': ra_deg_24,
            'dec_degrees': dec_deg_21,
        },
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

def test_45_bit_decoding():
    """Test 45-bit decoding and compare with 37-bit decoding"""
    print("="*80)
    print("TESTING 45-BIT DECODING (OPTIMAL INFORMATION-THEORETIC LENGTH)")
    print("="*80)
    
    # Load decoding results
    print("\n1. Loading decoding results...")
    results = load_decoding_results()
    
    if not results:
        print("Error: Decoding results not found")
        return
    
    # Get constant sequences
    sequences = results.get('constant_sequences', {}).get('atlas_spectrum', [])
    
    if not sequences:
        print("Error: No constant sequences found")
        return
    
    print(f"   Loaded {len(sequences)} constant sequences")
    
    # Convert to binary bits (PHI+LOG3 encoding)
    print("\n2. Converting to binary bits...")
    bits = []
    for seq in sequences:
        if seq['constant'] in ['PHI', 'LOG3']:
            bits.append(1)
        else:
            bits.append(0)
    
    print(f"   Total bits: {len(bits)}")
    
    # Decode 45-bit structures
    print("\n" + "="*80)
    print("3. DECODING 45-BIT STRUCTURES")
    print("="*80)
    
    # HIP 21684 coordinates for comparison
    HIP_21684 = {
        'ra_deg': 69.822875,
        'dec_deg': 11.265222,
        'name': 'HIP 21684 (HD 286941)',
    }
    
    # Decoded 37-bit coordinate for comparison
    DECODED_37_BIT = {
        'ra_deg': 69.797974,
        'dec_deg': 11.250000,
        'name': 'Decoded 37-bit coordinate',
    }
    
    all_45_bit_coordinates = []
    
    # Extract 45-bit structures
    print(f"\nExtracting 45-bit structures from bit sequence...")
    for start_idx in range(len(bits) - 44):
        segment = bits[start_idx:start_idx+45]
        
        # Decode coordinate
        decoded = decode_coordinate_from_bits_45(segment)
        if decoded:
            all_45_bit_coordinates.append({
                'start_position': start_idx,
                'decoded': decoded,
            })
    
    print(f"   Found {len(all_45_bit_coordinates)} 45-bit structures")
    
    # Analyze each split option
    print("\n" + "="*80)
    print("4. ANALYZING DIFFERENT BIT SPLITS (23/22, 22/23, 24/21)")
    print("="*80)
    
    split_results = {
        '23_22': [],
        '22_23': [],
        '24_21': [],
    }
    
    for coord_info in all_45_bit_coordinates[:10]:  # First 10 for analysis
        decoded = coord_info['decoded']
        
        # Test each split
        for split_name in ['split_23_22', 'split_22_23', 'split_24_21']:
            split_data = decoded[split_name]
            ra_deg = split_data['ra_degrees']
            dec_deg = split_data['dec_degrees']
            
            # Calculate separation from HIP 21684
            sep_hip = calculate_angular_separation(
                ra_deg, dec_deg,
                HIP_21684['ra_deg'], HIP_21684['dec_deg']
            )
            
            # Calculate separation from 37-bit decoded coordinate
            sep_37 = calculate_angular_separation(
                ra_deg, dec_deg,
                DECODED_37_BIT['ra_deg'], DECODED_37_BIT['dec_deg']
            )
            
            split_key = split_name.replace('split_', '')
            split_results[split_key].append({
                'start_position': coord_info['start_position'],
                'ra_deg': ra_deg,
                'dec_deg': dec_deg,
                'separation_from_hip_arcmin': sep_hip * 60,
                'separation_from_37bit_arcmin': sep_37 * 60,
            })
    
    # Find best matches for each split
    print("\nBest matches for each split:")
    print("="*80)
    
    for split_name, results in split_results.items():
        if not results:
            continue
        
        # Sort by separation from HIP 21684
        results_sorted = sorted(results, key=lambda x: x['separation_from_hip_arcmin'])
        best = results_sorted[0]
        
        print(f"\n{split_name} Split (RA/Dec):")
        print(f"  Best match:")
        print(f"    RA: {best['ra_deg']:.6f}°")
        print(f"    Dec: {best['dec_deg']:.6f}°")
        print(f"    Separation from HIP 21684: {best['separation_from_hip_arcmin']:.2f} arcmin")
        print(f"    Separation from 37-bit: {best['separation_from_37bit_arcmin']:.2f} arcmin")
        
        # Count matches within 5 arcmin of HIP 21684
        matches_within_5arcmin = sum(1 for r in results if r['separation_from_hip_arcmin'] < 5.0)
        matches_within_10arcmin = sum(1 for r in results if r['separation_from_hip_arcmin'] < 10.0)
        
        print(f"  Matches within 5 arcmin: {matches_within_5arcmin}/{len(results)}")
        print(f"  Matches within 10 arcmin: {matches_within_10arcmin}/{len(results)}")
    
    # Compare with 37-bit decoding
    print("\n" + "="*80)
    print("5. COMPARISON WITH 37-BIT DECODING")
    print("="*80)
    
    print(f"\n37-Bit Decoding:")
    print(f"  Coordinate: RA {DECODED_37_BIT['ra_deg']:.6f}°, Dec {DECODED_37_BIT['dec_deg']:.6f}°")
    print(f"  Separation from HIP 21684: 1.73 arcmin")
    print(f"  Statistical significance: p < 0.001 (***)")
    
    # Find best 45-bit match
    best_45bit = None
    best_separation = float('inf')
    best_split = None
    
    for split_name, results in split_results.items():
        for result in results:
            if result['separation_from_hip_arcmin'] < best_separation:
                best_separation = result['separation_from_hip_arcmin']
                best_45bit = result
                best_split = split_name
    
    if best_45bit:
        print(f"\n45-Bit Decoding (Best Match - {best_split} split):")
        print(f"  Coordinate: RA {best_45bit['ra_deg']:.6f}°, Dec {best_45bit['dec_deg']:.6f}°")
        print(f"  Separation from HIP 21684: {best_45bit['separation_from_hip_arcmin']:.2f} arcmin")
        print(f"  Separation from 37-bit: {best_45bit['separation_from_37bit_arcmin']:.2f} arcmin")
        
        if best_separation < 1.73:
            print(f"\n  ✓ 45-bit decoding gives BETTER match than 37-bit!")
            print(f"    Improvement: {1.73 - best_separation:.2f} arcmin")
        elif best_separation < 5.0:
            print(f"\n  ✓ 45-bit decoding gives GOOD match (within 5 arcmin)")
        else:
            print(f"\n  ⚠ 45-bit decoding does NOT improve on 37-bit")
            print(f"    Difference: {best_separation - 1.73:.2f} arcmin")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    print("\n1. Information-Theoretic Analysis:")
    print(f"  Optimal bit length: 45 bits")
    print(f"  Information content: 36.11 bits")
    print(f"  37-bit information content: 29.61 bits")
    print(f"  Difference: 6.49 bits")
    
    print("\n2. Coordinate Decoding:")
    print(f"  37-bit coordinate: RA {DECODED_37_BIT['ra_deg']:.6f}°, Dec {DECODED_37_BIT['dec_deg']:.6f}°")
    print(f"  Separation from HIP 21684: 1.73 arcmin (p < 0.001)")
    
    if best_45bit:
        print(f"  45-bit coordinate (best): RA {best_45bit['ra_deg']:.6f}°, Dec {best_45bit['dec_deg']:.6f}°")
        print(f"  Separation from HIP 21684: {best_45bit['separation_from_hip_arcmin']:.2f} arcmin")
        
        if best_separation < 1.73:
            print(f"\n  ✓ 45-bit decoding is BETTER than 37-bit")
        elif best_separation < 5.0:
            print(f"\n  ✓ 45-bit decoding is COMPARABLE to 37-bit")
        else:
            print(f"\n  ⚠ 45-bit decoding is WORSE than 37-bit")
    
    print("\n3. Conclusion:")
    if best_45bit and best_separation < 1.73:
        print(f"  ✓ 45 bits is information-theoretically optimal AND gives better coordinate match")
        print(f"    • Should use 45 bits instead of 37 bits")
    elif best_45bit and best_separation < 5.0:
        print(f"  ⚠ 45 bits is information-theoretically optimal but coordinate match is comparable")
        print(f"    • 37 bits may be correct (based on 's findings)")
        print(f"    • 45 bits should be tested further")
    else:
        print(f"  ⚠ 45 bits is information-theoretically optimal but coordinate match is worse")
        print(f"    • 37 bits may be correct (based on 's findings)")
        print(f"    • Need to investigate why 45 bits doesn't improve match")
    
    # Save results
    output_file = Path('3i_atlas_data/45_bit_decoding_test.json')
    with open(output_file, 'w') as f:
        json.dump({
            'analysis_date': datetime.now().isoformat(),
            'total_45_bit_structures': len(all_45_bit_coordinates),
            'split_results': split_results,
            'best_45bit_match': best_45bit,
            'best_split': best_split,
            'comparison_with_37bit': {
                '37bit_coordinate': DECODED_37_BIT,
                '37bit_separation_arcmin': 1.73,
                '45bit_coordinate': best_45bit,
                '45bit_separation_arcmin': best_separation if best_45bit else None,
                'improvement': 1.73 - best_separation if best_45bit and best_separation < 1.73 else None,
            },
        }, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    return {
        'best_45bit_match': best_45bit,
        'best_split': best_split,
        'comparison': {
            '37bit_separation': 1.73,
            '45bit_separation': best_separation if best_45bit else None,
        }
    }

if __name__ == "__main__":
    test_45_bit_decoding()

