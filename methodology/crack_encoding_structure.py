#!/usr/bin/env python3
"""
Crack the Encoding Structure - Structural Analysis
Analyze the structure of constant sequences to determine the actual encoding scheme
Similar to how we analyzed Lehmer sequences and Diophantine equations
"""

import json
import numpy as np
import math
from pathlib import Path
from datetime import datetime
from collections import Counter, defaultdict
from scipy import stats
from scipy.stats import entropy

# Mathematical constants
CONSTANTS = {
    'PHI': (1 + math.sqrt(5)) / 2,
    'E': math.e,
    'PI': math.pi,
    'TAU': 2 * math.pi,
    'SQRT2': math.sqrt(2),
    'SQRT3': math.sqrt(3),
    'SQRT5': math.sqrt(5),
    'GAMMA': 0.5772156649015329,
    'LOG2': math.log(2),
    'LOG3': math.log(3),
    'LOG_PHI': math.log((1 + math.sqrt(5)) / 2),
}

def load_decoding_results():
    """Load message decoding results"""
    results_file = Path('3i_atlas_data/message_decoding_analysis.json')
    if not results_file.exists():
        return None
    
    with open(results_file) as f:
        return json.load(f)

def analyze_constant_structure(sequences):
    """Analyze the structure of constant sequences"""
    if not sequences:
        return None
    
    # Group sequences by constant
    constant_groups = defaultdict(list)
    for seq in sequences:
        constant_groups[seq['constant']].append(seq)
    
    # Analyze patterns
    structure_analysis = {
        'constant_distribution': {k: len(v) for k, v in constant_groups.items()},
        'constant_positions': {k: [s['position'] for s in v] for k, v in constant_groups.items()},
        'constant_intervals': {},
        'constant_patterns': {},
    }
    
    # Calculate intervals between constants
    for const_name, positions in structure_analysis['constant_positions'].items():
        if len(positions) > 1:
            intervals = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
            structure_analysis['constant_intervals'][const_name] = {
                'intervals': intervals,
                'mean': float(np.mean(intervals)) if intervals else None,
                'std': float(np.std(intervals)) if intervals else None,
                'min': int(min(intervals)) if intervals else None,
                'max': int(max(intervals)) if intervals else None,
            }
    
    # Look for patterns in constant sequences
    # Pattern 1: Alternating patterns
    # Pattern 2: Repeating patterns
    # Pattern 3: Fibonacci-like patterns
    # Pattern 4: Diophantine structures
    
    # Analyze constant transitions
    constant_transitions = []
    for i in range(len(sequences) - 1):
        const1 = sequences[i]['constant']
        const2 = sequences[i+1]['constant']
        constant_transitions.append((const1, const2))
    
    transition_counts = Counter(constant_transitions)
    structure_analysis['constant_transitions'] = {f"{k[0]}_{k[1]}": v for k, v in transition_counts.items()}
    
    # Find most common transitions
    most_common_transitions = transition_counts.most_common(10)
    structure_analysis['most_common_transitions'] = [{'from': t[0][0], 'to': t[0][1], 'count': t[1]} for t in most_common_transitions]
    
    return structure_analysis

def analyze_diophantine_structure(sequences):
    """Analyze Diophantine structure in constant sequences"""
    if not sequences:
        return None
    
    # Convert constants to (2,3) energy framework
    energy_structures = []
    for seq in sequences:
        const_name = seq['constant']
        const_value = CONSTANTS.get(const_name, 0)
        
        if const_value > 0:
            # Express constant in (2,3) framework
            # E(n) = 2^v₂(n) × 3^v₃(n)
            # For constants, we need to find v₂ and v₃
            
            # Approximate: log₂(E) = v₂ + v₃ * log₂(3)
            # log₃(E) = v₂ * log₃(2) + v₃
            
            log2_val = math.log2(const_value)
            log3_val = math.log(const_value) / math.log(3)
            
            # Solve for v₂ and v₃
            # v₂ = (log2_val - log3_val * log2(3)) / (1 - log2(3) * log3(2))
            # v₃ = (log3_val - log2_val * log3(2)) / (1 - log2(3) * log3(2))
            
            # Simplified: round to nearest integers
            v2 = round(log2_val)
            v3 = round(log3_val)
            
            energy = (2 ** v2) * (3 ** v3)
            
            energy_structures.append({
                'position': seq['position'],
                'constant': const_name,
                'value': const_value,
                'v2': int(v2),
                'v3': int(v3),
                'energy_23': float(energy),
                'log2': float(log2_val),
                'log3': float(log3_val),
            })
    
    # Analyze energy patterns
    energy_analysis = {
        'energy_structures': energy_structures,
        'v2_distribution': dict(Counter([e['v2'] for e in energy_structures])),
        'v3_distribution': dict(Counter([e['v3'] for e in energy_structures])),
        'energy_patterns': {},
    }
    
    # Look for patterns in v2 and v3
    v2_sequence = [e['v2'] for e in energy_structures]
    v3_sequence = [e['v3'] for e in energy_structures]
    
    # Analyze v2 mod 2 and v3 mod 2 (bit encoding)
    v2_mod2 = [v % 2 for v in v2_sequence]
    v3_mod2 = [v % 2 for v in v3_sequence]
    
    energy_analysis['v2_mod2'] = v2_mod2
    energy_analysis['v3_mod2'] = v3_mod2
    energy_analysis['v2_mod2_distribution'] = dict(Counter(v2_mod2))
    energy_analysis['v3_mod2_distribution'] = dict(Counter(v3_mod2))
    
    # Combined bit encoding
    combined_bits = [v2_mod2[i] * 2 + v3_mod2[i] for i in range(len(v2_mod2))]
    energy_analysis['combined_bits'] = combined_bits
    energy_analysis['combined_bits_distribution'] = dict(Counter(combined_bits))
    
    return energy_analysis

def analyze_lehmer_like_structure(sequences):
    """Analyze Lehmer-like structure in constant sequences"""
    if not sequences:
        return None
    
    # Convert to ratio sequence (like Lehmer sequences)
    ratio_sequence = []
    for i in range(len(sequences) - 1):
        const1 = sequences[i]['constant']
        const2 = sequences[i+1]['constant']
        
        val1 = CONSTANTS.get(const1, 0)
        val2 = CONSTANTS.get(const2, 0)
        
        if val1 > 0:
            ratio = val2 / val1
            ratio_sequence.append({
                'position': i,
                'ratio': float(ratio),
                'const1': const1,
                'const2': const2,
            })
    
    # Analyze ratio patterns
    ratio_analysis = {
        'ratio_sequence': ratio_sequence,
        'ratio_values': [r['ratio'] for r in ratio_sequence],
        'ratio_statistics': {
            'mean': float(np.mean([r['ratio'] for r in ratio_sequence])) if ratio_sequence else None,
            'std': float(np.std([r['ratio'] for r in ratio_sequence])) if ratio_sequence else None,
            'min': float(min([r['ratio'] for r in ratio_sequence])) if ratio_sequence else None,
            'max': float(max([r['ratio'] for r in ratio_sequence])) if ratio_sequence else None,
        },
    }
    
    # Look for convergence to constants (like Lehmer sequences)
    # Check if ratios converge to PHI, LOG3, etc.
    convergence_analysis = {}
    for const_name, const_value in CONSTANTS.items():
        if const_value > 0:
            # Check how many ratios are close to this constant
            matches = [r for r in ratio_sequence if abs(r['ratio'] - const_value) < 0.1]
            convergence_analysis[const_name] = {
                'matches': len(matches),
                'percentage': len(matches) / len(ratio_sequence) * 100 if ratio_sequence else 0,
            }
    
    ratio_analysis['convergence_analysis'] = convergence_analysis
    
    return ratio_analysis

def test_encoding_schemes(sequences):
    """Test different encoding schemes to find the correct one"""
    if not sequences:
        return None
    
    # Different encoding schemes to test
    encoding_schemes = {
        'phi_log3': lambda s: 1 if s['constant'] in ['PHI', 'LOG3'] else 0,
        'phi_only': lambda s: 1 if s['constant'] == 'PHI' else 0,
        'log3_only': lambda s: 1 if s['constant'] == 'LOG3' else 0,
        'log2_log3': lambda s: 1 if s['constant'] in ['LOG2', 'LOG3'] else 0,
        'phi_log2': lambda s: 1 if s['constant'] in ['PHI', 'LOG2'] else 0,
        'all_logs': lambda s: 1 if s['constant'] in ['LOG2', 'LOG3', 'LOG_PHI'] else 0,
        'diophantine_v2': lambda s: 1 if (CONSTANTS.get(s['constant'], 0) > 0 and (round(math.log2(CONSTANTS.get(s['constant'], 1))) % 2 == 1)) else 0,
        'diophantine_v3': lambda s: 1 if (CONSTANTS.get(s['constant'], 0) > 0 and (round(math.log(CONSTANTS.get(s['constant'], 1)) / math.log(3)) % 2 == 1)) else 0,
    }
    
    # Test each encoding scheme
    scheme_results = {}
    
    for scheme_name, encoding_func in encoding_schemes.items():
        bits = []
        for seq in sequences:
            try:
                bit = encoding_func(seq)
                bits.append(bit)
            except:
                bits.append(0)
        
        # Analyze bit sequence
        bit_counts = Counter(bits)
        entropy_val = entropy([bit_counts[0]/len(bits), bit_counts[1]/len(bits)], base=2) if len(bits) > 0 else 0
        
        # Look for 37-bit structures
        bit_structures_37 = []
        for start_idx in range(len(bits) - 36):
            segment = bits[start_idx:start_idx+37]
            bit_string = ''.join(map(str, segment))
            coordinate_value = int(bit_string, 2)
            bit_structures_37.append({
                'start_position': start_idx,
                'bit_string': bit_string,
                'coordinate_value': coordinate_value,
            })
        
        scheme_results[scheme_name] = {
            'bits': bits,
            'bit_counts': dict(bit_counts),
            'entropy': float(entropy_val),
            'bit_structures_37': bit_structures_37,
            'num_structures': len(bit_structures_37),
        }
    
    return scheme_results

def crack_encoding_structure():
    """Crack the encoding structure through structural analysis"""
    print("="*80)
    print("CRACKING ENCODING STRUCTURE - STRUCTURAL ANALYSIS")
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
    
    # Analyze constant structure
    print("\n" + "="*80)
    print("2. ANALYZING CONSTANT STRUCTURE")
    print("="*80)
    
    structure_analysis = analyze_constant_structure(sequences)
    
    if structure_analysis:
        print(f"\nConstant Distribution:")
        for const_name, count in sorted(structure_analysis['constant_distribution'].items(), key=lambda x: x[1], reverse=True):
            print(f"  {const_name}: {count} occurrences")
        
        print(f"\nMost Common Transitions:")
        for transition in structure_analysis['most_common_transitions'][:10]:
            print(f"  {transition['from']} → {transition['to']}: {transition['count']} times")
        
        print(f"\nConstant Intervals:")
        for const_name, interval_data in structure_analysis['constant_intervals'].items():
            if interval_data['mean']:
                print(f"  {const_name}: mean = {interval_data['mean']:.2f}, std = {interval_data['std']:.2f}")
    
    # Analyze Diophantine structure
    print("\n" + "="*80)
    print("3. ANALYZING DIOPHANTINE STRUCTURE")
    print("="*80)
    
    diophantine_analysis = analyze_diophantine_structure(sequences)
    
    if diophantine_analysis:
        print(f"\nDiophantine Energy Structure:")
        print(f"  Total energy structures: {len(diophantine_analysis['energy_structures'])}")
        
        print(f"\nV2 Distribution (mod 2):")
        for bit, count in sorted(diophantine_analysis['v2_mod2_distribution'].items()):
            print(f"  {bit}: {count} occurrences")
        
        print(f"\nV3 Distribution (mod 2):")
        for bit, count in sorted(diophantine_analysis['v3_mod2_distribution'].items()):
            print(f"  {bit}: {count} occurrences")
        
        print(f"\nCombined Bits Distribution:")
        for bit, count in sorted(diophantine_analysis['combined_bits_distribution'].items()):
            print(f"  {bit}: {count} occurrences")
        
        # Test if Diophantine encoding gives better results
        print(f"\nTesting Diophantine Encoding (v2 mod 2, v3 mod 2)...")
        diophantine_bits = diophantine_analysis['v2_mod2']
        print(f"  Total bits: {len(diophantine_bits)}")
        print(f"  Bit distribution: {Counter(diophantine_bits)}")
    
    # Analyze Lehmer-like structure
    print("\n" + "="*80)
    print("4. ANALYZING LEHMER-LIKE STRUCTURE")
    print("="*80)
    
    lehmer_analysis = analyze_lehmer_like_structure(sequences)
    
    if lehmer_analysis:
        print(f"\nRatio Sequence Analysis:")
        stats = lehmer_analysis['ratio_statistics']
        print(f"  Mean ratio: {stats['mean']:.6f}")
        print(f"  Std ratio: {stats['std']:.6f}")
        print(f"  Min ratio: {stats['min']:.6f}")
        print(f"  Max ratio: {stats['max']:.6f}")
        
        print(f"\nConvergence Analysis:")
        for const_name, conv_data in sorted(lehmer_analysis['convergence_analysis'].items(), key=lambda x: x[1]['matches'], reverse=True):
            if conv_data['matches'] > 0:
                print(f"  {const_name}: {conv_data['matches']} matches ({conv_data['percentage']:.2f}%)")
    
    # Test different encoding schemes
    print("\n" + "="*80)
    print("5. TESTING DIFFERENT ENCODING SCHEMES")
    print("="*80)
    
    scheme_results = test_encoding_schemes(sequences)
    
    if scheme_results:
        print(f"\nEncoding Scheme Comparison:")
        print(f"{'Scheme':<20} {'Entropy':<12} {'Structures':<12} {'Bit Dist':<20}")
        print("-" * 70)
        
        for scheme_name, result in scheme_results.items():
            bit_dist = f"{result['bit_counts'].get(0, 0)}/{result['bit_counts'].get(1, 0)}"
            print(f"{scheme_name:<20} {result['entropy']:<12.6f} {result['num_structures']:<12} {bit_dist:<20}")
        
        # Find best encoding scheme (highest entropy, most structures)
        best_scheme = max(scheme_results.items(), key=lambda x: x[1]['entropy'] * x[1]['num_structures'])
        print(f"\nBest Encoding Scheme: {best_scheme[0]}")
        print(f"  Entropy: {best_scheme[1]['entropy']:.6f}")
        print(f"  Structures: {best_scheme[1]['num_structures']}")
    
    # Summary
    print("\n" + "="*80)
    print("STRUCTURAL ANALYSIS SUMMARY")
    print("="*80)
    
    print("\n1. Constant Structure:")
    if structure_analysis:
        print(f"  • {len(structure_analysis['constant_distribution'])} different constants")
        print(f"  • Most common: {max(structure_analysis['constant_distribution'].items(), key=lambda x: x[1])}")
        if structure_analysis['most_common_transitions']:
            first_transition = structure_analysis['most_common_transitions'][0]
            print(f"  • Most common transition: {first_transition['from']} → {first_transition['to']} ({first_transition['count']} times)")
        else:
            print(f"  • Most common transition: None")
    
    print("\n2. Diophantine Structure:")
    if diophantine_analysis:
        print(f"  • Energy structures: {len(diophantine_analysis['energy_structures'])}")
        print(f"  • V2 mod 2 distribution: {dict(diophantine_analysis['v2_mod2_distribution'])}")
        print(f"  • V3 mod 2 distribution: {dict(diophantine_analysis['v3_mod2_distribution'])}")
        print(f"  • Combined bits distribution: {dict(diophantine_analysis['combined_bits_distribution'])}")
    
    print("\n3. Lehmer-Like Structure:")
    if lehmer_analysis:
        print(f"  • Ratio sequence length: {len(lehmer_analysis['ratio_sequence'])}")
        print(f"  • Mean ratio: {lehmer_analysis['ratio_statistics']['mean']:.6f}")
        best_conv = max(lehmer_analysis['convergence_analysis'].items(), key=lambda x: x[1]['matches'])
        print(f"  • Best convergence: {best_conv[0]} ({best_conv[1]['matches']} matches, {best_conv[1]['percentage']:.2f}%)")
    
    print("\n4. Encoding Schemes:")
    if scheme_results:
        print(f"  • Tested {len(scheme_results)} encoding schemes")
        print(f"  • Best scheme: {best_scheme[0]}")
        print(f"  • Current scheme (phi_log3): {scheme_results.get('phi_log3', {}).get('entropy', 0):.6f}")
    
    # Save results
    output_file = Path('3i_atlas_data/encoding_structure_analysis.json')
    with open(output_file, 'w') as f:
        json.dump({
            'analysis_date': datetime.now().isoformat(),
            'constant_structure': structure_analysis,
            'diophantine_structure': diophantine_analysis,
            'lehmer_like_structure': lehmer_analysis,
            'encoding_schemes': {k: {**v, 'bits': v['bits'][:100]} for k, v in scheme_results.items()},  # Save first 100 bits only
        }, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    return {
        'constant_structure': structure_analysis,
        'diophantine_structure': diophantine_analysis,
        'lehmer_like_structure': lehmer_analysis,
        'encoding_schemes': scheme_results,
    }

if __name__ == "__main__":
    crack_encoding_structure()

