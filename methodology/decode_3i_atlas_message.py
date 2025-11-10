#!/usr/bin/env python3
"""
Decode 3I/ATLAS Message from Constant Patterns
Analyze WHERE and HOW constants appear to decode underlying message
"""

import json
import numpy as np
import math
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# Mathematical constants
CONSTANTS = {
    'PHI': (1 + math.sqrt(5)) / 2,  # Golden ratio ≈ 1.618
    'E': math.e,  # Euler's number ≈ 2.718
    'PI': math.pi,  # Pi ≈ 3.142
    'TAU': 2 * math.pi,  # Tau ≈ 6.283
    'SQRT2': math.sqrt(2),  # ≈ 1.414
    'SQRT3': math.sqrt(3),  # ≈ 1.732
    'SQRT5': math.sqrt(5),  # ≈ 2.236
    'GAMMA': 0.5772156649015329,  # Euler-Mascheroni constant
    'LOG2': math.log(2),  # ≈ 0.693
    'LOG3': math.log(3),  # ≈ 1.099
    'LOG_PHI': math.log((1 + math.sqrt(5)) / 2),  # ≈ 0.481
}

def find_closest_constant(value, tolerance=0.1):
    """Find the closest known constant to a value"""
    closest = None
    min_diff = float('inf')
    for name, const_val in CONSTANTS.items():
        diff = abs(value - const_val)
        if diff < min_diff:
            min_diff = diff
            closest = (name, const_val, diff)
    if closest and closest[2] < tolerance:
        return closest
    return None

def analyze_constant_positions(data, data_name):
    """Analyze WHERE constants appear in the data"""
    if data is None:
        return None
    
    results = {
        'data_name': data_name,
        'constant_positions': defaultdict(list),
        'constant_sequences': [],
        'message_patterns': [],
    }
    
    # Flatten data for analysis
    if hasattr(data, 'flatten'):
        flat_data = data.flatten()
    else:
        flat_data = np.array(data).flatten()
    
    flat_data = flat_data[np.isfinite(flat_data)]
    
    # Analyze each position
    for i, value in enumerate(flat_data):
        const_match = find_closest_constant(value, tolerance=0.1)
        if const_match:
            name, const_val, diff = const_match
            results['constant_positions'][name].append({
                'position': int(i),
                'value': float(value),
                'constant': float(const_val),
                'diff': float(diff)
            })
    
    # Find sequences of constants
    constant_sequence = []
    for i in range(len(flat_data) - 1):
        if flat_data[i] != 0 and np.isfinite(flat_data[i+1]):
            ratio = flat_data[i+1] / flat_data[i]
            if np.isfinite(ratio):
                const_match = find_closest_constant(ratio, tolerance=0.1)
                if const_match:
                    constant_sequence.append({
                        'position': i,
                        'ratio': float(ratio),
                        'constant': const_match[0],
                        'value1': float(flat_data[i]),
                        'value2': float(flat_data[i+1])
                    })
    
    results['constant_sequences'] = constant_sequence
    
    # Look for message patterns
    # Pattern 1: Constant sequences that might encode bits
    if len(constant_sequence) > 0:
        # Group consecutive constants
        message_bits = []
        current_constant = None
        current_count = 0
        
        for seq in constant_sequence:
            if seq['constant'] == current_constant:
                current_count += 1
            else:
                if current_constant:
                    message_bits.append({
                        'constant': current_constant,
                        'count': current_count,
                        'bit_value': 1 if current_constant in ['PHI', 'LOG3', 'LOG2'] else 0
                    })
                current_constant = seq['constant']
                current_count = 1
        
        if current_constant:
            message_bits.append({
                'constant': current_constant,
                'count': current_count,
                'bit_value': 1 if current_constant in ['PHI', 'LOG3', 'LOG2'] else 0
            })
        
        results['message_patterns'] = message_bits
    
    return results

def decode_binary_message(constant_sequence, encoding_scheme='phi_log3'):
    """Decode binary message from constant sequence"""
    if not constant_sequence:
        return None
    
    # Encoding schemes
    if encoding_scheme == 'phi_log3':
        # PHI = 1, LOG3 = 1, others = 0
        bits = []
        for seq in constant_sequence:
            if seq['constant'] in ['PHI', 'LOG3']:
                bits.append(1)
            else:
                bits.append(0)
    elif encoding_scheme == 'phi_only':
        # PHI = 1, others = 0
        bits = []
        for seq in constant_sequence:
            if seq['constant'] == 'PHI':
                bits.append(1)
            else:
                bits.append(0)
    elif encoding_scheme == 'log3_only':
        # LOG3 = 1, others = 0
        bits = []
        for seq in constant_sequence:
            if seq['constant'] == 'LOG3':
                bits.append(1)
            else:
                bits.append(0)
    else:
        return None
    
    return bits

def analyze_37_bit_structure(constant_sequence):
    """Analyze for 37-bit coordinate structure"""
    if len(constant_sequence) < 37:
        return None
    
    # Look for 37-bit patterns
    for start_idx in range(len(constant_sequence) - 36):
        segment = constant_sequence[start_idx:start_idx+37]
        
        # Decode as binary
        bits = []
        for seq in segment:
            if seq['constant'] in ['PHI', 'LOG3', 'LOG2']:
                bits.append(1)
            else:
                bits.append(0)
        
        # Convert to integer
        bit_string = ''.join(map(str, bits))
        coordinate_value = int(bit_string, 2)
        
        # Check if it's a valid coordinate range
        # 37 bits = 2^37 ≈ 1.37×10^11 possible values
        # For RA/Dec: RA in [0, 24h] = [0, 86400] arcseconds
        # Dec in [-90, +90] degrees = [-324000, +324000] arcseconds
        
        yield {
            'start_position': start_idx,
            'bits': bits,
            'bit_string': bit_string,
            'coordinate_value': coordinate_value,
            'segment': segment
        }

def analyze_energy_23_message(data):
    """Analyze (2,3) energy framework for message encoding"""
    if data is None:
        return None
    
    results = {
        'energy_values': [],
        'energy_sequence': [],
        'message_bits': [],
    }
    
    # Flatten data
    if hasattr(data, 'flatten'):
        flat_data = data.flatten()
    else:
        flat_data = np.array(data).flatten()
    
    flat_data = flat_data[np.isfinite(flat_data)]
    flat_data = flat_data[flat_data > 0]
    
    # Calculate energy for each value
    for i, value in enumerate(flat_data[:1000]):  # Sample for performance
        log2_val = math.log2(value)
        log3_val = math.log(value) / math.log(3)
        
        v2 = round(log2_val)
        v3 = round(log3_val)
        
        energy = (2 ** v2) * (3 ** v3)
        
        results['energy_values'].append({
            'position': i,
            'value': float(value),
            'v2': int(v2),
            'v3': int(v3),
            'energy_23': float(energy)
        })
        
        # Encode as bits: v2 mod 2, v3 mod 2
        bit1 = v2 % 2
        bit2 = v3 % 2
        results['message_bits'].append({
            'position': i,
            'bit1': bit1,
            'bit2': bit2,
            'combined_bit': bit1 * 2 + bit2  # 0, 1, 2, or 3
        })
    
    return results

def decode_coordinate_from_bits(bits):
    """Decode coordinate from bit sequence"""
    if len(bits) < 37:
        return None
    
    # Take first 37 bits
    coordinate_bits = bits[:37]
    bit_string = ''.join(map(str, coordinate_bits))
    coordinate_value = int(bit_string, 2)
    
    # Interpret as RA/Dec
    # 37 bits = 2^37 ≈ 1.37×10^11
    # Split into RA (18 bits) and Dec (19 bits) or similar
    
    # Option 1: RA = first 18 bits, Dec = next 19 bits
    ra_bits = coordinate_bits[:18]
    dec_bits = coordinate_bits[18:37]
    
    ra_value = int(''.join(map(str, ra_bits)), 2)
    dec_value = int(''.join(map(str, dec_bits)), 2)
    
    # Convert to degrees/arcseconds
    # RA: 18 bits = 2^18 = 262,144 possible values
    # If interpreted as arcseconds: 262,144 / 3600 = 72.8 hours (needs scaling)
    # Dec: 19 bits = 2^19 = 524,288 possible values
    # If interpreted as arcseconds: 524,288 / 3600 = 145.6 degrees (needs scaling)
    
    # Alternative: Interpret as fractional coordinates
    ra_deg = (ra_value / (2**18)) * 360  # 0 to 360 degrees
    dec_deg = (dec_value / (2**19)) * 180 - 90  # -90 to +90 degrees
    
    return {
        'bit_string': bit_string,
        'coordinate_value': coordinate_value,
        'ra_bits': ra_bits,
        'dec_bits': dec_bits,
        'ra_value': ra_value,
        'dec_value': dec_value,
        'ra_degrees': ra_deg,
        'dec_degrees': dec_deg,
        'ra_hours': ra_deg / 15,  # Convert to hours
        'dec_arcseconds': dec_deg * 3600,
    }

def analyze_all_data_for_message():
    """Analyze all 3I/ATLAS data for decodable message"""
    print("="*80)
    print("DECODING 3I/ATLAS MESSAGE FROM CONSTANT PATTERNS")
    print("="*80)
    
    data_dir = Path('3i_atlas_data')
    if not data_dir.exists():
        print("Error: Data directory not found")
        return
    
    # Load data files
    data_files = {
        'atlas_spectrum': data_dir / '3I-Reflectivity_JN_edit.csv',
        'atlas_lightcurve': data_dir / 'ATLAS.npy',
        'trappist_lightcurve': data_dir / 'TRAPPIST.npy',
    }
    
    all_results = {
        'analysis_date': datetime.now().isoformat(),
        'constant_positions': {},
        'constant_sequences': {},
        'message_patterns': {},
        '37_bit_structures': {},
        'energy_23_messages': {},
        'decoded_messages': {},
    }
    
    for key, filepath in data_files.items():
        if not filepath.exists():
            continue
        
        print(f"\n{'='*80}")
        print(f"ANALYZING: {filepath.name}")
        print(f"{'='*80}")
        
        # Load data
        if filepath.suffix == '.npy':
            data = np.load(filepath)
        elif filepath.suffix == '.csv':
            data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
        else:
            continue
        
        # Analyze constant positions
        print(f"\n1. Constant Position Analysis:")
        position_analysis = analyze_constant_positions(data, key)
        if position_analysis:
            all_results['constant_positions'][key] = position_analysis['constant_positions']
            all_results['constant_sequences'][key] = position_analysis['constant_sequences']
            
            print(f"   Total constant matches: {sum(len(v) for v in position_analysis['constant_positions'].values())}")
            for const_name, positions in position_analysis['constant_positions'].items():
                if len(positions) > 0:
                    print(f"   {const_name}: {len(positions)} matches")
                    # Show first few positions
                    print(f"     Positions: {[p['position'] for p in positions[:5]]}")
            
            # Analyze sequences
            sequences = position_analysis['constant_sequences']
            print(f"\n2. Constant Sequence Analysis:")
            print(f"   Total sequences: {len(sequences)}")
            
            if len(sequences) > 0:
                # Group by constant
                sequence_by_const = defaultdict(list)
                for seq in sequences:
                    sequence_by_const[seq['constant']].append(seq)
                
                print(f"   Sequences by constant:")
                for const_name, seqs in sorted(sequence_by_const.items(), key=lambda x: len(x[1]), reverse=True):
                    print(f"     {const_name}: {len(seqs)} sequences")
                
                # Look for 37-bit structures
                print(f"\n3. 37-Bit Structure Analysis:")
                bit_structures = list(analyze_37_bit_structure(sequences))
                if bit_structures:
                    all_results['37_bit_structures'][key] = []
                    for i, structure in enumerate(bit_structures[:5]):  # First 5
                        print(f"   Structure {i+1}:")
                        print(f"     Start position: {structure['start_position']}")
                        print(f"     Bit string: {structure['bit_string'][:50]}...")
                        print(f"     Coordinate value: {structure['coordinate_value']:,}")
                        
                        # Decode coordinate
                        decoded = decode_coordinate_from_bits(structure['bits'])
                        if decoded:
                            print(f"     RA: {decoded['ra_degrees']:.6f}° ({decoded['ra_hours']:.6f} hours)")
                            print(f"     Dec: {decoded['dec_degrees']:.6f}° ({decoded['dec_arcseconds']:.2f} arcsec)")
                            all_results['37_bit_structures'][key].append(decoded)
                else:
                    print(f"   No 37-bit structures found")
                
                # Decode binary messages
                print(f"\n4. Binary Message Decoding:")
                for encoding_scheme in ['phi_log3', 'phi_only', 'log3_only']:
                    bits = decode_binary_message(sequences, encoding_scheme)
                    if bits and len(bits) >= 37:
                        print(f"   Encoding: {encoding_scheme}")
                        print(f"     Total bits: {len(bits)}")
                        print(f"     First 37 bits: {''.join(map(str, bits[:37]))}")
                        
                        # Decode coordinate
                        decoded = decode_coordinate_from_bits(bits)
                        if decoded:
                            print(f"     RA: {decoded['ra_degrees']:.6f}°")
                            print(f"     Dec: {decoded['dec_degrees']:.6f}°")
                            all_results['decoded_messages'][f"{key}_{encoding_scheme}"] = decoded
        
        # Analyze (2,3) energy framework
        print(f"\n5. Energy (2,3) Framework Message:")
        energy_message = analyze_energy_23_message(data)
        if energy_message:
            all_results['energy_23_messages'][key] = energy_message
            print(f"   Energy values analyzed: {len(energy_message['energy_values'])}")
            print(f"   Message bits: {len(energy_message['message_bits'])}")
            
            # Show first few bits
            if energy_message['message_bits']:
                first_bits = [b['combined_bit'] for b in energy_message['message_bits'][:37]]
                print(f"   First 37 bits (combined): {first_bits}")
    
    # Save results
    output_file = data_dir / 'message_decoding_analysis.json'
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("MESSAGE DECODING ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    # Summary
    print(f"\n{'='*80}")
    print("MESSAGE DECODING SUMMARY")
    print(f"{'='*80}")
    
    if all_results['37_bit_structures']:
        print(f"\n37-Bit Structures Found:")
        for key, structures in all_results['37_bit_structures'].items():
            print(f"  {key}: {len(structures)} structures")
            for i, struct in enumerate(structures):
                print(f"    Structure {i+1}:")
                print(f"      RA: {struct['ra_degrees']:.6f}° ({struct['ra_hours']:.6f} hours)")
                print(f"      Dec: {struct['dec_degrees']:.6f}°")
    
    if all_results['decoded_messages']:
        print(f"\nDecoded Messages:")
        for key, message in all_results['decoded_messages'].items():
            print(f"  {key}:")
            print(f"    RA: {message['ra_degrees']:.6f}°")
            print(f"    Dec: {message['dec_degrees']:.6f}°")
    
    print(f"\n{'='*80}")
    
    return all_results

if __name__ == "__main__":
    analyze_all_data_for_message()

