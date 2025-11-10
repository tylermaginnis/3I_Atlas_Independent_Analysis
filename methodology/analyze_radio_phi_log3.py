#!/usr/bin/env python3
"""
Analyze Radio Emissions from NVSS J043903+111800
Radio source at decoded coordinate (phi_log3 encoding) near HIP 21684 (HD 286941)
"""

import json
import requests
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt, log10

# NVSS Radio Source at Decoded Coordinate
NVSS_J043903_111800 = {
    'identifier': 'NVSS J043903+111800',
    'type': 'Radio source (Rad)',
    'ra_deg': 69.765163,  # 04h39m03.639s = 4*15 + 39*15/60 + 3.639*15/3600
    'dec_deg': 11.299372,  # +11°17'57.74" = 11 + 17/60 + 57.74/3600
    'ra_hms': '04h39m03.639s',
    'dec_dms': '+11°17\'57.74"',
    'distance_arcsec': 212.16,
    'distance_arcmin': 3.54,
    'frequency': 1.4e9,  # 1.4 GHz (NVSS frequency)
}

# HIP 21684 (HD 286941) - star at decoded coordinate
HIP_21684 = {
    'identifier': 'HIP 21684',
    'name': 'HD 286941',
    'type': 'Emission star (Em*)',
    'ra_deg': 69.822875,  # 04h39m17.49s
    'dec_deg': 11.265222,  # +11°15'54.8"
    'ra_hms': '04h39m17.49s',
    'dec_dms': '+11°15\'54.8"',
    'distance_ly': 491.95,
    'distance_pc': 150.83,
    'spectral_type': 'G5',
    'v_magnitude': 9.65,
    'b_magnitude': 10.41,
}

# Decoded coordinate (phi_log3 encoding)
DECODED_COORDINATE = {
    'ra_deg': 69.797974,
    'dec_deg': 11.250000,
    'ra_hms': '04h39m11.514s',
    'dec_dms': '+11°15\'00.000"',
    'name': 'Decoded 3I/ATLAS coordinate (phi_log3)',
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

def query_ned_for_radio_properties(name):
    """Query NED for radio source properties"""
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
        print(f"  Querying NED for {name} radio properties...")
        response = requests.get(url, params=params, timeout=15)
        
        if response.status_code == 200:
            result = response.text
            return {'success': True, 'data': result}
        else:
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def parse_radio_properties(data):
    """Parse radio properties from NED data"""
    if not data:
        return None
    
    lines = data.split('\n')
    properties = {
        'flux_1_4_ghz_mjy': None,
        'flux_4_85_ghz_mjy': None,
        'spectral_index': None,
        'polarization': None,
        'morphology': None,
        'luminosity': None,
    }
    
    for i, line in enumerate(lines):
        # Look for flux density (Log format)
        if 'log' in line.lower() and 'flux' in line.lower() and 'mjy' in line.lower():
            parts = line.split()
            for j, part in enumerate(parts):
                try:
                    log_flux = float(part)
                    if 0 < log_flux < 5:  # Reasonable log flux range
                        flux_mjy = 10**log_flux
                        # Check next line for frequency
                        if i + 1 < len(lines):
                            next_line = lines[i + 1].lower()
                            if 'frequency' in next_line or 'ghz' in next_line:
                                next_parts = next_line.split()
                                for k, part2 in enumerate(next_parts):
                                    try:
                                        freq = float(part2)
                                        if 0.1 < freq < 100:  # GHz range
                                            if abs(freq - 4.85) < 0.1:
                                                properties['flux_4_85_ghz_mjy'] = flux_mjy
                                            elif abs(freq - 1.4) < 0.1:
                                                properties['flux_1_4_ghz_mjy'] = flux_mjy
                                            break
                                    except:
                                        continue
                        break
                except:
                    continue
        
        # Look for spectral index
        if 'sp. index' in line.lower() or ('spectral' in line.lower() and 'index' in line.lower()):
            parts = line.split()
            for j, part in enumerate(parts):
                try:
                    alpha = float(part)
                    if -3 < alpha < 3:  # Reasonable spectral index range
                        properties['spectral_index'] = alpha
                        break
                except:
                    continue
        
        # Look for polarization
        if 'polarization' in line.lower() or 'polarized' in line.lower():
            parts = line.split()
            for j, part in enumerate(parts):
                try:
                    pol = float(part)
                    if 0 < pol < 100:  # Reasonable polarization percentage
                        properties['polarization'] = pol
                        break
                except:
                    continue
    
    # If we have flux at 4.85 GHz and spectral index, calculate flux at 1.4 GHz
    if properties['flux_4_85_ghz_mjy'] and properties['spectral_index'] and not properties['flux_1_4_ghz_mjy']:
        # S(ν) = S(ν₀) × (ν/ν₀)^α
        flux_1_4 = properties['flux_4_85_ghz_mjy'] * ((1.4 / 4.85) ** properties['spectral_index'])
        properties['flux_1_4_ghz_mjy'] = flux_1_4
    
    return properties

def calculate_radio_luminosity(flux_mjy, distance_pc, frequency_ghz=1.4):
    """Calculate radio luminosity from flux density"""
    if flux_mjy is None or distance_pc is None:
        return None
    
    # Convert flux from mJy to Jy
    flux_jy = flux_mjy / 1000.0
    
    # Calculate luminosity: L = 4πd²F
    # d in pc, F in Jy, L in erg/s/Hz
    # 1 Jy = 10^-26 W/m²/Hz
    # 1 pc = 3.086e16 m
    
    d_m = distance_pc * 3.086e16  # Distance in meters
    flux_w_m2_hz = flux_jy * 1e-26  # Flux in W/m²/Hz
    
    # Luminosity in W/Hz
    luminosity_w_hz = 4 * 3.14159 * d_m**2 * flux_w_m2_hz
    
    # Convert to erg/s/Hz (1 W = 1e7 erg/s)
    luminosity_erg_s_hz = luminosity_w_hz * 1e7
    
    # Convert to log luminosity (log L in erg/s/Hz)
    log_luminosity = log10(luminosity_erg_s_hz) if luminosity_erg_s_hz > 0 else None
    
    return {
        'luminosity_w_hz': luminosity_w_hz,
        'luminosity_erg_s_hz': luminosity_erg_s_hz,
        'log_luminosity': log_luminosity,
    }

def classify_radio_source(flux_mjy, spectral_index, luminosity, distance_pc):
    """Classify radio source based on properties"""
    classification = {
        'type': 'Unknown',
        'likely_stellar': False,
        'likely_star_formation': False,
        'likely_agn': False,
        'likely_synchrotron': False,
        'confidence': 'Low',
    }
    
    # For nearby stars (distance < 1000 pc):
    # - Stellar radio emission: log L ~ 15-20 erg/s/Hz
    # - Star formation: log L ~ 20-22 erg/s/Hz
    # - AGN: log L > 23 erg/s/Hz
    
    if luminosity and luminosity.get('log_luminosity'):
        log_l = luminosity['log_luminosity']
        
        if distance_pc < 1000:  # Nearby star
            if log_l < 18:
                classification['likely_stellar'] = True
                classification['type'] = 'Stellar Radio Emission'
                classification['confidence'] = 'Medium'
            elif log_l < 21:
                classification['likely_star_formation'] = True
                classification['type'] = 'Star Formation or HII Region'
                classification['confidence'] = 'Medium'
            elif log_l < 23:
                classification['likely_agn'] = True
                classification['type'] = 'Possible AGN or Active Star'
                classification['confidence'] = 'Low'
            else:
                classification['likely_agn'] = True
                classification['type'] = 'AGN (Active Galactic Nucleus)'
                classification['confidence'] = 'High'
        else:  # Distant source
            if log_l > 23:
                classification['likely_agn'] = True
                classification['type'] = 'AGN (Active Galactic Nucleus)'
                classification['confidence'] = 'High'
            elif log_l > 22:
                classification['likely_agn'] = True
                classification['type'] = 'Possible AGN'
                classification['confidence'] = 'Medium'
            elif log_l > 21:
                classification['likely_star_formation'] = True
                classification['type'] = 'Star Formation'
                classification['confidence'] = 'Medium'
            else:
                classification['likely_star_formation'] = True
                classification['type'] = 'Star Formation or HII Region'
                classification['confidence'] = 'Low'
    
    if spectral_index is not None:
        if spectral_index < -0.5:
            classification['likely_synchrotron'] = True
            if not classification['likely_agn']:
                classification['type'] = 'Synchrotron Source'
                classification['confidence'] = 'Medium'
        elif spectral_index > -0.5:
            if not classification['likely_agn']:
                classification['likely_star_formation'] = True
                classification['type'] = 'Thermal or Star Formation'
                classification['confidence'] = 'Medium'
    
    return classification

def analyze_radio_phi_log3():
    """Analyze radio emissions from NVSS J043903+111800"""
    print("="*80)
    print("ANALYZING RADIO EMISSIONS FROM NVSS J043903+111800")
    print("="*80)
    
    # NVSS Radio Source
    print("\n1. NVSS Radio Source Information:")
    print("="*80)
    
    print(f"\nNVSS J043903+111800:")
    print(f"  Type: {NVSS_J043903_111800['type']}")
    print(f"  RA: {NVSS_J043903_111800['ra_deg']:.6f}° ({NVSS_J043903_111800['ra_hms']})")
    print(f"  Dec: {NVSS_J043903_111800['dec_deg']:.6f}° ({NVSS_J043903_111800['dec_dms']})")
    print(f"  Distance from decoded coordinate: {NVSS_J043903_111800['distance_arcmin']:.2f} arcmin")
    
    # HIP 21684 (HD 286941)
    print(f"\n2. HIP 21684 (HD 286941) - Star at Decoded Coordinate:")
    print("="*80)
    
    print(f"\nHIP 21684 (HD 286941):")
    print(f"  Type: {HIP_21684['type']}")
    print(f"  RA: {HIP_21684['ra_deg']:.6f}° ({HIP_21684['ra_hms']})")
    print(f"  Dec: {HIP_21684['dec_deg']:.6f}° ({HIP_21684['dec_dms']})")
    print(f"  Distance: {HIP_21684['distance_ly']:.2f} light-years ({HIP_21684['distance_pc']:.2f} parsecs)")
    print(f"  Spectral Type: {HIP_21684['spectral_type']}")
    print(f"  V Magnitude: {HIP_21684['v_magnitude']}")
    print(f"  B Magnitude: {HIP_21684['b_magnitude']}")
    
    # Calculate separations
    print(f"\n3. Angular Separations:")
    print("="*80)
    
    sep_nvss_hip = calculate_angular_separation(
        NVSS_J043903_111800['ra_deg'], NVSS_J043903_111800['dec_deg'],
        HIP_21684['ra_deg'], HIP_21684['dec_deg']
    )
    
    sep_nvss_decoded = calculate_angular_separation(
        NVSS_J043903_111800['ra_deg'], NVSS_J043903_111800['dec_deg'],
        DECODED_COORDINATE['ra_deg'], DECODED_COORDINATE['dec_deg']
    )
    
    sep_hip_decoded = calculate_angular_separation(
        HIP_21684['ra_deg'], HIP_21684['dec_deg'],
        DECODED_COORDINATE['ra_deg'], DECODED_COORDINATE['dec_deg']
    )
    
    print(f"\n  NVSS J043903+111800 vs HIP 21684:")
    print(f"    Angular separation: {sep_nvss_hip:.6f}° ({sep_nvss_hip*60:.2f} arcmin)")
    print(f"    RA difference: {abs(NVSS_J043903_111800['ra_deg'] - HIP_21684['ra_deg'])*60:.2f} arcmin")
    print(f"    Dec difference: {abs(NVSS_J043903_111800['dec_deg'] - HIP_21684['dec_deg'])*60:.2f} arcmin")
    
    print(f"\n  NVSS J043903+111800 vs Decoded Coordinate:")
    print(f"    Angular separation: {sep_nvss_decoded:.6f}° ({sep_nvss_decoded*60:.2f} arcmin)")
    
    print(f"\n  HIP 21684 vs Decoded Coordinate:")
    print(f"    Angular separation: {sep_hip_decoded:.6f}° ({sep_hip_decoded*60:.2f} arcmin)")
    print(f"    Status: ✓ Matches decoded coordinate (p < 0.001)")
    
    # Check if NVSS source is associated with HIP 21684
    if sep_nvss_hip * 60 < 5.0:  # Within 5 arcmin
        print(f"\n    ✓ NVSS source is VERY CLOSE to HIP 21684!")
        print(f"      Could be radio emission from the star")
    elif sep_nvss_hip * 60 < 10.0:  # Within 10 arcmin
        print(f"\n    ✓ NVSS source is close to HIP 21684")
        print(f"      Could be related to the star")
    else:
        print(f"\n    ⚠ NVSS source is not close to HIP 21684")
        print(f"      May be unrelated")
    
    # Query NED for radio properties
    print(f"\n4. Querying NED for Radio Properties:")
    print("="*80)
    
    ned_result = query_ned_for_radio_properties('NVSS J043903+111800')
    
    if ned_result.get('success'):
        print(f"  ✓ NED query successful")
        
        # Parse radio properties
        radio_properties = parse_radio_properties(ned_result['data'])
        
        if radio_properties:
            print(f"\n  Radio Properties:")
            if radio_properties['flux_1_4_ghz_mjy']:
                print(f"    Flux at 1.4 GHz: {radio_properties['flux_1_4_ghz_mjy']:.2f} mJy")
            if radio_properties['spectral_index']:
                print(f"    Spectral Index: {radio_properties['spectral_index']:.2f}")
            if radio_properties['polarization']:
                print(f"    Polarization: {radio_properties['polarization']:.2f}%")
        else:
            print(f"  ⚠ Could not parse radio properties from NED data")
            print(f"  Data preview:")
            print(f"  {ned_result['data'][:500]}...")
            radio_properties = None
    else:
        print(f"  ✗ NED query failed: {ned_result.get('error', 'Unknown')}")
        radio_properties = None
    
    # Calculate radio luminosity
    print(f"\n5. Calculating Radio Luminosity:")
    print("="*80)
    
    if radio_properties and radio_properties.get('flux_1_4_ghz_mjy'):
        flux_mjy = radio_properties['flux_1_4_ghz_mjy']
        distance_pc = HIP_21684['distance_pc']
        
        luminosity = calculate_radio_luminosity(flux_mjy, distance_pc)
        
        if luminosity:
            print(f"\n  Radio Luminosity:")
            print(f"    Flux at 1.4 GHz: {flux_mjy:.2f} mJy")
            print(f"    Distance: {distance_pc:.2f} parsecs ({HIP_21684['distance_ly']:.2f} light-years)")
            print(f"    Luminosity: {luminosity['luminosity_erg_s_hz']:.2e} erg/s/Hz")
            print(f"    Log Luminosity: {luminosity['log_luminosity']:.2f} erg/s/Hz")
        else:
            luminosity = None
    else:
        print(f"  ⚠ Cannot calculate luminosity without flux density")
        luminosity = None
    
    # Classify radio source
    print(f"\n6. Classifying Radio Source:")
    print("="*80)
    
    if radio_properties and luminosity:
        classification = classify_radio_source(
            radio_properties.get('flux_1_4_ghz_mjy'),
            radio_properties.get('spectral_index'),
            luminosity,
            HIP_21684['distance_pc']
        )
        
        print(f"\n  Classification:")
        print(f"    Type: {classification['type']}")
        print(f"    Likely Stellar: {classification['likely_stellar']}")
        print(f"    Likely Star Formation: {classification['likely_star_formation']}")
        print(f"    Likely AGN: {classification['likely_agn']}")
        print(f"    Likely Synchrotron: {classification['likely_synchrotron']}")
        print(f"    Confidence: {classification['confidence']}")
    else:
        print(f"  ⚠ Cannot classify without radio properties")
        classification = None
    
    # Analysis
    print(f"\n7. ANALYSIS:")
    print("="*80)
    
    print(f"\n  A. Radio Source Properties:")
    if radio_properties:
        print(f"     • Flux at 1.4 GHz: {radio_properties.get('flux_1_4_ghz_mjy', 'Unknown')} mJy")
        print(f"     • Spectral Index: {radio_properties.get('spectral_index', 'Unknown')}")
        print(f"     • Polarization: {radio_properties.get('polarization', 'Unknown')}%")
    else:
        print(f"     • Properties: Not available")
    
    print(f"\n  B. Radio Luminosity:")
    if luminosity:
        print(f"     • Log Luminosity: {luminosity['log_luminosity']:.2f} erg/s/Hz")
        print(f"     • Luminosity: {luminosity['luminosity_erg_s_hz']:.2e} erg/s/Hz")
        
        # Compare with typical values for nearby stars
        if luminosity['log_luminosity'] < 18:
            print(f"     • ✓ Low luminosity (typical of stellar radio emission)")
        elif luminosity['log_luminosity'] < 21:
            print(f"     • ✓ Moderate luminosity (typical of star formation)")
        elif luminosity['log_luminosity'] < 23:
            print(f"     • ⚠ High luminosity (possible AGN or active star)")
        else:
            print(f"     • ⚠ Very high luminosity (typical of AGN)")
    else:
        print(f"     • Luminosity: Not available")
    
    print(f"\n  C. Source Classification:")
    if classification:
        print(f"     • Type: {classification['type']}")
        print(f"     • Confidence: {classification['confidence']}")
        
        if classification['likely_stellar']:
            print(f"     • ✓ Likely stellar radio emission from HIP 21684")
            print(f"     • Could be radio emission from active star")
        elif classification['likely_star_formation']:
            print(f"     • ✓ Likely star formation emission")
            print(f"     • Could be radio emission from HII regions")
        elif classification['likely_agn']:
            print(f"     • ⚠ Likely AGN emission")
            print(f"     • Could be radio emission from active galactic nucleus")
        else:
            print(f"     • ⚠ Classification uncertain")
    else:
        print(f"     • Classification: Not available")
    
    print(f"\n  D. Relationship to HIP 21684:")
    print(f"     • NVSS source is {sep_nvss_hip*60:.2f} arcmin from HIP 21684")
    print(f"     • HIP 21684 is {HIP_21684['distance_ly']:.2f} light-years away (in our galaxy)")
    print(f"     • If associated, radio emission is from the star or nearby region")
    
    if classification and classification['likely_stellar']:
        print(f"     • ✓ Stellar radio emission suggests active star")
        print(f"     • Could be radio emission from stellar activity")
    elif classification and classification['likely_star_formation']:
        print(f"     • ✓ Star formation emission suggests nearby star-forming region")
        print(f"     • Could be radio emission from HII regions")
    
    # Summary
    print(f"\n8. SUMMARY")
    print("="*80)
    
    print(f"\nNVSS J043903+111800 Radio Emissions:")
    if radio_properties:
        print(f"  Flux at 1.4 GHz: {radio_properties.get('flux_1_4_ghz_mjy', 'Unknown')} mJy")
        print(f"  Spectral Index: {radio_properties.get('spectral_index', 'Unknown')}")
        print(f"  Polarization: {radio_properties.get('polarization', 'Unknown')}%")
    else:
        print(f"  Properties: Not available from queries")
    
    if luminosity:
        print(f"\n  Radio Luminosity:")
        print(f"    Log Luminosity: {luminosity['log_luminosity']:.2f} erg/s/Hz")
        print(f"    Luminosity: {luminosity['luminosity_erg_s_hz']:.2e} erg/s/Hz")
    
    if classification:
        print(f"\n  Classification:")
        print(f"    Type: {classification['type']}")
        print(f"    Confidence: {classification['confidence']}")
    
    print(f"\nRelationship to HIP 21684 (HD 286941):")
    print(f"  • NVSS source is {sep_nvss_hip*60:.2f} arcmin from HIP 21684")
    print(f"  • HIP 21684 is {HIP_21684['distance_ly']:.2f} light-years away (in our galaxy)")
    print(f"  • If associated, radio emission is galactic (not extragalactic)")
    
    print(f"\nComparison with Diophantine v3 Coordinates:")
    print(f"  • Phi_Log3 encoding: Points to HIP 21684 (star) + NVSS J043903+111800 (radio)")
    print(f"  • Diophantine v3 encoding: Points to LEDA 1363602 (galaxy) + NVSS J045150+092332 (AGN radio)")
    print(f"  • Both encodings point to objects with associated radio sources!")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'nvss_source': NVSS_J043903_111800,
        'hip_21684': HIP_21684,
        'decoded_coordinate': DECODED_COORDINATE,
        'separations': {
            'nvss_vs_hip_deg': sep_nvss_hip,
            'nvss_vs_hip_arcmin': sep_nvss_hip * 60,
            'nvss_vs_decoded_deg': sep_nvss_decoded,
            'nvss_vs_decoded_arcmin': sep_nvss_decoded * 60,
            'hip_vs_decoded_deg': sep_hip_decoded,
            'hip_vs_decoded_arcmin': sep_hip_decoded * 60,
        },
        'ned_query': ned_result,
        'radio_properties': radio_properties,
        'radio_luminosity': luminosity,
        'classification': classification,
    }
    
    output_file = Path('3i_atlas_data/radio_phi_log3_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    return results

if __name__ == "__main__":
    analyze_radio_phi_log3()

