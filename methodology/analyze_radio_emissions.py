#!/usr/bin/env python3
"""
Analyze Radio Emissions from NVSS J045150+092332
Query NVSS catalog for flux density, spectral index, polarization, and morphology
"""

import json
import requests
from pathlib import Path
from datetime import datetime
from math import radians, degrees, sin, cos, acos, atan2, sqrt, log10

# NVSS Radio Source
NVSS_J045150_092332 = {
    'identifier': 'NVSS J045150+092332',
    'type': 'Radio source (Rad)',
    'ra_deg': 72.961750,  # 04h51m50.82s
    'dec_deg': 9.392222,  # +09°23'32.0"
    'ra_hms': '04h51m50.82s',
    'dec_dms': '+09°23\'32.0"',
    'frequency': 1.4e9,  # 1.4 GHz (NVSS frequency)
}

# LEDA 1363602 (galaxy)
LEDA_1363602 = {
    'identifier': 'LEDA 1363602',
    'type': 'Galaxy (G)',
    'ra_deg': 72.939583,
    'dec_deg': 9.334444,
    'distance_ly': 720e6,  # 720 million light-years
    'distance_mpc': 220.8,
    'redshift': 0.053915,
}

def query_nvss_catalog(ra_deg, dec_deg, radius_arcmin=10.0):
    """Query NVSS catalog for sources near coordinates"""
    # NVSS catalog browser URL
    url = "http://www.cv.nrao.edu/nvss/cgi-bin/NVSSlist.cgi"
    
    # Convert RA/Dec to NVSS format
    ra_h = int(ra_deg / 15.0)
    ra_m = int((ra_deg / 15.0 - ra_h) * 60)
    ra_s = ((ra_deg / 15.0 - ra_h) * 60 - ra_m) * 60
    
    dec_sign = '+' if dec_deg >= 0 else '-'
    dec_abs = abs(dec_deg)
    dec_d = int(dec_abs)
    dec_m = int((dec_abs - dec_d) * 60)
    dec_s = ((dec_abs - dec_d) * 60 - dec_m) * 60
    
    params = {
        'RA': f"{ra_h:02d}{ra_m:02d}{int(ra_s):02d}",
        'DEC': f"{dec_sign}{dec_d:02d}{dec_m:02d}{int(dec_s):02d}",
        'Radius': str(radius_arcmin / 60.0),  # Convert to degrees
        'RadiusUnits': 'deg',
    }
    
    try:
        print(f"  Querying NVSS catalog for RA {ra_deg:.6f}°, Dec {dec_deg:.6f}°...")
        response = requests.get(url, params=params, timeout=15)
        
        if response.status_code == 200:
            result = response.text
            return {'success': True, 'data': result}
        else:
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def parse_nvss_catalog(data):
    """Parse NVSS catalog data"""
    if not data:
        return None
    
    lines = data.split('\n')
    sources = []
    
    # Look for catalog table
    in_table = False
    for i, line in enumerate(lines):
        # Look for table header
        if 'RA' in line and 'DEC' in line and 'Flux' in line:
            in_table = True
            continue
        
        if in_table:
            # Parse source lines
            if line.strip() and not line.startswith('-') and not line.startswith('='):
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        # Try to extract source information
                        # Format: RA DEC Flux(mJy) ...
                        ra_str = parts[0] if len(parts) > 0 else ''
                        dec_str = parts[1] if len(parts) > 1 else ''
                        flux_str = parts[2] if len(parts) > 2 else ''
                        
                        # Try to parse flux
                        flux_mjy = None
                        try:
                            flux_mjy = float(flux_str)
                        except:
                            pass
                        
                        sources.append({
                            'ra': ra_str,
                            'dec': dec_str,
                            'flux_mjy': flux_mjy,
                            'raw_line': line,
                        })
                    except:
                        continue
    
    return sources

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
        'spectral_index': None,
        'polarization': None,
        'morphology': None,
        'luminosity': None,
    }
    
    for line in lines:
        # Look for flux density
        if '1.4' in line.lower() and ('ghz' in line.lower() or 'mhz' in line.lower()):
            parts = line.split()
            for i, part in enumerate(parts):
                try:
                    flux = float(part)
                    if 0.1 < flux < 10000:  # Reasonable flux range (mJy)
                        properties['flux_1_4_ghz_mjy'] = flux
                        break
                except:
                    continue
        
        # Look for spectral index
        if 'spectral' in line.lower() and 'index' in line.lower():
            parts = line.split()
            for i, part in enumerate(parts):
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
            for i, part in enumerate(parts):
                try:
                    pol = float(part)
                    if 0 < pol < 100:  # Reasonable polarization percentage
                        properties['polarization'] = pol
                        break
                except:
                    continue
    
    return properties

def calculate_radio_luminosity(flux_mjy, distance_mpc, frequency_ghz=1.4):
    """Calculate radio luminosity from flux density"""
    if flux_mjy is None or distance_mpc is None:
        return None
    
    # Convert flux from mJy to Jy
    flux_jy = flux_mjy / 1000.0
    
    # Calculate luminosity: L = 4πd²F
    # d in Mpc, F in Jy, L in W/Hz
    # 1 Jy = 10^-26 W/m²/Hz
    # 1 Mpc = 3.086e22 m
    
    d_m = distance_mpc * 3.086e22  # Distance in meters
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

def classify_radio_source(flux_mjy, spectral_index, luminosity):
    """Classify radio source based on properties"""
    classification = {
        'type': 'Unknown',
        'likely_agn': False,
        'likely_star_formation': False,
        'likely_synchrotron': False,
        'confidence': 'Low',
    }
    
    # AGN typically have:
    # - High luminosity (log L > 23 erg/s/Hz at 1.4 GHz)
    # - Steep spectral index (α < -0.5)
    # - High flux density
    
    # Star formation typically has:
    # - Lower luminosity (log L < 23 erg/s/Hz at 1.4 GHz)
    # - Flatter spectral index (α > -0.5)
    # - Lower flux density
    
    if luminosity and luminosity.get('log_luminosity'):
        log_l = luminosity['log_luminosity']
        
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

def analyze_radio_emissions():
    """Analyze radio emissions from NVSS J045150+092332"""
    print("="*80)
    print("ANALYZING RADIO EMISSIONS FROM NVSS J045150+092332")
    print("="*80)
    
    # Query NVSS catalog
    print("\n1. Querying NVSS Catalog:")
    print("="*80)
    
    nvss_result = query_nvss_catalog(
        NVSS_J045150_092332['ra_deg'],
        NVSS_J045150_092332['dec_deg'],
        radius_arcmin=10.0
    )
    
    if nvss_result.get('success'):
        print(f"  ✓ NVSS catalog query successful")
        
        # Parse catalog
        sources = parse_nvss_catalog(nvss_result['data'])
        
        if sources:
            print(f"  ✓ Found {len(sources)} sources in NVSS catalog")
            print(f"\n  Sources found:")
            for i, source in enumerate(sources[:5], 1):  # Show first 5
                print(f"    {i}. RA: {source['ra']}, Dec: {source['dec']}, Flux: {source['flux_mjy']} mJy")
        else:
            print(f"  ⚠ Could not parse NVSS catalog data")
            print(f"  Data preview:")
            print(f"  {nvss_result['data'][:500]}...")
    else:
        print(f"  ✗ NVSS catalog query failed: {nvss_result.get('error', 'Unknown')}")
        sources = None
    
    # Query NED for radio properties
    print(f"\n2. Querying NED for Radio Properties:")
    print("="*80)
    
    ned_result = query_ned_for_radio_properties('NVSS J045150+092332')
    
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
    print(f"\n3. Calculating Radio Luminosity:")
    print("="*80)
    
    if radio_properties and radio_properties.get('flux_1_4_ghz_mjy'):
        flux_mjy = radio_properties['flux_1_4_ghz_mjy']
        distance_mpc = LEDA_1363602['distance_mpc']
        
        luminosity = calculate_radio_luminosity(flux_mjy, distance_mpc)
        
        if luminosity:
            print(f"\n  Radio Luminosity:")
            print(f"    Flux at 1.4 GHz: {flux_mjy:.2f} mJy")
            print(f"    Distance: {distance_mpc:.1f} Mpc")
            print(f"    Luminosity: {luminosity['luminosity_erg_s_hz']:.2e} erg/s/Hz")
            print(f"    Log Luminosity: {luminosity['log_luminosity']:.2f} erg/s/Hz")
        else:
            luminosity = None
    else:
        print(f"  ⚠ Cannot calculate luminosity without flux density")
        luminosity = None
    
    # Classify radio source
    print(f"\n4. Classifying Radio Source:")
    print("="*80)
    
    if radio_properties and luminosity:
        classification = classify_radio_source(
            radio_properties.get('flux_1_4_ghz_mjy'),
            radio_properties.get('spectral_index'),
            luminosity
        )
        
        print(f"\n  Classification:")
        print(f"    Type: {classification['type']}")
        print(f"    Likely AGN: {classification['likely_agn']}")
        print(f"    Likely Star Formation: {classification['likely_star_formation']}")
        print(f"    Likely Synchrotron: {classification['likely_synchrotron']}")
        print(f"    Confidence: {classification['confidence']}")
    else:
        print(f"  ⚠ Cannot classify without radio properties")
        classification = None
    
    # Analysis
    print(f"\n5. ANALYSIS:")
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
        
        # Compare with typical values
        if luminosity['log_luminosity'] > 23:
            print(f"     • ⚠ Very high luminosity (typical of AGN)")
        elif luminosity['log_luminosity'] > 22:
            print(f"     • ⚠ High luminosity (possible AGN)")
        elif luminosity['log_luminosity'] > 21:
            print(f"     • ✓ Moderate luminosity (typical of star formation)")
        else:
            print(f"     • ✓ Low luminosity (typical of star formation or HII regions)")
    else:
        print(f"     • Luminosity: Not available")
    
    print(f"\n  C. Source Classification:")
    if classification:
        print(f"     • Type: {classification['type']}")
        print(f"     • Confidence: {classification['confidence']}")
        
        if classification['likely_agn']:
            print(f"     • ⚠ Likely AGN emission from LEDA 1363602")
            print(f"     • Could be radio emission from active galactic nucleus")
        elif classification['likely_star_formation']:
            print(f"     • ✓ Likely star formation emission from LEDA 1363602")
            print(f"     • Could be radio emission from star-forming regions")
        else:
            print(f"     • ⚠ Classification uncertain")
    else:
        print(f"     • Classification: Not available")
    
    print(f"\n  D. Relationship to LEDA 1363602:")
    print(f"     • NVSS source is 3.71 arcmin from LEDA 1363602")
    print(f"     • LEDA 1363602 is 720 million light-years away")
    print(f"     • If associated, radio emission is from the galaxy")
    
    if classification and classification['likely_agn']:
        print(f"     • ⚠ AGN emission suggests active galactic nucleus")
        print(f"     • Could be supermassive black hole activity")
    elif classification and classification['likely_star_formation']:
        print(f"     • ✓ Star formation emission suggests active star formation")
        print(f"     • Could be radio emission from HII regions")
    
    # Summary
    print(f"\n6. SUMMARY")
    print("="*80)
    
    print(f"\nNVSS J045150+092332 Radio Emissions:")
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
    
    print(f"\nRelationship to LEDA 1363602:")
    print(f"  • NVSS source is 3.71 arcmin from LEDA 1363602")
    print(f"  • LEDA 1363602 is 720 million light-years away")
    print(f"  • If associated, radio emission is extragalactic")
    
    # Save results
    results = {
        'analysis_date': datetime.now().isoformat(),
        'nvss_source': NVSS_J045150_092332,
        'leda_1363602': LEDA_1363602,
        'nvss_catalog_query': nvss_result,
        'nvss_sources_found': sources,
        'ned_query': ned_result,
        'radio_properties': radio_properties,
        'radio_luminosity': luminosity,
        'classification': classification,
    }
    
    output_file = Path('3i_atlas_data/radio_emissions_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")
    print(f"\nResults saved to: {output_file}")
    
    return results

if __name__ == "__main__":
    analyze_radio_emissions()

