#!/usr/bin/env python3
"""
Comprehensive Analysis of HD 286941
- Verify distance (Gaia parallax)
- Check for exoplanets
- Analyze stellar activity
- Check for binary system
- Cross-reference with 3I/ATLAS trajectory
"""

import json
import requests
import time
from pathlib import Path
from datetime import datetime
from urllib.parse import urlencode

def query_gaia_parallax(ra_deg, dec_deg, radius_arcsec=60):
    """
    Query Gaia Archive for parallax and distance data
    """
    print("="*80)
    print("QUERYING GAIA FOR HD 286941 PARALLAX AND DISTANCE")
    print("="*80)
    
    # Gaia Archive ADQL query
    query = f"""
    SELECT TOP 10
        source_id, ra, dec, parallax, parallax_error,
        pmra, pmdec, phot_g_mean_mag, bp_rp,
        distance_gspphot, distance_gspphot_lower, distance_gspphot_upper,
        ruwe, radial_velocity, radial_velocity_error
    FROM gaiadr3.gaia_source
    WHERE 1=CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {ra_deg:.6f}, {dec_deg:.6f}, {radius_arcsec/3600:.6f})
    )
    AND parallax > 0
    ORDER BY SQRT(POWER(ra - {ra_deg:.6f}, 2) + POWER(dec - {dec_deg:.6f}, 2))
    """
    
    url = "https://gea.esac.esa.int/tap-server/tap/sync"
    params = {
        'REQUEST': 'doQuery',
        'LANG': 'ADQL',
        'FORMAT': 'json',
        'QUERY': query,
    }
    
    try:
        print(f"\nQuerying Gaia Archive...")
        print(f"Coordinate: RA {ra_deg:.6f}°, Dec {dec_deg:.6f}°")
        print(f"Search radius: {radius_arcsec} arcsec")
        
        response = requests.get(url, params=params, timeout=20)
        
        if response.status_code == 200:
            data = response.json()
            
            if 'data' in data and len(data['data']) > 0:
                print(f"\n✓ Found {len(data['data'])} objects in Gaia")
                
                # Find closest match
                closest = data['data'][0]
                
                print(f"\nClosest match (likely HD 286941):")
                print(f"  Gaia Source ID: {closest.get('source_id', 'N/A')}")
                print(f"  RA: {closest.get('ra', 'N/A'):.6f}°")
                print(f"  Dec: {closest.get('dec', 'N/A'):.6f}°")
                
                # Parallax and distance
                parallax = closest.get('parallax')
                parallax_error = closest.get('parallax_error')
                
                if parallax:
                    print(f"\nParallax:")
                    print(f"  Parallax: {parallax:.3f} ± {parallax_error:.3f} mas")
                    
                    # Calculate distance from parallax
                    if parallax > 0:
                        distance_pc = 1000.0 / parallax  # parsecs
                        distance_ly = distance_pc * 3.26156  # light-years
                        
                        print(f"\nDistance (from parallax):")
                        print(f"  Distance: {distance_pc:.2f} parsecs")
                        print(f"  Distance: {distance_ly:.2f} light-years")
                
                # Gaia distance estimate
                distance_gspphot = closest.get('distance_gspphot')
                if distance_gspphot:
                    distance_gspphot_ly = distance_gspphot * 3.26156
                    print(f"\nDistance (Gaia GSP-Phot):")
                    print(f"  Distance: {distance_gspphot:.2f} parsecs")
                    print(f"  Distance: {distance_gspphot_ly:.2f} light-years")
                
                # Proper motion
                pmra = closest.get('pmra')
                pmdec = closest.get('pmdec')
                if pmra and pmdec:
                    print(f"\nProper Motion:")
                    print(f"  PM RA: {pmra:.3f} mas/yr")
                    print(f"  PM Dec: {pmdec:.3f} mas/yr")
                
                # Radial velocity
                rv = closest.get('radial_velocity')
                rv_error = closest.get('radial_velocity_error')
                if rv:
                    print(f"\nRadial Velocity:")
                    print(f"  RV: {rv:.2f} ± {rv_error:.2f} km/s")
                
                # Photometry
                g_mag = closest.get('phot_g_mean_mag')
                bp_rp = closest.get('bp_rp')
                if g_mag:
                    print(f"\nPhotometry:")
                    print(f"  G magnitude: {g_mag:.2f}")
                    if bp_rp:
                        print(f"  BP-RP color: {bp_rp:.2f}")
                
                # RUWE (Renormalised Unit Weight Error)
                ruwe = closest.get('ruwe')
                if ruwe:
                    print(f"\nRUWE: {ruwe:.3f}")
                    if ruwe > 1.4:
                        print(f"  ⚠ High RUWE - may indicate binary or unresolved companion")
                
                return {
                    'success': True,
                    'data': closest,
                    'parallax': parallax,
                    'parallax_error': parallax_error,
                    'distance_pc': 1000.0 / parallax if parallax and parallax > 0 else None,
                    'distance_ly': (1000.0 / parallax * 3.26156) if parallax and parallax > 0 else None,
                    'distance_gspphot_ly': distance_gspphot * 3.26156 if distance_gspphot else None,
                }
            else:
                print("\n✗ No objects found in Gaia Archive")
                return {'success': False, 'error': 'No objects found'}
        else:
            print(f"\n✗ Gaia query failed: HTTP {response.status_code}")
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        print(f"\n✗ Gaia query error: {e}")
        return {'success': False, 'error': str(e)}

def query_exoplanet_archive(ra_deg, dec_deg, radius_arcsec=60):
    """
    Query NASA Exoplanet Archive for exoplanets
    """
    print(f"\n{'='*80}")
    print("QUERYING NASA EXOPLANET ARCHIVE FOR HD 286941")
    print("="*80)
    
    # NASA Exoplanet Archive doesn't have direct coordinate search
    # We'll check by HD number
    print("\nNote: NASA Exoplanet Archive requires star name (HD number)")
    print("Checking for HD 286941...")
    
    # Try to query by HD number
    url = "https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI"
    params = {
        'table': 'exoplanets',
        'format': 'json',
        'where': 'hostname like "%HD 286941%"',
    }
    
    try:
        response = requests.get(url, params=params, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if len(data) > 0:
                print(f"\n✓ Found {len(data)} exoplanet(s) for HD 286941")
                for planet in data:
                    print(f"\n  Planet: {planet.get('pl_name', 'N/A')}")
                    print(f"    Period: {planet.get('pl_orbper', 'N/A')} days")
                    print(f"    Mass: {planet.get('pl_bmassj', 'N/A')} Jupiter masses")
                    print(f"    Radius: {planet.get('pl_radj', 'N/A')} Jupiter radii")
                return {'success': True, 'exoplanets': data}
            else:
                print("\n✗ No exoplanets found for HD 286941")
                return {'success': True, 'exoplanets': []}
        else:
            print(f"\n✗ Exoplanet Archive query failed: HTTP {response.status_code}")
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        print(f"\n✗ Exoplanet Archive query error: {e}")
        return {'success': False, 'error': str(e)}

def query_simbad_detailed(identifier='HD 286941'):
    """
    Query SIMBAD for detailed information about HD 286941
    """
    print(f"\n{'='*80}")
    print(f"QUERYING SIMBAD FOR DETAILED HD 286941 INFORMATION")
    print("="*80)
    
    url = "https://simbad.cds.unistra.fr/simbad/sim-id"
    params = {
        'Ident': identifier,
        'output.format': 'ASCII',
    }
    
    try:
        print(f"\nQuerying SIMBAD for {identifier}...")
        response = requests.get(url, params=params, timeout=10)
        
        if response.status_code == 200:
            result = response.text
            
            # Parse key information
            print("\n✓ SIMBAD query successful")
            print("\nKey Information:")
            
            # Extract spectral type
            if 'Spectral type:' in result:
                spectral_line = [l for l in result.split('\n') if 'Spectral type:' in l]
                if spectral_line:
                    print(f"  {spectral_line[0].strip()}")
            
            # Extract magnitude
            if 'V=' in result:
                mag_line = [l for l in result.split('\n') if 'V=' in l and 'mag' in l]
                if mag_line:
                    print(f"  {mag_line[0].strip()}")
            
            # Extract coordinates
            if 'Coordinates' in result:
                coord_line = [l for l in result.split('\n') if 'Coordinates' in l]
                if coord_line:
                    print(f"  {coord_line[0].strip()}")
            
            # Check for binary
            if 'Binary' in result or 'binary' in result:
                print("\n  ⚠ Binary system detected!")
            
            # Check for emission
            if 'Em*' in result or 'emission' in result.lower():
                print("\n  ⚠ Emission star - active stellar activity")
            
            return {'success': True, 'data': result}
        else:
            print(f"\n✗ SIMBAD query failed: HTTP {response.status_code}")
            return {'success': False, 'error': f'HTTP {response.status_code}'}
    except Exception as e:
        print(f"\n✗ SIMBAD query error: {e}")
        return {'success': False, 'error': str(e)}

def cross_reference_3i_atlas_trajectory(ra_deg, dec_deg):
    """
    Cross-reference HD 286941 with 3I/ATLAS trajectory
    """
    print(f"\n{'='*80}")
    print("CROSS-REFERENCING HD 286941 WITH 3I/ATLAS TRAJECTORY")
    print("="*80)
    
    # Load 3I/ATLAS trajectory data if available
    trajectory_file = Path('3i_atlas_data/trajectory_data.json')
    
    if trajectory_file.exists():
        with open(trajectory_file) as f:
            trajectory = json.load(f)
        
        print("\n3I/ATLAS Trajectory Data Found")
        print("Analyzing proximity to HD 286941...")
        # Add trajectory analysis here
    else:
        print("\n3I/ATLAS Trajectory Data Not Available")
        print("Note: Need 3I/ATLAS trajectory data to cross-reference")
        print(f"HD 286941 location: RA {ra_deg:.6f}°, Dec {dec_deg:.6f}°")
    
    return {'note': 'Trajectory cross-reference requires 3I/ATLAS trajectory data'}

def analyze_hd286941():
    """Comprehensive analysis of HD 286941"""
    print("="*80)
    print("COMPREHENSIVE ANALYSIS OF HD 286941")
    print("="*80)
    
    # HD 286941 coordinates (from SIMBAD lookup)
    ra_deg = 69.797974  # Approximate from decoded coordinate
    dec_deg = 11.250000
    
    # More precise coordinates from SIMBAD result
    # RA 04h 39m 17.537s = 69.823071°
    # Dec +11° 15' 54.264" = 11.265073°
    ra_precise = 69.823071
    dec_precise = 11.265073
    
    results = {
        'analysis_date': datetime.now().isoformat(),
        'object': 'HD 286941',
        'coordinates': {
            'ra_deg': ra_precise,
            'dec_deg': dec_precise,
        },
        'gaia_analysis': None,
        'exoplanet_analysis': None,
        'simbad_analysis': None,
        'trajectory_analysis': None,
    }
    
    # 1. Query Gaia for parallax and distance
    gaia_result = query_gaia_parallax(ra_precise, dec_precise, radius_arcsec=60)
    results['gaia_analysis'] = gaia_result
    
    # 2. Query exoplanet archive
    exoplanet_result = query_exoplanet_archive(ra_precise, dec_precise, radius_arcsec=60)
    results['exoplanet_analysis'] = exoplanet_result
    
    # 3. Query SIMBAD for detailed info
    simbad_result = query_simbad_detailed('HD 286941')
    results['simbad_analysis'] = simbad_result
    
    # 4. Cross-reference with 3I/ATLAS trajectory
    trajectory_result = cross_reference_3i_atlas_trajectory(ra_precise, dec_precise)
    results['trajectory_analysis'] = trajectory_result
    
    # Save results
    output_file = Path('3i_atlas_data/hd286941_analysis.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_file}")
    
    # Summary
    print(f"\n{'='*80}")
    print("ANALYSIS SUMMARY")
    print("="*80)
    
    if gaia_result.get('success'):
        distance_ly = gaia_result.get('distance_ly')
        if distance_ly:
            print(f"\nHD 286941 Distance:")
            print(f"  Distance: {distance_ly:.2f} light-years (from Gaia catalog)")
    
    if exoplanet_result.get('success'):
        exoplanets = exoplanet_result.get('exoplanets', [])
        if exoplanets:
            print(f"\n✓ Found {len(exoplanets)} exoplanet(s) for HD 286941")
        else:
            print("\n✗ No exoplanets found for HD 286941")
    
    print(f"\n{'='*80}")
    
    return results

if __name__ == "__main__":
    analyze_hd286941()

