#!/usr/bin/env python3
"""
Download and Process 3I/ATLAS Data
Downloads data from official sources and analyzes for Diophantine patterns
"""

import os
import json
import numpy as np
import requests
from pathlib import Path
from datetime import datetime
import zipfile
import tarfile

# Data sources from https://3i-atlas.github.io/data.html
DATA_SOURCES = {
    'tess_photometry': {
        'name': 'TESS Photometry',
        'format': '.np',
        'source': 'Zenodo',
        'url': 'https://zenodo.org/records/16943786',
        'reference': 'Feinstein, Noonan, & Seligman (2025)',
        'file_pattern': '*.np'
    },
    'soar_spectra': {
        'name': 'SOAR/GOODMAN Reflectance Spectra',
        'format': '.dat, .txt',
        'source': 'Zenodo',
        'url': 'https://zenodo.org/records/15881487',
        'reference': 'Puzia et al. (2025)',
        'file_pattern': '*.dat, *.txt'
    },
    'palomar_apo': {
        'name': 'Palomar P200 and APO Imaging and Spectroscopy',
        'format': '.fits',
        'source': 'CalTech Data Records',
        'url': 'https://data.caltech.edu/records/qdce4-pvm83',
        'reference': 'Bolin et al. under review; Belyakov et al. submitted',
        'file_pattern': '*.fits'
    },
    'atlas_spectrum': {
        'name': 'ATLAS Spectrum',
        'format': '.csv',
        'source': 'GitHub',
        'url': 'https://github.com/3I-ATLAS/discovery-paper/blob/main/src/data/3I-Reflectivity_JN_edit.csv',
        'reference': 'Seligman et al. (2025)',
        'file_pattern': '*.csv',
        'raw_url': 'https://raw.githubusercontent.com/3I-ATLAS/discovery-paper/main/src/data/3I-Reflectivity_JN_edit.csv'
    },
    'atlas_lightcurve': {
        'name': 'ATLAS Light Curve',
        'format': '.npy',
        'source': 'GitHub',
        'url': 'https://github.com/3I-ATLAS/discovery-paper/blob/main/src/data/ATLAS.npy',
        'reference': 'Seligman et al. (2025)',
        'file_pattern': '*.npy',
        'raw_url': 'https://raw.githubusercontent.com/3I-ATLAS/discovery-paper/main/src/data/ATLAS.npy'
    },
    'trappist_lightcurve': {
        'name': 'TRAPPIST Light Curve',
        'format': '.npy',
        'source': 'GitHub',
        'url': 'https://github.com/3I-ATLAS/discovery-paper/blob/main/src/data/TRAPPIST.npy',
        'reference': 'Seligman et al. (2025)',
        'file_pattern': '*.npy',
        'raw_url': 'https://raw.githubusercontent.com/3I-ATLAS/discovery-paper/main/src/data/TRAPPIST.npy'
    },
    'faulkes_lightcurve': {
        'name': 'Faulkes Telescope Light Curve',
        'format': '.npy',
        'source': 'GitHub',
        'url': 'https://github.com/3I-ATLAS/discovery-paper/blob/main/src/data/TJO.npy',
        'reference': 'Seligman et al. (2025)',
        'file_pattern': '*.npy',
        'raw_url': 'https://raw.githubusercontent.com/3I-ATLAS/discovery-paper/main/src/data/TJO.npy'
    },
    'lco_lightcurves': {
        'name': 'Las Cumbres Observatory Light Curves',
        'format': '.npy',
        'source': 'GitHub',
        'url': 'https://github.com/3I-ATLAS/discovery-paper/blob/main/src/data/LCO_1.npy',
        'reference': 'Seligman et al. (2025)',
        'file_pattern': '*.npy',
        'raw_urls': [
            'https://raw.githubusercontent.com/3I-ATLAS/discovery-paper/main/src/data/LCO_1.npy',
            'https://raw.githubusercontent.com/3I-ATLAS/discovery-paper/main/src/data/LCO_2.npy'
        ]
    },
    'atlas_volume_limits': {
        'name': 'ATLAS Volume Detection Limits',
        'format': '.txt',
        'source': 'GitHub',
        'url': 'https://github.com/3I-ATLAS/discovery-paper/blob/main/src/data/3I_xc.txt',
        'reference': 'Seligman et al. (2025)',
        'file_pattern': '*.txt',
        'raw_url': 'https://raw.githubusercontent.com/3I-ATLAS/discovery-paper/main/src/data/3I_xc.txt'
    }
}

def download_file(url, output_path, chunk_size=8192):
    """Download a file from URL"""
    print(f"  Downloading: {url}")
    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"    Progress: {percent:.1f}%", end='\r')
        
        print(f"    Downloaded: {output_path} ({downloaded:,} bytes)")
        return True
    except Exception as e:
        print(f"    Error downloading {url}: {e}")
        return False

def download_github_file(raw_url, output_path):
    """Download file from GitHub raw URL"""
    return download_file(raw_url, output_path)

def download_zenodo_record(record_id, output_dir):
    """Download data from Zenodo record"""
    print(f"  Downloading Zenodo record: {record_id}")
    # Note: Zenodo requires API access or manual download
    # For now, provide instructions
    print(f"    Please download manually from: https://zenodo.org/records/{record_id}")
    print(f"    Save files to: {output_dir}")
    return False

def download_caltech_record(record_url, output_dir):
    """Download data from CalTech Data Records"""
    print(f"  Downloading CalTech record: {record_url}")
    # Note: CalTech records may require authentication
    print(f"    Please download manually from: {record_url}")
    print(f"    Save files to: {output_dir}")
    return False

def download_all_data(output_dir='3i_atlas_data'):
    """Download all 3I/ATLAS data from official sources"""
    print("="*80)
    print("DOWNLOADING 3I/ATLAS DATA")
    print("="*80)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    downloaded_files = {}
    
    for key, source in DATA_SOURCES.items():
        print(f"\n{source['name']} ({source['format']}):")
        print(f"  Source: {source['source']}")
        print(f"  Reference: {source['reference']}")
        
        if source['source'] == 'GitHub':
            # Download from GitHub raw URLs
            if 'raw_url' in source:
                filename = Path(source['raw_url']).name
                output_file = output_path / filename
                if download_github_file(source['raw_url'], output_file):
                    downloaded_files[key] = str(output_file)
            elif 'raw_urls' in source:
                files = []
                for i, raw_url in enumerate(source['raw_urls']):
                    filename = Path(raw_url).name
                    output_file = output_path / filename
                    if download_github_file(raw_url, output_file):
                        files.append(str(output_file))
                if files:
                    downloaded_files[key] = files
        elif source['source'] == 'Zenodo':
            # Zenodo requires manual download or API access
            record_id = source['url'].split('/')[-1]
            zenodo_dir = output_path / f"zenodo_{record_id}"
            zenodo_dir.mkdir(exist_ok=True)
            download_zenodo_record(record_id, zenodo_dir)
        elif source['source'] == 'CalTech Data Records':
            # CalTech may require authentication
            caltech_dir = output_path / "caltech_data"
            caltech_dir.mkdir(exist_ok=True)
            download_caltech_record(source['url'], caltech_dir)
    
    # Save download manifest
    manifest = {
        'download_date': datetime.now().isoformat(),
        'sources': DATA_SOURCES,
        'downloaded_files': downloaded_files
    }
    
    manifest_file = output_path / 'download_manifest.json'
    with open(manifest_file, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    print("\n" + "="*80)
    print("DOWNLOAD SUMMARY")
    print("="*80)
    print(f"Downloaded files: {len(downloaded_files)}")
    print(f"Manifest saved to: {manifest_file}")
    print(f"\nNote: Some sources (Zenodo, CalTech) require manual download")
    print(f"Please download those files and place them in: {output_path}")
    
    return downloaded_files, manifest

def load_npy_file(filepath):
    """Load NumPy array file"""
    try:
        data = np.load(filepath)
        return data
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def load_csv_file(filepath):
    """Load CSV file"""
    try:
        data = np.genfromtxt(filepath, delimiter=',', skip_header=1)
        return data
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def load_txt_file(filepath):
    """Load text file"""
    try:
        data = np.loadtxt(filepath)
        return data
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def analyze_data_for_constants(data, data_name):
    """Analyze data for mathematical constants"""
    if data is None:
        return None
    
    results = {
        'data_name': data_name,
        'shape': data.shape if hasattr(data, 'shape') else str(type(data)),
        'dtype': str(data.dtype) if hasattr(data, 'dtype') else str(type(data)),
    }
    
    # Flatten data for analysis
    if hasattr(data, 'flatten'):
        flat_data = data.flatten()
    else:
        flat_data = np.array(data).flatten()
    
    # Remove NaN and infinite values
    flat_data = flat_data[np.isfinite(flat_data)]
    
    if len(flat_data) == 0:
        return results
    
    # Calculate statistics
    results['statistics'] = {
        'count': len(flat_data),
        'mean': float(np.mean(flat_data)),
        'median': float(np.median(flat_data)),
        'std': float(np.std(flat_data)),
        'min': float(np.min(flat_data)),
        'max': float(np.max(flat_data)),
    }
    
    # Check for mathematical constants
    constants = {
        'PHI': (1 + np.sqrt(5)) / 2,
        'PI': np.pi,
        'E': np.e,
        'SQRT2': np.sqrt(2),
        'SQRT3': np.sqrt(3),
        'SQRT5': np.sqrt(5),
        'LOG2': np.log(2),
        'LOG3': np.log(3),
    }
    
    constant_matches = {}
    for name, const_val in constants.items():
        # Check if constant appears in data
        matches = np.abs(flat_data - const_val) < 0.01
        match_count = np.sum(matches)
        if match_count > 0:
            constant_matches[name] = {
                'constant': float(const_val),
                'matches': int(match_count),
                'match_values': flat_data[matches][:10].tolist()  # First 10 matches
            }
    
    results['constant_matches'] = constant_matches
    
    return results

def process_downloaded_data(downloaded_files, output_dir='3i_atlas_data'):
    """Process downloaded data files"""
    print("\n" + "="*80)
    print("PROCESSING DOWNLOADED DATA")
    print("="*80)
    
    output_path = Path(output_dir)
    analysis_results = {}
    
    for key, filepath in downloaded_files.items():
        if isinstance(filepath, list):
            # Multiple files
            for fp in filepath:
                process_file(fp, key, analysis_results)
        else:
            process_file(filepath, key, analysis_results)
    
    # Save analysis results
    results_file = output_path / 'data_analysis.json'
    with open(results_file, 'w') as f:
        json.dump(analysis_results, f, indent=2)
    
    print(f"\nAnalysis results saved to: {results_file}")
    return analysis_results

def process_file(filepath, key, analysis_results):
    """Process a single file"""
    filepath = Path(filepath)
    if not filepath.exists():
        print(f"  File not found: {filepath}")
        return
    
    print(f"\nProcessing: {filepath.name}")
    
    # Load based on file extension
    if filepath.suffix == '.npy':
        data = load_npy_file(filepath)
    elif filepath.suffix == '.csv':
        data = load_csv_file(filepath)
    elif filepath.suffix == '.txt':
        data = load_txt_file(filepath)
    else:
        print(f"  Unsupported format: {filepath.suffix}")
        return
    
    # Analyze for constants
    analysis = analyze_data_for_constants(data, filepath.name)
    if analysis:
        analysis_results[key] = analysis
        print(f"  Shape: {analysis.get('shape', 'N/A')}")
        if 'statistics' in analysis:
            stats = analysis['statistics']
            print(f"  Mean: {stats['mean']:.6f}, Std: {stats['std']:.6f}")
        if 'constant_matches' in analysis and analysis['constant_matches']:
            print(f"  Constant matches: {list(analysis['constant_matches'].keys())}")

def main():
    """Main function"""
    print("="*80)
    print("3I/ATLAS DATA DOWNLOAD AND ANALYSIS")
    print("="*80)
    print("\nData sources from: https://3i-atlas.github.io/data.html")
    
    # Download data
    downloaded_files, manifest = download_all_data()
    
    # Process downloaded files
    if downloaded_files:
        analysis_results = process_downloaded_data(downloaded_files)
        
        print("\n" + "="*80)
        print("ANALYSIS SUMMARY")
        print("="*80)
        for key, analysis in analysis_results.items():
            if analysis and 'constant_matches' in analysis:
                matches = analysis['constant_matches']
                if matches:
                    print(f"\n{key}:")
                    for const_name, match_data in matches.items():
                        print(f"  {const_name}: {match_data['matches']} matches")
    
    print("\n" + "="*80)
    print("NEXT STEPS")
    print("="*80)
    print("1. Manually download Zenodo and CalTech data")
    print("2. Place files in 3i_atlas_data/ directory")
    print("3. Run analysis script to process all data")
    print("4. Search for Diophantine patterns and mathematical constants")

if __name__ == "__main__":
    main()

