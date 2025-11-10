#!/usr/bin/env python3
"""
Generate visualizations for Diophantine v3 coordinates pointing to LEDA 1363602 galaxy
Shows extragalactic coordinate encoding and dual-scale analysis
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patches as mpatches

def load_data():
    """Load Diophantine v3 and dual-scale encoding data"""
    data = {}
    
    # Load dual-scale encoding analysis
    dual_scale_file = Path('results/dual_scale_encoding_analysis.json')
    if dual_scale_file.exists():
        with open(dual_scale_file) as f:
            data['dual_scale'] = json.load(f)
    
    # Load Diophantine v3 objects
    diophantine_file = Path('results/diophantine_v3_objects.json')
    if diophantine_file.exists():
        with open(diophantine_file) as f:
            data['diophantine'] = json.load(f)
    
    # Load decoded coordinates for comparison
    decoded_file = Path('results/decoded_coordinates.json')
    if decoded_file.exists():
        with open(decoded_file) as f:
            data['decoded'] = json.load(f)
    
    return data

def create_extragalactic_sky_map(data):
    """Create sky map showing Diophantine v3 coordinates and LEDA 1363602"""
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(111, projection='aitoff')
    
    # Diophantine v3 coordinates
    diophantine_ra = 72.888794
    diophantine_dec = 9.364471
    
    # LEDA 1363602 galaxy
    leda_ra = 72.939583
    leda_dec = 9.334444
    
    # Convert RA to radians (for aitoff projection)
    def ra_to_rad(ra_deg):
        return np.radians(ra_deg - 180)
    
    # Plot Diophantine v3 coordinate
    ax.scatter(ra_to_rad(diophantine_ra), np.radians(diophantine_dec), 
              c='purple', s=200, marker='D', label='Diophantine v3 Coordinate', 
              edgecolors='darkviolet', linewidths=2, alpha=0.8, zorder=10)
    ax.text(ra_to_rad(diophantine_ra), np.radians(diophantine_dec), 
           '  Diophantine v3', fontsize=9, color='purple', weight='bold', ha='left', zorder=11)
    
    # Plot LEDA 1363602 galaxy
    ax.scatter(ra_to_rad(leda_ra), np.radians(leda_dec), 
              c='red', s=300, marker='*', label='LEDA 1363602 (Galaxy)', 
              edgecolors='darkred', linewidths=2, alpha=0.9, zorder=11)
    ax.text(ra_to_rad(leda_ra), np.radians(leda_dec), 
           '  LEDA 1363602', fontsize=9, color='red', weight='bold', ha='left', zorder=12)
    
    # Draw line connecting Diophantine v3 to LEDA 1363602
    ax.plot([ra_to_rad(diophantine_ra), ra_to_rad(leda_ra)], 
           [np.radians(diophantine_dec), np.radians(leda_dec)],
           'g--', linewidth=2, alpha=0.5, label='3.51 arcmin separation', zorder=5)
    
    # Grid and labels
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.set_xlabel('Right Ascension', fontsize=12, labelpad=20)
    ax.set_ylabel('Declination', fontsize=12, labelpad=20)
    ax.set_title('Diophantine v3 Extragalactic Coordinates - Sky Map\n(Points to LEDA 1363602 Galaxy)', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Legend
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    
    plt.tight_layout()
    return fig

def create_dual_scale_comparison(data):
    """Create comparison plot showing dual-scale encoding"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Diophantine v3 Extragalactic Encoding - Dual-Scale Analysis', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    dual_scale = data.get('dual_scale', {})
    encoded_objects = dual_scale.get('encoded_objects', {})
    
    # Extract data
    phi_log3 = encoded_objects.get('phi_log3', {})
    diophantine_v3 = encoded_objects.get('diophantine_v3', {})
    
    phi_obj = phi_log3.get('objects', [{}])[0] if phi_log3.get('objects') else {}
    dioph_obj = diophantine_v3.get('objects', [{}])[0] if diophantine_v3.get('objects') else {}
    
    # Plot 1: Coordinate positions
    ax1 = axes[0, 0]
    
    # Diophantine v3 and LEDA 1363602
    diophantine_ra = 72.888794
    diophantine_dec = 9.364471
    leda_ra = 72.939583
    leda_dec = 9.334444
    
    ax1.scatter(diophantine_ra, diophantine_dec, c='purple', s=200, marker='D', 
               label='Diophantine v3', edgecolors='darkviolet', linewidths=2, alpha=0.8)
    ax1.scatter(leda_ra, leda_dec, c='red', s=300, marker='*', 
               label='LEDA 1363602 (Galaxy)', edgecolors='darkred', linewidths=2, alpha=0.9)
    
    # Draw connection line
    ax1.plot([diophantine_ra, leda_ra], [diophantine_dec, leda_dec], 
            'g--', linewidth=2, alpha=0.5, label='3.51 arcmin')
    
    ax1.set_xlabel('Right Ascension (degrees)', fontsize=11)
    ax1.set_ylabel('Declination (degrees)', fontsize=11)
    ax1.set_title('Coordinate Positions', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)
    
    # Plot 2: Distance comparison (log scale)
    ax2 = axes[0, 1]
    
    if phi_obj.get('distance_ly') and dioph_obj.get('distance_ly'):
        distances = [phi_obj['distance_ly'], dioph_obj['distance_ly']]
        labels = ['Galactic\n(HIP 21684)', 'Extragalactic\n(LEDA 1363602)']
        colors = ['blue', 'red']
        
        bars = ax2.bar(labels, distances, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax2.set_yscale('log')
        ax2.set_ylabel('Distance (light-years, log scale)', fontsize=11)
        ax2.set_title('Distance Comparison', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, dist in zip(bars, distances):
            height = bar.get_height()
            if dist < 1e6:
                label = f'{dist:.1f} ly'
            else:
                label = f'{dist/1e6:.0f}M ly'
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    label, ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Plot 3: Radio luminosity comparison
    ax3 = axes[1, 0]
    
    if phi_obj.get('radio_source') and dioph_obj.get('radio_source'):
        galactic_lum = phi_obj['radio_source'].get('log_luminosity', 0)
        extragalactic_lum = dioph_obj['radio_source'].get('log_luminosity', 0)
        
        luminosities = [galactic_lum, extragalactic_lum]
        labels = ['Galactic\n(Star Formation)', 'Extragalactic\n(AGN)']
        colors = ['blue', 'red']
        
        bars = ax3.bar(labels, luminosities, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax3.set_ylabel('Log Radio Luminosity (erg/s/Hz)', fontsize=11)
        ax3.set_title('Radio Luminosity Comparison', fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, lum in zip(bars, luminosities):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height,
                    f'{lum:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Plot 4: Encoding entropy comparison
    ax4 = axes[1, 1]
    
    if phi_log3.get('entropy') and diophantine_v3.get('entropy'):
        entropies = [phi_log3['entropy'], diophantine_v3['entropy']]
        labels = ['Phi+Log3\n(Galactic)', 'Diophantine v3\n(Extragalactic)']
        colors = ['blue', 'purple']
        
        bars = ax4.bar(labels, entropies, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
        ax4.set_ylabel('Encoding Entropy', fontsize=11)
        ax4.set_title('Encoding Entropy Comparison', fontsize=12, fontweight='bold')
        ax4.set_ylim([0, 1.1])
        ax4.grid(True, alpha=0.3, axis='y')
        
        # Add value labels
        for bar, ent in zip(bars, entropies):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height,
                    f'{ent:.4f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    return fig

def create_objects_found_plot(data):
    """Create plot showing all objects found near Diophantine v3 coordinates"""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    diophantine = data.get('diophantine', {})
    objects = diophantine.get('objects_found', [])
    
    if not objects:
        ax.text(0.5, 0.5, 'No object data available', 
               ha='center', va='center', fontsize=14, transform=ax.transAxes)
        ax.set_title('Objects Found Near Diophantine v3 Coordinates', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        return fig
    
    # Extract data
    identifiers = []
    distances = []
    types = []
    
    for obj in objects:
        identifiers.append(obj.get('identifier', 'Unknown'))
        dist_arcsec = float(obj.get('distance_arcsec', 0))
        distances.append(dist_arcsec / 60.0)  # Convert to arcmin
        obj_type = obj.get('object_type', 'Unknown')
        types.append(obj_type)
    
    # Color coding by type
    color_map = {
        'G': 'red',      # Galaxy
        'Em*': 'blue',   # Emission star
        'Rad': 'orange', # Radio source
        '*': 'green',   # Star
        'Unknown': 'gray'
    }
    colors = [color_map.get(t, 'gray') for t in types]
    
    # Create bar plot
    bars = ax.barh(identifiers, distances, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    
    # Add value labels
    for bar, dist in zip(bars, distances):
        width = bar.get_width()
        ax.text(width, bar.get_y() + bar.get_height()/2.,
               f'{dist:.2f} arcmin', ha='left', va='center', fontsize=9, fontweight='bold')
    
    ax.set_xlabel('Distance from Diophantine v3 Coordinate (arcmin)', fontsize=11)
    ax.set_ylabel('Object Identifier', fontsize=11)
    ax.set_title('Objects Found Near Diophantine v3 Coordinates\n(LEDA 1363602 is closest galaxy)', 
                fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='x')
    
    # Add legend for object types
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color_map.get(t, 'gray'), label=t) 
                      for t in set(types)]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
    plt.tight_layout()
    return fig

def main():
    """Main function to generate all Diophantine v3 visualizations"""
    print("="*80)
    print("GENERATING DIOPHANTINE V3 EXTRAGALACTIC VISUALIZATIONS")
    print("="*80)
    
    # Load data
    print("\nLoading data...")
    data = load_data()
    print(f"  Loaded dual-scale encoding data: {bool(data.get('dual_scale'))}")
    print(f"  Loaded Diophantine v3 objects: {bool(data.get('diophantine'))}")
    print(f"  Loaded decoded coordinates: {bool(data.get('decoded'))}")
    
    # Create output directory
    output_dir = Path('visualizations')
    output_dir.mkdir(exist_ok=True)
    
    # Generate sky map
    print("\nGenerating extragalactic sky map...")
    fig1 = create_extragalactic_sky_map(data)
    sky_map_path = output_dir / 'diophantine_v3_extragalactic_sky_map.png'
    fig1.savefig(sky_map_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {sky_map_path}")
    plt.close(fig1)
    
    # Generate dual-scale comparison
    print("\nGenerating dual-scale comparison plot...")
    fig2 = create_dual_scale_comparison(data)
    comparison_path = output_dir / 'diophantine_v3_dual_scale_comparison.png'
    fig2.savefig(comparison_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {comparison_path}")
    plt.close(fig2)
    
    # Generate objects found plot
    print("\nGenerating objects found plot...")
    fig3 = create_objects_found_plot(data)
    objects_path = output_dir / 'diophantine_v3_objects_found.png'
    fig3.savefig(objects_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {objects_path}")
    plt.close(fig3)
    
    print("\n" + "="*80)
    print("VISUALIZATION GENERATION COMPLETE")
    print("="*80)
    print(f"\nGenerated files:")
    print(f"  1. {sky_map_path}")
    print(f"  2. {comparison_path}")
    print(f"  3. {objects_path}")

if __name__ == "__main__":
    main()

