#!/usr/bin/env python3
"""
Generate visualizations for 3I/ATLAS decoded coordinates
Creates sky map and comparison plots
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.patches as mpatches

def load_coordinates():
    """Load decoded coordinates from JSON file"""
    results_file = Path('results/decoded_coordinates.json')
    if not results_file.exists():
        raise FileNotFoundError(f"Results file not found: {results_file}")
    
    with open(results_file) as f:
        return json.load(f)

def ra_to_plot_angle(ra_deg):
    """Convert RA (degrees) to plot angle (0-360, increasing left to right)"""
    # RA increases from right to left in standard plots
    # Convert so 0h is at right, 12h is at left
    return (360 - ra_deg) % 360

def create_sky_map(data):
    """Create sky map visualization showing all decoded coordinates"""
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(111, projection='aitoff')
    
    # Extract all coordinates
    all_coords = data.get('all_coordinates', [])
    primary = data.get('primary_coordinate', {})
    convergence_groups = data.get('convergence_groups', [])
    
    # Collect coordinates by type
    structure_coords = []
    encoding_coords = []
    primary_coord = None
    
    for coord in all_coords:
        ra = coord.get('ra_deg', 0)
        dec = coord.get('dec_deg', 0)
        
        if 'structure_number' in coord:
            structure_coords.append((ra, dec, coord.get('structure_number', 0)))
        elif 'encoding_scheme' in coord:
            encoding_coords.append((ra, dec, coord.get('encoding_scheme', '')))
    
    if primary:
        primary_coord = (primary.get('ra_deg', 0), primary.get('dec_deg', 0))
    
    # Convert RA to radians (for aitoff projection, RA needs to be in range -180 to 180)
    def ra_to_rad(ra_deg):
        return np.radians(ra_deg - 180)  # Shift so 0h is at center
    
    # Plot structure coordinates
    if structure_coords:
        ra_vals = [ra_to_rad(ra) for ra, dec, _ in structure_coords]
        dec_vals = [np.radians(dec) for ra, dec, _ in structure_coords]
        ax.scatter(ra_vals, dec_vals, c='blue', s=100, marker='o', 
                  label='37-bit Structures', alpha=0.7, edgecolors='darkblue', linewidths=1.5)
        
        # Add labels for structure numbers
        for ra, dec, num in structure_coords:
            ax.text(ra_to_rad(ra), np.radians(dec), f'  S{num}', 
                   fontsize=8, color='blue', ha='left')
    
    # Plot encoding scheme coordinates
    if encoding_coords:
        ra_vals = [ra_to_rad(ra) for ra, dec, _ in encoding_coords]
        dec_vals = [np.radians(dec) for ra, dec, _ in encoding_coords]
        ax.scatter(ra_vals, dec_vals, c='red', s=100, marker='s', 
                  label='Encoding Schemes', alpha=0.7, edgecolors='darkred', linewidths=1.5)
        
        # Add labels for encoding schemes
        scheme_labels = {
            'atlas_spectrum_phi_log3': 'PHI+LOG3',
            'atlas_spectrum_phi_only': 'PHI',
            'atlas_spectrum_log3_only': 'LOG3'
        }
        for ra, dec, scheme in encoding_coords:
            label = scheme_labels.get(scheme, scheme.replace('_', ' ').upper())
            ax.text(ra_to_rad(ra), np.radians(dec), f'  {label}', 
                   fontsize=8, color='red', ha='left')
    
    # Plot primary coordinate
    if primary_coord:
        ra_rad = ra_to_rad(primary_coord[0])
        dec_rad = np.radians(primary_coord[1])
        ax.scatter(ra_rad, dec_rad, c='gold', s=300, marker='*', 
                  label='Primary Coordinate', edgecolors='orange', linewidths=2, zorder=10)
        ax.text(ra_rad, dec_rad, '  PRIMARY', fontsize=10, color='gold', 
               weight='bold', ha='left', zorder=11)
    
    # Plot convergence groups
    for i, group in enumerate(convergence_groups):
        if group.get('num_coordinates', 0) > 1:
            center_ra = group.get('center_ra', 0)
            center_dec = group.get('center_dec', 0)
            ra_rad = ra_to_rad(center_ra)
            dec_rad = np.radians(center_dec)
            
            # Draw circle around convergence region (simplified - just mark with larger point)
            ax.scatter(ra_rad, dec_rad, c='green', s=200, marker='o', 
                      edgecolors='darkgreen', linewidths=2, alpha=0.3, zorder=5)
    
    # Grid and labels
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.set_xlabel('Right Ascension', fontsize=12, labelpad=20)
    ax.set_ylabel('Declination', fontsize=12, labelpad=20)
    ax.set_title('3I/ATLAS Decoded Coordinates - Sky Map', fontsize=14, fontweight='bold', pad=20)
    
    # Legend
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    
    plt.tight_layout()
    return fig

def create_comparison_plot(data):
    """Create comparison plot showing coordinate convergence"""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('3I/ATLAS Decoded Coordinates - Comparison Analysis', 
                 fontsize=16, fontweight='bold', y=0.995)
    
    all_coords = data.get('all_coordinates', [])
    primary = data.get('primary_coordinate', {})
    convergence_groups = data.get('convergence_groups', [])
    
    # Extract coordinates by type
    structure_ras = []
    structure_decs = []
    structure_nums = []
    encoding_ras = []
    encoding_decs = []
    encoding_schemes = []
    
    for coord in all_coords:
        ra = coord.get('ra_deg', 0)
        dec = coord.get('dec_deg', 0)
        
        if 'structure_number' in coord:
            structure_ras.append(ra)
            structure_decs.append(dec)
            structure_nums.append(coord.get('structure_number', 0))
        elif 'encoding_scheme' in coord:
            encoding_ras.append(ra)
            encoding_decs.append(dec)
            encoding_schemes.append(coord.get('encoding_scheme', ''))
    
    # Plot 1: RA vs Dec scatter
    ax1 = axes[0, 0]
    if structure_ras:
        ax1.scatter(structure_ras, structure_decs, c='blue', s=150, marker='o', 
                   label='37-bit Structures', alpha=0.7, edgecolors='darkblue', linewidths=1.5)
        for i, num in enumerate(structure_nums):
            ax1.annotate(f'S{num}', (structure_ras[i], structure_decs[i]), 
                        fontsize=9, color='blue', xytext=(5, 5), textcoords='offset points')
    
    if encoding_ras:
        ax1.scatter(encoding_ras, encoding_decs, c='red', s=150, marker='s', 
                   label='Encoding Schemes', alpha=0.7, edgecolors='darkred', linewidths=1.5)
        scheme_labels = {
            'atlas_spectrum_phi_log3': 'PHI+LOG3',
            'atlas_spectrum_phi_only': 'PHI',
            'atlas_spectrum_log3_only': 'LOG3'
        }
        for i, scheme in enumerate(encoding_schemes):
            label = scheme_labels.get(scheme, scheme.replace('_', ' ').upper())
            ax1.annotate(label, (encoding_ras[i], encoding_decs[i]), 
                        fontsize=9, color='red', xytext=(5, 5), textcoords='offset points')
    
    if primary:
        ax1.scatter(primary.get('ra_deg', 0), primary.get('dec_deg', 0), 
                   c='gold', s=400, marker='*', label='Primary Coordinate', 
                   edgecolors='orange', linewidths=2, zorder=10)
        ax1.annotate('PRIMARY', (primary.get('ra_deg', 0), primary.get('dec_deg', 0)), 
                    fontsize=11, color='gold', weight='bold', 
                    xytext=(10, 10), textcoords='offset points')
    
    ax1.set_xlabel('Right Ascension (degrees)', fontsize=11)
    ax1.set_ylabel('Declination (degrees)', fontsize=11)
    ax1.set_title('Coordinate Distribution', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)
    
    # Plot 2: RA distribution
    ax2 = axes[0, 1]
    all_ras = structure_ras + encoding_ras
    if all_ras:
        ax2.hist(all_ras, bins=20, color='steelblue', alpha=0.7, edgecolor='black')
        if primary:
            ax2.axvline(primary.get('ra_deg', 0), color='gold', linestyle='--', 
                        linewidth=2, label='Primary RA')
        ax2.set_xlabel('Right Ascension (degrees)', fontsize=11)
        ax2.set_ylabel('Frequency', fontsize=11)
        ax2.set_title('RA Distribution', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
        if primary:
            ax2.legend(fontsize=9)
    
    # Plot 3: Dec distribution
    ax3 = axes[1, 0]
    all_decs = structure_decs + encoding_decs
    if all_decs:
        ax3.hist(all_decs, bins=20, color='coral', alpha=0.7, edgecolor='black')
        if primary:
            ax3.axvline(primary.get('dec_deg', 0), color='gold', linestyle='--', 
                        linewidth=2, label='Primary Dec')
        ax3.set_xlabel('Declination (degrees)', fontsize=11)
        ax3.set_ylabel('Frequency', fontsize=11)
        ax3.set_title('Declination Distribution', fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3, axis='y')
        if primary:
            ax3.legend(fontsize=9)
    
    # Plot 4: Convergence groups
    ax4 = axes[1, 1]
    if convergence_groups:
        group_sizes = [g.get('num_coordinates', 0) for g in convergence_groups]
        group_labels = [f"Group {i+1}" for i in range(len(convergence_groups))]
        colors = plt.cm.viridis(np.linspace(0, 1, len(convergence_groups)))
        
        bars = ax4.bar(range(len(convergence_groups)), group_sizes, color=colors, 
                      alpha=0.7, edgecolor='black')
        ax4.set_xlabel('Convergence Group', fontsize=11)
        ax4.set_ylabel('Number of Coordinates', fontsize=11)
        ax4.set_title('Convergence Group Sizes', fontsize=12, fontweight='bold')
        ax4.set_xticks(range(len(convergence_groups)))
        ax4.set_xticklabels(group_labels, rotation=45, ha='right')
        ax4.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for i, (bar, size) in enumerate(zip(bars, group_sizes)):
            height = bar.get_height()
            ax4.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(size)}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    return fig

def main():
    """Main function to generate all visualizations"""
    print("="*80)
    print("GENERATING 3I/ATLAS VISUALIZATIONS")
    print("="*80)
    
    # Load data
    print("\nLoading coordinate data...")
    data = load_coordinates()
    print(f"  Loaded {len(data.get('all_coordinates', []))} coordinates")
    
    # Create output directory
    output_dir = Path('visualizations')
    output_dir.mkdir(exist_ok=True)
    
    # Generate sky map
    print("\nGenerating sky map visualization...")
    fig1 = create_sky_map(data)
    sky_map_path = output_dir / 'decoded_coordinates_sky_map.png'
    fig1.savefig(sky_map_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {sky_map_path}")
    plt.close(fig1)
    
    # Generate comparison plot
    print("\nGenerating comparison plot...")
    fig2 = create_comparison_plot(data)
    comparison_path = output_dir / 'decoded_coordinates_comparison.png'
    fig2.savefig(comparison_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {comparison_path}")
    plt.close(fig2)
    
    print("\n" + "="*80)
    print("VISUALIZATION GENERATION COMPLETE")
    print("="*80)
    print(f"\nGenerated files:")
    print(f"  1. {sky_map_path}")
    print(f"  2. {comparison_path}")

if __name__ == "__main__":
    main()

