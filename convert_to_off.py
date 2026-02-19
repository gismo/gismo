#!/usr/bin/env python3
"""
Convert Val3Fct1.xml geometries to .off format.
Each geometry has 25 points in a 5x5 layout.
Extract the 4 corner points (indices 0, 4, 20, 24) from each patch.
"""

import xml.etree.ElementTree as ET
from collections import OrderedDict

def parse_xml(filename):
    """Parse the XML file and extract geometries."""
    tree = ET.parse(filename)
    root = tree.getroot()
    
    geometries = []
    for geom in root.findall('Geometry'):
        coefs_elem = geom.find('.//coefs')
        if coefs_elem is not None:
            # Parse the coefficients (3D points)
            coefs_text = coefs_elem.text.strip()
            values = [float(x) for x in coefs_text.split()]
            
            # Group into 3D points
            points = []
            for i in range(0, len(values), 3):
                points.append((values[i], values[i+1], values[i+2]))
            
            geometries.append(points)
    
    return geometries

def extract_corners(geometries):
    """Extract the 4 corner points from each 5x5 patch."""
    all_corners = []
    
    for geom in geometries:
        # In a 5x5 grid (indices 0-24), corners are at:
        # 0 (top-left), 4 (top-right), 20 (bottom-left), 24 (bottom-right)
        corners = [
            geom[0],   # (0, 0)
            geom[4],   # (0, 4)
            geom[20],  # (4, 0)
            geom[24]   # (4, 4)
        ]
        all_corners.append(corners)
    
    return all_corners

def build_mesh(patch_corners):
    """Build mesh by identifying unique vertices and creating faces."""
    # Collect all unique vertices
    vertex_map = OrderedDict()
    vertex_index = 0
    
    # For each patch, map its corners to vertex indices
    patch_indices = []
    
    for corners in patch_corners:
        indices = []
        for corner in corners:
            # Round to avoid floating point comparison issues
            key = tuple(round(x, 10) for x in corner)
            
            if key not in vertex_map:
                vertex_map[key] = vertex_index
                vertex_index += 1
            
            indices.append(vertex_map[key])
        
        patch_indices.append(indices)
    
    # Convert vertex map to list
    vertices = list(vertex_map.keys())
    
    return vertices, patch_indices

def write_off(filename, vertices, faces):
    """Write the mesh to an OFF file."""
    with open(filename, 'w') as f:
        f.write("OFF\n")
        f.write(f"{len(vertices)} {len(faces)} 0\n")
        
        # Write vertices
        for v in vertices:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")
        
        # Write faces (quads)
        for face in faces:
            f.write(f"4 {face[0]} {face[1]} {face[3]} {face[2]}\n")

def main():
    import sys
    
    if len(sys.argv) > 1:
        # Use command-line argument
        input_file = sys.argv[1]
        # Extract base name for output
        import os
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = f"{base_name}_corners.off"
    else:
        # Default
        input_file = "filedata/freeformSubdivision/fitting_functions/Val3Fct1.xml"
        output_file = "Val3Fct1_corners.off"
    
    print(f"Reading {input_file}...")
    geometries = parse_xml(input_file)
    print(f"Found {len(geometries)} geometries")
    
    print("Extracting corner points...")
    patch_corners = extract_corners(geometries)
    
    print("Building mesh and identifying shared vertices...")
    vertices, faces = build_mesh(patch_corners)
    
    print(f"Mesh has {len(vertices)} unique vertices and {len(faces)} faces")
    
    print(f"Writing to {output_file}...")
    write_off(output_file, vertices, faces)
    
    print("Done!")

if __name__ == "__main__":
    main()
