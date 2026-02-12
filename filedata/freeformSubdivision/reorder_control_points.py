#!/usr/bin/env python3
"""
Reorder control points in 5x5 Bezier patches so they radiate outward from shared center.
Each of the 4 patches shares a common center corner and needs different rotation.
"""

import xml.etree.ElementTree as ET
import sys

def rotate_180(points):
    """Rotate 180° - center at bottom-right (24) moves to top-left (0)"""
    if len(points) != 25:
        raise ValueError(f"Expected 25 points, got {len(points)}")
    # Simply reverse the list
    return list(reversed(points))

def rotate_90_ccw(points):
    """Rotate 90° counter-clockwise - center at bottom-left (20) moves to top-left (0)"""
    if len(points) != 25:
        raise ValueError(f"Expected 25 points, got {len(points)}")
    # Grid position (i,j) -> (j, 4-i) in rotated grid
    # Index conversion: old[i*5+j] -> new[j*5+(4-i)]
    result = [None] * 25
    for i in range(5):
        for j in range(5):
            old_idx = i * 5 + j
            new_i = j
            new_j = 4 - i
            new_idx = new_i * 5 + new_j
            result[new_idx] = points[old_idx]
    return result

def rotate_90_cw(points):
    """Rotate 90° clockwise - center at top-right (4) moves to top-left (0)"""
    if len(points) != 25:
        raise ValueError(f"Expected 25 points, got {len(points)}")
    # Grid position (i,j) -> (4-j, i) in rotated grid
    # Index conversion: old[i*5+j] -> new[(4-j)*5+i]
    result = [None] * 25
    for i in range(5):
        for j in range(5):
            old_idx = i * 5 + j
            new_i = 4 - j
            new_j = i
            new_idx = new_i * 5 + new_j
            result[new_idx] = points[old_idx]
    return result

def no_rotation(points):
    """No rotation - center already at top-left (0)"""
    return points

def process_file(input_file, output_file, rotation_func, rotation_name):
    """Process an XML file containing a 25x2 matrix of control points."""
    tree = ET.parse(input_file)
    root = tree.getroot()
    
    # Find the Matrix element
    matrix = root.find('Matrix')
    if matrix is None:
        raise ValueError("No Matrix element found in XML")
    
    # Get matrix dimensions
    rows = int(matrix.get('rows'))
    cols = int(matrix.get('cols'))
    
    if rows != 25:
        raise ValueError(f"Expected 25 rows, got {rows}")
    
    # Parse the matrix data
    lines = matrix.text.strip().split('\n')
    points = []
    for line in lines:
        line = line.strip()
        if line:
            points.append(line)
    
    if len(points) != 25:
        raise ValueError(f"Expected 25 data lines, got {len(points)}")
    
    # Reorder the points using specified rotation
    reordered_points = rotation_func(points)
    
    # Create new matrix text
    new_matrix_text = '\n    ' + '\n    '.join(reordered_points) + '\n  '
    matrix.text = new_matrix_text
    
    # Write output
    tree.write(output_file, encoding='UTF-8', xml_declaration=True)
    print(f"Processed: {input_file} ({rotation_name}) -> {output_file}")

if __name__ == '__main__':
    # Each file needs a different rotation to bring shared center to position 0
    # Patch 1: center at (4,4) bottom-right -> rotate 180°
    # Patch 2: center at (4,0) bottom-left -> rotate 90° CCW
    # Patch 3: center at (0,0) top-left -> no rotation
    # Patch 4: center at (0,4) top-right -> rotate 90° CW
    
    # Process only valence 3
    valences = [3]
    
    rotations = [
        (1, rotate_180, 'rotate 180°'),
        (2, rotate_90_ccw, 'rotate 90° CCW'),
        (3, no_rotation, 'no rotation'),
        (4, rotate_90_cw, 'rotate 90° CW'),
    ]
    
    base_path = '/home/linus/Development/gismo/filedata/freeformSubdivision/'
    
    for valence in valences:
        print(f"\n=== Processing valence {valence} ===")
        for patch_num, rotation_func, rotation_name in rotations:
            filename = f'control_net_d5_v{valence}_fine_{patch_num}.xml'
            input_path = base_path + filename
            output_path = base_path + filename
            try:
                process_file(input_path, output_path, rotation_func, rotation_name)
            except Exception as e:
                print(f"Error processing {filename}: {e}", file=sys.stderr)
                sys.exit(1)
    
    print("\n=== All files processed successfully! ===")
