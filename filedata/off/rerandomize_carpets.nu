#!/usr/bin/env nu

# Re-randomize Z coordinates in all carpet OFF files

def process_carpet_file [file: string] {
    print $"Processing ($file)..."
    
    let lines = open $file | lines
    
    let new_lines = $lines | enumerate | each { |it|
        let i = $it.index
        let line = $it.item
        
        if $i == 0 or $i == 1 {
            # Keep header lines as-is
            $line
        } else {
            let stripped = $line | str trim
            
            if ($stripped | is-empty) {
                $line
            } else {
                let parts = $stripped | split row ' '
                
                # Check if this is a vertex line (3 coordinates)
                if ($parts | length) == 3 {
                    try {
                        let x = $parts | get 0 | into float
                        let y = $parts | get 1 | into float
                        # Generate new random Z between 0 and 1
                        let z = (random float 0..1)
                        $"($x) ($y) ($z)"
                    } catch {
                        # Not a vertex line, keep as-is
                        $line
                    }
                } else {
                    # Face line or other, keep as-is
                    $line
                }
            }
        }
    }
    
    # Write back to file
    $new_lines | str join (char newline) | save -f $file
    print $"  Updated ($file)"
}

# Process all carpet files
let carpet_files = glob ev_valence*_carpet.off
print $"Found ($carpet_files | length) carpet files"

for file in $carpet_files {
    process_carpet_file $file
}

print "Done re-randomizing all carpet files!"
