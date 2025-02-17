#!/bin/bash

#Make It Executable : chmod +x ../examples/run_file.sh
 
# Delete existing error analysis file if it exists
rm -r ../build/error_analysis.txt

# Build r_refinement_square before running
#make r_refinement_square -j 15
make r_refinement_ComplexGeometry -j 15

# Path to the executable
#EXECUTABLE="./bin/r_refinement_square"
EXECUTABLE="./bin/r_refinement_ComplexGeometry"

# Parameters (tags) to run the executable with
TAGS=(
    "-u 4 -f 10. -l 5 -a 0.3 -c 1 -p 0"
    "-u 4 -f 10. -l 6 -a 0.7 -c 2 -p 0"
)

# Run the executable with each set of parameters
for TAG in "${TAGS[@]}"; do
    echo "Running $EXECUTABLE with parameters: $TAG"
    $EXECUTABLE --errorsave $TAG
    echo "-------------------------------------------------"
done