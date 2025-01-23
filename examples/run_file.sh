#!/bin/bash

#Make It Executable : chmod +x ../examples/run_r_refinement_square.sh
# Path to the executable
EXECUTABLE="./bin/r_refinement_square"
#EXECUTABLE="./bin/r_refinement_ComplexGeometry"

# Parameters (tags) to run the executable with
TAGS=(
    "-u 3 -f 6. -l 5 -a 0."
    "-u 3 -f 6. -l 5 -a 0.25"
    "-u 3 -f 6. -l 5 -a 0.5"
    "-u 3 -f 6. -l 6 -a 0.7"
    "-u 3 -f 6. -l 10 -a 0.9"
    "-u 3 -f 6. -l 5 -a 0. -p 0.2"
)

# Run the executable with each set of parameters
for TAG in "${TAGS[@]}"; do
    echo "Running $EXECUTABLE with parameters: $TAG"
    $EXECUTABLE $TAG
    echo "-------------------------------------------------"
done