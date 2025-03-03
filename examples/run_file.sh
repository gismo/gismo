#!/bin/bash

#Make It Executable : chmod +x ../examples/run_file.sh
 
# Delete existing error analysis file if it exists
rm -r ../build/error_analysis.txt

# Build r_refinement_square before running
make rh_refinement_example -j 15
#make r_refinement_ComplexGeometry -j 15

# Path to the executable
EXECUTABLE="./bin/rh_refinement_example"
#EXECUTABLE="./bin/r_refinement_ComplexGeometry"

# Parameters (tags) to run the executable with -r 1: GARU, 2: PUCA, 3: BULK, 4: PBULK
TAGS=(
    "-r 2 -u 4 -f 13. -l 2 -a 0.0 -c 0 -p 0 -e 1"
    #"-r 2 -u 4 -f 10. -l 5 -a 0.7 -c 3 -p 2 -e 0"
)

# Run the executable with each set of parameters
for TAG in "${TAGS[@]}"; do
    echo "Running $EXECUTABLE with parameters: $TAG"
    $EXECUTABLE --errorsave --plot $TAG
    echo "-------------------------------------------------"
done