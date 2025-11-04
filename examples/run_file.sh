#!/bin/bash

#Make It Executable : chmod +x ../examples/run_file.sh
 
# Delete existing error analysis file if it exists
rm -r ../build/error_analysis.txt

# Build r_refinement_square before running
make rh_refinement_example -j 2
# make rh_adaptiveAdvectiondiffusion -j 15
# make rh_elasticity_example -j 15

# Path to the executable
EXECUTABLE="./bin/rh_refinement_example"
# EXECUTABLE="./bin/rh_adaptiveAdvectiondiffusion"
# EXECUTABLE="./bin/rh_elasticity_example"

# Parameters (tags) to run the executable with -r 1: GARU, 2: PUCA, 3: BULK, 4: PBULK
TAGS=(
    # "-l -1 -u 3 -f 12. -q 1 -c 1 -v 0 -e 0"
    # .. Poisson 2D ...
    # "-r 2 -u 4 -f  0.  -l 5 -a 0.0 -c 0 -e 0"
    # "-r 2 -u 4 -f  0.  -l 5 -a 0.5 -c 1 -e 0"
    # "-r 2 -u 4 -f  0.  -l 5 -a 0.7 -c 1 -e 0"
    # "-r 2 -u 4 -f  12. -l 5 -a 0.0 -c 0 -e 0"
    # "-r 2 -u 4 -f  12. -l 5 -a 0.5 -c 0 -e 0"
    # "-r 2 -u 4 -f  12. -l 2 -a 0.7 -c 1 -e 0"
    # 3 dimensions case
    # "-r 2 -u 0 -f  0.  -l 5 -a 0.0 -c 0 -e 0 -d "volumes/GshapedVolume.xml""
    # "-r 2 -u 0 -f  0.  -l 6 -a 0.5 -c 0 -e 0 -d "volumes/GshapedVolume.xml""
    # "-r 2 -u 0 -f  0.  -l 6 -a 0.7 -c 1 -e 0 -d "volumes/GshapedVolume.xml""
    # "-r 2 -u 0 -f  12. -l 4 -a 0.0 -c 0 -e 0 -p 0.5 -d "volumes/GshapedVolume.xml""
    "-r 2 -u 0 -f  12. -l 5 -a 0.5 -c 0 -e 0 -d "volumes/GshapedVolume.xml""
    "-r 2 -u 0 -f  12. -l 5 -a 0.7 -c 1 -e 0 -d "volumes/GshapedVolume.xml""
    #.. Advection diffusion ... -f 0. : without r-refinement
    # "-r 2 -u 2  -f   0. -l 6 -a 0.7 -c 1 -e 1"
    # "-r 2 -u 2  -f  12. -l 6 -a 0.7 -c 1 -e 1"
    # .. Elasticity 2D ...
    #"-r 2 -u 2 -f  0.  -l 6 -a 0.0 -e 0"
    #"-r 2 -u 2 -f  0.  -l 6 -a 0.7 -e 0"
    # "-r 2 -u 2 -f  12. -l 3 -a 0.0 -e 0"
    # "-r 2 -u 2 -f  12. -l 6 -a 0.7  -e 0"
)

# Run the executable with each set of parameters
for TAG in "${TAGS[@]}"; do
    echo "Running $EXECUTABLE with parameters: $TAG"
    $EXECUTABLE --errorsave  $TAG
    echo "-------------------------------------------------"
done