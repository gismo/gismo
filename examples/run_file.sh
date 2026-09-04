#!/bin/bash

#Make It Executable : chmod +x ../examples/run_file.sh
 
# Delete existing error analysis file if it exists
rm -r ../build/error_analysis.txt

# Build r_refinement_square before running
make rh_poisson_example -j 2
# make rh_adaptiveAdvectiondiffusion -j 15
# make rh_elasticity_example -j 15

# Path to the executable
EXECUTABLE="./bin/rh_poisson_example"
# EXECUTABLE="./bin/rh_adaptiveAdvectiondiffusion"

# Parameters (tags) to run the executable with -r 1: GARU, 2: PUCA, 3: BULK, 4: PBULK
TAGS=(
    "-l -1 -u 3 -f 12. -q 1 -c 1 -v 0 -e 0"
    # .. Poisson 2D ...
    # "-r 2 -u 5 -f  0.  -l 2 -a 0.0 -c 0 -e 0"
    # "-r 2 -u 5 -f  0.  -l 2 -a 0.5 -c 0 -e 0"
    # "-r 2 -u 5 -f  0.  -l 2 -a 0.7 -c 1 -e 0"
    # "-r 2 -u 5 -f  10. -l 5 -a 0.0 -c 0 -e 0"
    # "-r 2 -u 5 -f  10. -l 6 -a 0.5 -c 0 -e 0"
    # "-r 2 -u 5 -f  10. -l 6 -a 0.7 -c 1 -e 0"
    # 3 dimensions case
    # "-r 2 -u 1 -f   0.  -l 0 -a 0.0 -c 0 -d "pde/example3D.xml""
    # "-r 2 -u 1 -f   0.  -l 4 -a 0.5 -c 0 -d "pde/example3D.xml""
    # "-r 2 -u 1 -f   0.  -l 4 -a 0.7 -c 1 -d "pde/example3D.xml""
    # "-r 2 -u 1 -f    9.  -l 3 -a 0.0 -c 0 -d "pde/example3D.xml""
    # "-r 2 -u 1 -f  10.  -l 4 -a 0.5 -c 0 -d "pde/example3D.xml""
    # "-r 2 -u 1 -f  10.  -l 4 -a 0.7 -c 1 -d "pde/example3D.xml""
    #.. Advection diffusion ... -f 0. : without r-refinement
    # "-r 2 -u 2  -f   0. -l 6 -a 0.7 -c 1 -e 1"
    # "-r 2 -u 2  -f  12. -l 6 -a 0.7 -c 1 -e 1"
)


# Run the executable with each set of parameters
for TAG in "${TAGS[@]}"; do
    echo "Running $EXECUTABLE with parameters: $TAG"
    $EXECUTABLE --errorsave $TAG
    echo "-------------------------------------------------"
done