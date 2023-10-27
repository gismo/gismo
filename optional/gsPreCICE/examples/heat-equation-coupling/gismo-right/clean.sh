#!/bin/sh
set -e -u

# From precice/tutorials/tools/cleaning-tools.sh
clean_precice_logs() {
    (
        set -e -u
        cd "$1"
        echo "---- Cleaning up preCICE logs in $(pwd)"
        rm -fv ./precice-*-iterations.log \
            ./precice-*-convergence.log \
            ./precice-*-events.json \
            ./precice-*-events-summary.log \
            ./precice-postProcessingInfo.log \
            ./precice-*-watchpoint-*.log \
            ./precice-*-watchintegral-*.log \
            ./core
    )
}

clean_gismo() {
    (
        set -e -u
        cd "$1"
        echo "--- Cleaning up Gismo case in $(pwd)"
        rm -fv ./*.vtp
        rm -fv ./*.vts
        rm -fv ./*.pvd
        rm -rfv ./preCICE-output/
        clean_precice_logs .
    )
}

clean_gismo .
