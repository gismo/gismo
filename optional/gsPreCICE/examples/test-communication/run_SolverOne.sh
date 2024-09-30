#!/bin/bash
set -e -u

../../../../build/bin/test-communication -n "SolverOne" -c ./precice-config.xml
