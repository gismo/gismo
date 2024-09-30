#!/bin/bash
set -e -u

../../../../build/bin/test-communication -n "SolverTwo" -c ./precice-config.xml
