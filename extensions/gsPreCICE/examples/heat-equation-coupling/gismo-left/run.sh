#!/bin/bash
set -e -u

../../../../../build/bin/heat-equation-coupling -c ../precice-config.xml --plot -s 0
