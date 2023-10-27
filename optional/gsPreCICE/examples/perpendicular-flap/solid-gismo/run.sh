#!/bin/bash
set -e -u

../../../../../build_deb/bin/perpendicular-flap -c ../precice-config.xml --plot -r 2 -e 1 -m 10
