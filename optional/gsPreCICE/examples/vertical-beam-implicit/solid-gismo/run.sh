#!/bin/bash
set -e -u

../../../../../build_deb/bin/vertical-beam-implicit -c ../precice-config.xml --plot -r 1 -e 2 -m 10
