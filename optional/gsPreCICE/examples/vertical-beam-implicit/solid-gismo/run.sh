#!/bin/bash
set -e -u

../../../../../build/bin/vertical-beam-implicit -c ../precice-config.xml --plot -r 1 -e 2 -m 1
