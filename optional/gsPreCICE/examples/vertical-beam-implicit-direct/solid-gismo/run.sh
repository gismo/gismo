#!/bin/bash
set -e -u

../../../../../build/bin/vertical-beam-implicit-direct-solid -c ../precice-config.xml -r 1 -e 2
