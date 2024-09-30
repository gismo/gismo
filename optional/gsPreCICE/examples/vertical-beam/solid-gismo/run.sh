#!/bin/bash
set -e -u

../../../../../build_deb/bin/vertical-beam -c ../precice-config.xml --plot -r 1 -e 1 -m 10
