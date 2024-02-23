#!/bin/bash
set -e -u

../../../../../build/bin/vertical-beam-implicit-direct-fluid -c ../precice-config.xml --plot -m 1
