#!/bin/bash
set -e -u

../../../../../../../../../build/bin/test-communication-IGA-IGA-solid -c ../precice_config.xml -r 1 -e 2
