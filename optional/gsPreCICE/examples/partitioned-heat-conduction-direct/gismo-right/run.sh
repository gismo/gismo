#!/bin/bash
set -e -u

../../../../../build/bin/partitioned-heat-conduction -c ../precice-config.xml --plot -s 1
