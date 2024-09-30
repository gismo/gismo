#!/bin/sh
set -e -u

. ../../tools/cleaning-tools.sh

DIR="../precice-run"

if [ -d "$DIR" ]; then
    rm -rf "$DIR"
    echo "Deleted"
fi
