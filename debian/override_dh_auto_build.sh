#!/bin/bash
set -e
wmake -j -a src
wmake -j -a applications
wmake -j unittest
