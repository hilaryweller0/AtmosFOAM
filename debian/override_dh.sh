#!/bin/bash
ARGS=$@
shift $#
source $FOAM_BASHRC
set -e
export FOAM_USER_APPBIN=$DH_ROOT_DIR/build/bin
export FOAM_USER_LIBBIN=$DH_ROOT_DIR/build/lib
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH
$ARGS
