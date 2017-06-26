#!/bin/bash
set -e

display_usage() {
	echo -e "Usage: singularity.bootstrap.sh <version> <codename>\n"
}

if [ $# -le 1 ]
then
	display_usage
	exit 1
fi

VERSION=$1
CODENAME=$2

VERSION=$VERSION envsubst < Singularity > Singularity.$CODENAME
sudo singularity create --size 2048 $CODENAME.img
sudo singularity bootstrap $CODENAME.img Singularity.$CODENAME
