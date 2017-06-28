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

./singularity.bootstrap.sh $VERSION $CODENAME
./dist.sh $CODENAME
