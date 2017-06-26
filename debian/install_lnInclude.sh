#!/bin/bash
set -e

display_usage() {
	echo -e "Usage: install_lnInclude <target_dir>\n"
}

if [ $# -lt 1 ]
then
	display_usage
	exit 1
fi

for dir in $(find . -name lnInclude -type d); do
	install -m 644 -D $dir/* $1
done


