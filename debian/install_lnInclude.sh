#!/bin/bash
set -e

display_usage() {
	echo -e "Usage: install_lnInclude <source_dir> <target_dir>\n"
}

if [ $# -le 1 ]
then
	display_usage
	exit 1
fi

src_dir=$1
target_dir=$2

for dir in $(find $src_dir -name lnInclude -type d); do
	install -m 644 -D $dir/* $target_dir
done


