#!/bin/bash

nargs=$#
if ! [ $nargs -eq 2 ]; then
echo "usage: proj_all.sh dir(cat*.fits) config.par -> dir/cls_cat*"
exit 1
fi

dir=$1
conf=$2

ls $dir/cat*.fits
n=`ls -1  $dir/cat*.fits | wc -l`

echo "$n files : OK? [y/n]"
read answer

if [ $answer != 'y' ] ; then
exit 1
fi


if ! [ $? -eq 0 ] ; then
exit
fi

for filein in $dir/cat*.fits ; do

grep -v file $conf > tmp.par
echo "filein=$filein" >> tmp.par
dir=$(dirname $filein)
f=$(basename $filein)
fileout="!"$dir/cls_${f}
echo "fileout=$fileout" >> tmp.par
HEAD/$CMTCONFIG/proj tmp.par

done
