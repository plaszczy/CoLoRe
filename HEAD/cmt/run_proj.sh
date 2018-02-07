grep -v file tophat4.par > tmp.par
filein=$1
echo "filein=$filein" >> tmp.par
dir=$(dirname $filein)
f=$(basename $filein)
fileout="!"$dir/cls_${f}
echo "fileout=$fileout" >> tmp.par
../$CMTCONFIG/proj tmp.par
