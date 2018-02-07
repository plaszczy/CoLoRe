dir=$1

for filein in $dir/*.fits ; do

grep -v file HEAD/cmt/tophat4.par > tmp.par
echo "filein=$filein" >> tmp.par
dir=$(dirname $filein)
f=$(basename $filein)
fileout="!"$dir/cls_${f}
echo "fileout=$fileout" >> tmp.par
HEAD/$CMTCONFIG/proj tmp.par

done
