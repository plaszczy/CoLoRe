
file=$1

grep -v file tophat4.par > tmp.par

grep -v file tophat4.par > tmp.par
echo "filein=$file" >> tmp.par
dir=$(dirname $file)
f=$(basename $file)
fileout="!"$dir/cls_${f}
echo "fileout=$fileout" >> tmp.par
../$CMTCONFIG/proj tmp.par
