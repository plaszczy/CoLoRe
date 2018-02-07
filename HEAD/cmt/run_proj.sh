grep -v file config.par > tmp.par
filein=$1
echo "filein=$filein" >> tmp.par
fileout=cls_${filein}
echo "fileout=$fileout" >> tmp.par
../$CMTCONFIG/proj tmp.par
