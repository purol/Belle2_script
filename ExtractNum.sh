N_B2BIIFixMdst=0
N_B2BIIConvertMdst=0
N_RootOutput=0
for file in $(find "$1" -name "*.log")
do
    temp=`cat $file | grep "B2BIIFixMdst   " | cut -d"|" -f2 | tr -d " "`
    temp2=`cat $file | grep "B2BIIConvertMdst   " | cut -d"|" -f2 | tr -d " "`
    temp3=`cat $file | grep "RootOutput   " | cut -d"|" -f2 | tr -d " "`
    N_B2BIIFixMdst=`expr $N_B2BIIFixMdst + $temp`
    N_B2BIIConvertMdst=`expr $N_B2BIIConvertMdst + $temp2`
    N_RootOutput=`expr $N_RootOutput + $temp3`
done
echo "the number of total event: $N_B2BIIFixMdst"
echo "the number of total event after FixMdst: $N_B2BIIConvertMdst"
echo "the number of total event after Skimming: $N_RootOutput"
