#20240321
if [ $# -lt 1 ]; then
   echo ""
   echo "USAGE: `basename $0` com.pdb"
   echo ""
   exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/" && pwd)"
cdir=$(pwd)
tmpdir=$(mktemp -d)
cd $tmpdir

#ln -s $cdir/$1 A.pdb    # only work for relative path
ln -s $1 A.pdb          # only work for absolute path
#input_file=$(readlink -f "$1")
#ln -s $input_file A.pdb

num=0
for i in $(cat A.pdb |awk '{if(substr($0,1,4)=="ATOM")print substr($0,22,1)}' |sort -u)
do
   num=$((num+1))
   cat A.pdb |awk -v ch=$i '{if(substr($0,22,1)==ch)print $0}' > $num.pdb 
   echo $i >> tem
done

tscore=0.0
num2=$((num-1))
for i in $(seq 2 $num)
do
  num3=$((i-1))
  for j in $(seq 1 $num3)
  do
     score=`${SCRIPT_DIR}/ppscore $j.pdb $i.pdb -nomin`
#     score=`2pscore $j.pdb $i.pdb -nomin`
#     score=`ppscore $j.pdb $i.pdb`
     ch1=`head -n $i tem |tail -n 1`
     ch2=`head -n $j tem |tail -n 1`
#     echo $ch1 $ch2 $score
     tscore=$(echo "$tscore + $score" | bc)
  done 
done
#    echo "total score=  $tscore"
echo $tscore
cd - >/dev/null
rm -rf $tmpdir


