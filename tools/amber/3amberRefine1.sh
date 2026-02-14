#!/bin/bash 


# check input parameters
if [ $# -lt 3 ]; then
  echo ""
  echo "USAGE: `basename $0`  com.pdb  nmax  output.pdb"
  echo ""
  exit 1
fi

date
tim1=$(date +%s)

source /share/amber14/amber.sh


# parse input parameters
input_pdb=$1
nmax=$2
output_pdb=$3

echo " "
echo "nstep=  $nmax"
echo " "


# create temporary directory
cdir=$(pwd)
tmpdir=$(mktemp -d)
cd $tmpdir


# link the input pdb file
ln -s $cdir/$input_pdb

id=$$

rec=.temp${id}rec
complex=.temp${id}complex
tleapin=.temp${id}leap.in
minin=.temp${id}min.in
temppdb=.temp${id}.pdb
newpdb=.temp${id}new.pdb

fixamberinput.sh $input_pdb > $rec
cat $rec > $complex.pdb


# set tleap parameters
cat <<EOF > $tleapin
source leaprc.ff14SB
source leaprc.phosaa10
mol = loadpdb $complex.pdb
saveamberparm mol $complex.prmtop $complex.prmcrd
quit
EOF

tleap -f $tleapin


# set energy minimization parameters
cat <<EOF > $minin
Minimization
 &cntrl
  imin   = 1,
  maxcyc = $nmax,
  ncyc   = 500,
  ntb    = 0,
  igb    = 0,
  cut    = 9,
  ntr=1,
  restraint_wt=1.0,
  restraintmask=':1-9999@CA,C,N',
 /
EOF

# run sander for energy minimization
sander -O -i $minin -o min.out -p $complex.prmtop -c $complex.prmcrd -ref $complex.prmcrd -r $complex.restrt1

# use ambpdb to generate new pdb file and remove hydrogen atoms
ambpdb -p $complex.prmtop < $complex.restrt1 | delh.awk > $temppdb

# merge and clean pdb files
cat $rec | awk '{if(FILENAME=="-"){if($1=="ATOM"||substr($1,1,6)=="HETATM"){s=substr($0,18,10);if(s!=s0){n++;rnum[n]=substr($0,23,5);ch[n]=substr($0,22,1);s0=s}}}else{if($1=="ATOM"||substr($1,1,6)=="HETATM"){t=substr($0,18,10);if(t!=t0){k++;t0=t};printf"%s%s%s%s\n",substr($0,1,21),ch[k],rnum[k],substr($0,28)}else{print}}}' - $temppdb | sed 's/ HIE / HIS /g;s/ HID / HIS /g' > $newpdb

nrecres=`awk '{if($1=="ATOM"||substr($1,1,6)=="HETATM"){s=substr($0,18,10);if(s!=s0){n++;s0=s}}}END{print n}' $rec`

cleanpdb.sh $newpdb | awk -v nrec=$nrecres '{if($1=="ATOM"||substr($1,1,6)=="HETATM"){s=substr($0,18,10);if(s!=s0){n++;s0=s}}; if(n<=nrec){print $0}}' > $output_pdb

mv -f $output_pdb $cdir 2>/dev/null


# record time and clean up temporary directory
tim2=$(date +%s)
t_tim=$((tim2-tim1))
echo "time for minization is $t_tim s"
date

cd - >/dev/null
rm -rf $tmpdir
