#!/bin/bash

if [ $# -lt 1 ]; then
        echo ""
        echo "USAGE: `basename $0`  file1.pdb"
        echo ""
        exit 1
fi

protreslib="/share/amber14/dat/reslib/leap/protein.amberua.prepin"
nareslib="/share/amber14/dat/reslib/leap/na.amberua.prepin"


cat $protreslib $nareslib | \
awk '{if($2=="INT"){resn=$1;fw=1};if($0=="")fw=0;if(fw>0&&NF==11&&$2!="DUMM"&&substr($2,1,1)!="H")print resn,$2 }' | sort -u | \
awk '{\
	if(FILENAME=="-"){\
		s[$1,$2]=1\
	}else{\
	 	if($1=="ATOM"||substr($0,1,6)=="HETATM"){\
			atmn=substr($0,13,4);\
			gsub(/ /,"",atmn);\
			resn=substr($0,18,3);\
			gsub(/ /,"",resn);\
			if(resn=="HIS")resn="HIE"; \
			if(resn=="ADE")resn="DA"; \
			if(resn=="GUA")resn="DG"; \
			if(resn=="CYT")resn="DC"; \
			if(resn=="THY")resn="DT"; \
			if(resn=="URA")resn="RU"; \
			if(resn=="A")resn="RA"; \
			if(resn=="G")resn="RG"; \
			if(resn=="C")resn="RC"; \
			if(resn=="T"||resn=="RT")resn="DT"; \
			if(resn=="U"||resn=="DU")resn="RU"; \
			if(s[resn,atmn]>0)printf"ATOM  %s %-3s %-3s %s\n",substr($0,7,6),atmn,resn,substr($0,22)\
		}else if($1=="TER"){\
			print \
		}\
	}\
}' - $1 | delalt.awk | delh.awk | \
awk '{if($1=="ATOM"){s=substr($0,18,10);if(s!=s0){nres++;n=substr($0,23,4)*1;ch=substr($0,22,1);if(nres>1&&(ch!=ch0||n-n0<0||n-n0>1))print"TER";n0=n;ch0=ch;s0=s}; print substr($0,1,54)}else if($1=="TER"){print $1} }' | uniq | \
awk '{if($1=="TER")nres=0;if($1=="ATOM"||substr($1,1,6)=="HETATM"){s=substr($0,18,10);if(s!=s0){nres++;s0=s}; if(!(nres==1&&(index($0," P ")||index($0," O1P ")||index($0," O2P ")||index($0," O3P ")))>0)print}else{print}}' | \
awk '{s=substr($0,18,10);if(s!=s0){i++;s0=s};n[i]++;a[i,n[i]]=$0;atmn=substr($0,13,4);gsub(/ /,"",atmn);if(atmn=="O2'\''")t[i]=1}END{nres=i;for(i=1;i<=nres;i++){for(k=1;k<=n[i];k++){s=a[i,k];at=substr(s,1,6);if(at=="ATOM  "||at=="HETATM"){resn=substr(s,18,3);gsub(/ /,"",resn);ns=length(resn);if(ns<3){if(t[i]==1){resnx="R"substr(resn,ns,1)}else{resnx="D"substr(resn,ns,1)}}else{resnx=resn};if(resnx=="DU")resnx="RU";if(resnx=="RT")resnx="DT";printf"%s%-3s%s\n",substr(s,1,17),resnx,substr(s,21)}else{print s} }}}'
