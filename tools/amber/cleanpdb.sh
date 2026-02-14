#!/bin/bash

if [ $# -lt 1 ];then
        echo "Usage: `basename $0`  *.pdb"
        exit 1
fi


egrep -v "^CONECT" $1 | awk '{if($1!="ANISOU")print;if($1=="ENDMDL")exit}' | cleanpdb.awk
