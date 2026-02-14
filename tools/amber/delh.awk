#!/bin/awk -f
{
	if(substr($1,1,4)=="ATOM"||substr($1,1,6)=="HETATM"){
		atmn=substr($0,13,4);
		gsub(/ /,"",atmn);
		if(substr(atmn,1,1)=="H"||substr($0,14,1)=="H") next
		print $0;
	}else{
		print $0 
	}
}
