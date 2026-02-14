#!/bin/awk -f
BEGIN{
	nmodel=0
	nter=0
}
{
	if($0=="")next

	if($1=="MODEL")nmodel++
	if($1=="TER")nter++

	if(/^ANISOU/||/^CONECT/)next  # added on 2021/05/29

        if(substr($0,18,3)=="HOH"||substr($0,18,3)=="WAT")next

	if(substr($1,1,4)=="ATOM"||substr($1,1,6)=="HETATM"){

		s=substr($0,18,10)
		if(s!=s0){
			nres++
			res0=substr(s0,5,6)
			s0=s
		}

                atmn=substr($0,13,4);
                gsub(/ /,"",atmn);
                if(substr(atmn,1,1)=="H"||substr($0,14,1)=="H") next

#               fix the altloc issue for small ligands on 2021/05/29

		altloc=substr($0,17,1);
		atom=substr($0,13,4)" "substr($0,18,10);
		if(altloc==" "||arec[nmodel,nter,atom]==0){
#			print $0;
			res=substr($0,22,6)
#			print substr($0,1,16)" "substr($0,18);
 			if(res!=res0)print substr($0,1,16)" "substr($0,18,50);
			if (altloc!=" ") {
				arec[nmodel,nter,atom]=1;
			}
		}
	}else{
		print $0 
	}
}
