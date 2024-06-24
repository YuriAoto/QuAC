#!/usr/bin/awk
BEGIN{
    if(ARGC<3){
	print "Give the name of potential";
	exit 1;
    }
    pot_name = ARGV[2];
    points=0;
    orig=0;
    if(ARGC>3){
	if(ARGV[3] == "points"){
	    points=1;
	}else{
	    print "Unknown argument";
	    exit 1;
	}
	
    }

    p=0;
}
{
    if(p<2)
    {
	if($0 ~  "Potential " pot_name){p=p+1};
	if(p==1 && /^Distance/){p=p+1};
	if(p==1 && points==1 && /R               V/){p=p+1};
    }
    else
    {
	if(/---/){exit};
	print $0
    }
}
