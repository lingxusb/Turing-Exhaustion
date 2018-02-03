#include"header_file.h"

void generate_net(){

	int i;
	long unsigned int n,m;

	for ( n = 0 ; n < NetN ; n++ )
	{
		m=n;
		for ( i = 1 ; i <= EdgeN ; i++ )
		{
			Net[n].topology[i] = m % NodeN -1 ;
			m /= NodeN ;
		}
		Net[n].validity=1;
	}

	net_validity_isolateIF();
//	net_noactivateIF();
	net_validity_ABCsymmetryIF();
	initialize_net();
	text_validnet();

}

void net_validity_isolateIF(){

	long unsigned int n;
	for(n=0;n<NetN;n++)
	{
		if (Net[n].validity==-1)
			continue;
		if ((Net[n].topology[2]==0&&Net[n].topology[3]==0)||(Net[n].topology[4]==0&&Net[n].topology[7]==0))
			Net[n].validity=-1;
		else if ((Net[n].topology[4]==0&&Net[n].topology[6]==0)||(Net[n].topology[2]==0&&Net[n].topology[8]==0))
			Net[n].validity=-1;
		else if ((Net[n].topology[7]==0&&Net[n].topology[8]==0)||(Net[n].topology[3]==0&&Net[n].topology[6]==0))
			Net[n].validity=-1;
	}
}


void net_noactivateIF(){

	long unsigned int n;
	int i;
	for(n=0;n<NetN;n++){
		if(Net[n].validity==-1)
			continue;
		for(i=0;i<NodeN;i++){
			if(Net[n].topology[i*NodeN+1]!=1&&Net[n].topology[i*NodeN+2]!=1&&Net[n].topology[i*NodeN+3]!=1)
			{
				Net[n].validity=-1;
				break;
			}
		}
	}
}

void net_validity_ABCsymmetryIF(){


	long unsigned int i,n;
	for(n=0;n<NetN;n++){

		if(Net[n].validity==-1)
			continue;

    	for(i=0;i<NetN;i++){
    		if(Net[i].validity!=-1 && i!=n)
	    	{
    			if(Net[n].topology[1]==Net[i].topology[5] && Net[n].topology[2]==Net[i].topology[4]
         	  &&Net[n].topology[3]==Net[i].topology[6] && Net[n].topology[4]==Net[i].topology[2]
     		  &&Net[n].topology[5]==Net[i].topology[1] && Net[n].topology[6]==Net[i].topology[3]
    		 && Net[n].topology[7]==Net[i].topology[8] && Net[n].topology[8]==Net[i].topology[7]
	    	 && Net[n].topology[9]==Net[i].topology[9]){
    				Net[i].validity=-1;
	    			//in_valid_one(n,i);
	    			break;
    			}
    			else if(Net[n].topology[1]==Net[i].topology[1] && Net[n].topology[2]==Net[i].topology[3]
         	  &&Net[n].topology[3]==Net[i].topology[2] && Net[n].topology[4]==Net[i].topology[7]
     		  &&Net[n].topology[5]==Net[i].topology[9] && Net[n].topology[6]==Net[i].topology[8]
    		 && Net[n].topology[7]==Net[i].topology[4] && Net[n].topology[8]==Net[i].topology[6]
	    	 && Net[n].topology[9]==Net[i].topology[5]){
    				Net[i].validity=-1;
	    			//in_valid_one(n,i);
	    			break;
    			}
    			else if(Net[n].topology[1]==Net[i].topology[9] && Net[n].topology[2]==Net[i].topology[7]
         	  &&Net[n].topology[3]==Net[i].topology[8] && Net[n].topology[4]==Net[i].topology[3]
     		  &&Net[n].topology[5]==Net[i].topology[1] && Net[n].topology[6]==Net[i].topology[2]
    		 && Net[n].topology[7]==Net[i].topology[6] && Net[n].topology[8]==Net[i].topology[4]
	    	 && Net[n].topology[9]==Net[i].topology[5]){
    				Net[i].validity=-1;
	    			//in_valid_one(n,i);
	    			break;
    			}
    			else if(Net[n].topology[1]==Net[i].topology[5] && Net[n].topology[2]==Net[i].topology[6]
         	  &&Net[n].topology[3]==Net[i].topology[4] && Net[n].topology[4]==Net[i].topology[8]
     		  &&Net[n].topology[5]==Net[i].topology[9] && Net[n].topology[6]==Net[i].topology[7]
    		 && Net[n].topology[7]==Net[i].topology[2] && Net[n].topology[8]==Net[i].topology[3]
	    	 && Net[n].topology[9]==Net[i].topology[1]){
    				Net[i].validity=-1;
	    			//in_valid_one(n,i);
	    			break;
    			}
    			else if(Net[n].topology[1]==Net[i].topology[9] && Net[n].topology[2]==Net[i].topology[8]
         	  &&Net[n].topology[3]==Net[i].topology[7] && Net[n].topology[4]==Net[i].topology[6]
     		  &&Net[n].topology[5]==Net[i].topology[5] && Net[n].topology[6]==Net[i].topology[4]
    		 && Net[n].topology[7]==Net[i].topology[3] && Net[n].topology[8]==Net[i].topology[2]
	    	 && Net[n].topology[9]==Net[i].topology[1]){
    				Net[i].validity=-1;
	    			//in_valid_one(n,i);
	    			break;
    			}
	    	}
    	}
	}
}

void in_valid_one(long unsigned int n, long unsigned int i){

	if(Net[n].topology[5]<Net[i].topology[5]){
		Net[i].validity=-1;return;
	}
	else if(Net[n].topology[5]>Net[i].topology[5]){
		Net[n].validity=-1;return;
	}
	else{
		if(Net[n].topology[1]>Net[i].topology[1]){
			Net[i].validity=-1;return;
		}
		else if(Net[n].topology[1]<Net[i].topology[1]){
			Net[n].validity=-1;return;
		}
		else{
			if(Net[n].topology[2]<Net[i].topology[2]){
				Net[i].validity=-1;return;
			}
			else if(Net[n].topology[2]>Net[i].topology[2]){
				Net[n].validity=-1;return;
			}
			else{
				if(Net[n].topology[4]>Net[i].topology[4]){
					Net[i].validity=-1;return;
				}
				else if(Net[n].topology[4]<Net[i].topology[4]){
					Net[n].validity=-1;return;
				}
				else{
					if(Net[n].topology[1]+Net[n].topology[2]+Net[n].topology[3]<Net[i].topology[4]+Net[i].topology[5]+Net[i].topology[6]){
						Net[n].validity=-1;return;
					}
					else{
						Net[i].validity=-1;return;
					}
				}
			}
		}
	}
}


void initialize_net(){

	long unsigned int n;
	int i;

	for ( n = 0 ; n < NetN ; n++ ){

		Net[n].test_if=0;
		Net[n].q_value=0;

		Net[n].edge=0;
		for ( i = 1 ; i <= EdgeN ; i++ )
			if(Net[ n ].topology[ i ] != 0)
				++Net[n].edge;
		for ( i = 0; i<6;i++)
			Net[n].role[i]=0;
	}

	PosiQSum=0;

	ValiditySum=0;
	for ( n = 0 ; n < NetN ; n++ )
	{
		if (Net[n].validity==1)
			ValiditySum++;
	}
	printf("\nRatio of valid net:%5f\n",(double)ValiditySum/(long unsigned int)NetN);

}

