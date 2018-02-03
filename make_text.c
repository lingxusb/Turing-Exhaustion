#include"header_file.h"

void text_state(){

	text_validnet();
	text_params();
	text_sample();

}


void text_results(){

	text_valid_matrix();
	text_sample_matrix();
//	text_robust();
	text_sampleWork();

}


void text_params(){

	char path[PATH_LEN];
	sprintf(path,"%s/State/Parameters",Catalog_Path);

	FILE*fp=NULL;

	if((fp=fopen(path,"r"))!=NULL){
		fclose(fp);
		fp=NULL;
		int n=0;
		do{
			++n;
			char path2 [PATH_LEN];
			sprintf(path2,"%s-(%d)",path,n);
			if((fp=fopen(path,"r"))==NULL)
			{
				strcpy(path,path2);
				break;
			}
			else{
				fclose(fp);
				fp=NULL;
			}
		}while(n<100);
	}

	if((fp=fopen(path,"w"))==NULL)
	{
		printf("\nfail to create a Params file!\n");
		return;
	}

	fprintf(fp,"\nPath: %s\n",Catalog_Path);
	fprintf(fp,"\n\n\nConst Parameters:\n");

	fprintf(fp,"\n\nGlobal Const:\n\n");
	fprintf(fp,"\nSample Sum:%-6lu   Random Sum:%-6lu   Parameter Sum:%-5d\n",
			(long unsigned int)SampleN,(long unsigned int)RandomN,ParamsN);

	fprintf(fp,"\n\nProgram Setting:\n\n");
	fprintf(fp,"\nMax Reset Time for Solve Root:%d   Max Iteration Time for Solve Root:%d   "
			"Stop Residual for Solve Root:%-8.7f\nInitial Variable Scale for Solve Root:"
			"Scale of Initial Variables for Solve Root-- A0:%5f~%-5f   B0:%5f~%-5f"
			"   C0:%5f~%-5f\nScale of EigenValue Scan:%5f~%-5f   ",SOLVEROOT_RESET_MAXTIME,
			SOLVEROOT_ITERATION_MAXTIME,SOLVEROOT_STOP_RESIDUAL,SOLVEROOT_MIN_INITIAL_A,SOLVEROOT_MAX_INITIAL_A,
			SOLVEROOT_MIN_INITIAL_B,SOLVEROOT_MAX_INITIAL_B,SOLVEROOT_MIN_INITIAL_C,SOLVEROOT_MAX_INITIAL_C,
			EIGENVALUESCAN_MINKAPPA,EIGENVALUESCAN_MAXKAPPA);

	fprintf(fp,"\n\nDynamical Coefficient:\n\n");
	fprintf(fp,"\nn1:%2d  n2:%2d  n3:%2d  n4:%2d  n5:%2d  n6:%2d  n7:%2d  n8:%2d  n9:%2d \t",
			N1,N2,N3,N4,N5,N6,N7,N8,N9);
	fprintf(fp,"\nBa:%7.4f~%-7.4f    Bb:%7.4f~%-7.4f    Bc:%7.4f~%-7.4f\n",MIN_Ba,MAX_Ba,MIN_Bb,MAX_Bb,MIN_Bc,MAX_Bc);
	fprintf(fp,"Da:%7.4f\t",_Da);
	fprintf(fp,"\tDb:%7.4f~%-7.4f\t",MIN_Db,MAX_Db);
	fprintf(fp,"\tDc:%7.4f~%-7.4f\t",MIN_Dc,MAX_Dc);
	fprintf(fp,"\nha:%7.4f~%-7.4f    hb:%7.4f~%-7.4f    hc:%7.4f~%-7.4f\n",MIN_ha,MAX_ha,MIN_hb,MAX_hb,MIN_hc,MAX_hc);
	fprintf(fp,"ga:%7.4f~%-7.4f    gb:%7.4f~%-7.4f    gc:%7.4f~%-7.4f\n",MIN_ga,MAX_gb,MIN_gc,MAX_ga,MIN_gb,MAX_gc);
	fprintf(fp,"\nK1:   1.000           K2:%7.4f~%-7.4f   K3:%7.4f~%-7.4f\n",MIN_K2,MAX_K2,MIN_K3,MAX_K3);
	fprintf(fp,"K4:%7.4f~%-7.4f   K5:   1.000           K6:%7.4f~%-7.4f\n",MIN_K4,MAX_K4,MIN_K6,MAX_K6);
	fprintf(fp,"K7:%7.4f~%-7.4f   K8:%7.4f~%-7.4f   K9:   1.000\n",MIN_K7,MAX_K7,MIN_K8,MAX_K8);

	fprintf(fp,"\n\n\nNet Sum:%-5lu   Valid Net Sum:%-5lu\n",(long unsigned int)NetN,(long unsigned int)ValiditySum);

	fclose(fp);
	printf("\nParameters File created.\n");

}


void text_validnet()
{
	char path[PATH_LEN];
	sprintf(path,"%s/State/ValidNet",Catalog_Path);

	FILE*fp=NULL;

	if((fp=fopen(path,"r"))!=NULL){
		fclose(fp);
		fp=NULL;
		int k=0;
		do{
			++k;
			char path2 [PATH_LEN];
			sprintf(path2,"%s-(%d)",path,k);
			if((fp=fopen(path,"r"))==NULL)
			{
				strcpy(path,path2);
				break;
			}
			else{
				fclose(fp);
				fp=NULL;
			}
		}while(k<100);
	}

	if((fp=fopen(path,"w"))==NULL)
	{
		printf("\nfail to create a ValidNet file!\n");
		return;
	}

	long unsigned int n,i=0,j;
	char topo[20];

	fprintf(fp,"\nSum of valid net:%5lu\n",(long unsigned int)ValiditySum);
	fprintf(fp,"\nRatio of valid net:%5f\n",(double)ValiditySum/(long unsigned int)NetN);
	fprintf(fp,"\nSum of positive-Q net:%5d\n",PosiQSum);
	fprintf(fp,"\n\nValid Net:\n\n");

	for (n=0;n<NetN;n++){
		if (Net[n].validity==1)
		{
			if(n<=Find_End&&n>=Find_Start)
				fprintf(fp,"*No.%-5lu  ",++i);
			else
				fprintf(fp," No.%-5lu  ",++i);
			fprintf(fp,"Seq: %-5lu    ",n);
			seq2topo_divide(n,topo);
			fprintf(fp,"Topo: %s  ",topo);
			fprintf(fp," EdgeN: %-5d ",Net[n].edge);
			if(Net[n].test_if==0){
				fprintf(fp,"TestIf: Not been tested yet\n");
			}
			else{
				fprintf(fp,"Q-value: %-5d ",Net[n].q_value);
				fprintf(fp,"Q-Ratio: %-4f  ", (double)Net[n].q_value/(long unsigned int)SampleN);
				fprintf(fp,"MultiRT: %-3d ", Net[n].multiInterval);
				fprintf(fp,"ImYes: %-3d ", Net[n].image[0]);
				fprintf(fp,"ImAll: %-3d ", Net[n].image[1]);
				fprintf(fp,"Roles: ");
				for(j=0;j<6;j++)
					fprintf(fp,"%-3d ",Net[n].role[j]);
				fprintf(fp,"\n");
			}
		}
	}
	fclose(fp);
	printf("\nValidNet file created\n");
}

void text_paraWork(	struct RootSet RS[SampleWorkN][MultiRootN],long unsigned int net)
{
	char topo[20],path[PATH_LEN];
	seq2topo_divide(net,topo);
	sprintf(path,"%s/Topo/%s",Catalog_Path,topo);
	if(Net[net].image[0]>0)
		strcat(path,"*");  // image part exists
	if(Net[net].multiInterval>0)
		strcat(path,"#"); // multiroot exists


	FILE*fp;
	if((fp=fopen(path,"w"))==NULL)
	{
		printf("\nfail to create a paraWork file!\n");
		return;
	}

	int sam;
	int rn,i,j;

	for(sam=0;sam<Net[net].q_value;++sam)
	{
		fprintf(fp,"\nNo.%-3d\tSeq.%-6lu\t Topo: ",sam+1,RS[sam][0].sam_No);
	    for(i=1;i<EdgeN;i++) fprintf(fp,"%2d,",Net[net].topology[i]);
	    fprintf(fp,"%2d\t N: ",Net[net].topology[i]);
	    for(i=1;i<EdgeN;i++) fprintf(fp,"%d,",N[i]);
	    fprintf(fp,"%d\t ",N[i]);
	    fprintf(fp,"B: %-11.8f, %-11.8f, %-11.8f\n",Sample[RS[sam][0].sam_No].B[0],Sample[RS[sam][0].sam_No].B[1],Sample[RS[sam][0].sam_No].B[2]);
	    fprintf(fp,"D: %-11.8f, %-11.8f, %-11.8f\t h: %-11.8f, %-11.8f, %-11.8f\t g: %-11.8f, %-11.8f, %-11.8f\n",
	    		Sample[RS[sam][0].sam_No].D[0],Sample[RS[sam][0].sam_No].D[1],Sample[RS[sam][0].sam_No].D[2],Sample[RS[sam][0].sam_No].h[0],Sample[RS[sam][0].sam_No].h[1],Sample[RS[sam][0].sam_No].h[2],Sample[RS[sam][0].sam_No].g[0],Sample[RS[sam][0].sam_No].g[1],Sample[RS[sam][0].sam_No].g[2]);
	    fprintf(fp,"K: %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f\n",
	    		Sample[RS[sam][0].sam_No].K[1],Sample[RS[sam][0].sam_No].K[2],Sample[RS[sam][0].sam_No].K[3],Sample[RS[sam][0].sam_No].K[4],Sample[RS[sam][0].sam_No].K[5],Sample[RS[sam][0].sam_No].K[6],Sample[RS[sam][0].sam_No].K[7],Sample[RS[sam][0].sam_No].K[8],Sample[RS[sam][0].sam_No].K[9]);
		for(rn=0;rn<RS[sam][0].root_Sum;++rn)
		{
			fprintf(fp,"Steady.%-2d:  (%-11.8f, %-11.8f, %-11.8f)\n",rn+1,RS[sam][rn].steady[0],RS[sam][rn].steady[1],RS[sam][rn].steady[2]);
			fprintf(fp,"\tJacobi: ");
			for(i=0;i<EdgeN-1;i++) fprintf(fp,"%-11.8f,  ",RS[sam][rn].partial[i]);
			fprintf(fp,"%-11.8f\n",RS[sam][rn].partial[i]);
			if(RS[sam][rn].intervalReal==0){
				fprintf(fp,"\tkc-Scale: no instable k found within %7.4f~%-7.4f for root.%d\n",KC_MIN,KC_MAX,rn+1);
			}
			else{
				for(i=0;i<RS[sam][rn].intervalReal;i++){
					if(i>=RealIntervalN){
						fprintf(fp,"\tSum of RealInterval overflow!\n");
						break;
					}
					fprintf(fp,"\tWindow.%1d:\t",i+1);
					if(RS[sam][rn].RE[i].k_min<0) fprintf(fp,"kc: <%7.4f~",KC_MIN);
					else fprintf(fp,"kc: %7.4f~",RS[sam][rn].RE[i].k_min);
					if(RS[sam][rn].RE[i].k_max<0) fprintf(fp,">%-7.4f\t",KC_MAX);
					else fprintf(fp,"%-7.4f\t",RS[sam][rn].RE[i].k_max);
					if(RS[sam][rn].RE[i].k_mid<-5.0 && RS[sam][rn].RE[i].k_mid>-10.0) fprintf(fp,"kc-max: 0\tLc: infinite\n");
					else if(RS[sam][rn].RE[i].k_mid<0 && RS[sam][rn].RE[i].k_mid>-5.0) fprintf(fp,"kc-max: >%-7.4f\tLc: <%-7.4f\n",KC_MAX,2*PI/KC_MAX);
					else fprintf(fp,"kc: %-7.4f\tLc: %-7.4f\n",RS[sam][rn].RE[i].k_mid,2*PI/RS[sam][rn].RE[i].k_mid);
					fprintf(fp,"\tkc-eigen: E1=%10.7f+%10.7fi    E2=%10.7f+%10.7fi    E3=%10.7f+%10.7fi\n",RS[sam][rn].RE[i].eigen[0][0],RS[sam][rn].RE[i].eigen[0][1],
							RS[sam][rn].RE[i].eigen[1][0],RS[sam][rn].RE[i].eigen[1][1],RS[sam][rn].RE[i].eigen[2][0],RS[sam][rn].RE[i].eigen[2][1]);
					if(RS[sam][rn].RE[i].intervalImage>0){
						for(j=0;j<RS[sam][rn].RE[i].intervalImage;j++){
							if(j>=ImageIntervalN){
								fprintf(fp,"\tSum of ImageIntervalN overflow!\n");
								break;
							}
							fprintf(fp,"\t\tImageInterval.%1d:\t",j+1);
							if(RS[sam][rn].RE[i].IM[j].k_minIM<0) fprintf(fp,"kIm: <%7.4f~",RS[sam][rn].RE[i].k_min);
							else fprintf(fp,"kIm: %7.4f~",RS[sam][rn].RE[i].IM[j].k_minIM);
							if(RS[sam][rn].RE[i].IM[j].k_maxIM<0) fprintf(fp,">%-7.4f\t",RS[sam][rn].RE[i].k_max);
							else fprintf(fp,"%-7.4f\t",RS[sam][rn].RE[i].IM[j].k_maxIM);
							if(RS[sam][rn].RE[i].IM[j].k_midIM<-5.0 && RS[sam][rn].RE[i].IM[j].k_midIM>-10.0) fprintf(fp,"kIm-max = kc-left\n");
							else if(RS[sam][rn].RE[i].IM[j].k_midIM<0.0 && RS[sam][rn].RE[i].IM[j].k_midIM>-5.0) fprintf(fp, "kIm-max = kc-right\n");
							else fprintf(fp,"kIm-max: %-7.4f\tLc: %-7.4f\n",RS[sam][rn].RE[i].IM[j].k_midIM,2*PI/RS[sam][rn].RE[i].IM[j].k_midIM);
							fprintf(fp,"\t\tkIm-eigen: E1=%10.7f+%10.7fi    E2=%10.7f+%10.7fi    E3=%10.7f+%10.7fi\n",RS[sam][rn].RE[i].IM[j].eigenIM[0][0],RS[sam][rn].RE[i].IM[j].eigenIM[0][1],
									RS[sam][rn].RE[i].IM[j].eigenIM[1][0],RS[sam][rn].RE[i].IM[j].eigenIM[1][1],RS[sam][rn].RE[i].IM[j].eigenIM[2][0],RS[sam][rn].RE[i].IM[j].eigenIM[2][1]);
						}
					}
				}
			}
		}
	}
	fclose(fp);
	printf("paraWork file created\n\n");
}

void text_sampleWork(){

	char path[PATH_LEN];
	sprintf(path,"%s/SampleWork",Catalog_Path);
	FILE*fp;
	if((fp=fopen(path,"w"))==NULL)
	{
		printf("\nfail to create a SampleWork file!\n");
		return;
	}

	long unsigned int i, j;
	char topo[20];
	for(i=0;i<SampleN;i++)
	{
		if(!Sample[Sam_Order[i]].work>0)
			break;
		fprintf(fp,"\n\nNo.%-6lu\tSample.%-6lu\t\t",i,Sam_Order[i]);
	    fprintf(fp,"\tB:%-11.8f, %-11.8f, %-11.8f          ",Sample[Sam_Order[i]].B[0],Sample[Sam_Order[i]].B[1],Sample[Sam_Order[i]].B[2]);
	    fprintf(fp,"N:%2d %2d %2d %2d %2d %2d %2d %2d %2d\n",N[1],N[2],N[3],N[4],N[5],N[6],N[7],N[8],N[9]);
	    fprintf(fp,"D: %-11.8f, %-11.8f, %-11.8f\t h: %-11.8f, %-11.8f, %-11.8f\t g: %-11.8f, %-11.8f, %-11.8f\n",
	    		Sample[Sam_Order[i]].D[0],Sample[Sam_Order[i]].D[1],Sample[Sam_Order[i]].D[2],Sample[Sam_Order[i]].h[0],Sample[Sam_Order[i]].h[1],Sample[Sam_Order[i]].h[2],Sample[Sam_Order[i]].g[0],Sample[Sam_Order[i]].g[1],Sample[Sam_Order[i]].g[2]);
	    fprintf(fp,"K: %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f, %-11.8f\n",
	    		Sample[Sam_Order[i]].K[1],Sample[Sam_Order[i]].K[2],Sample[Sam_Order[i]].K[3],Sample[Sam_Order[i]].K[4],Sample[Sam_Order[i]].K[5],Sample[Sam_Order[i]].K[6],Sample[Sam_Order[i]].K[7],Sample[Sam_Order[i]].K[8],Sample[Sam_Order[i]].K[9]);
	    fprintf(fp,"Topo-Work: \t");
	    for(j=0;j<Sample[Sam_Order[i]].work;j++){
	    	if(j%5==0&&j!=0) fprintf(fp,"\n           ");
	    	seq2topo_divide(Sample[Sam_Order[i]].workspace[j],topo); fprintf(fp,"%s\t",topo);
	    }
	    if(Sample[Sam_Order[i]].image[0]>0){
	    	fprintf(fp,"\nTopo-Image: \t");
	    	for(j=0;j<Sample[Sam_Order[i]].image[0];j++){
	    		seq2topo_divide(Sample[Sam_Order[i]].imagespace[0][j],topo); fprintf(fp,"%s\t",topo);
	    		if(j%5==0) fprintf(fp,"\n            ");
	    	}
	    }
	    if(Sample[Sam_Order[i]].image[1]>0){
	    	fprintf(fp,"\nTopo-ImageALLWindow: \t");
	    	for(j=0;j<Sample[Sam_Order[i]].image[1];j++){
	    		if(j%5==0) fprintf(fp,"\n                     ");
	    		seq2topo_divide(Sample[Sam_Order[i]].imagespace[1][j],topo); fprintf(fp,"%s\t",topo);
	    	}
	    }
	    if(Sample[Sam_Order[i]].multiInterval>0){
	    	fprintf(fp,"\nTopo-MultiWindow: \t");
	    	for(j=0;j<Sample[Sam_Order[i]].multiInterval;j++){
	    		if(j%5==0) fprintf(fp,"\n                  ");
	    		seq2topo_divide(Sample[Sam_Order[i]].multispace[j],topo); fprintf(fp,"%s\t",topo);
	    	}
	    }
	}
	fclose(fp);
}

void text_sample(){

	char path[PATH_LEN];
	sprintf(path,"%s/State/Sample-%lu",Catalog_Path,(long unsigned int)SampleN);

	FILE*fp=NULL;

	if((fp=fopen(path,"r"))!=NULL){
		fclose(fp);
		fp=NULL;
		int n=0;
		do{
			++n;
			char path2 [PATH_LEN];
			sprintf(path2,"%s-(%d)",path,n);
			if((fp=fopen(path,"r"))==NULL)
			{
				strcpy(path,path2);
				break;
			}
			else{
				fclose(fp);
				fp=NULL;
			}
		}while(n<100);
	}

	if((fp=fopen(path,"w"))==NULL)
	{
		printf("\nfail to create a Sample file!\n");
		return;
	}

	    fprintf(fp,"     Db        Dc       ha        hb        hc        ga        gb        gc        K2        K3        K4        K6        K7        K8\n");
    
    long unsigned int n;
	for(n=0;n<SampleN;n++)
	{
		fprintf(fp," %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n",
				Sample[n].D[1],Sample[n].D[2],Sample[n].h[0],Sample[n].h[1],Sample[n].h[2],Sample[n].g[0],Sample[n].g[1],Sample[n].g[2],
				Sample[n].K[2],Sample[n].K[3],Sample[n].K[4],Sample[n].K[6],Sample[n].K[7],Sample[n].K[8]);
	}
	fclose(fp);
	printf("\nSample file created\n");
}


void text_random(){

	char path[PATH_LEN];
	sprintf(path,"%s/State/Random-%lu",Catalog_Path,(unsigned long int)RandomN);

	long unsigned int n;
	FILE*fp;

	if((fp=fopen(path,"w"))==NULL)
	{
		printf("\nfail to create a Random file!\n");
		return;
	}

	for(n=0;n<RandomN;n++)
	{
		if(n%10==0)
			fprintf(fp,"\n");
		fprintf(fp,"%11.8f ",Random_Number[n]);
	}

	fclose(fp);
	printf("\nRandom file created\n");
}

void text_robust(){

	int i=0;
	for(i=1;i<=ROBUST_POOL_SUM;i++){

		text_topN((PosiQSum*i+ROBUST_POOL_SUM-1)/ROBUST_POOL_SUM);

	}
}

void text_topN(int n){

	int daughter[n][n];
	int mother[n][n];
	int skeleton[n];

	int i,j,k;

	for(i=0;i<n;i++)
	{
		skeleton[i]=1;
		for(j=0;j<n;j++)
		{
			if(j==i){
				daughter[i][j]=0;
				mother[i][j]=0;
				continue;
			}
			if(Net[ Q_Order[i] ].edge<Net[ Q_Order[j] ].edge )
			{
				for(k=1;k<=EdgeN;k++)
					if(Net[ Q_Order[i] ].topology[ k ]!=0)
						if(  Net[ Q_Order[i] ].topology[ k ]!=Net[ Q_Order[j] ].topology[ k ]  )
							break;
				if(k>EdgeN)
					daughter[i][j]=1;
				else
					daughter[i][j]=0;
			}
			else
			{
				for(k=1;k<=EdgeN;k++)
					if(Net[ Q_Order[j] ].topology[ k ]!=0)
						if(  Net[ Q_Order[i] ].topology[ k ]!=Net[ Q_Order[j] ].topology[ k ]  )
							break;
				if(k>EdgeN)
				{
					mother[i][j]=1;
					skeleton[i]=0;
				}
				else
					mother[i][j]=0;
			}
		}
	}

	char path[PATH_LEN];
	sprintf(path,"%s/Cluster/top-%d",Catalog_Path,n);

	FILE* fp;
	fp=fopen(path,"w");
	fprintf(fp,"\nThreshold Q：%d\n\n", Net[Q_Order[n-1]].q_value);

	for(i=0;i<n;i++)
	{
		char topo[10];
		seq2topo(Q_Order[i],topo);
		if(skeleton[i]==1)
			fprintf(fp,"\n*");
		else
			fprintf(fp,"\n ");
		fprintf(fp,"N.o:%-5d \tSeq:%-6lu\tTopo:%s \tEdgeN:%-5d \tGoodParaN:%-5d \tQ.Value:%-5f \t",
				i+1, Q_Order[i], topo, Net[Q_Order[i]].edge, Net[Q_Order[i]].q_value, (double)Net[Q_Order[i]].q_value/(long unsigned int)SampleN);
		for(i=0;i<6;i++)
			fprintf(fp," %3d ",Net[n].role[i]);
		fprintf(fp,"\n");
		fprintf(fp,"Daughter: ");
		for(j=0;j<n;j++)
			if(daughter[i][j]==1)
				fprintf(fp,"%d, ",j);
		if(skeleton[i]==0){
		    fprintf(fp,"\n  Mother: ");
		    for(j=0;j<n;j++)
		    	if(mother[i][j]==1)
		    		fprintf(fp,"%d, ",j);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void text_sample_matrix(){
    
    char path[PATH_LEN];
    FILE*fp;
    
    sprintf(path,"%s/SampleWorkMatrix",Catalog_Path);
    
    if((fp=fopen(path,"w"))==NULL)
    {
        printf("\nfail to create a SampleWorkMatrix file!\n");
        return;
    }
    
    int i, j;
    
    fprintf(fp,"\n N.o    Seq    Q-value    Image    ImageAll   MultiInterval      Family\n");
    for(i=0;i<SampleN;i++)
    {
        if(!Sample[Sam_Order[i]].work>0)
            break;
        fprintf(fp,"%3d   %6lu ",i+1,Sam_Order[i]);
        fprintf(fp," %5d    ",Sample[Sam_Order[i]].work);
        fprintf(fp," %5d    ",Sample[Sam_Order[i]].image[0]);
        fprintf(fp," %5d        ",Sample[Sam_Order[i]].image[1]);
        fprintf(fp," %5d         ",Sample[Sam_Order[i]].multiInterval);
        for(j=0;j<Sample[Sam_Order[i]].work;j++){
            fprintf(fp,"%6lu",Sample[Sam_Order[i]].workspace[j]);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}


void text_valid_matrix(){

	char path[PATH_LEN];
	int n,i;
	FILE*fp;

	sprintf(path,"%s/ValidMatrix",Catalog_Path);

	if((fp=fopen(path,"w"))==NULL)
	{
		printf("\nfail to create a ValidMatrix file!\n");
		return;
	}

	fprintf(fp,"\n N.o    Seq             Topo                Edge    Q-value    Q-ratio    MultiRoot   Image   ImageAll              Roles\n");
	for (n=0;n<PosiQSum;n++){

		fprintf(fp,"%3d    %6lu ",n+1,Q_Order[n]);
		for(i=1;i<=EdgeN;i++)
			fprintf(fp,"%2d ",Net[Q_Order[n]].topology[i]);
        fprintf(fp," %5d   ",Net[Q_Order[n]].edge);
		fprintf(fp," %5d     ",Net[Q_Order[n]].q_value);
		fprintf(fp," %5f     ",(double)Net[Q_Order[n]].q_value/(long unsigned int)SampleN);
		fprintf(fp," %3d      ",Net[Q_Order[n]].multiInterval);
		fprintf(fp," %3d     ",Net[Q_Order[n]].image[0]);
		fprintf(fp," %3d     ",Net[Q_Order[n]].image[1]);
		for(i=0;i<6;i++)
			fprintf(fp," %3d ",Net[Q_Order[n]].role[i]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}
