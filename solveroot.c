#include"header_file.h"

int solveroot( struct RootSet RS[SampleWorkN][MultiRootN], long unsigned int net, long unsigned int sam, long unsigned int samwork) {

	int iter,i,j,temp,max_i,rn=0;
	int status;
	double step,start;
	long unsigned int para[2]={net,sam};

	const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridsj; // hybridsj : using derivatives;
	gsl_multiroot_fdfsolver *s;
	gsl_vector *x;
	gsl_multiroot_function_fdf f = {&dynfunc_f,&dynfunc_df,&dynfunc_fdf,NodeN,para};

	max_i=SOLVEROOT_RESET_MAXTIME;

	for (i=0;i<max_i;i++)
	{
    	iter=0;
    	s = gsl_multiroot_fdfsolver_alloc (T, NodeN);
    	x = gsl_vector_alloc( NodeN );
    	if (s==NULL||x==NULL) continue;

    	if(rn==0){
    		srand((int)time(0));
        	for(j=0;j<NodeN;j++)
        		gsl_vector_set( x,j,minInitialSR[j]*pow(maxInitialSR[j]/minInitialSR[j],rand()/(double)RAND_MAX));
    	}
        else{
    		temp=i;
    		for(j=0;j<NodeN;j++)
    		{
    			gsl_vector_set( x,j,RS[samwork][0].steady[j]/SOLVEROOT_MULTI_SCALERATIO*pow(SOLVEROOT_MULTI_SCALERATIO*SOLVEROOT_MULTI_SCALERATIO,start+(temp%SOLVEROOT_MULTI_AXISPOINT)*step));
    			temp=(int)temp/SOLVEROOT_MULTI_AXISPOINT;
    		}
    	}

    	gsl_multiroot_fdfsolver_set (s, &f, x);

    	do{
	    	status = gsl_multiroot_fdfsolver_iterate (s);
    		if (status==GSL_EBADFUNC || status==GSL_ENOPROG)  break;
    		status = gsl_multiroot_test_residual (s->f, SOLVEROOT_STOP_RESIDUAL);
    	}while (status == GSL_CONTINUE && ++iter < SOLVEROOT_ITERATION_MAXTIME);

    	if (status==GSL_SUCCESS && steadyIF(s,net,sam)==GSL_SUCCESS)
    	{
    		if(rn==0)
    		{
    		    for(j=0;j<NodeN;j++)
    		    	RS[samwork][0].steady[j]=gsl_vector_get( s->x,j);
    		    rn=1;
    		    max_i=indeX(SOLVEROOT_MULTI_AXISPOINT,NodeN);
    		    step=1.0/SOLVEROOT_MULTI_AXISPOINT;
    		    start=step*0.5;
    	    }
    	    else if(differRootIF(RS,sam,samwork,rn,s)==GSL_SUCCESS)
         	{
         		if(rn>=MultiRootN)
         		{
         			RS[samwork][0].root_Sum=MultiRootN+1;
         			gsl_multiroot_fdfsolver_free (s);
         			gsl_vector_free (x);
         			return GSL_CONTINUE;
       			}
         		else{
         			for(j=0;j<NodeN;j++)
         				RS[samwork][rn].steady[j]=gsl_vector_get( s->x,j);
         			++rn;
         		}
         	}
    	}
    	gsl_multiroot_fdfsolver_free (s);
    	gsl_vector_free (x);
	}
	RS[samwork][0].root_Sum=rn;
	if(rn==0)
		return GSL_CONTINUE;
	else
		return GSL_SUCCESS;
}

int differRootIF(struct RootSet RS[SampleWorkN][MultiRootN],long unsigned int sam,long unsigned int samwork,int rn,gsl_multiroot_fdfsolver *s){

	int k,j,tag;
	for(k=0;k<rn;k++){
		tag=0;
 		for(j=0;j<NodeN;j++){
 			if(fabs(RS[samwork][k].steady[j]-gsl_vector_get(s->x,j))/(RS[samwork][k].steady[j]+gsl_vector_get(s->x,j))>SOLVEROOT_RELATIVE_DIFFER)
 				tag++;
 		}
 		if(tag>1)
 			continue;
 		else
 			return GSL_CONTINUE;
	}
	return GSL_SUCCESS;

}


int steadyIF(gsl_multiroot_fdfsolver *s, long unsigned int net, long unsigned int sam){

	if(gsl_vector_get(s->x,0)<SOLVEROOT_STEADY_NONZERO||
	   gsl_vector_get(s->x,1)<SOLVEROOT_STEADY_NONZERO||
	   gsl_vector_get(s->x,2)<SOLVEROOT_STEADY_NONZERO)
		return GSL_CONTINUE;

	double partial[EdgeN];
	double posi,nega,regu,dyn;
	int i,j;

	for(i=1;i<=EdgeN;i++)
	{
		if(Net[net].topology[i]!=0)
		{
			nega=1.0;
			for(j=((int)(i-1)/NodeN)*NodeN+1;j<=((int)(i-1)/NodeN)*NodeN+NodeN;j++){
				if(i==j)continue;
				if(Net[net].topology[j]==-1)
					nega/=(1.0+indeX(gsl_vector_get(s->x,(j-1)%NodeN)/Sample[sam].K[j],N[j]));
			}
			if(Net[net].topology[i]==1){
				regu=indeX(gsl_vector_get(s->x,(i-1)%NodeN)/Sample[sam].K[i],N[i]);
				dyn=Sample[sam].h[(i-1)/NodeN]*N[i]*regu/(indeX(1.0+regu,2)*gsl_vector_get(s->x,(i-1)%NodeN))*nega;
			}
			else{
				regu=indeX(Sample[sam].K[i]/gsl_vector_get(s->x,(i-1)%NodeN),N[i]);
				posi=0.0;
				for(j=((int)(i-1)/NodeN)*NodeN+1;j<=((int)(i-1)/NodeN)*NodeN+NodeN;j++){
					if(Net[net].topology[j]==1)
						posi+=1.0/(1.0+indeX(Sample[sam].K[j]/gsl_vector_get(s->x,(j-1)%NodeN),N[j]));
				}
				dyn=-Sample[sam].h[(i-1)/NodeN]*N[i]*regu/(indeX(1.0+regu,2)*gsl_vector_get(s->x,(i-1)%NodeN))*posi*nega;
			}
		}
		else
			dyn=0.0;
		if((i-1)/NodeN==(i-1)%NodeN)
			dyn-=Sample[sam].g[(i-1)/NodeN];
		partial[i-1]=dyn;
	}

	gsl_matrix_view m;
	gsl_vector_complex * eval;
	gsl_eigen_nonsymm_workspace * w;

	m = gsl_matrix_view_array (partial, NodeN, NodeN);
    eval = gsl_vector_complex_alloc (NodeN);
    w = gsl_eigen_nonsymm_alloc (NodeN);
   	gsl_eigen_nonsymm (&m.matrix, eval, w);

	if(GSL_REAL( gsl_vector_complex_get( eval, 0 ) ) > 0 ||
	   GSL_REAL( gsl_vector_complex_get( eval, 1 ) ) > 0 ||
	   GSL_REAL( gsl_vector_complex_get( eval, 2 ) ) > 0)
	{
		gsl_vector_complex_free(eval);
		gsl_eigen_nonsymm_free (w);
		return GSL_CONTINUE;
	}
	else{
		gsl_vector_complex_free(eval);
		gsl_eigen_nonsymm_free (w);
		return GSL_SUCCESS;
	}
}
