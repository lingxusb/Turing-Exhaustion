#include"header_file.h"

int dynfunc_f( const gsl_vector * x, long unsigned int para[2], gsl_vector * f ) {

	double posi,nega;
	int i,j;
	long unsigned int net=para[0],sam=para[1];

	for(i=0;i<NodeN;i++)
	{
		posi=0.0;
    	for(j=i*NodeN+1;j<=(i+1)*NodeN;j++){
	    	if(Net[net].topology[j]==1)
	    		posi+=1.0/(1.0+indeX(Sample[sam].K[j]/gsl_vector_get(x,(j-1)%NodeN),N[j]));
	    	else if(Net[net].topology[j]==-1)
	    		nega/=(1.0+indeX(gsl_vector_get(x,(j-1)%NodeN)/Sample[sam].K[j],N[j]));
	    }
    	gsl_vector_set( f, i, Sample[sam].B[i]+Sample[sam].h[i]*posi*nega-Sample[sam].g[i]*gsl_vector_get(x,i));
	}
	return GSL_SUCCESS;
}


int dynfunc_df( const gsl_vector * x, long unsigned int para[2], gsl_matrix * mm ) {

	double posi,nega,regu,dyn;
	int i,j;
	long unsigned int net=para[0],sam=para[1];

	for(i=1;i<=EdgeN;i++)
	{
		if(Net[net].topology[i]!=0)
		{
			nega=1.0;
			for(j=((int)(i-1)/NodeN)*NodeN+1;j<=((int)(i-1)/NodeN)*NodeN+NodeN;j++){
				if(i==j)continue;
				if(Net[net].topology[j]==-1)
					nega/=(1.0+indeX(gsl_vector_get(x,(j-1)%NodeN)/Sample[sam].K[j],N[j]));
			}
			if(Net[net].topology[i]==1){
				regu=indeX(Sample[sam].K[i]/gsl_vector_get(x,(i-1)%NodeN),N[i]);
				dyn=Sample[sam].h[(i-1)/NodeN]*N[i]*regu/(indeX(1.0+regu,2)*gsl_vector_get(x,(i-1)%NodeN))*nega;
			}
			else{
				regu=indeX(gsl_vector_get(x,(i-1)%NodeN)/Sample[sam].K[i],N[i]);
				posi=0.0;
				for(j=((int)(i-1)/NodeN)*NodeN+1;j<=((int)(i-1)/NodeN)*NodeN+NodeN;j++){
					if(Net[net].topology[j]==1)
						posi+=1.0/(1.0+indeX(Sample[sam].K[j]/gsl_vector_get(x,(j-1)%NodeN),N[j]));
				}
				dyn=-Sample[sam].h[(i-1)/NodeN]*N[i]*regu/(indeX(1.0+regu,2)*gsl_vector_get(x,(i-1)%NodeN))*posi*nega;
			}
		}
		else
			dyn=0.0;
		if((i-1)/NodeN==(i-1)%NodeN)
			dyn-=Sample[sam].g[(i-1)/NodeN];
		gsl_matrix_set( mm, (i-1)/NodeN, (i-1)%NodeN, dyn );
	}
	return GSL_SUCCESS;
}

int dynfunc_fdf( const gsl_vector * x, long unsigned int para[2], gsl_vector * f, gsl_matrix * J ) {

	dynfunc_f( x, para, f );
	dynfunc_df( x, para, J );

	return GSL_SUCCESS;
}

