#include"header_file.h"

int turing_criteria (struct RootSet RS[SampleWorkN][MultiRootN], long unsigned int net, long unsigned int sam, long unsigned int samwork,int rn){

	double partial[EdgeN];
	double posi,nega,regu,dyn;
	int i,j;

	for(i=1;i<=EdgeN;i++)
	{
		if(Net[net].topology[i]!=0)
		{
			nega=1.0;
			for(j=((int)(i-1)/NodeN)*NodeN+1;j<=((int)(i-1)/NodeN)*NodeN+NodeN;j++){
				if(i==j) continue;
				else if(Net[net].topology[j]==-1)
					nega/=(1.0+indeX(RS[samwork][rn].steady[(j-1)%NodeN]/Sample[sam].K[j],N[j]));
			}
			if(Net[net].topology[i]==1){
				regu=indeX(Sample[sam].K[i]/RS[samwork][rn].steady[(i-1)%NodeN],N[i]);
				dyn=Sample[sam].h[(i-1)/NodeN]*N[i]*regu/(indeX(1.0+regu,2)*RS[samwork][rn].steady[(i-1)%NodeN])*nega;
			}
			else{
				regu=indeX(RS[samwork][rn].steady[(i-1)%NodeN]/Sample[sam].K[i],N[i]);
				posi=0.0;
				for(j=((int)(i-1)/NodeN)*NodeN+1;j<=((int)(i-1)/NodeN)*NodeN+NodeN;j++){
					if(Net[net].topology[j]==1)
						posi+=1.0/(1.0+indeX(Sample[sam].K[j]/RS[samwork][rn].steady[(j-1)%NodeN],N[j]));
				}
				dyn=-Sample[sam].h[(i-1)/NodeN]*N[i]*regu/(indeX(1.0+regu,2)*RS[samwork][rn].steady[(i-1)%NodeN])*posi*nega;
			}
		}
		else
			dyn=0.0;
		if((i-1)/NodeN==(i-1)%NodeN)
			dyn-=Sample[sam].g[(i-1)/NodeN];
		partial[i-1]=dyn;
		RS[samwork][rn].partial[i-1]=dyn;
	}

	double point=0.5*EIGENSCAN_STEP;
	double max_eigen=0,tem_eigen,max_point=0,left_point=0.0,right_point=1.0;
	int on_off=0,n_interval=0,kappawork=0;

	do{
		if(eigenscan_step(partial,sam,point,&tem_eigen)==GSL_SUCCESS)
		{
			if(max_eigen<tem_eigen)
			{
				max_eigen=tem_eigen;
				max_point=point;
			}
			if(on_off==0)
			{
				if(n_interval>=RealIntervalN){
					n_interval++;
					break;
				}
				RS[samwork][rn].RE[n_interval].k_min=sqrt( EIGENVALUESCAN_MINKAPPA*pow(KAPPA_SCALE_RATIO, point) );
				left_point=point;
				on_off=1;
			}
			else if(1.0-point<EIGENSCAN_STEP)
			{
				right_point=1.0;
				RS[samwork][rn].RE[n_interval].k_max=-1;
				if(max_point-left_point<EIGENSCAN_STEP)
					RS[samwork][rn].RE[n_interval].k_mid=-7.0;  // -10~-5--> kc=0
				else if(right_point-max_point<EIGENSCAN_STEP)
					RS[samwork][rn].RE[n_interval].k_mid=-3.0; // -5~0--> kc more than max_k
				else
					RS[samwork][rn].RE[n_interval].k_mid=sqrt(EIGENVALUESCAN_MINKAPPA*pow(KAPPA_SCALE_RATIO, max_point));

				eigen_set(partial,RS,max_point,sam,samwork,rn,n_interval);
				image_find(RS, partial, net, sam, samwork, rn, n_interval, left_point, right_point);
				n_interval++;
				break;
			}
		}
		else if(on_off==1)
		{
			right_point=point;
			RS[samwork][rn].RE[n_interval].k_max=sqrt( EIGENVALUESCAN_MINKAPPA*pow(KAPPA_SCALE_RATIO, point) );
			kappawork++;

			if(max_point-left_point<EIGENSCAN_STEP){
				RS[samwork][rn].RE[n_interval].k_mid=-7.0;  // -10~-5--> kc=0
			}
			else if(right_point-max_point<EIGENSCAN_STEP){
				RS[samwork][rn].RE[n_interval].k_mid=-3.0; // -5~0--> kc more than max_k
			}
			else{
				RS[samwork][rn].RE[n_interval].k_mid=sqrt(EIGENVALUESCAN_MINKAPPA*pow(KAPPA_SCALE_RATIO, max_point));
			}
			eigen_set(partial,RS,max_point,sam,samwork,rn,n_interval);
			image_find(RS, partial, net, sam, samwork, rn, n_interval, left_point, right_point);

			max_eigen=0;
			on_off=0;
			n_interval++;
		}
		point+=EIGENSCAN_STEP;
	}while(point<1.0);

	RS[samwork][rn].intervalReal=n_interval;

	if(n_interval==0 || kappawork==0)
		return GSL_CONTINUE;
	else
		return GSL_SUCCESS;
}


void eigen_set(double partial[EdgeN], struct RootSet RS[SampleN][MultiRootN], double point, long unsigned int sam, long unsigned int samwork,int rn, int realInterval){

	double kappa;
	double data[EdgeN];
	int i;

	kappa = EIGENVALUESCAN_MINKAPPA * pow(EIGENVALUESCAN_MAXKAPPA/EIGENVALUESCAN_MINKAPPA, point);
	for (i=0;i<EdgeN;i++){
		data[i]=partial[i];
		if( (int)i/NodeN == i%NodeN ) data[i]-=kappa*Sample[sam].D[(int)i/NodeN];
	}

	gsl_matrix_view m;
	gsl_vector_complex * eval;
	gsl_eigen_nonsymm_workspace * w;

	m = gsl_matrix_view_array (data, NodeN, NodeN);
    eval = gsl_vector_complex_alloc (NodeN);
    w = gsl_eigen_nonsymm_alloc (NodeN);
   	gsl_eigen_nonsymm (&m.matrix, eval, w);

	for ( i = 0; i < NodeN; i ++ )
	{
		RS[samwork][rn].RE[realInterval].eigen[i][0]=GSL_REAL( gsl_vector_complex_get( eval, i ));
		RS[samwork][rn].RE[realInterval].eigen[i][1]=GSL_IMAG( gsl_vector_complex_get( eval, i ));
	}
	gsl_vector_complex_free(eval);
	gsl_eigen_nonsymm_free (w);

}


int eigenscan_step(double partial[EdgeN],long unsigned int sam,double point,double*eigen){

	double kappa,max_eigen=0;
	double data[EdgeN];
	int i,sign=0;

	kappa = EIGENVALUESCAN_MINKAPPA * pow(EIGENVALUESCAN_MAXKAPPA/EIGENVALUESCAN_MINKAPPA, point);
	for (i=0;i<EdgeN;i++){
		data[i]=partial[i];
		if( (int)i/NodeN == i%NodeN ) data[i]-=kappa*Sample[sam].D[(int)i/NodeN];
	}

	gsl_matrix_view m;
	gsl_vector_complex * eval;
	gsl_eigen_nonsymm_workspace * w;

	m = gsl_matrix_view_array (data, NodeN, NodeN);
    eval = gsl_vector_complex_alloc (NodeN);
    w = gsl_eigen_nonsymm_alloc (NodeN);
   	gsl_eigen_nonsymm (&m.matrix, eval, w);

	for ( i = 0; i < NodeN; i ++ )
		if ( fabs(GSL_IMAG( gsl_vector_complex_get( eval, i )))<1e-15 &&
				GSL_REAL( gsl_vector_complex_get( eval, i ) ) > 0)
		{
			sign = 1;
		    if(max_eigen<GSL_REAL( gsl_vector_complex_get( eval, i )))
		    	max_eigen=GSL_REAL( gsl_vector_complex_get( eval, i ));
		}

		gsl_vector_complex_free(eval);
		gsl_eigen_nonsymm_free (w);

	if (sign==1)
	{
	*eigen=max_eigen;
		return GSL_SUCCESS;
	}
	else
	{
		*eigen=-1;
		return GSL_CONTINUE;
	}
}


void image_find(struct RootSet RS[SampleN][MultiRootN], double partial[EdgeN], long unsigned int net, long unsigned int sam, long unsigned int samwork,int rn, int n_interval, double left_point, double right_point){

	double point=left_point+0.5*EIGENSCAN_STEP;
	double max_eigen=0,tem_eigen,max_point,left_pointIM,right_pointIM;
	int on_off=0;
	int n_image=0;

	do{
		if(eigenscan_stepIM(partial,sam,point,&tem_eigen)==GSL_SUCCESS)
		{
			if(max_eigen<tem_eigen)
			{
				max_eigen=tem_eigen;
				max_point=point;
			}
			if(on_off==0)
			{
				if(n_image>=ImageIntervalN){
					n_image++;
					break;
				}
				if(point-left_point<EIGENSCAN_STEP)
					RS[samwork][rn].RE[n_interval].IM[n_image].k_minIM=-1;
				else
					RS[samwork][rn].RE[n_interval].IM[n_image].k_minIM=sqrt( EIGENVALUESCAN_MINKAPPA*pow(KAPPA_SCALE_RATIO, point) );
				left_pointIM=point;
				on_off=1;
			}
		}
		else if(on_off==1)
		{
			right_pointIM=point;
			if(right_point-point<EIGENSCAN_STEP)
				RS[samwork][rn].RE[n_interval].IM[n_image].k_maxIM=-1;
			else
				RS[samwork][rn].RE[n_interval].IM[n_image].k_maxIM=sqrt( EIGENVALUESCAN_MINKAPPA*pow(KAPPA_SCALE_RATIO, point) );
			if(max_point-left_pointIM<EIGENSCAN_STEP){
				RS[samwork][rn].RE[n_interval].IM[n_image].k_midIM=-7.0;  // -10~-5--> kc=k_min
			}
			else if(right_pointIM-max_point<EIGENSCAN_STEP){
				RS[samwork][rn].RE[n_interval].IM[n_image].k_midIM=-3.0; // -5~0--> kc=k_max
			}
			else{
				RS[samwork][rn].RE[n_interval].IM[n_image].k_midIM=sqrt(EIGENVALUESCAN_MINKAPPA*pow(KAPPA_SCALE_RATIO, max_point));
			}
			eigen_setIM(partial,RS,max_point,sam,samwork,rn,n_interval,n_image);

			on_off=0;
			max_eigen=0;
			n_image++;
		}
		point+=EIGENSCAN_STEP;
	}while(point<right_point);

	if(on_off==1)
	{
		RS[samwork][rn].RE[n_interval].IM[n_image].k_maxIM=-1.0;
		if(max_point-left_pointIM<EIGENSCAN_STEP){
			RS[samwork][rn].RE[n_interval].IM[n_image].k_midIM=-7.0;  // -10~-5--> kc=k_min
		}
		else if(right_pointIM-max_point<EIGENSCAN_STEP){
			RS[samwork][rn].RE[n_interval].IM[n_image].k_midIM=-3.0; // -5~0--> kc=k_max
		}
		else{
			RS[samwork][rn].RE[n_interval].IM[n_image].k_midIM=sqrt(EIGENVALUESCAN_MINKAPPA*pow(KAPPA_SCALE_RATIO, max_point));
		}
		eigen_setIM(partial,RS,max_point,sam,samwork,rn,n_interval,n_image);
		n_image++;
	}

	RS[samwork][rn].RE[n_interval].intervalImage=n_image;
}


int eigenscan_stepIM(double partial[EdgeN],long unsigned int sam,double point,double*eigen){

	double kappa;
	double data[EdgeN];
	int i;

	kappa = EIGENVALUESCAN_MINKAPPA * pow(EIGENVALUESCAN_MAXKAPPA/EIGENVALUESCAN_MINKAPPA, point);
	for (i=0;i<EdgeN;i++){
		data[i]=partial[i];
		if( (int)i/NodeN == i%NodeN ) data[i]-=kappa*Sample[sam].D[(int)i/NodeN];
	}

	gsl_matrix_view m;
	gsl_vector_complex * eval;
	gsl_eigen_nonsymm_workspace * w;

	m = gsl_matrix_view_array (data, NodeN, NodeN);
    eval = gsl_vector_complex_alloc (NodeN);
    w = gsl_eigen_nonsymm_alloc (NodeN);
   	gsl_eigen_nonsymm (&m.matrix, eval, w);

	for ( i = 0; i < NodeN; i ++ ){
		if ( fabs(GSL_IMAG( gsl_vector_complex_get( eval, i ))) > 1e-15 &&
				GSL_REAL( gsl_vector_complex_get( eval, i ) ) > 0 )
		{
			*eigen=GSL_REAL( gsl_vector_complex_get( eval, i ));
			gsl_vector_complex_free(eval);
			gsl_eigen_nonsymm_free (w);
			return GSL_SUCCESS;
		}
	}

	gsl_vector_complex_free(eval);
	gsl_eigen_nonsymm_free (w);
	return GSL_CONTINUE;

}


void eigen_setIM(double partial[EdgeN],struct RootSet RS[SampleN][MultiRootN],double point,long unsigned int sam,long unsigned int samwork,int rn, int realInterval, int imageInterval){

	double kappa;
	double data[EdgeN];
	int i;

	kappa = EIGENVALUESCAN_MINKAPPA * pow(EIGENVALUESCAN_MAXKAPPA/EIGENVALUESCAN_MINKAPPA, point);
	for (i=0;i<EdgeN;i++){
		data[i]=partial[i];
		if( (int)i/NodeN == i%NodeN ) data[i]-=kappa*Sample[sam].D[(int)i/NodeN];
	}

	gsl_matrix_view m;
	gsl_vector_complex * eval;
	gsl_eigen_nonsymm_workspace * w;

	m = gsl_matrix_view_array (data, NodeN, NodeN);
    eval = gsl_vector_complex_alloc (NodeN);
    w = gsl_eigen_nonsymm_alloc (NodeN);
   	gsl_eigen_nonsymm (&m.matrix, eval, w);

	for ( i = 0; i < NodeN; i ++ )
	{
		RS[samwork][rn].RE[realInterval].IM[imageInterval].eigenIM[i][0]=GSL_REAL( gsl_vector_complex_get( eval, i ));
		RS[samwork][rn].RE[realInterval].IM[imageInterval].eigenIM[i][1]=GSL_IMAG( gsl_vector_complex_get( eval, i ));
	}
	gsl_vector_complex_free(eval);
	gsl_eigen_nonsymm_free (w);

}
