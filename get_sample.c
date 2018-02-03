#include"header_file.h"

void get_sample(){

	if (!NewSamIF)
	{
		if(Sample_Read()!=GSL_SUCCESS)
			exit(0);
		else
			return;
	}

	long unsigned int n;
	srand((int)time(0));
    long unsigned int temp=rand()%(long unsigned int)RandomN;
    get_random();

	for (n=0;n<SampleN;n++)
	{
		Sample[n].B[0] = MIN_Ba;
		Sample[n].B[1] = MIN_Bb;
		Sample[n].B[2] = MIN_Bc;
//		Sample[n].B[0] = MIN_Ba * pow(MAX_Ba/MIN_Ba,Random_Number[temp++%RandomN]);
//		Sample[n].B[1] = MIN_Bb * pow(MAX_Bb/MIN_Bb,Random_Number[temp++%RandomN]);
//		Sample[n].B[2] = MIN_Bc * pow(MAX_Bc/MIN_Bc,Random_Number[temp++%RandomN]);
		Sample[n].g[0] = MIN_ga * pow(MAX_ga/MIN_ga,Random_Number[temp++%RandomN]);
		Sample[n].g[1] = MIN_gb * pow(MAX_gb/MIN_gb,Random_Number[temp++%RandomN]);
		Sample[n].g[2] = MIN_gc * pow(MAX_gc /MIN_gc,Random_Number[temp++%RandomN]);
		Sample[n].h[0] = MIN_ha * pow(MAX_ha/MIN_ha,Random_Number[temp++%RandomN]);
		Sample[n].h[1] = MIN_hb * pow(MAX_hb/MIN_hb,Random_Number[temp++%RandomN]);
		Sample[n].h[2] = MIN_hc * pow(MAX_hc /MIN_hc,Random_Number[temp++%RandomN]);
		Sample[n].K[1] =_K1;
		Sample[n].K[2] = MIN_K2 * pow(MAX_K2/MIN_K2,Random_Number[temp++%RandomN]);
		Sample[n].K[3] = MIN_K3 * pow(MAX_K3/MIN_K3,Random_Number[temp++%RandomN]);
		Sample[n].K[4] = MIN_K4 * pow(MAX_K4/MIN_K4,Random_Number[temp++%RandomN]);
		Sample[n].K[5] =_K5;
		Sample[n].K[6] = MIN_K6 * pow(MAX_K6/MIN_K6,Random_Number[temp++%RandomN]);
		Sample[n].K[7] = MIN_K7 * pow(MAX_K7/MIN_K7,Random_Number[temp++%RandomN]);
		Sample[n].K[8] = MIN_K8 * pow(MAX_K8/MIN_K8,Random_Number[temp++%RandomN]);
		Sample[n].K[9] =_K9;
		Sample[n].D[0] = _Da;
		Sample[n].D[1] = MIN_Db * pow(MAX_Db/MIN_Db,Random_Number[temp++%RandomN]);
		Sample[n].D[2] = MIN_Dc * pow(MAX_Dc/MIN_Dc,Random_Number[temp++%RandomN]);

		Sample[n].work=0;
		Sample[n].image[0]=0;
		Sample[n].image[1]=0;
		Sample[n].multiInterval=0;
	}
	Sample_Create();
}


void get_random(){

	const gsl_rng_type * T;
	gsl_rng * r;

	long unsigned int i;
	srand((int)time(0));
	gsl_rng_env_setup();

	T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    printf ("\ngenerator type: %s\n", gsl_rng_name (r));
    printf ("seed = %lu \n", gsl_rng_default_seed);
    printf ("first value =  %lu\n", gsl_rng_get (r));

    for (i = 0; i < RandomN; i++)
    	Random_Number[i] = (double) gsl_rng_uniform (r);

    gsl_rng_free (r);
    text_random();

}

int Sample_Create(){

	FILE* fp;

	if ((fp = fopen( Catalog_SamplePath , "w" ))==NULL)
	{
		printf("\nfail to create a Sample file!\n");
		return GSL_CONTINUE;
	}
	long unsigned int n;
	for (n=0;n<SampleN;n++){
		fprintf(fp,"%11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",
				Sample[n].D[1],Sample[n].D[2],Sample[n].h[0],Sample[n].h[1],Sample[n].h[2],Sample[n].g[0],Sample[n].g[1],Sample[n].g[2],
				Sample[n].K[2],Sample[n].K[3],Sample[n].K[4],Sample[n].K[6],Sample[n].K[7],Sample[n].K[8]);
	}
	fclose(fp);

	printf("\nnew Sample file created!\n");
	return GSL_SUCCESS;
}

int Sample_Read(){

	FILE* fp;

	if ((fp = fopen( Catalog_SamplePath , "r" ))==NULL)
	{
		printf("\nfail to open Sample file!\n");
		return GSL_CONTINUE;
	}
	long unsigned int n;
	for (n=0;n<SampleN;n++){
		fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
				&Sample[n].D[1],&Sample[n].D[2],&Sample[n].h[0],&Sample[n].h[1],&Sample[n].h[2],&Sample[n].g[0],&Sample[n].g[1],&Sample[n].g[2],
				&Sample[n].K[2],&Sample[n].K[3],&Sample[n].K[4],&Sample[n].K[6],&Sample[n].K[7],&Sample[n].K[8]);
	}
	fclose(fp);

	for (n=0;n<SampleN;n++){
		Sample[n].B[0] = MIN_Ba;
		Sample[n].B[1] = MIN_Bb;
		Sample[n].B[2] = MIN_Bc;
		Sample[n].K[1] =_K1;
		Sample[n].K[5] =_K5;
		Sample[n].K[9] =_K9;
		Sample[n].D[0] = _Da;
	}

	printf("\nSample file imported!\n");
	return GSL_SUCCESS;
}
