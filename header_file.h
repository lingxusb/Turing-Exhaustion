
//head file:

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include<time.h>
#include<sys/types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include "parameters.h"

long unsigned int Q_Order[NetN];
long unsigned int Sam_Order[NetN];
int N[EdgeN+1];

int ValiditySum;
int PosiQSum;

char Catalog_Path [PATH_LEN];
char Catalog_SamplePath [PATH_LEN];
double Random_Number [ RandomN ];

struct ParaSample{

	double B[NodeN];
	double D[NodeN];
	double h[NodeN];
	double g[NodeN];
	double K[EdgeN+1];

	int work;
	long unsigned int workspace[TopoWorkN];
	int image[2];
	int imagespace[2][TopoImageN];
	int multiInterval;
	int multispace[TopoMutirootN];

}Sample[SampleN];

struct NetCharacter{
	int topology[EdgeN+1];
	int validity;
	int test_if;
	long unsigned int order;
	int edge;
	int role[6]; //A<B<C, A<C<B, B<A<C, B<C<A， C<A<B, C<B<A

	int q_value;
	int image[2]; // image[0]->exist, image[1]->whole window
	int multiInterval;

}Net[NetN];

struct RootSet{

	long unsigned int sam_No;
	int root_Sum;

	double steady[NodeN];
	double partial[EdgeN];

	int intervalReal;
	struct{
		double eigen[NodeN][2];
		double k_min;
		double k_mid;
		double k_max;

		int intervalImage;
		struct{
			double eigenIM[NodeN][2];
			double k_minIM;
			double k_midIM;
			double k_maxIM;
		}IM[ImageIntervalN];

	}RE[RealIntervalN];
};

double minInitialSR[NodeN];
double maxInitialSR[NodeN];

double KC_MAX;
double KC_MIN;
double KAPPA_SCALE_RATIO;

// functions declare:

// main.c
int get_command();
int get_net_number();
int free_all();
void prepare();
void view_params();

// dynfuc.c:
int dynfunc_f( const gsl_vector * x, long unsigned int paras[2], gsl_vector * f );
int dynfunc_df( const gsl_vector * x, long unsigned int paras[2], gsl_matrix * mm );
int dynfunc_fdf( const gsl_vector * x, long unsigned int paras[2], gsl_vector * f, gsl_matrix * J );

// generate_net.c
void generate_net();
void net_validity_isolateIF();
void net_noactivateIF();
void net_validity_ABCsymmetryIF();
void initialize_net();
void in_valid_one(long unsigned int n, long unsigned int i);

// find_robust_topology.c
void test_nets();
void get_generate_params();
void test_one_net(long unsigned int n);

// get_sample.c
void get_sample();
void get_random();
int Sample_Read();
int Sample_Create();
int generate_random_number_file();

// open_files.c
void open_files();

// solveroot.c
int solveroot( struct RootSet RS[SampleWorkN][MultiRootN], long unsigned int net, long unsigned int sam, long unsigned int samwork);
int steadyIF(gsl_multiroot_fdfsolver *s, long unsigned int net, long unsigned int sam);
int differRootIF(struct RootSet RS[SampleWorkN][MultiRootN],long unsigned int sam,long unsigned int samwork,int rn,gsl_multiroot_fdfsolver *s);


// test_specific_topology.c
void test_specific_topology(long unsigned int n);
int* create_robustparams_set();

// tools.c
void create_folder();
void seq2topo(long unsigned int n,char topo[]);
void topo2seq(long unsigned int*n,char topo[]);
void set_robusetness();
void q_sort();
double indeX(double x, int y);
void seq2topo_divide(long unsigned int n, char chip[]);
void q_sortSam();

// turing_criteria.c
int turing_criteria (struct RootSet RS[SampleWorkN][MultiRootN], long unsigned int net, long unsigned int sam, long unsigned int samwork,int rn);
void eigen_set_zero(double partial[EdgeN], struct RootSet RS[SampleN][MultiRootN],long unsigned int sam,int rn);
void eigen_set(double partial[EdgeN], struct RootSet RS[SampleN][MultiRootN], double point, long unsigned int sam, long unsigned int samwork,int rn, int realInterval);
int eigenscan_step(double partial[EdgeN],long unsigned int sam,double point,double*eigen);
int eigenscan_step_zero(double partial[EdgeN],int sam,double*eigen);
int eigenscan_stepIM(double partial[EdgeN],long unsigned int sam,double point,double*eigen);
void image_find(struct RootSet RS[SampleN][MultiRootN], double partial[EdgeN], long unsigned int net, long unsigned int sam, long unsigned int samwork,int rn, int n_interval, double left_point, double right_point);
void eigen_setIM(double partial[EdgeN],struct RootSet RS[SampleN][MultiRootN],double point,long unsigned int sam,long unsigned int samwork,int rn, int realInterval, int imageInterval);
// check.c
void show_robust_topology();
void show_valid_topology();
void show_specific_topology(int n);
void view_present_params();
void view_generate_params();

// make_text.c
void text_state();
void text_robust();
void text_random();
void text_sample();
void text_topN(int n);
void text_results();
void text_params();
void text_paraWork(	struct RootSet RS[SampleWorkN][MultiRootN],long unsigned int net);
void text_valid_matrix();
void text_validnet();
void text_sampleWork();
void text_sample_matrix();
