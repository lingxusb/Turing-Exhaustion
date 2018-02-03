
//program setting:

#define NewSamIF 	0
#define FolderName  "Turing-3mor-2N_"
//#define SamPath     ../0_code/Sample

#define SampleN  	50000
#define RandomN 	50000

#define Find_Start  	0
#define Find_End   	19682

#define PATH_LEN 	120
#define SAVE_STEP  	100

#define SOLVEROOT_RESET_MAXTIME        	5
#define SOLVEROOT_ITERATION_MAXTIME  	3000
#define SOLVEROOT_STOP_RESIDUAL 	 	1e-8
#define SOLVEROOT_STEADY_NONZERO 		1e-5
#define SOLVEROOT_RELATIVE_DIFFER	   	1e-2
#define SOLVEROOT_MIN_INITIAL_A    		1e-2
#define SOLVEROOT_MIN_INITIAL_B 		1e-2
#define SOLVEROOT_MIN_INITIAL_C 		1e-2
#define SOLVEROOT_MAX_INITIAL_A 		1e+2
#define SOLVEROOT_MAX_INITIAL_B			1e+2
#define SOLVEROOT_MAX_INITIAL_C 		1e+2
#define SOLVEROOT_MULTI_SCALERATIO  	1e+4
#define SOLVEROOT_MULTI_AXISPOINT   	2
#define EIGENSCAN_STEP 					1e-3
#define EIGENVALUESCAN_MINKAPPA 		1e-4
#define EIGENVALUESCAN_MAXKAPPA 		1e+4

#define ROBUST_POOL_SUM 	20

#define SampleWorkN  		1000
#define MultiRootN      	3
#define RealIntervalN    	3
#define ImageIntervalN  	3
#define TopoWorkN       	200
#define TopoImageN      	50
#define TopoMutirootN  		50

//global const:
#define NodeN 3
#define EdgeN 9 // (NodeN*NodeN)
#define NetN 19683
#define ParamsN  14
#define PI 3.14159265358979324

//dynamical function:

#define _Da        	1.0
#define MIN_Db    	0.01  /* Def:Da=1,Dc=0 */
#define MAX_Db 		100.0
#define MIN_Dc      0.01
#define MAX_Dc      100.0

#define N0   2
#define N1   2
#define N2   2
#define N3   2
#define N4   2
#define N5   2
#define N6   2
#define N7   2
#define N8   2
#define N9   2

#define MIN_Ba   0.1
#define MIN_Bb   0.1
#define MIN_Bc   0.1

#define MAX_Ba   0.1
#define MAX_Bb   0.1
#define MAX_Bc   0.1

#define MIN_ha   0.01
#define MIN_hb   0.01
#define MIN_hc   0.01

#define MAX_ha   100.0
#define MAX_hb   100.0
#define MAX_hc   100.0

#define MIN_ga   0.01
#define MIN_gb   0.01
#define MIN_gc   0.01

#define MAX_ga   100.0
#define MAX_gb   100.0
#define MAX_gc   100.0

#define   _K1    1.0
#define   _K5    1.0
#define   _K9    1.0

#define MIN_K2   0.01  /* Def:K1=K5=K9=1 */
#define MIN_K3   0.01
#define MIN_K4   0.01
#define MIN_K6   0.01
#define MIN_K7   0.01
#define MIN_K8 	 0.01

#define MAX_K2  100.0
#define MAX_K3  100.0
#define MAX_K4  100.0
#define MAX_K6  100.0
#define MAX_K7  100.0
#define MAX_K8  100.0


