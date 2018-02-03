#include"header_file.h"


int main() {
	
	prepare();
	create_folder();
	generate_net();
	get_sample();
	text_state();
	test_nets();
	text_results();
}


void prepare(){

	N[0]=N0;
	N[1]=N1;
	N[2]=N2;
	N[3]=N3;
	N[4]=N4;
	N[5]=N5;
	N[6]=N6;
	N[7]=N7;
	N[8]=N8;
	N[9]=N9;

	minInitialSR[0]=SOLVEROOT_MIN_INITIAL_A;
	minInitialSR[1]=SOLVEROOT_MIN_INITIAL_B;
	minInitialSR[2]=SOLVEROOT_MIN_INITIAL_C;

	maxInitialSR[0]=SOLVEROOT_MAX_INITIAL_A;
	maxInitialSR[1]=SOLVEROOT_MAX_INITIAL_B;
	maxInitialSR[2]=SOLVEROOT_MAX_INITIAL_C;

	KC_MAX=sqrt(EIGENVALUESCAN_MAXKAPPA);
	KC_MIN=sqrt(EIGENVALUESCAN_MINKAPPA);
	KAPPA_SCALE_RATIO=EIGENVALUESCAN_MAXKAPPA/EIGENVALUESCAN_MINKAPPA;
}
