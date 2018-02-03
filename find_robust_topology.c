#include"header_file.h"

void test_nets(){

	long unsigned int n;
	int Test_i=0;

	for (n=Find_Start;n<=Find_End;n++)
	{
		if ( Net[n].validity == -1 || Net[n].test_if==1)
			continue;
		test_one_net(n);
		++Test_i;
	}

	q_sort();
	q_sortSam();
}

void test_one_net(long unsigned int net){

	printf("seq:%-5lu\n",net);//    topo:%s    thread:%d \n",n,topo,omp_get_thread_num());

	struct RootSet RS[SampleWorkN][MultiRootN];
	long unsigned int sam;
	int rn,i,root_work,samwork=0,image,image_whole,mutiInterval;
	Net[net].q_value=0;Net[net].image[0]=0;Net[net].image[1]=0;Net[net].multiInterval=0;

	for ( sam=0; sam < SampleN; sam++ ){

		RS[samwork][0].sam_No=sam;

		if ( solveroot( RS, net, sam ,samwork) != GSL_SUCCESS)
			continue;

		root_work=0;image=0;image_whole=0;mutiInterval=0;
		for(rn=0;rn<RS[samwork][0].root_Sum;rn++){
			if(turing_criteria(RS,net,sam,samwork,rn)==GSL_SUCCESS)
			{
				Sample[sam].workspace[Sample[sam].work++]=net;root_work++; // exists steady that no instabe kc exists
				if(RS[samwork][rn].intervalReal>1){
					Sample[sam].multispace[Sample[sam].multiInterval++]=net;mutiInterval=1;
				}
				for(i=0;i<RS[samwork][rn].intervalReal;i++)
					if(RS[samwork][rn].RE[i].intervalImage>0){
						Sample[sam].imagespace[0][Sample[sam].image[0]++]=net;image=1;
						if(RS[samwork][rn].RE[i].intervalImage==1 && RS[samwork][rn].RE[i].IM[0].k_minIM<0 && RS[samwork][rn].RE[i].IM[0].k_maxIM<0){
							Sample[sam].imagespace[1][Sample[sam].image[1]++]=net;image_whole=1;
						}
					}
			}
		}
		if(root_work>0){
			samwork++;
			Net[net].q_value++;
			Net[net].image[0]+=image;
			Net[net].image[1]+=image_whole;
			Net[net].multiInterval+=mutiInterval;
		}
	}
	if(samwork>0){
		PosiQSum++;
		for(sam=0;sam<samwork;++sam)
		{
			if(Sample[RS[sam][0].sam_No].D[0]<Sample[RS[sam][0].sam_No].D[1] && Sample[RS[sam][0].sam_No].D[1]<Sample[RS[sam][0].sam_No].D[2]) Net[net].role[0]++;
			else if(Sample[RS[sam][0].sam_No].D[0]<Sample[RS[sam][0].sam_No].D[2] && Sample[RS[sam][0].sam_No].D[2]<Sample[RS[sam][0].sam_No].D[1]) Net[net].role[1]++;
			else if(Sample[RS[sam][0].sam_No].D[1]<Sample[RS[sam][0].sam_No].D[0] && Sample[RS[sam][0].sam_No].D[0]<Sample[RS[sam][0].sam_No].D[2]) Net[net].role[2]++;
			else if(Sample[RS[sam][0].sam_No].D[1]<Sample[RS[sam][0].sam_No].D[2] && Sample[RS[sam][0].sam_No].D[2]<Sample[RS[sam][0].sam_No].D[0]) Net[net].role[3]++;
			else if(Sample[RS[sam][0].sam_No].D[2]<Sample[RS[sam][0].sam_No].D[0] && Sample[RS[sam][0].sam_No].D[0]<Sample[RS[sam][0].sam_No].D[1]) Net[net].role[4]++;
			else if(Sample[RS[sam][0].sam_No].D[2]<Sample[RS[sam][0].sam_No].D[1] && Sample[RS[sam][0].sam_No].D[1]<Sample[RS[sam][0].sam_No].D[0]) Net[net].role[5]++;
		}
		Net[net].q_value=samwork;
		text_paraWork(RS,net);
	}
	Net[net].test_if=1;
}

