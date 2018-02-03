#include"header_file.h"


void create_folder(){

	int makeif=0;

	while(makeif==0)
	{

        char path[PATH_LEN];
        char cluster[PATH_LEN];
        char topo[PATH_LEN];
        char state[PATH_LEN];
        char foldername[PATH_LEN]=FolderName;

        getcwd(path,sizeof(path));
        sprintf(Catalog_SamplePath,"%s/Sample.txt",path);

        sprintf(path,"%s/%s_%lu-%lu",path,foldername,(long unsigned int)Find_Start,(long unsigned int)Find_End);

        if(access(path,0)==-1){//access函数是查看文件是不是存在

		mkdir(path,0777);
		sprintf(cluster,"%s/Cluster",path);
		sprintf(topo,"%s/Topo",path);
		sprintf(state,"%s/State",path);
		mkdir(cluster,0777);
		mkdir(topo,0777);
		mkdir(state,0777);
		strcpy(Catalog_Path,path);
		makeif=1;
	}
        else
            exit(1);
	}
}

void seq2topo(long unsigned int n,char topo[]){

	int i;

	for ( i = 1 ; i <= EdgeN ; i++ )
	{
		if( Net[ n ].topology[ i ]==1)
			topo[i-1]='a';
		else if(Net[ n ].topology[ i ]==-1)
			topo[i-1]='i';
		else if(Net[ n ].topology[ i ]==0)
			topo[i-1]='n';
	}
	topo[EdgeN]='\0';
}


void topo2seq(long unsigned int*n,char topo[]){

	int i,x,y=0;
	for(i=EdgeN-1;i>=0;i--)
	{
		if (topo[i]=='a')
			x=2;
		else if (topo[i]=='n')
			x=1;
		else if (topo[i]=='i')
			x=0;
		else
		{
			printf("\n invalid topology input! \n");
			return;
		}
		y=y*3+x;
	}
	*n=y;
}



void seq2topo_divide(long unsigned int n, char chip[]){

	int i,j;

	for (i=0,j=0;j<EdgeN;i++,j++)
	{
		if (j%NodeN==0&&j!=0)
			chip[i++]=',';
		if(Net[n].topology[j+1]==1)
			chip[i]='a';
		else if(Net[n].topology[j+1]==-1)
			chip[i]='i';
		else
			chip[i]='n';
	}
	chip[i]='\0';
}

void q_sort(){

	long unsigned int i,j,temp;

	for(i=0;i<NetN;i++)
		Q_Order[i]=i;

	for(i=0;i<NetN;i++){
		for(j=0;j<NetN-i-1;j++)
			if(Net[Q_Order[j]].q_value<Net[Q_Order[j+1]].q_value)
			{
				temp=Q_Order[j];
				Q_Order[j]=Q_Order[j+1];
				Q_Order[j+1]=temp;
			}
	}
	for(i=0;i<NetN;i++)
		Net[Q_Order[i]].order=i;
}

void q_sortSam(){

	long unsigned int i,j,temp;

	for(i=0;i<SampleN;i++)
		Sam_Order[i]=i;

	for(i=0;i<SampleN;i++){
		for(j=0;j<SampleN-i-1;j++)
			if(Sample[Sam_Order[j]].work<Sample[Sam_Order[j+1]].work)
			{
				temp=Sam_Order[j];
				Sam_Order[j]=Sam_Order[j+1];
				Sam_Order[j+1]=temp;
			}
	}
}


double indeX(double x, int y){

	double result=1.0;

	if (y>0)
		do{
			result*=x;
		}while(--y>0);

	else if (y<0)
		do{
			result/=x;
		}while(++y<0);

	return result;
}
