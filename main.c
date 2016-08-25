#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <time.h>
#include <ctype.h>
#include "lib.h"

//-----------------------------------------
int get_opt(char flag[],char* argl[],int argc,char args[])
//returns 0 if flag found 1 if not and -1 if overflow
{
    int i=1,ret=1;
    while(ret==1 && i<argc)
    {
		if(ainb(flag,argl[i])==0)
		{
			if(i<(argc-1)){strcpy(args,argl[i+1]);ret=0;}
			else{ret=-1;}
			break;
		}
		i++;
	}
	return ret;
}

int main(int argc,char* argl[])
{
    int i,N,het,co,n=0;
    int j,tot,opt,iters=0;
    double sdist;
    char path[40],cmd[80],nums[10],rns[10],pathi[90];
    char *arglist[8];

    opt=get_opt("-N",argl,argc,cmd);
    if(opt==0)
    {
        N=(int)strtod(cmd,NULL);
        printf("\n%d non-ligand nodes\n",N);
        if(N<=2){printf("\nNode number must be greater than 2\n");exit(0);}
    }
    else if(opt==-1){printf("\nNode number must be greater than 2\n");exit(0);}
    else
    {
        printf("\n-N flag not found.\nSTOP\n");
        exit(0);
    }

    opt=get_opt("-sd",argl,argc,cmd);
    if(opt==0)
    {
        sdist=strtod(cmd,NULL);
        printf("\nBinding site seperation = %lf\n",sdist);
        if(sdist<=8){printf("\nBinding site seperation must be greater than 8\n");exit(0);}
    }
    else if(opt==-1){printf("\nBinding site seperation must be greater than 8\n");exit(0);}
    else
    {
        printf("\n-sd flag not found.\nSTOP\n");
        exit(0);
    }

    if(get_opt("-var",argl,argc,cmd)==0)
    {
        printf("\nVariable binding site position\n");
        het=1;
    }
    else
    {
        printf("\n-var flag not found. Using fixed site default\n");
        het=0;
    }

    if(get_opt("-aco",argl,argc,cmd)==0)
    {
        printf("\nOptimising for anti-cooperativity\n");
        co=1;
    }
    else
    {
        printf("\n-aco flag not found. Optimising for cooperativity\n");
        co=-1;
    }

    opt=get_opt("-tot",argl,argc,cmd);
    if(opt==0)
    {
        tot=(int)strtod(cmd,NULL);
        printf("\n%d runs to be completed\n",tot);
    }
    else
    {
        tot=1;
    }
    if(tot<1){tot=1;}

    if(get_opt("-sam",argl,argc,cmd)==0)
    {
		n=(int)strtod(cmd,NULL);
        printf("\nCreating new sample of %d files\n",n);
		sprintf(path,"sample_N%d_h%d/",N,het);
        if(n>0){samples(path,n,N,sdist,het);}
        else{printf("\nSample size not found\nSTOP\n");exit(0);}
    }
    else
    {
        printf("\nUsing existing sample file\n");
    }



    FILE* fp;
	sprintf(pathi,"sample_N%d_h%d/",N,het);
	sprintf(path,"%ssampleruns.txt",pathi);
    fp=fopen(path,"r");//check existence
    if(fp==NULL)
    {
		
        printf("\nSample file could not be opened\nUse [-sam samplesize] option to create sample\n\n");
        exit(0);
    }
    else{fclose(fp);}

    char sol[60], path2[60],fname[60];
	sprintf(path,"%d_h%d_t%d/",co+1,het,tot);
    for(i=1;i<=tot;i++)
    //runs to times
    {
        mkdir(path,S_IRWXU);
        j = monte(N+2,sdist,het,co,&iters,pathi);
        //strcat(path,"/");

        sprintf(sol,"%ssol_%d_h%d_t%d_iter%d.pdb",path,co+1,het,i,iters);
        sprintf(fname,"runs/%d.pdb",iters);
        copy(fname,sol);

        sprintf(sol,"%siter%d.vmd",path,iters);
		sprintf(fname,"runs/%d.vmd",iters);
        copy(fname,sol);

		sprintf(path2,"%srun%d",path,i);
        cpdir("runs",path2);
		sprintf(path2,"%srun%d.txt",path,i);
        mv("runs.txt",path2);
		rmvdir("runs");
        //DIAGNOSTICS
    }

}

