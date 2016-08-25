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


//depenedencies: all

void make_seed(FILE *runs,char dir[],int N,double sdist,double G[5],node *atoms,int **indexs)
//generate random structure and save to dir as d.pd,
{
    printf("\n----------------------\nGenerating seed Structure\n");

    double ac=0.5*sdist/sqrt(3);//stores site co-ords
    char path[50];
    int i=0,d=0,cyc=0,count=0;
    char fname[20]="h2.pdb";
    int s=time(NULL);
	static int seed;
	time_t t;
	seed = -(int)time(&t);
	printf("Writing ligands\n");
// Writing initial details
    //first ligand
    atoms[i].atnum = N+1;
    atoms[i].resnum = N;
	/*cpystr(atoms[i].atm,"HETATM",6);
    cpystr(atoms[i].res,"LIG",3);
    cpystr(atoms[i].name,"CA",2);
    cpystr(atoms[i].chain,"A",1);
    cpystr(atoms[i].elem,"C",1);*/
	strncpy(atoms[i].atm,"HETATM",6);
    strncpy(atoms[i].res,"LIG",3);
    strncpy(atoms[i].name,"CA",2);
    strncpy(atoms[i].chain,"A",1);
    strncpy(atoms[i].elem,"C",1);
    atoms[i].occ = 1.00;
    atoms[i].bfac = urand(&seed,15,20);
    atoms[i].x = -ac;
	atoms[i].y = -ac;
	atoms[i].z = -ac;
    i++;
    //second ligand
    atoms[i].atnum = N+2;
    atoms[i].resnum = N+1;
	/*cpystr(atoms[i].atm,"HETATM",6);
    cpystr(atoms[i].res,"LIG",3);
    cpystr(atoms[i].name,"CA",2);
    cpystr(atoms[i].chain,"A",1);
    cpystr(atoms[i].elem,"C",1);*/
	strncpy(atoms[i].atm,"HETATM",6);
    strncpy(atoms[i].res,"LIG",3);
    strncpy(atoms[i].name,"CA",2);
    strncpy(atoms[i].chain,"A",1);
    strncpy(atoms[i].elem,"C",1);
    atoms[i].occ = 1.00;
    atoms[i].bfac = urand(&seed,15,20);
    atoms[i].x = ac;
	atoms[i].y = ac;
	atoms[i].z = ac;
	printf("Ligands written\n");
	i++;
    
    do{
		atoms[i].bfac = rbfac(&seed);
    	atoms[i].atnum = 1;
    	atoms[i].resnum = 1;
		strncpy(atoms[i].res,"TOY",3);
		strncpy(atoms[i].atm,"ATOM",4);
    	strncpy(atoms[i].name,"CA",2);
    	strncpy(atoms[i].chain,"A",1);
    	strncpy(atoms[i].elem,"C",1);
    	atoms[i].occ = 1.00;
		ratmpos(i,0,atoms,&seed);
		i++;
		atoms[i].bfac = rbfac(&seed);
    	atoms[i].atnum = 2;
    	atoms[i].resnum = 2;
		strncpy(atoms[i].res,"TOY",3);
		strncpy(atoms[i].atm,"ATOM",4);
    	strncpy(atoms[i].name,"CA",2);
    	strncpy(atoms[i].chain,"A",1);
    	strncpy(atoms[i].elem,"C",1);
    	atoms[i].occ = 1.00;
		ratmpos(i,1,atoms,&seed);
		for(i=4;i<(N);i++)//allocates details
        {
            ratmpos(i,-1,atoms,&seed);
            //gen_data(atoms,i,seed);
			atoms[i].bfac = rbfac(&seed);
    		atoms[i].atnum = atoms[i-1].atnum + 1;
    		atoms[i].resnum = atoms[i-1].resnum + 1;
			atoms[i].occ = 1.00;

			strncpy(atoms[i].res,"TOY",3);
			strncpy(atoms[i].atm,"ATOM",4);
    		strncpy(atoms[i].name,"CA",2);
    		strncpy(atoms[i].chain,"A",1);
    		strncpy(atoms[i].elem,"C",1);
			//printf("%s %s %d ",atoms[i].res,atoms[i].atm,atoms[i].atnum);
			//printf("%s %s\n",atoms[i].name,atoms[i].chain);
        }
        connections(indexs,N,atoms);
        cyc = cycle(indexs,N,atoms);
		printf("\nconnectivity check %d\n",cyc);
		i = 2; count++;
    }while(cyc!=1);

    write_pdb(fname,G,N,atoms);
    sprintf(path,"%s/0.pdb",dir);
    copy("h2.pdb",path);

    fprintf(runs,"%-6d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",d,G[0],G[1],G[2],G[3],G[4]);//calculate allosteric free energy and save in seperate file
    printf("\nSeed structure complete | %d attempts required\n----------------\n",count);

}

void samples(char acc[],int n,int N,double sdist,int het)
{
    double xyzo[3],G0,G[5];
    int d=1,k,i=0,pa;
    char patha[60],pathb[99];
    time_t t0,t1;
    int d_t,tip=420;
	static int seed;
	seed = -(int)time(&t0);

    FILE *runsa;
    mkdir(acc,S_IRWXU); //make directory for storing sample
	sprintf(pathb,"%ssampleruns.txt",acc);
    runsa=fopen(pathb,"w");
    if(runsa==NULL)//checking for success
    {
        printf("\nRun data storage files could not be created\nSTOP");
        exit(0);
    }
    fprintf(runsa,"#         G0          G11          G12          G2           DDG\n");

	sprintf(pathb,"%sstructs/",acc); //make directory for storing structs in sample
	mkdir(pathb,S_IRWXU);
	
    t0=time(NULL);


    int **indexs;

    //make dynamic arrays
    indexs=malloc(N * sizeof(int*));
    node *atoms;
    atoms=calloc(N,sizeof(node));
    if(atoms==NULL || indexs==NULL)
    {
        printf("\nMemory could not be allocated for line data array\nSTOP");
        exit(0);
    }

    for(k=0;k<N;k++)
    {
        indexs[k]=malloc((N+1)*sizeof(int));
        if(indexs[k]==NULL)
        {
            printf("\nMemory could not be allocated for line data arrays\nSTOP");
            exit(0);
        }
    }


    make_seed(runsa,pathb,N,sdist,G,atoms,indexs);

    k=0;
    connections(indexs,N,atoms);

    while(d<n)
    {
        k=HandOfGod(atoms,N,het,xyzo,indexs,&seed);
        if(k>N){exit(0);}
        sprintf(patha,"%s%d.pdb",pathb,d);
        write_pdb(patha,G,N,atoms);
        fprintf(runsa,"%-6d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",d,G[0],G[1],G[2],G[3],G[4]);
        /*pa = aprob(G0,G[4],1,co,seed);//check for acceptance
        if(pa!=0)
        {
            atoms[k].x = xyzo[0];
            atoms[k].y = xyzo[1];
            atoms[k].z = xyzo[2];
        }*/
        d+=1;
    }

    printf("\n-------------------------------------------\nsample Complete\n");
    fclose(runsa);
    free(atoms);
    for(i=0;i<N;i++)
    {
        free(indexs[i]);
    }
    free(indexs);
}

int monte(int N,double sdist,int het,int co,int *iters,char sdir[])
{
    double T,cool=0.95,sig,cdist,dist,msqrt;
    double xyzo[3],G0,G[5],Nm[3],Gmax,cdistbest;
    int d=1,io=0,k,j=0,pa,count=0,i=0,best;
    char patha[60],acc[20]="runs";
    time_t t0,t1;
    int d_t,tip=420;
	static int seed;
 	seed = -(int)time(&t0);

    FILE *runsa;
    mkdir(acc,S_IRWXU);
    runsa = fopen("runs.txt","w");
    if(runsa==NULL)//checking for success
    {
        printf("\nRun data storage files could not be created\nSTOP");
        exit(0);
    }
    fprintf(runsa,"#         G0          G11          G12          G2           DDG         T        Msqr    C.O.M   LigDist\n");

    t0=time(NULL);

    int **indexs;

    //make dynamic arrays
    indexs = malloc(N * sizeof(*indexs));

    //if(x==NULL || y==NULL || z==NULL  || occ==NULL || bfac==NULL || atnum==NULL ||  resnum==NULL ||  atm==NULL ||  name==NULL ||  res==NULL ||  elem==NULL ||  chain==NULL || indexs==NULL)
    if(indexs==NULL)
    {
        printf("\nMemory could not be allocated for line data array\nSTOP");
        exit(0);
    }

    for(k=0;k<N;k++)
    {
        indexs[k] = malloc((N+1) * sizeof(*indexs[k]));
        //if(atm[k]==NULL || name[k]==NULL || res[k]==NULL  || elem[k]==NULL  || chain[k]==NULL || indexs[k]==NULL)
        if(indexs[k]==NULL)
        {
            printf("\nMemory could not be allocated for line data arrays\nSTOP");
            exit(0);
        }
    }

    node *atoms;
    atoms = calloc(N,sizeof(*atoms));
    //OPTION FOR SEED FROM SAMPLE
    make_seed(runsa,"runs",N,sdist,G,atoms,indexs);
    enms(j);j++;
    G0 = G[4];
    Gmax = G0;
	sprintf(patha,"%ssampleruns.txt",sdir);
    T = Temp(patha);
    Nm[0]=0;Nm[1]=0;Nm[2]=0;
    connections(indexs,N,atoms);

    while(T>fabs(Gmax*5e-4))
    {
        io = HandOfGod(atoms,N,het,xyzo,indexs,&seed);
        write_pdb("h2.pdb",G,N,atoms);
        pa = aprob(G0,G[4],T,co,&seed);//check for acceptance

        if(pa!=0)
        {
            atoms[io].x = xyzo[0];
            atoms[io].y = xyzo[1];
            atoms[io].z = xyzo[2];
            if(count>50){printf("\nTerminated after %d unsuccessful move attempts\n",count);break;}
            count++;
        }
        else
        {
            sprintf(patha,"runs/%d.pdb",j);//check for folder
            copy("h2.pdb",patha);
            cdist = m_centre(atoms,N);
            dist = lig_dist(atoms,N);
            fprintf(runsa,"%-6d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %8.3lf %8.3lf\n",j,G[0],G[1],G[2],G[3],G[4],T,msqrt,cdist,dist);
            if((co*G[4])>(co*Gmax))
            {
                Gmax = G[4];
                *iters = j;
            }
            enms(j); j++;

            Nm[0] += co*(G[4]-G0);
            Nm[1] += SQR(G[4]-G0);
            Nm[2] += 1;
            sig = sqrt((Nm[1]-(SQR(Nm[0])/Nm[2]))/(Nm[2]-1));//variance calculation
            msqrt = sqrt(Nm[1])/Nm[2];
            if(Nm[2]>=10 && msqrt<T)//when 10 or more accepted
            {
                Nm[0]=0;
                Nm[1]=0;
                Nm[2]=0;
                T*=cool;
                if(sig==0){break;}
            }
            G0=G[4];
            count=0;
        }
    }


    printf("\n-------------------------------------------\nRun Complete\n");
    fclose(runsa);
    free(atoms);
    for(i=0;i<N;i++)
    {
        free(indexs[i]);
    }
    free(indexs);
    return j;
}


