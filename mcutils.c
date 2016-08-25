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

// dependencies: sysutils,sest,lib,rwutils

//--------------------------------------------------------
int aprob(double G0,double G1,double T,int i,int *seed)
//move=G0->G1:returns 0 if accepted 1 if rejected
//i is -1 if negative DDG and 1 for positive
{
    int ret=0;
    double eg,x;
    eg=exp(i*(G1-G0)/T);//maximise DDG
    if(eg<1)
    {
        x = urand(seed,0,1);
        ret = floor(x/eg);//ret=1 if x>eg 0 otherwise
    }
    return ret;
}

double Temp(char filename[])
//calculates temperature for metropolis criterion from data in sample runs file fp
{
    int i=0,N=0,j=0,c;
    char lign[72],hold[13],*spit;
    double G,Nmu2=0,Nmu=0,sig;
    FILE* fp;
    fp = fopen(filename,"r");//check existence
    if(fp==NULL)
    {
        printf("\nSample file: %s \nDoes not exist\nSTOP\n",filename);
        exit(0);
    }
    printf("\n--------------------------\nCalculating Temperature\nReading Sample file\n");
    while(fgets(lign,sizeof lign,fp)!=NULL)
    {
        fgets(lign,sizeof lign,fp);
        if(j==0){j++;continue;}
        for(c=0;c<12;c++)
        {
            hold[c] = lign[59+c];
        }
        G = strtod(hold,&spit);
        N++; Nmu+=G; Nmu2+=SQR(G);
        printf("\nN=%d\tNmu= %lf\tNmu2= %lf",N,Nmu,Nmu2);
    }
    fclose(fp);

    printf("\n");
    sig = sqrt((Nmu2-(SQR(Nmu)/N))/(N-1));
    printf("\nTemperature is: %lf\n--------------------------\n",sig);
    return sig;
}


//Moves//

int HandOfGod(node *atoms,int N,int het,double xyzo[3],int **indexs,int *seed)
//N is number of nodes, calculates montecarlo move.
//het is 0 if fixed ligda site distance & 1 for variable
//moves a random node
{
    //if(con>2.5e4){return con;}
    int cont=0,con=0;
    int rn, an;
    int i=1,hetc=0,j;
    double rdist;
    rn=floor(urand(seed,0,N));
    while(/*con<=(1e5) &&*/cont==0)
    {
        if((con/50)>=1){
			con=0;
			rn = floor(urand(seed,0,N));
		}
        hetc = 0;
        xyzo[0] = atoms[rn].x;
        xyzo[1] = atoms[rn].y;
        xyzo[2] = atoms[rn].z;
        i = ainb("HETATM",atoms[rn].atm);
        if(i==0 && het==0){rn = floor(urand(seed,0,N)); continue;}//if fixed ligands and ligand selected for move
        an=floor(urand(seed,0,N));
        if(rn!=an)//selected nodes okay
        {
            rdist=urand(seed,2,8);
            if(i==1){
                ratmpos(rn,an,atoms,seed);
            }
            else if(i==0){
                ratmpos(rn,rn,atoms,seed);
            }

            if(spache(atoms,rn,N)==1)//test physical space (x,y,z) unoccupied
            {
                cont = 0;
                atoms[rn].x = xyzo[0];
				atoms[rn].y = xyzo[1];
				atoms[rn].z = xyzo[2];
                printf("initial | x= %5.2lf | y= %5.2lf | z= %5.2lf\n",atoms[rn].x,atoms[rn].y,atoms[rn].z);
                continue;
            }

            connections(indexs,N,atoms);

            if(cycle(indexs,N,atoms)!=1)
            //if move disconnects
            {
                //printf("\nDisconnected... Retry...");
                cont = 0;
                con++;
                atoms[rn].x = xyzo[0];
				atoms[rn].y = xyzo[1];
				atoms[rn].z = xyzo[2];
                printf("\ntrying cyc %d\n",con);
                rn = floor(urand(seed,0,N));
                continue;
            }
            else{
                printf("\nmove found %d\n",con);
				cont = 1;
            }
        }
        else
        //selected nodes not okay
        {
            cont = 0;
        }
    }
    if(cont==0){printf("\nNo move found\n");return con;}
    return rn;
}

//--------------------------
