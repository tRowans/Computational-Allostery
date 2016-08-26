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

//--------------------------------

//seed utils
//dependencies: lib, rwutils

int spache_seed(node *atoms, int i)
//check occupation
{
    int j=0,bools=0,k=0,l=0;//0 is empty 1 is filled
    double dist;
    for(j=0;j<i;j++)
    {
        //if(x[j]!=x[j]){continue;}
        dist = SQR(atoms[i].x - atoms[j].x) + SQR(atoms[i].y - atoms[j].y) + SQR(atoms[i].z - atoms[j].z);
        if(sqrt(dist)<A_CUT)//filled if within sqrt(3) units of another residue
        {
			//printf("sep = %5.3g, j = %d\n",dist,j);
            bools=1;
            break;
        }
    }
    //if(bools==1){printf("\nSpace occupied...\tRetrying...");}
    return bools;
}


void ratmpos(int i,int j,node *atoms,int *seed)
//generate random atom position
{
    int u=1,k=0;
    
    double rdist;
    double theta,phi;
	//printf("j = %d, rad = %5.3g\n",j,rdist);
    do
    {
		if(j==-1){j=floor(urand(seed,0,i));}
		rdist = urand(seed,A_CUT,M_CUT);
    	phi = urand(seed,0,2*M_PI);
    	theta = acos(1- 2*ran2(seed));
    	atoms[i].x = atoms[j].x + rdist*sin(theta)*cos(phi);
    	//nrdist(rdist,atoms[i].x);

    	atoms[i].y = atoms[j].y + rdist*sin(theta)*sin(phi);
    	//nrdist(rdist,atoms[i].y);

    	atoms[i].z = atoms[j].z + rdist*cos(theta);

   		u = spache_seed(atoms,i);
		j = -1;
    } while(u==1);
	//printf("u = %d\n",u);
}

double rbfac(int *seed)
//random bfactor
{
    return urand(seed,20,70);
}

void gen_data(node *atoms,int i,int *seed)
//generates data for the new atom
{
    //getres(atoms[i].res);
    
    atoms[i].bfac = rbfac(seed);
    atoms[i].atnum = atoms[i-1].atnum + 1;
    atoms[i].resnum = atoms[i-1].resnum + 1;
	/*cpystr(atoms[i].res,"TOY",3);
	cpystr(&(atoms[i].atm),"ATOM",4);
    cpystr(&(atoms[i].name),"CA",2);
    cpystr(&(atoms[i].chain),"A",1);
    cpystr(&(atoms[i].elem),"C",1);*/
    atoms[i].occ = 1.00;

	strncpy(atoms[i].res,"TOY",3);
	strncpy(atoms[i].atm,"ATOM",4);
    strncpy(atoms[i].name,"CA",2);
    strncpy(atoms[i].chain,"A",1);
    strncpy(atoms[i].elem,"C",1);
	

}

//----------------------------

//constraint utils
//dependeancies: random utils, read utils

int spache(node *atoms, int i,int N)
//check occupation of x,y,z
{
    int j,bools=0,het;//0 is empty 1 is filled
    double dist;
    if(ainb("HETATM",atoms[i].atm)==0)
    {
        het = 0;
    }
    for(j=0;j<N;j++)
    {
        if(i==j){continue;}
        //dist=SQR(x[i]-x[j])+SQR(y[i]-y[j])+SQR(z[i]-z[j]);
        dist = rad(atoms[i].x-atoms[j].x,atoms[i].y-atoms[j].y,atoms[i].z-atoms[j].z);
        if(het==0 && ainb("HETATM",atoms[j].atm)==0 && dist<=M_CUT)
        {
            bools = 1;
            printf("\nLigands too close\n");
            break;
        }
        if(dist<A_CUT)//filled if within cut units of another residue
        {
            bools = 1;
            printf("\nspache %5lf\n",dist);
            printf("j= %4d | x= %5.2lf | y= %5.2lf | z= %5.2lf\n",j,atoms[j].x,atoms[j].y,atoms[j].z);
            printf("rn= %4d | x= %5.2lf | y= %5.2lf | z= %5.2lf\n",i,atoms[i].x,atoms[i].y,atoms[i].z);
            break;
        }
    }
    //if(bools==1){printf("\nSpace occupied...\tRetrying...");}
    return bools;
}

void connections(int **indexs,int N,node *atoms)
//Loops through all nodes finding which nodes are within the cutoff for connection and stores the results in indexs
//Does not connect backbone nodes
{
    int i=0,count,a=0;
    int size,j=0;
    double dist;
    for(a=0;a<N;a++)
    {
        count = 1;
        for(i=0;i<N;i++)
        {
            if(i==a){continue;}  //Stops connection to self
			if (a > 2 && fabs(i - a) == 1) { continue; } //Stops connection to neighbouring backbone nodes
			if (a == 2 && (i - a) == 1) { continue; } //Stops connection to neighbouring backbone node for first node in chain
            dist = SQR(atoms[i].x - atoms[a].x) + SQR(atoms[i].y - atoms[a].y) + SQR(atoms[i].z - atoms[a].z);
            if(sqrt(dist)<M_CUT)
            {
                indexs[a][count]=i;
                count++;
            }
        }
        indexs[a][0]=count;
		indexs[a] = realloc(indexs[a],(count)*sizeof(int));
    }
}

int cycle(int **index,int N,node *atoms)
// returns 1 if structure is connected
{
	int **visits, *cycs;
	int i,j,k,cmax,c,ev,ab,ret;
	int ni,nj,nk,count;
	ev = 0;
	c = 0;
	cmax = 0;
	i = 0;
	ret = 1;
	count = 0;
	visits = malloc(2*sizeof(*visits));
	visits[0] = calloc(N,sizeof(int));
	visits[1] = malloc(N*sizeof(int));
	cycs = calloc(N,sizeof(int));
	do{
		//i lig check
		if(visits[0][i]==0){
			ab = ainb("HETATM",atoms[i].atm);
			if(ab==0){visits[0][i]=1; continue;}
			cycs[i] = c;
			visits[1][ev] = i;
			visits[0][i] = 1;
			ev++;
			for(nj=1;nj<index[i][0];nj++){
				j = index[i][nj];
				ab = ainb("HETATM",atoms[j].atm);
				if(ab==0){visits[0][j]=1; continue;}
				//j lig check
				if(visits[0][j]==0){
					i = j;
					break;
				}
				else{
					if(cycs[j]<cycs[i]){
						for(nk=0;nk<ev;nk++){
							k = visits[1][nk];
							if(cycs[k]==cycs[i]){
								cycs[k] = cycs[j];
							}
						}
						if(c>cmax){
							cmax = c;
						}
						c = cycs[j];
					}
					else if(cycs[j]>cycs[i]){
						for(nk=0;nk<ev;nk++){
							k = visits[1][nk];
							if(cycs[k]==cycs[j]){
								cycs[k] = cycs[i];
							}
						}
					}
				}
			}
			count++;
		}
		else{
			for(i=0;i<N;i++){
				if(visits[0][i]==0){
					if(count!=0){c = cmax + 1;}
					break;
				}
			}
		}
	}while(ev<(N-2));

	for(i=0;i<N;i++){
		//i lig check
		ab = ainb("HETATM",atoms[i].atm);
		if(ab==0 && index[i][0]<=1){ret=0; break;}
		else if(ab!=0 && cycs[i]!=0){
			ret=0; break;
		}
	}

	free(cycs);
	free(visits[0]);
	free(visits[1]);
	free(visits);
	return ret;
}


//------------------------------------

int ep_alloc(ep *p,int n)
{
    p->atoms=calloc(n,sizeof(node));
    if(p->atoms==NULL){return 1;}
    p->N=n;
    return 0;
}

int ep_free(ep *p,int n)
{
	if(p->atoms != NULL){	
		free(p->atoms);
	}
	return 0;
}

int ep_copy(ep *p1,ep *p2)
//copies p2 into p1 (overwrites & resizes where necessary)
{
    if(p2->N!=p1->N){p1->atoms=realloc(p1->atoms,sizeof(*(p2->atoms)));}
    if(p1->atoms==NULL){return 1;}
    memcpy(p1->atoms,p2->atoms,sizeof(*(p2->atoms)));

    p1->N = p2->N;p1->iters=p2->iters;//p1->T=p2->T;
	memcpy(p1->G,p2->G,sizeof(p2->G));
    return 0;
}

int ep_trans_o(ep *p,double u[3])
{
    for(int i=0;i<p->N;i++)
    {
        p->atoms[i].x-=u[0];
        p->atoms[i].y-=u[1];
        p->atoms[i].z-=u[2];
    }
    return 0;
}

int ep_trans_com(ep *p)
{
    double com[3]={0};
    for(int i=0;i<p->N;i++)
    {
        com[0]+=(p->atoms[i].x)/(double)p->N;
        com[1]+=(p->atoms[i].y)/(double)p->N;
        com[2]+=(p->atoms[i].z)/(double)p->N;
    }
    ep_trans_o(p,com);
    return 0;
}

