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

int spache_seed(node *atoms, int i)
//check occupation
{
	int j = 0, bools = 0, k = 0, l = 0;//0 is empty 1 is filled
	double dist;
	for (j = 0; j<i; j++)
	{
		//if(x[j]!=x[j]){continue;}
		dist = SQR(atoms[i].x - atoms[j].x) + SQR(atoms[i].y - atoms[j].y) + SQR(atoms[i].z - atoms[j].z);
		if (sqrt(dist)<A_CUT)//filled if within sqrt(3) units of another residue
		{
			//printf("sep = %5.3g, j = %d\n",dist,j);
			bools = 1;
			break;
		}
	}
	//if(bools==1){printf("\nSpace occupied...\tRetrying...");}
	return bools;
}





void ratmpos_seed(int i, int j, node *atoms, int *seed)
//Puts atom i at a random position a distance rdist from atom j
{
	int u = 1;
	double rdist;
	double theta, phi;

	do
	{
		rdist = urand(seed, A_CUT, M_CUT);
		phi = urand(seed, 0, 2 * M_PI);
		theta = acos(1 - 2 * ran2(seed));
		atoms[i].x = atoms[j].x + rdist*sin(theta)*cos(phi);

		atoms[i].y = atoms[j].y + rdist*sin(theta)*sin(phi);

		atoms[i].z = atoms[j].z + rdist*cos(theta);

		u = spache_seed(atoms, i);
	} while (u == 1);
}




void connections(int **indexs, int N, node *atoms)
//stores their indices of connected nodes in index
//no connections to ligand sites recorded
//assumes ligands at end
{
	int i = 0, count, a = 0;
	int size, j = 0;
	double dist;
	for (a = 0; a<N; a++)
	{
		count = 1;
		//if(ainb("HETATM",atm[a])==0){continue;}//stops connection to ligands
		//printf("\n%i\t",a);
		indexs[a] = realloc(indexs[a], N * sizeof(int));
		for (i = 0; i<N; i++)
		{

			//if(i==a || ainb("HETATM",atm[a])==0){continue;}//stops connection to self or ligands
			if (i == a) { continue; }
			dist = SQR(atoms[i].x - atoms[a].x) + SQR(atoms[i].y - atoms[a].y) + SQR(atoms[i].z - atoms[a].z);
			if (sqrt(dist)<M_CUT)
			{
				indexs[a][count] = i;
				//printf("%i\t",indexs[a][count]);
				count++;
			}
		}
		indexs[a][0] = count;
		indexs[a] = realloc(indexs[a], (count) * sizeof(int));
		//printf("count= %d\n",count);
	}
}





void make_seed_b(FILE *runs, char dir[], int N, double sdist, double G[5], node *atoms, int **indexs)
//Generates a model protein formed from N nodes joined by a backbone.
{
    printf("\n----------------------\nGenerating Seed Structure\n");

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
	atoms[i].atnum = N+3;
	atoms[i].resnum = N+2;
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
	atoms[i].atnum = N+4;
	atoms[i].resnum = N+3;
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
    //writing atoms
	for (i = 2; i < N + 2; i++)
	{
		atoms[i].atnum = i - 1;
		atoms[i].resnum = i - 1;
		strncpy(atoms[i].res, "TOY", 3);
		strncpy(atoms[i].atm, "ATOM", 4);
		strncpy(atoms[i].name, "CA", 2);
		strncpy(atoms[i].chain, "A", 1);
		strncpy(atoms[i].elem, "C", 1);
		atoms[i].occ = 1.00;
		atoms[i].bfac = rbfac(&seed);
		ratmpos(i, i - 1, atoms, &seed);
	}

	connections(indexs, N, atoms);

	write_pdb(fname,G,N,atoms);
	sprintf(path,"%s/0.pdb",dir);
	copy("h2.pdb",path);

	fprintf(runs,"%-6d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",d,G[0],G[1],G[2],G[3],G[4]);//calculate allosteric free energy and save in seperate file
	printf("\nSeed structure complete | %d attempts required\n----------------\n",count);
}




int get_opt(char flag[], char* argl[], int argc, char args[])
//returns 0 if flag found 1 if not and -1 if overflow
{
	int i = 1, ret = 1;
	while (ret == 1 && i<argc)
	{
		if (ainb(flag, argl[i]) == 0)
		{
			if (i<(argc - 1)) { strcpy(args, argl[i + 1]); ret = 0; }
			else { ret = -1; }
			break;
		}
		i++;
	}
	return ret;
}





int main(int argc, char *argl[])
{
	int N, **indexs;
	double sdist, G[5];
	int opt, k, i;
	char cmd[80], path[40];

	//Getting number of nodes
	opt = get_opt("-N", argl, argc, cmd);
	if (opt == 0)
	{
		N = (int)strtod(cmd, NULL);
	}

	if (opt == 1) 
	{
		printf("No -N flag found");
	}

	//Getting ligand separation distance
	opt = get_opt("-sdist", argl, argc, cmd);
	if (opt == 0) 
	{
		sdist = (int)strtod(cmd, NULL);
	}
	if (opt == 1)
	{
		printf("No -sdist flag found");
	}

	FILE *fp;  //File pointer to DDG.txt file
	sprintf(path, "test_seed");
	mkdir(path, S_IRWXU);  //Making directory for storing pdb and DDG.txt
	fp = fopen("test_seed/DDG.txt", 'w');  //Creating DDG.txt file
	fprintf(fp, "#         G0          G11          G12          G2           DDG\n");

	//make dynamic arrays
	indexs = malloc((N+2) * sizeof(int*));
	node *atoms;
	atoms = calloc((N+2), sizeof(node));
	if (atoms == NULL || indexs == NULL)
	{
		printf("\nMemory could not be allocated for line data array\nSTOP");
		exit(0);
	}

	for (k = 0; k<N+2; k++)
	{
		indexs[k] = malloc((N + 3) * sizeof(int));
		if (indexs[k] == NULL)
		{
			printf("\nMemory could not be allocated for line data arrays\nSTOP");
			exit(0);
		}
	}

	make_seed_b(fp, path, N, sdist, G, atoms, indexs);

	fclose(fp);
	free(atoms);
	for (i = 0; i<N; i++)
	{
		free(indexs[i]);
	}
	free(indexs);
}