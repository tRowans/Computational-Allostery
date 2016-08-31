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
//check file read/write and DDG calculation
//rotation and similarity

int main()
{
	static int seed;
	time_t t;
	int j,n=29,N=27;
	FILE *out;
	out = fopen("mtst.txt","w");
	ep p;
	seed = -(int)time(&t);
	/*node *atoms;
	atoms = calloc(n,sizeof(*atoms));*/
	/*p.N = N;
	p.atoms = calloc(N,sizeof(node));*/
	ep_alloc(&p,N);

	int **indexs;

    //make dynamic arrays
    indexs=malloc(N * sizeof *indexs);
    
    for(int k=0;k<N;k++)
    {
        indexs[k]=malloc((N+1) * sizeof *indexs[k]);
        if(indexs[k]==NULL)
        {
            printf("\nMemory could not be allocated for line data arrays\nSTOP");
            exit(0);
        }
    }
	printf("trialing seeding");

	make_seed(out,"tst",N,11.0,p.G,p.atoms,indexs);
	for(j=0;j<N;j++){free(indexs[j]);}
	free(indexs);
	ep_free(&p,N);
	fclose(out);
	return 0;
}
