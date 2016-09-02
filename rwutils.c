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

// File Read and writing
// dependencies: sysutils

void cpystr(char *dest[],char a[],int len)
//copies string a into pointer dest
//len=length of string a to be copied
//assumes dest is large enough
{
    //dest=realloc(dest,len * sizeof *dest);
    int i;
    for(i=0;i<len;i++)
    {
        (*dest)[i] = a[i];
    }
}

int ainb(char a[],char b[])
//returns 0 if str b contains a returns 1 otherwise
{
    int i=0,j0=-1,j1=0,count=0;
	int lena = strlen(a);
	int lenb = strlen(b);
    if(lenb<lena){return 1;}

    for(i=0;i<lenb;i++){
        if((b[i]==a[j1])&&(j1==j0+1))
        {
            j0=j1;
            j1++;
        }
        else
        {
            j0=-2;
            j1=0;
        }

        if((j1+1)==lena)//escape early if string is found
        {break;}
    }
    if((j1+1)==lena)//string found
    {return 0;}

    else{return 1;}
}

void extract(char lign[],int a,int b,char *dep[])
//extract string length b from lign and put in dep[]
//a is starting point in lign. len(dep)>=b
{
    int c;
    for(c=0;c<b;c++)
    {
        (*dep)[c] = lign[a+c];
    }
}


//write utils
//dependencies: system utils

double ddg(double k[])
//calculate allosteric free energy
//expects k to have at least 4 elements, first four are G0,G1a,G1b,G2
{
    double dg;
    dg=(k[3]-k[2])-(k[1]-k[0]);
    return dg;
}


int countline_pdb(FILE* fp)
//returns number of lines in open file fp
{
    int i=0,j=1,k;
	char buff[8],lign[99];
	printf("hi\n");
	//rewind(fp);
    while(fgets(buff,7,fp)!= NULL)
    {
		//printf("reading line %d\n",j); j++;      
		//fgets(buff,7,fp);
		
        if((ainb("ATOM",buff)==0)||(ainb("HETATM",buff)==0))
		{
			//if(i==0){k=j;}
			//printf("%s\n",buff);
			i++;
		}
		fgets(lign,sizeof lign,fp);//move to next line
    }
    rewind(fp); printf("\n%d Atoms found\n",i);
    return i;
}

void read_pdb(FILE* fp,node *atoms)
//reads file lines and stores in arrays
//make ligs first in list
{
    printf("\nReading pdb data\n");
    int i=0,j=1,fcheck;
    char *buff,*b1,lign[10],dump[99],altlock[2];
    fpos_t position;
    do
    {        
		fgetpos(fp,&position);buff=fgets(lign,7,fp);
		//fsetpos(fp,&position);b1=fgets(dump,99,fp);
        if((ainb("ATOM",lign)==0)||(ainb("HETATM",lign)==0))
        {
            fsetpos(fp,&position);
            fcheck=fscanf(fp,"%6s%4d%4s %3s %1s %4d %8lf %8lf %8lf %6lf %6lf %2s \n",atoms[i].atm,&atoms[i].atnum,atoms[i].name,atoms[i].res,atoms[i].chain,&atoms[i].resnum,&atoms[i].x,&atoms[i].y,&atoms[i].z,&atoms[i].occ,&atoms[i].bfac,atoms[i].elem);
            if(fcheck!=12)
            {
				printf("%-6s\n%-4d\n%-4s\n%-5s\n%-1s\n%-4d\n%-8lf\n%-8lf\n%-8lf\n%-6lf\n%-6lf\n%-2s\n",atoms[i].atm,atoms[i].atnum,atoms[i].name,atoms[i].res,atoms[i].chain,atoms[i].resnum,atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].occ,atoms[i].bfac,atoms[i].elem);
				printf("\nNot all data read from line %d\nSTOP\nfcheck=%d\n",j,fcheck);
				printf("\n%s\n",dump);
                exit(0);
            }
	
            i++;
        }
		else{fgets(dump,99,fp);}//move to next line
		j++;
    }while(buff!= NULL);

    printf("\n%d Atoms read\n",i);
}

double readen()
//read free energy from mode.energy file
{
    char lign[60],*test,*check;
    int k=0,c=0;
    double mod=0,buff;
    test=calloc(11,sizeof *test);
    FILE* fre;
    fre=fopen("mode.energy","r");
    if(fre==NULL)
    {
        printf("\nmode.energy file could not be read\nSTOP");
        exit(0);
    }
    printf("\nReading \"mode.energy\" file\n");
    do
    {
        check=fgets(lign,sizeof lign,fre);
        if(k<1)//skip first line
        {
            //printf("\ncheck line skip\n");
            k+=1;
            continue;
        }
        for(c=0;c<10;c++)
        {
            test[c]=lign[19+c];
        }
        //printf("\nk=%d\n",k);
        buff=strtod(test,NULL);
        mod+=buff;
        k+=1;
    }while(check!=NULL);
    fclose(fre);
    free(test);
    return mod;
}

void readen_init(char *path, double G[5])
//Reads G0, G11, G12, G2 and DDG for initial seed states from specified file (sampleruns.txt)
{
	char line[72], *buffer;
	int i=0, j=0, k=0, d=0;
	FILE* fp;

	fp = fopen(path, "r");
	if (fp == NULL)
	{
		printf("\nSample file: %s \nDoes not exist\nSTOP\n", path);
		exit(0);
	}

	buffer = (char *)calloc(50, sizeof(char));

	for (i = 0; i < 2; i++)
	{
		fgets(line, 72, fp);
		if (i == 0) { continue; }
		for (j = 0; j < 5; j++)
		{
			for (k = 0; k < 12; k++)
			{
				buffer[k] = line[7 + 13 * j + k];
			}
			G[j] = strtod(buffer, NULL);
		}
	}

	fclose(fp);
	
}

//-----------------------------------

void wlpdb(FILE* fname, node atom)
//write line data to pdb file
{
    fprintf(fname,"%-6s%5d  %-4s%3s %1s%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s\n",atom.atm,atom.atnum,atom.name,atom.res,atom.chain,atom.resnum,atom.x,atom.y,atom.z,atom.occ,atom.bfac,atom.elem);
}

void terch(FILE* fname,char chain[1],char res[3],int atnum, int resnum)
//write chain termination
{
    fprintf(fname,"TER   %5d      %3s %1s%4d\n",atnum,res,chain,resnum);
}

void write_apo(double G[5],int N,node *atoms)
//writes apo data to file fpath and returns G for that file
{
    int i=0,j=0,tmp;
    //double k[3];
    FILE *fp;
    fp=fopen("apo.pdb","w");

    if(fp==NULL)
    {
        printf("\npdb file could not be created\nSTOP");
        exit(0);
    }

    fprintf(fp,"HEADER    Autogenerated structure\n");
    //printf("\ntest\n");
    while(i<N)
    {
        if(ainb("HETATM",atoms[i].atm)!=0)//if not ligand
        {
			//if(ainb(atoms[i-1].chain,atoms[i].chain)!=0){terch(fp,atoms[j].chain,atoms[j].res,atoms[j].atnum +1,atoms[j].resnum);}
            wlpdb(fp,atoms[i]);
            j=i;
        }
        i++;
    }
    terch(fp,atoms[j].chain,atoms[j].res,atoms[j].atnum +1,atoms[j].resnum);
    fclose(fp);
    ddpt("apo.pdb");
    G[0]=readen();
}

//writes holo1 files and polices the binding order
//returns 1 if h11 used and 2 if h12 used
void write_h1(double G[5],int N,node *atoms)
{
    double Gh11,Gh12;
    int i=0,j=0,k=0;
    FILE *h11,*h12;
    copy("apo.pdb","h11.pdb");
    copy("apo.pdb","h12.pdb");
    h11=fopen("h11.pdb","a");
    h12=fopen("h12.pdb","a");
    if(h11==NULL || h12==NULL)//checking for success
    {
        printf("\npdb files could not be created\nSTOP");
        exit(0);
    }

    while(i<N)
    {
		k=ainb("HETATM",atoms[i].atm);
		if(k==0 && j==0)
        {
            wlpdb(h11,atoms[i]);
            //terch(h11,chain[i],res[i],atnum[i]+3,resnum[i]);
            fclose(h11);
            ddpt("h11.pdb");
            Gh11=readen();
            //printf("\ntest\n");
            j++;
			i++;
			continue;
        }
        else if(k==0 && j==1)
        {
            wlpdb(h12,atoms[i]);
            //terch(h12,chain[i],res[i],atnum[i]+3,resnum[i]);
            fclose(h12);
            ddpt("h12.pdb");
            Gh12=readen();
            //printf("\ntest\n");
            break;
        }
        i++;
    }

    G[1]=Gh11;
    G[2]=Gh12;
    //copy("h11.pdb","h1.pdb");
}

void write_pdb(char fname[],double G[5],int N,node *atoms)
//adds second binding site and writes output pdb file as fname
//calculates DDG and stores in G[4]
//assumes ligand entries are consecutive
{
    FILE *h2;
    printf("\n\nWriting structure to file\n--------------------\n");
    int i=0,j=0;
    write_apo(G,N,atoms);
    write_h1(G,N,atoms);
    copy("h11.pdb",fname);
    h2=fopen(fname,"a");
    if(h2==NULL)
    {
        printf("\npdb files could not be created\nSTOP");
        exit(0);
    }

    do
    {
        if(ainb("HETATM",atoms[i].atm)==0 && j<2)
        {			
			if(j==0){i++;j++;continue;}            
			wlpdb(h2,atoms[i]);
            terch(h2,atoms[i].chain,atoms[i].res,atoms[i].atnum +1,atoms[i].resnum);
            fclose(h2);
            ddpt(fname);
            G[3]=readen();
			

            break;
        }
        i++;
    }while(i<N);

    G[4]=ddg(G);
    printf("\nStructure Complete\n");
    printf("\nAllosteric free energy is: %lf\n--------------------\n",G[4]);
}

//--------------------------------
