#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <time.h>
#include <ctype.h>

#ifndef lib_h__
#define lib_h__

#define M_CUT 8.0
#define A_CUT 2.0
#define M_PI 3.14159265359

typedef struct
{
	int atnum;
    char atm[7];
    char name[6];
    char res[4];
    char chain[2];
	char elem[3];
    /*char *atm;
    char *name;
    char *res;
    char *chain;
	char *elem;*/
    int resnum;
    double x;
    double y;
    double z;
    double occ;
    double bfac;

}node;//holds node data

typedef struct
{
    node *atoms;
    int N;//no. of nodes
    double G[5];//ddpt outputs: {G0,G11,G12,G2,DDG}
    //double T;//temperature
    int iters;//no. of accepected steps to reach
    int ligs[2];
    //Diagnostics//
    int con[2];//lig connectivity: {l1,l2}
    double Is[3][3];//principal axes in columns of Is
    double lig_com[3];//lig COM
    double lig_ax[3];//lig axis

}ep;//stores system of nodes and it's diagnostics

double SQR(double x);
double rad(double x,double y,double z);

double urand(int *seed,double a, double b);
int rsign();
void nrdist(double rdist, double cord);
double ran2(int *idum);
//------------------------------------

//SYS utilities
int spawn(char* program,char** arg_list);
void copy(char path[],char dpath[]);
void mv(char path[],char dpath[]);
void ddpt(char fname[]);
void enms(int d);
void cpdir(char src[],char dest[]);
void rmvdir(char src[]);
//-------------

//RW utilities
void cpystr(char *dest[],char a[],int len);
int ainb(char a[],char b[]);
void extract(char lign[],int a,int b,char *dep[]);

double ddg(double k[]);

int countline_pdb(FILE* fp);
void read_pdb(FILE* fp,node *atoms);
void getres(char *resi);
double readen();
void readen_init(char *path, double G[5]);

void wlpdb(FILE* fname, node atom);
void terch(FILE* fname,char chain[1],char res[3],int atnum, int resnum);
void write_apo(double G[5],int N,node *atoms);
void write_h1(double G[5],int N,node *atoms);
void write_pdb(char fname[],double G[5],int N,node *atoms);
//------------------------

//Seed&Struct utilities
int spache_seed(node *atoms, int i);
void ratmpos(int i,int j,node *atoms,int *seed);
double rbfac(int *seed);
void gen_data(node *atoms,int i,int *seed);

int spache(node *atoms, int i,int N);
double connections(int **indexs,double **lens,int N,node *atoms);
int cycle(int **indexs,int N,node *atoms);

int ep_alloc(ep *p,int n);
int ep_free(ep *p,int n);
int ep_copy(ep *p1,ep *p2);

int ep_trans_o(ep *p,double u[3]);
int ep_trans_com(ep *p);
//-------------------------

//MC utilities
int aprob(double G0,double G1,double T,int i,int *seed);
double Temp(char filename[99]);

//moves
int HandOfGod(node *atoms,int N,int het,double xyzo[3],int **indexs,int *see);
//----------------------

//avo
void make_seed_b(FILE *runs, FILE *bonds,char dir[],int N,double sdist,double G[5],node *atoms,int **indexs);
void samples(char acc[],int n,int N,double sdist,int het);
int monte(int N,double sdist,int het,int co,int *iters,char sdir[]);
//---------------------

//Diagnostic utilities
double lig_dist(node *atoms,int N);
double m_centre(node *atoms,int N);
int lig_id(ep *p);
int count_connections(ep *p);
int lig_com(ep *p);
int count_diff(int N, node *atoms, node *old);

void dsytrd3(double A[3][3], double Q[3][3], double d[3], double e[2]);
int dsyevq3(double A[3][3], double Q[3][3], double w[3]);
void inertia_t(ep *p);//calculate axes

//-----------------------------------


#endif
