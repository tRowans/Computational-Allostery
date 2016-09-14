#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <time.h>
#include <complex.h>
#include <float.h>
#include "lib.h"

//Measures & Diagnostics
//dependencies: random utils, read utils, constraint utils, system utils

double lig_dist(node *atoms,int N)
//returns distance between binding sites
{
    int i,j=-1;
    double Rl[3]={0,0,0};
    for(i=0;i<N;i++)
    {
        if(ainb("HETATM",atoms[i].atm)==0)
        {
            Rl[0]+=j*atoms[i].x;
            Rl[1]+=j*atoms[i].x;
            Rl[2]+=j*atoms[i].z;
            j+=2;
        }
        if(j>1){break;}
    }

    return rad(Rl[0],Rl[1],Rl[2]);
}

double m_centre(node *atoms,int N)
//returns normalized distance between structure centre of mass and binding site centre of mass
{
    int i,j;
    double siz,maxsiz=0;
    double cdist;
    double R[3]={0,0,0},Rl[3]={0,0,0};
    for(i=0;i<N;i++)
    {
        R[0]+=atoms[i].x/(double)N;
        R[1]+=atoms[i].y/(double)N;
        R[2]+=atoms[i].z/(double)N;
        if(ainb("HETATM",atoms[i].atm)==0)
        {
            Rl[0]+=atoms[i].x*0.5;
            Rl[1]+=atoms[i].x*0.5;
            Rl[2]+=atoms[i].x*0.5;
        }
        for(j=i;j<N;j++)
        {
            if(i==j){continue;}
            siz = rad(atoms[i].x-atoms[j].x,atoms[i].y-atoms[j].y,atoms[i].z-atoms[j].z);
            if(siz>maxsiz){maxsiz=siz;}
        }
    }
    //Rl[0]*=0.5;Rl[1]*=0.5;Rl[2]*=0.5;
    //R[0]=R[0]/(double)N;R[1]=R[1]/(double)N;R[2]=R[2]/(double)N;
    cdist=rad(R[0]-Rl[0],R[1]-Rl[1],R[2]-Rl[2]);
    printf("\nmaxsize= %lf\tcdist= %lf\n",maxsiz,cdist);
    cdist=cdist/maxsiz;
    return cdist;
}

int lig_com(ep *p)
//determines lig axes
{
    int lig1=p->ligs[0], lig2=p->ligs[1];
    p->lig_com[0]=(p->atoms[lig1].x + p->atoms[lig2].x)/2;
    p->lig_com[1]=(p->atoms[lig1].y + p->atoms[lig2].y)/2;
    p->lig_com[2]=(p->atoms[lig1].z + p->atoms[lig2].z)/2; //printf("\ninside com\n%lf,%lf,%lf\n",p->lig_com[0],p->lig_com[1],p->lig_com[2]);

    p->lig_ax[0]=(p->atoms[lig1].x - p->atoms[lig2].x);
    p->lig_ax[1]=(p->atoms[lig1].y - p->atoms[lig2].y);
    p->lig_ax[2]=(p->atoms[lig1].z - p->atoms[lig2].z);
	//printf("\ninside1\n%lf,%lf,%lf\n",(p->atoms[lig1].x),(p->atoms[lig1].y),(p->atoms[lig1].z));
	//printf("inside2\n%lf,%lf,%lf\n",(p->atoms[lig2].x),(p->atoms[lig2].y),(p->atoms[lig2].z));
    /*double r=rad(p->lig_ax[0],p->lig_ax[1],p->lig_ax[2]);
    p->lig_ax[0]*=1/r;
    p->lig_ax[1]*=1/r;
    p->lig_ax[2]*=1/r;*/
	
    return 0;
}

int lig_id(ep *p)
{
    int lig1=-1, lig2=-1;
    for(int i=0;i<p->N;i++)
    {
        if(ainb("HETATM",p->atoms[i].atm)==0)
        {
            if(lig1==-1){lig1=i;}
            else if(lig2==-1){lig2=i;}
        }
        if(lig1!=-1 && lig2!=-1){break;}
    }
    p->ligs[0]=lig1; p->ligs[1]=lig2;
	//printf("\nlig1=%d\tlig2=%d\n",lig1,lig2);
    return 0;
}

int count_connections(ep *p)
{
    int **indi;
    indi=malloc(p->N *sizeof(*indi));
	for(int i=0;i<p->N;i++){indi[i]=calloc(p->N,sizeof(**indi));}
    connections(indi,p->N,p->atoms);
    int lig1=p->ligs[0], lig2=p->ligs[1];
    int l1=0,l2=0;
    for(int i=1;i<indi[lig1][0];i++)
    {
		//if(l0==-2){store[k]=indi[lig1][i]; k++;}
		for(int j=1;j<indi[lig2][0];j++)
        {
            if(indi[lig1][i]==indi[lig2][j]){l2++;break;}
			//else if(l0==-1){store[k]=indi[lig2][j]; k++;}
        }
    }
    l1=indi[lig2][0]+indi[lig1][0]-2*l2;
    p->con[0]=l1;p->con[1]=l2;
    //p->ligs[0]=lig1; p->ligs[1]=lig2;
    //ep_trans_com(p);


    for(int i=0;i<p->N;i++){free(indi[i]);}
    free(indi);
    return 0;
}


int count_diff(int N, node *atoms, char *path)
//Compares current state of *atoms to previously saved version and determines the number of atoms in different positions
//Returns this number
{
	node *old;
	int i=0, j=0;
	FILE *fp;

	printf("\nReading from %s\n", path);

	old = calloc(N, sizeof(node));

	fp = fopen(path, "r");
	read_pdb(fp, old);
	for (i = 0; i < N; i++)
	{
		if (atoms[i].x != old[i].x || atoms[i].y != old[i].y || atoms[i].z != old[i].z)
		{
			printf("\nNode %d is different\n", i);
			j++;
		}
	}

	fclose(fp);
	return(j);
}

//Diagonalization routine courtesy of:
/*Joachim Kopp
Numerical diagonalization of hermitian 3x3 matrices
arXiv.org preprint: physics/0610206
Int. J. Mod. Phys. C19 (2008) 523-548
*/

// ----------------------------------------------------------------------------
void dsytrd3(double A[3][3], double Q[3][3], double d[3], double e[2])
// ----------------------------------------------------------------------------
// Reduces a symmetric 3x3 matrix to tridiagonal form by applying
// (unitary) Householder transformations:
//            [ d[0]  e[0]       ]
//    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
//            [       e[1]  d[2] ]
// The function accesses only the diagonal and upper triangular parts of
// A. The access is read-only.
// ---------------------------------------------------------------------------
{
  const int n = 3;
  double u[n], q[n];
  double omega, f;
  double K, h, g;

  // Initialize Q to the identitity matrix
#ifndef EVALS_ONLY
  for (int i=0; i < n; i++)
  {
    Q[i][i] = 1.0;
    for (int j=0; j < i; j++)
      Q[i][j] = Q[j][i] = 0.0;
  }
#endif

  // Bring first row and column to the desired form
  h = SQR(A[0][1]) + SQR(A[0][2]);
  if (A[0][1] > 0)
    g = -sqrt(h);
  else
    g = sqrt(h);
  e[0] = g;
  f    = g * A[0][1];
  u[1] = A[0][1] - g;
  u[2] = A[0][2];

  omega = h - f;
  if (omega > 0.0)
  {
    omega = 1.0 / omega;
    K     = 0.0;
    for (int i=1; i < n; i++)
    {
      f    = A[1][i] * u[1] + A[i][2] * u[2];
      q[i] = omega * f;                  // p
      K   += u[i] * f;                   // u* A u
    }
    K *= 0.5 * SQR(omega);

    for (int i=1; i < n; i++)
      q[i] = q[i] - K * u[i];

    d[0] = A[0][0];
    d[1] = A[1][1] - 2.0*q[1]*u[1];
    d[2] = A[2][2] - 2.0*q[2]*u[2];

    // Store inverse Householder transformation in Q
#ifndef EVALS_ONLY
    for (int j=1; j < n; j++)
    {
      f = omega * u[j];
      for (int i=1; i < n; i++)
        Q[i][j] = Q[i][j] - f*u[i];
    }
#endif

    // Calculate updated A[1][2] and store it in e[1]
    e[1] = A[1][2] - q[1]*u[2] - u[1]*q[2];
  }
  else
  {
    for (int i=0; i < n; i++)
      d[i] = A[i][i];
    e[1] = A[1][2];
  }
}

// ----------------------------------------------------------------------------
int dsyevq3(double A[3][3], double Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using the QL algorithm with implicit shifts, preceded by a
// Householder reduction to tridiagonal form.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error (no convergence)
// ----------------------------------------------------------------------------
// Dependencies:
//   dsytrd3()
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double e[3];                   // The third element is used only as temporary workspace
  double g, r, p, f, b, s, c, t; // Intermediate storage
  int nIter;
  int m;

  // Transform A to real tridiagonal form by the Householder method
  dsytrd3(A, Q, w, e);

  // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
  // with the QL method
  //
  // Loop over all off-diagonal elements
  for (int l=0; l < n-1; l++)
  {
    nIter = 0;
    while (1)
    {
      // Check for convergence and exit iteration loop if off-diagonal
      // element e(l) is zero
      for (m=l; m <= n-2; m++)
      {
        g = fabs(w[m])+fabs(w[m+1]);
        if (fabs(e[m]) + g == g)
          break;
      }
      if (m == l)
        break;

      if (nIter++ >= 30)
        return -1;

      // Calculate g = d_m - k
      g = (w[l+1] - w[l]) / (e[l] + e[l]);
      r = sqrt(SQR(g) + 1.0);
      if (g > 0)
        g = w[m] - w[l] + e[l]/(g + r);
      else
        g = w[m] - w[l] + e[l]/(g - r);

      s = c = 1.0;
      p = 0.0;
      for (int i=m-1; i >= l; i--)
      {
        f = s * e[i];
        b = c * e[i];
        if (fabs(f) > fabs(g))
        {
          c      = g / f;
          r      = sqrt(SQR(c) + 1.0);
          e[i+1] = f * r;
          c     *= (s = 1.0/r);
        }
        else
        {
          s      = f / g;
          r      = sqrt(SQR(s) + 1.0);
          e[i+1] = g * r;
          s     *= (c = 1.0/r);
        }

        g = w[i+1] - p;
        r = (w[i] - g)*s + 2.0*c*b;
        p = s * r;
        w[i+1] = g + p;
        g = c*r - b;

        // Form eigenvectors
#ifndef EVALS_ONLY
        for (int k=0; k < n; k++)
        {
          t = Q[k][i+1];
          Q[k][i+1] = s*Q[k][i] + c*t;
          Q[k][i]   = c*Q[k][i] - s*t;
        }
#endif
      }
      w[l] -= p;
      e[l]  = g;
      e[m]  = 0.0;
    }
  }

  return 0;
}


void inertia_t(ep *p)
//calculate principle axes
{

    double Is[3][3]={0},eigvecs[3][3],eigvals[3];
    //calculating inertia tensor
    for(int i=0;i<p->N;i++)
    {
        Is[0][0] += SQR(p->atoms[i].y) + SQR(p->atoms[i].z);//Ixx
        Is[0][1] -= (p->atoms[i].x)*(p->atoms[i].y);//Ixy
        Is[0][2] -= (p->atoms[i].x)*(p->atoms[i].z);//Ixz
        Is[1][1] += SQR(p->atoms[i].x) + SQR(p->atoms[i].z);//Iyy
        Is[1][2] -= (p->atoms[i].y)*(p->atoms[i].z);//Iyz
        Is[2][2] += SQR(p->atoms[i].x) + SQR(p->atoms[i].y);//Izz
    }
    Is[1][0] = Is[0][1];       //Iyx
    Is[2][0] = Is[0][2];       //Izx
    Is[2][1] = Is[1][2];       //Izy

    dsyevq3(Is,eigvecs,eigvals);

    memcpy(p->Is,eigvecs,sizeof(eigvecs));
    //return 0;
}
