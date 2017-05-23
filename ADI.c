#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.1415926535897932384

int main(){

	/* Open txt file to write data */
   	FILE *finalHeatDistribution;
   	finalHeatDistribution = fopen("finalHeatDistribution.txt","w");

    if(finalHeatDistribution == NULL)
  	{
      printf("Error! finalHeatDistribution");   
      exit(1);             
   	}

  /* Declare variables */ 
    int i,j,k,n;
    int N=10;
    /*int Ny=20;*/
    int Nt;
    float alpha, dx, dy, dt, rx, ry, Lx, Ly, tmax;
    float d[N+1],T[N+1][N+1], ax[N+1], bx[N+1], cx[N+1], ay[N+1], by[N+1], cy[N+1];
    float Ax[N+1][N+1],Ay[N+1][N+1];

    /* Space domain */
    Lx=10.0;
    dx=Lx/N;

    Ly=10.0;
    dy=Ly/N;

    /* Time step */ 
    Nt=10;
    tmax=10.0;
    dt=tmax/Nt;

    /* Parameters */ 
    alpha=0.05;
    rx=alpha*dt/(2*dx*dx);
    ry=alpha*dt/(2*dy*dy);

  /* Populating matrices and vectors with zeros or initial condition*/
    for (i=0;i<=N;i++){
      for(j=0;j<=N;j++){
        T[i][j]=3.0;
      }
      d[i]=0;
    }

  /*Populating Thomas algorithm coefficient vectors ax,bx and ay,by which remain constant */
    for (i=0;i<=N;i++){
        if (i==0){
          ax[i]=0;
          ay[i]=0;

          bx[i]=1+2*rx;
          by[i]=1+2*ry;
        } 
        else if (i==N){
          ax[i]=0;
          ay[i]=0;

          bx[i]=1+2*rx;
          by[i]=1+2*ry;
        }
        else {
        ax[i]=-rx;
        bx[i]=1+2*rx;

        ay[i]=-ry;
        by[i]=1+2*ry;
        }
      }

  /* Homogeneous boundary conditions */
    for (k=0;k<=N;k++){
        T[0][k]=0;
        T[k][0]=0;
        T[N][k]=0;
        T[k][N]=0;
      }

  /* Step through time */
    for (n=1;n<=Nt;n++){
      
      /* X SWEEP */
      /* Forming vector dij */
      for (j=1;j<=N-1;j++){
        for (i=0;i<=N;i++){
          d[i]=ry*T[i][j-1]+(1-2*ry)*T[i][j]+ry*T[i][j+1];
        }

        /* Thomas algorithm to solve A*T[:][j]=d[:] */
        /*Populating Thomas algorithm coefficient vector c */
        for (k=0;k<=N;k++){
          if (k==N) {
            cx[k]=0;
          }
          else if (k==0){
            cx[k]=0;
          }
          else {
            cx[k]=-rx;
          }
        }

        /* Forward substitution */
        for(k=1;k<=N;k++){
          cx[k]=cx[k]/(bx[k]-ax[k]*cx[k-1]);
          d[k]=(d[k]-ax[k]*d[k-1])/(bx[k]-ax[k]*cx[k-1]);
        }

        T[N][j]=d[N];

        /* Backwards substitution */
        for(k=N-1;k>=0;k-=1){
          T[k][j]=d[k]-cx[k]*T[k+1][j];
        }
      }/* end X SWEEP */

      /* Y SWEEP */
      /* Forming vector dij */

      for (i=1;i<=N-1;i++){
        for (j=0;j<=N;j++){
          d[j]=rx*T[i-1][j]+(1-2*rx)*T[i][j]+rx*T[i+1][j];
        }

        /* Thomas algorithm to solve A*T[:][j]=d[:] */
        /*Populating Thomas algorithm coefficient vector c */
        for (k=0;k<=N;k++){
          if (k==N) {
            cy[k]=0;
          }
          else if (k==0){
            cy[k]=0;
          }
          else {
            cy[k]=-ry;
          }
        }
        /* Forward substitution */
        for(k=1;k<N;k++){
          cy[k]=cy[k]/(by[k]-ay[k]*cy[k-1]);
          d[k]=(d[k]-ay[k]*d[k-1])/(by[k]-ay[k]*cy[k-1]);
        }

        T[i][N]=d[N];

        /* Backwards substitution */
        for(k=N-1;k>=0;k-=1){
          T[i][k]=d[k]-cy[k]*T[i][k+1];
        }
      }/* end Y SWEEP */
    }/* end time  */

  /* Writing data to file finalHeatDistribution.txt */
    for (j=0;j<=N;j++){
      for (i=0;i<=N;i++){
          fprintf(finalHeatDistribution,"%f %f %f\n",dx*i,dy*j,T[i][j]);     
      }
    }

    fclose(finalHeatDistribution);
    return 0;
}