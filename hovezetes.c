#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define eps 1e-10

// az egész A mátrix van tárolva a programban ezt még lehet csökkenteni
// , ha csak a sávban lévõket tároljuk, majd 
// A[x,y] helyett A[x,y-x+felsavszelesseg+1]-et írunk, ezzel a memóriaigény O(L1^2*L2^2*L3)-ra csökkenne

int main()

{

    int seconds=time(NULL);

    FILE* in;
    FILE* out;
    
    in=fopen("hovezetes_input.txt","r");
    out=fopen("hovezetes_output.txt","a+");
    
    int n;
    fscanf(in,"%d",&n);  // aktiv cellak szama

    int L1,L2,L3;
    double H1,H2,H3;
    fscanf(in,"%d %d %d",&L1,&L2,&L3);
    fscanf(in,"%lf %lf %lf",&H1,&H2,&H3);
    
    int d,h,i,j,k;
    int nb=L1*L2*L3;
    
    int L[4]={0,L1,L2,L3};
    double H[4]={0.0,H1,H2,H3};
    
    double K;
    
    double **P;
    P=(double**)  (malloc) ((n+1)*sizeof(double));
    for(i=1;i<=n;i++)  P[i]=(double*) (malloc) (8*sizeof(double));

    int **Q;
    Q=(int**) (malloc) ((n+1)*sizeof(int));
    for(i=1;i<=n;i++)  Q[i]=(int*) (malloc) (5*sizeof(int));
   
    for(i=1;i<=n;i++)  {
        for(j=1;j<=4;j++)  fscanf(in,"%d",&Q[i][j]);
        for(j=2;j<=7;j++)  fscanf(in,"%lf",&P[i][j]);
    }
    fclose(in);

    double **U;
    U=(double**)  (malloc) ((nb+1)*sizeof(double));
    for(i=1;i<=nb;i++)  U[i]=(double*) (malloc) (8*sizeof(double));

    int *tipus_U;
    tipus_U=(int*) (malloc) ((nb+1)*sizeof(int));
    
    for(i=1;i<=nb;i++)  tipus_U[i]=0;  // ha nincs definialva, akkor inaktiv a cella

    int Tip[8];

    int w[4][7]={{0,0,0,0,0,0,0},
                 {0,1,-1,0,0,0,0},
                 {0,0,0,L1,-L1,0,0},
                 {0,0,0,0,0,L1*L2,-L1*L2}};

    int R[7]={0,1,-1,L1,-L1,L1*L2,-L1*L2};

    double **A;
    A=(double**) (malloc) ((nb+1)*sizeof(double));
    for(i=1;i<=nb;i++)  A[i]=(double*) (malloc) ((nb+1)*sizeof(double));

    int felsavszelesseg=L1*L2;
    int minimum,maximum;

    for(i=1;i<=nb;i++)  {
        minimum=i-felsavszelesseg;
        if(minimum<1)  minimum=1;
        maximum=i+felsavszelesseg;
        if(maximum>nb) maximum=nb;
        for(j=minimum;j<=maximum;j++)  A[i][j]=0.0;
    }


    double *b;
    b=(double*) (malloc) ((nb+1)*sizeof(double));
    
    for(i=1;i<=nb;i++)  b[i]=0.0;

    int F[7]={0,0,1,0,1,0,1};

    int k1,k2,k3;
    int P0,P1; 
    int G[4];    
    double s[8],s2[8];
    int cellatipus,Tipus;
    double c[4];
    double fact[4];
    double V0;
    
    for(i=1;i<=n;i++)  {
         for(j=2;j<=7;j++)  s[j]=P[i][j];
         k1=Q[i][1];
         k2=Q[i][2];
         k3=Q[i][3];
         d=Q[i][4];
         P0=k1+L1*(k2-1)+L1*L2*(k3-1);
         tipus_U[P0]=d;
         for(j=2;j<=7;j++)  U[P0][j]=s[j];
     }

     for(i=1;i<=n;i++)  {
         for(j=2;j<=7;j++)  s[j]=P[i][j];
         k1=Q[i][1];
         k2=Q[i][2];
         k3=Q[i][3];
         G[1]=k1;
         G[2]=k2;
         G[3]=k3;
         P0=k1+L1*(k2-1)+L1*L2*(k3-1);
         cellatipus=tipus_U[P0];
         for(j=2;j<=7;j++)  s[j]=U[P0][j];
         if(cellatipus==0)   A[P0][P0]=1.0,b[P0]=0.0;
         if((cellatipus==1)&&((k1%L1<2)||(k2%L2<2)||(k3%L3<2)))  {
             printf("x=%d,y=%d,z=%d pozicioban hotermelo cella a modell szelen van. Ez nem lehet!\n",k1,k2,k3);
             exit(1);
         }
         if(cellatipus==1)  {
            for(h=1;h<=6;h++)  {
                d=(h+1)>>1;
                P1=P0+R[h];
                Tipus=tipus_U[P1];
                if(Tipus==0)  {
                   printf("Nem lezart a modell x=%d,y=%d,z=%d pozicioban!\n",k1+w[1][h],k2+w[2][h],k3+w[3][h]);
                   exit(1);
                }
             }
          }
          c[1]=H1;
          c[2]=H2;
          c[3]=H3;
          
           for(h=1;h<=6;h++)  {
               d=(h+1)>>1;
               P1=P0+R[h];
               if((G[d]%L[d]==F[h])||(tipus_U[P1]==0))  c[d]=0.5*H[d];
           }
           fact[1]=c[2]*c[3];
           fact[2]=c[1]*c[3];
           fact[3]=c[1]*c[2];
           if(cellatipus==0)  ;
           else if(cellatipus==3)  A[P0][P0]=1.0,b[P0]=s[4];
           else {
               for(h=1;h<=6;h++)  {
                   d=(h+1)>>1;
                   P1=P0+R[h];
                   if(G[d]%L[d]==F[h])  Tip[h]=0;
                   else                 Tip[h]=tipus_U[P1];
               }
               V0=c[1]*c[2]*c[3];
               b[P0]-=s[6]*V0;
               for(h=1;h<=6;h++)  {
                   d=(h+1)>>1;
                   P1=P0+R[h];
                   Tipus=Tip[h];
                   if(Tipus!=0)  {
                      for(j=2;j<=7;j++)  s2[j]=U[P1][j];
                      K=2.0/(1.0/s[2]+1.0/s2[2]);
                      A[P0][P0]-=fact[d]*K/H[d];
                      A[P0][P1]=fact[d]*K/H[d];  // += volt itt, de ez egyenloseg!
                   }
                   else if(cellatipus==4)  b[P0]-=fact[d]*s[3];
                   else if(cellatipus==5)  {
                           b[P0]-=fact[d]*s[5]*s[7];
                           A[P0][P0]-=fact[d]*s[5];
                   }
               }
          }
     }

   // A*x=b megoldása savos Gauss eliminációval, feltetelezve, hogy nem kell sort cserelni
   // felsavszelesseg=L1*L2 teljesül az A-ra.
   // nb=L1*L2*L3

   printf("Matrixot felepitette a program.\n");

   double u;
   double szorzo;
   int hatar;

   for(i=1;i<=nb;i++)  {
       u=A[i][i];

       for(j=i;(j<=nb)&&(j<=i+felsavszelesseg);j++)  A[i][j]/=u;
       b[i]/=u;
       
       for(j=i+1;(j<=nb)&&(j<=i+felsavszelesseg);j++)  {
           if((A[j][i]>eps)||(A[j][i]<-eps))  {
               szorzo=A[j][i];
               hatar=i+felsavszelesseg;
               if(hatar>nb)  hatar=nb;
               for(k=i;k<=hatar;k++)  A[j][k]-=szorzo*A[i][k];
               b[j]-=szorzo*b[i];
           }
       }
    }

    double *x;
    x=(double*) (malloc) ((nb+1)*sizeof(double));
    
    // felso haromszog alakra van hozva ezzel a matrix
    
    for(i=nb;i>=1;i--)  {
        u=0.0;
        j=i+felsavszelesseg;
        if(j>nb)  j=nb;
        for(k=j;k>i;k--)  u+=A[i][k]*x[k];
        x[i]=(b[i]-u)/A[i][i];
    }

   // eredmeny kiiratasa

   for(k=1;k<=L3;k++)  {
       for(j=1;j<=L2;j++)  {
           for(i=1;i<=L1;i++)  {
              if(tipus_U[n]!=0)  fprintf(out,"T[%d,%d,%d]=%lf\n",i,j,k,x[i+(j-1)*L1+(k-1)*L1*L2]);
           }
       }
   }

   printf("Hasznalt ido=%ld masodperc\n",time(NULL)-seconds);

   fclose(out);

   return 1;
}
