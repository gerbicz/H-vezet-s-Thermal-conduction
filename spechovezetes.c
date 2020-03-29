#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define eps 1e-10

// sávos Gauss elimináció és gazdaságos memóriahasználat a 2 spec. feladatra:
// tesztfeladat+termes
// O(n^5) memória, kb. 16*n^5 byte
// O(n^7) idõ

int main()

{

    int seconds;

    FILE* out;
    out=fopen("hovezetes_output.txt","a+");
    
    int a0,b0,c0,n,q;
    int eset,add;
    
    int d,h,i,j,k;
    int nb;  // nem inaktiv cellak szama
   
    printf("Hovezetes 3 dimenzioban\n");
    
    double **P;
    int **Q;

    int L1,L2,L3;
    double H1,H2,H3;

    while(1) {
           printf("Tesztfeladat vagy tanar es diak feladat. Valasszon esetet (1/2): ");
           scanf("%d",&eset);
           if((eset==1)||(eset==2))  break;
           else printf("Csak 1-t vagy 2-t lehet valasztani!\n");
    }

    switch(eset)  {
    case 1: 
           printf("A peremen adott homerseklet. H1=H2=H3=1 meter. Hovezetesi egyutthato==1\n");
           printf("nxnxn-es feladat. Adja meg n-et: ");
           scanf("%d",&n);
           printf("T[x,y,z]=a*x^3+b*y*z+c*z^2 az elert homerseklet ( a,b,c egeszek ), peremen is\n");
           printf("a=");
           scanf("%d",&a0);
           printf("b=");
           scanf("%d",&b0);
           printf("c=");
           scanf("%d",&c0);
           seconds=time(NULL);
           L1=n,L2=n,L3=n;
           nb=n*n*n;
           H1=1.0,H2=1.0,H3=1.0;
           
           Q=(int**) (malloc) ((nb+1)*sizeof(int));
           for(i=1;i<=nb;i++)  Q[i]=(int*) (malloc) (5*sizeof(int));

           P=(double**)  (malloc) ((nb+1)*sizeof(double));
           for(i=1;i<=nb;i++)  P[i]=(double*) (malloc) (8*sizeof(double));         
           
           for(i=1;i<=n;i++)  {
               for(j=1;j<=n;j++)  {
                   for(k=1;k<=n;k++)  {
                   if((i%n<2)||(j%n<2)||(k%n<2))  {
                       q=i+(j-1)*n+(k-1)*n*n;
                       Q[q][1]=i;
                       Q[q][2]=j;
                       Q[q][3]=k;
                       Q[q][4]=3;
                       P[q][2]=1.0;
                       P[q][3]=0.0;
                       P[q][4]=(double) a0*i*i*i+b0*j*k+c0*k*k;
                       P[q][5]=0.0;
                       P[q][6]=0.0;
                       P[q][7]=0.0;
                   }
                   else  {
                       q=i+(j-1)*n+(k-1)*n*n;
                       Q[q][1]=i;
                       Q[q][2]=j;
                       Q[q][3]=k;
                       Q[q][4]=1;
                       P[q][2]=1.0;
                       P[q][3]=0.0;
                       P[q][4]=0.0;
                       P[q][5]=0.0;
                       P[q][6]=(double) -6*a0*i-2*c0;
                       P[q][7]=0.0;
                   }
                   }
                }
             }
             printf("Felsavszelesseg=%d\n",n*n);
             n=n*n*n;
             printf("Ismeretlenek szama=%d\n",n);
             if(c0>0) printf("Hoforras a hotermelo cellakban: Q[x,y,z]=%d*x%d\n",-6*a0,-2*c0);
             else      printf("Hoforras a hotermelo cellakban: Q[x,y,z]=%d*x+%d\n",-6*a0,-2*c0);
             break;
    case 2:
           printf("Tanar es diak feladat\n");
           printf("nxnxn-es kockara. Adja meg n-et: \n");
           scanf("%d",&n);
           printf("z=1 sikon (talaj) turbulens hocsere a kulvilaggal  T_0=290 K, a=5 W/(m^2*K)\n");
           printf("x=1 vagy x=%d vagy y=1 vagy z=%d-re szigetelt perem.\n",n,n);
           printf("y=%d-re (ablak miatt) eloirt hoaram: j=1 W/m^2\n",n);
           printf("x+z=7 (sikon) pontokban vannak a diakok (hotermelo cella) : Q[x,y,z]=5 W/m^3\n");
           printf("tabla elott a tanar: Q[%d,%d,2]=15 W/m^3\n",n>>1,n-1);
           printf("futotest: Q[2,%d,2]=45 W/m^3\n",n>>1);
           seconds=time(NULL);
           L1=n,L2=n,L3=n;
           nb=n*n*n;
           H1=1.0,H2=1.0,H3=1.0;

           Q=(int**) (malloc) ((nb+1)*sizeof(int));
           for(i=1;i<=nb;i++)  Q[i]=(int*) (malloc) (5*sizeof(int));

           P=(double**)  (malloc) ((nb+1)*sizeof(double));
           for(i=1;i<=nb;i++)  P[i]=(double*) (malloc) (8*sizeof(double));         
           
           for(i=1;i<=n;i++)  {
               for(j=1;j<=n;j++)  {
                   for(k=1;k<=n;k++)  {
                       q=i+(j-1)*n+(k-1)*n*n;
                       Q[q][1]=i;
                       Q[q][2]=j;
                       Q[q][3]=k;
                   
                   if(k==1)  {
                       Q[q][4]=5;
                       P[q][2]=1.0;
                       P[q][3]=0.0;
                       P[q][4]=0.0;
                       P[q][5]=5.0;
                       P[q][6]=0.0;
                       P[q][7]=290.0;
                   }
                   else if((i==1)||(i==n)||(j==1)||(k==n)) {
                       Q[q][4]=2;
                       P[q][2]=1.0;
                       P[q][3]=0.0;
                       P[q][4]=0.0;
                       P[q][5]=0.0;
                       P[q][6]=0.0;
                       P[q][7]=0.0;
                   }
                   else if(j==n) {
                       Q[q][4]=4;
                       P[q][2]=1.0;
                       P[q][3]=1.0;
                       P[q][4]=0.0;
                       P[q][5]=0.0;
                       P[q][6]=0.0;
                       P[q][7]=0.0;
                   }
                   else if(j+k==7) {
                       Q[q][4]=1;
                       P[q][2]=1.0;
                       P[q][3]=0.0;
                       P[q][4]=0.0;
                       P[q][5]=0.0;
                       P[q][6]=5.0;
                       P[q][7]=0.0;
                   }
                   else if((i==(n>>1))&&(j==n-1)&&(k==2)) {
                       Q[q][4]=1;
                       P[q][2]=1.0;
                       P[q][3]=0.0;
                       P[q][4]=0.0;
                       P[q][5]=0.0;
                       P[q][6]=15.0;
                       P[q][7]=0.0;
                   }
                   else if((i==2)&&(j==(n>>1))&&(k==2)) {
                       Q[q][4]=1;
                       P[q][2]=1.0;
                       P[q][3]=0.0;
                       P[q][4]=0.0;
                       P[q][5]=0.0;
                       P[q][6]=45.0;
                       P[q][7]=0.0;
                   }
                   else {
                       Q[q][4]=1;
                       P[q][2]=1.0;
                       P[q][3]=0.0;
                       P[q][4]=0.0;
                       P[q][5]=0.0;
                       P[q][6]=0.0;
                       P[q][7]=0.0;
                   }
                   }
                }
             }
             printf("Felsavszelesseg=%d\n",n*n);
             n=n*n*n;
             printf("Ismeretlenek szama=%d\n",n);
             break;
    }

    int L[4]={0,L1,L2,L3};
    double H[4]={0.0,H1,H2,H3};
    
    double K;
    
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

    int felsavszelesseg=L1*L2;

    double **A;
    A=(double**) (malloc) ((nb+1)*sizeof(double));
    for(i=1;i<=nb;i++)  A[i]=(double*) (malloc) ((2*felsavszelesseg+2)*sizeof(double));

    int minimum,maximum;

    for(i=1;i<=nb;i++)
        for(j=0;j<=2*felsavszelesseg+1;j++)  A[i][j]=0.0;

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
         if(cellatipus==0)   A[P0][felsavszelesseg+1]=1.0,b[P0]=0.0;
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
           else if(cellatipus==3)  A[P0][felsavszelesseg+1]=1.0,b[P0]=s[4];
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
                      A[P0][felsavszelesseg+1]-=fact[d]*K/H[d];
                      A[P0][P1-P0+felsavszelesseg+1]=fact[d]*K/H[d];  // += volt itt, de ez egyenloseg!
                   }
                   else if(cellatipus==4)  b[P0]-=fact[d]*s[3];
                   else if(cellatipus==5)  {
                           b[P0]-=fact[d]*s[5]*s[7];
                           A[P0][felsavszelesseg+1]-=fact[d]*s[5];
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
       u=A[i][felsavszelesseg+1];

       for(j=i;(j<=nb)&&(j<=i+felsavszelesseg);j++)  A[i][j-i+felsavszelesseg+1]/=u;
       b[i]/=u;
       
       for(j=i+1;(j<=nb)&&(j<=i+felsavszelesseg);j++)  {
           if((A[j][i-j+felsavszelesseg+1]>eps)||(A[j][i-j+felsavszelesseg+1]<-eps))  {
               szorzo=A[j][i-j+felsavszelesseg+1];
               hatar=i+felsavszelesseg;
               if(hatar>nb)  hatar=nb;
               minimum=i-j+felsavszelesseg+1;
               maximum=hatar-j+felsavszelesseg+1;
               add=j-i;
               for(k=minimum;k<=maximum;k++)  A[j][k]-=szorzo*A[i][k+add];
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
        for(k=j;k>i;k--)  u+=A[i][k-i+felsavszelesseg+1]*x[k];
        x[i]=(b[i]-u)/A[i][felsavszelesseg+1];
    }

   // eredmeny kiiratasa

   for(k=1;k<=L3;k++)  {
       for(j=1;j<=L2;j++)  {
           for(i=1;i<=L1;i++)  {
              if(tipus_U[n]!=0)  fprintf(out,"T[%d,%d,%d]=%lf\n",i,j,k,x[i+(j-1)*L1+(k-1)*L1*L2]);
           }
       }
   }

   printf("Savos Gauss eliminacio kesz.\n");
   printf("Hasznalt ido=%ld masodperc\n",time(NULL)-seconds);
   printf("Eredmeny a hovezetes_output.txt fileban\n");

   if(eset==1)  {
      n=L1;  // L1=L2=L3=n volt
      double maxerror=0.0;
      for(i=1;i<=n;i++)  {
          for(j=1;j<=n;j++)  {
              for(k=1;k<=n;k++)  {
                  q=i+(j-1)*n+(k-1)*n*n;
                  d=a0*i*i*i+b0*j*k+c0*k*k;
                  if(x[q]>(double) d+maxerror)  maxerror=(double) x[q]-d;
                  if(x[q]<(double) d-maxerror)  maxerror=(double) d-x[q];
              }
          }
      }
      printf("Hibavektor maximum normaja=%lf\n",maxerror);
   }

   double mini=x[1];
   double maxi=x[1];
   int ax=1,ay=1,az=1;
   int bx=1,by=1,bz=1;
   
   n=L1;
   for(i=1;i<=n;i++)  {
       for(j=1;j<=n;j++)  {
           for(k=1;k<=n;k++)  {
               q=i+(j-1)*n+(k-1)*n*n;
               if(x[q]<mini)  mini=x[q],ax=i,ay=j,az=k;
               if(x[q]>maxi)  maxi=x[q],bx=i,by=j,bz=k;
           }
       }
   }
   
   printf("Kialakult legalacsonyabb homerseklet=%lf K [%d,%d,%d] helyen\n",mini,ax,ay,az);
   printf("Kialakult legmagasabb homerseklet=%lf K [%d,%d,%d] helyen\n",maxi,bx,by,bz);
   
   fclose(out);

   return 1;
}
