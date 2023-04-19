#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define pi atan(1)*4
#define Max(a,b) (((a)>(b))?(a):(b))
#define Min(a,b) (((a)<(b))?(a):(b))
#define Sgn(a) (((a)==0.0)?0.0:((a)/(fabs(a))))

// #define sta 1

// argument of (x[1],x[2]) in [0,2\pai) or -1 (the origin)

// double arg( double *x ){
//   double eps=1.0e-10, a; 
//   if ( x[1]*x[1]+x[2]*x[2] < eps ){
//     return -1.0;
//   }
//   if ( x[1] >= fabs(x[2]) ){
//     a = atan(x[2]/x[1]);
//     return (a>=0.0)?a:(a+2*pi);
//   }
//   if ( x[2] >= fabs(x[1]) ){
//     a = atan(-x[1]/x[2]);
//     return a + 0.5*pi;
//   }
//   if ( x[1] <= - fabs(x[2]) ){
//     a = atan(x[2]/x[1]);
//     return a+pi;
//   }
//   if ( x[2] <= - fabs(x[1]) ){
//     a = atan(-x[1]/x[2]);
//     return a + 1.5*pi;
//   }
//   exit(0);
// }


double *dvector(unsigned int i,unsigned  int j){
  double *a;
  int length=j-i+1;
  if ( (a=(double *)malloc( (length*sizeof(double)) ) ) == NULL ){
    printf("dvector: memory allocation is failed \n");
    exit(1);
  }
  //最初の位置が０であるようなベクトルをつくる。
 for(int n=0;n<j-i+1;n++){
    *(a+n)=0.0;
  } 
  
  return (a-i);
}

void free_dvector(double *a,unsigned int i,unsigned int j){
  if(i!=j){
    free((double*)(a+i));
  }else{
    printf("Non Length!\n");
    exit(1);
  }
  a=NULL;
}

int *ivector(unsigned int i,unsigned int j){
  int *a;
  if ( (a=(int *)malloc( ((j-i+1)*sizeof(int))) ) == NULL ){
    printf("ivector: memory allocation is failed \n");
    exit(1);
  }  
  for(int n=0;n<j-i+1;n++){
    *(a+n)=1;
  } 
  return (a-i);
}

void free_ivector(int *a,unsigned int i,unsigned int j){
  if(i!=j){
    free(a+i);
  }else{
    printf("Non Length!\n");
    exit(1);
  }
  a=NULL;
}

double **dmatrix(unsigned int nr1,unsigned int nr2,unsigned int nl1,unsigned int nl2){
  int i, nrow, ncol; 
  double **a; 
  nrow = nr2 - nr1 + 1 ; // number of row
  ncol = nl2 - nl1 + 1 ; // number of column
  if ( ( a=(double **)malloc( nrow*sizeof(double *) ) ) == NULL ){
    printf("dmatrix: memory allocation is failed \n");
    exit(1);
  }
  a = a - nr1; 
  for( i=nr1; i<=nr2; i++) a[i] = (double *)malloc(ncol*sizeof(double)); 
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1; 

  for(i=nr1;i<=nr2;i++){
    for(int j=nl2;j<=nl2;j++){
      a[i][j]=0.0;
    }    
  }        
  return a;
}

void free_dmatrix(double **a,unsigned int nr1,unsigned int nr2,unsigned int nl1,unsigned int nl2){
  for (int i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
  return;
}


int **imatrix(unsigned int nr1,unsigned int nr2,unsigned int nl1,unsigned int nl2){
  int i, nrow, ncol;
  int **a; 
  nrow = nr2 - nr1 + 1 ; /* number of row */
  ncol = nl2 - nl1 + 1 ; /* number of column */
  if ( ( a=(int **)malloc( nrow*sizeof(int *) ) ) == NULL ){
    printf("imatrix: memory allocation is failed \n");
    exit(1);
  }
  a = a - nr1; 
  for( i=nr1; i<=nr2; i++) a[i] = (int *)malloc(ncol*sizeof(int));
  for( i=nr1; i<=nr2; i++) a[i] = a[i]-nl1;

 for(i=nr1;i<=nr2;i++){
    for(int j=nl2;j<=nl2;j++){
      a[i][j]=0;
    }    
  }

  return (a);
}

void free_imatrix(int **a, unsigned int nr1,unsigned int nr2,unsigned int nl1,unsigned int nl2){
  int i=nl2;

  for ( i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
  free((void *)(a+nr1));
  return;
}

double *matrix_vector_product(double **a, double *b ,unsigned int n){
  int sta=1;
  int end=sta+n-1;
  double *c=dvector(1,n);

  double wk;
  int i, j;
  for ( i = sta; i <= end; i++){
    wk = 0.0;
    for ( j = sta; j <= end; j++ ){
      wk += a[i][j]*b[j];
    }
    c[i] = wk;
  }
  return c;
}

double inner_product( int m, int n, double *a, double *b){
  int i;
  double s = 0.0;
  for( i = m; i <= n; i++) s += a[i]*b[i];
  return s;
}
double vector_norm1( double *a, int m, int n ){
  int i; 
  double norm = 0.0;
  double num=0;
  for ( i = m; i <= n; i++ ){
    num+=1.0;
    norm += fabs(a[i]);
  }
  return norm/num; 
}
//ベクトルのノルム
double norm(double *X1,int m, int n){
  double sum=0.0;
  double num=0;
  for(int i=m;i<=n;i++){
    num+=1.0;
    sum+=pow(X1[i],2);
  }
  return sqrt(sum/num);
}


//２本のベクトルからなす角を計算する
double arg(double *x1,double *x2,int M){
  int sta=1;
  int end=sta+M;
  return acos(inner_product(0,M,x1,x2)/(norm(x1,sta,end)*norm(x2,sta,end)));
}

int factorial(int n){
  if(n<=0){
    return 1;
  }
    return n*factorial(n-1);
}

void pivod(double **A,double *b,int M){
  int sta=1;
  int end=sta+M;

  int j=sta;//sta列目に関してを考える
  while(j<end){
    /*=============j列目において最大絶対値となる行の探索====================*/
    int pivod=j;//交換する行の更新
    double pivod_base=fabs(A[j][j]);//対角成分を基準にとる
  for(int k=j+1;k<end;k++){//j+1~endまでの最大絶対値を探索
    if(fabs(A[k][j])>pivod_base){
      pivod_base=fabs(A[k][j]);
      pivod=k;
    }
  }
    // printf("%f\n",pivod_base);
    //k行目が最大絶対値になる
    /*================================================================*/

    if(pivod_base!=0.0){
      /*=================pivod列目とj列目の全要素の置換=====================*/
      for(int col=sta;col<end;col++){
        double swap=A[j][col];
        A[j][col]=A[pivod][col];
        A[pivod][col]=swap;
      }
      double swap_b=b[j];
      b[j]=b[pivod];
      b[pivod]=swap_b;
      /*================================================================*/
    }else{
      /*=========対角成分が0になる場合，同列の0でない部分を加える===============*/
      int non_zero=0;
      for(int k=sta;k<end;k++){
        if(A[k][j]!=0.0){//k行目をj行目に加える
          non_zero=k;
        }
      }
      for(int i=sta;i<end;i++){
        A[j][i]+=A[non_zero][i];
      }
      b[j]+=b[non_zero];
      /*================================================================*/
    }
    j++;
  }
}

void LU(double **Acoef,int M,double **L_r,double **U_r){//要素が行列のベクトルを返す関数
  printf("LU ");
  int sta=1;
  int end=sta+M;
  double **L=dmatrix(sta,end-1,sta,end-1);//alloc_matrix(M,M);
  double **U=dmatrix(sta,end-1,sta,end-1);//alloc_matrix(M,M);
  
  // // 対角成分に0があったらピボットの実行
  // double check;
  // for(int i=sta;i<end;i++){
  //   if(Acoef[i][i]==0.0){
  //     check=1;
  //   }
  // }
  
  // if(check==1){
    double *B=dvector(sta,end-1);
    pivod(Acoef,B,M);//ピボットがあるとうまくいかない
    free_dvector(B,sta,end-1);
  // }

  for(int i=sta;i<=end-1;i++){
    U[i][i]=1.0;
    for(int j=i;j<=end-1;j++){
      if(i==sta){
        L[j][i]=Acoef[j][i];
      }else{
        double sumU=0.0;
        for(int sumnumber=sta;sumnumber<i-1;sumnumber++){
          sumU+=L[i-1][sumnumber]*U[sumnumber][j];
        }
        U[i-1][j]=(Acoef[i-1][j]-sumU)/L[i-1][i-1];

        double sumL=0.0;
        for(int sumnumber=sta;sumnumber<i;sumnumber++){
          sumL+=L[j][sumnumber]*U[sumnumber][i];
        }
        L[j][i]=Acoef[j][i]-sumL;
      }
    }
  }
  
  printf("L\n");
  for(int i=sta;i<end;i++){
    for(int j=sta;j<end;j++){
      L_r[i][j]=L[i][j];
    }
  }

  printf("U\n");
  for(int i=sta;i<end;i++){
    for(int j=sta;j<end;j++){
      U_r[i][j]=U[i][j];
    }
  }

  free_dmatrix(L,sta,end-1,sta,end-1);
  free_dmatrix(U,sta,end-1,sta,end-1);
}


//LU分解
double *LU_Decomp(double **L,double **U,double *B,int M){
  printf("LU\n");
  int sta=1;
  int end=sta+M;
  double *W_out=dvector(sta,end-1);//alloc_vector(M);
  // double **L=dmatrix(sta,end-1,sta,end-1);//alloc_matrix(M,M);
  // double **U=dmatrix(sta,end-1,sta,end-1);//alloc_matrix(M,M);

  // // 対角成分に0があったらピボットの実行
  // double check;
  // for(int i=sta;i<end;i++){
  //   if(Acoef[i][i]==0.0){
  //     check=1;
  //   }
  // }
  
  // if(check==1){
  //   pivod(Acoef,B,M);//ピボットがあるとうまくいかない
  // }

  // for(int i=sta;i<=end-1;i++){
  //   U[i][i]=1.0;
  //   for(int j=i;j<=end-1;j++){
  //     if(i==sta){
  //       L[j][i]=Acoef[j][i];
  //     }else{
  //       double sumU=0.0;
  //       for(int sumnumber=sta;sumnumber<i-1;sumnumber++){
  //         sumU+=L[i-1][sumnumber]*U[sumnumber][j];
  //       }
  //       U[i-1][j]=(Acoef[i-1][j]-sumU)/L[i-1][i-1];

  //       double sumL=0.0;
  //       for(int sumnumber=sta;sumnumber<i;sumnumber++){
  //         sumL+=L[j][sumnumber]*U[sumnumber][i];
  //       }
  //       L[j][i]=Acoef[j][i]-sumL;
  //     }
  //   }
  // }

  double *by=dvector(sta,end-1);
  //Ly=bを解く
  for(int j=sta;j<=end-1;j++){
    double sumby=0.0;
    for(int i=sta;i<j;i++){
      sumby+=L[j][i]*by[i];
    }
    by[j]=(B[j]-sumby)/L[j][j];//asymptotic expression (漸化式) of computing Ly=b
  }
  
  //Ux=y
  for(int i=end-1;i>=sta;i--){
    double sumx=0;
    for(int k=i+1;k<=end-1;k++){
      sumx+=U[i][k]*W_out[k];
    }
    W_out[i]=by[i]-sumx;//asymptotic expression (漸化式) of computing Ux=y
  }
  

  // //誤差の表示
  // double sum=0.0;
  // double *Ans=matrix_vector_product(Acoef,W_out,M);
  // for(int i=sta;i<=end-1;i++){
  //   sum+=pow(Ans[i]-B[i],2);
  // }
  // printf("Err=%f\n",sum);

  free_dvector(by,sta,end-1);

  return W_out;
}

