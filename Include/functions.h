#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define pi atan(1)*4
#define Max(a,b) (((a)>(b))?(a):(b))
#define Min(a,b) (((a)<(b))?(a):(b))
#define Sgn(a) (((a)==0.0)?0.0:((a)/(fabs(a))))

// #define sta 1

//argument of (x[1],x[2]) in [0,2\pai) or -1 (the origin)
double argument( double *x ){
  double eps=1.0e-10, a; 
  if ( x[1]*x[1]+x[2]*x[2] < eps ){
    return -1.0;
  }
  if ( x[1] >= fabs(x[2]) ){
    a = atan(x[2]/x[1]);
    return (a>=0.0)?a:(a+2*pi);
  }
  if ( x[2] >= fabs(x[1]) ){
    a = atan(-x[1]/x[2]);
    return a + 0.5*pi;
  }
  if ( x[1] <= - fabs(x[2]) ){
    a = atan(x[2]/x[1]);
    return a+pi;
  }
  if ( x[2] <= - fabs(x[1]) ){
    a = atan(-x[1]/x[2]);
    return a + 1.5*pi;
  }
  exit(0);
}


double *dvector(int i, int j){
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

void free_dvector(double *a,int i,int j){
  if(i!=j){
    free((double*)(a+i));
  }else{
    printf("Non Length!\n");
    exit(1);
  }
  a=NULL;
}

int *ivector(int i,int j){
  int *a;
  if ( (a=(int *)malloc( ((j-i+1)*sizeof(int))) ) == NULL ){
    printf("ivector: memory allocation is failed \n");
    exit(1);
  }  
  for(int n=0;n<j-i+1;n++){
    *(a+n)=-1;
  } 
  return (a-i);
}

void free_ivector(int *a,int i,int j){
  if(i!=j){
    free(a+i);
  }else{
    printf("Non Length!\n");
    exit(1);
  }
  a=NULL;
}

double **dmatrix(int nr1,int nr2,int nl1,int nl2){
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

void free_dmatrix(double **a,int nr1,int nr2,int nl1,int nl2){
  if(nl2!=nl1 && nr2!=nr1){
    for (int i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
    free((void *)(a+nr1));
  }else{
    printf("Non Matrix dimension!\n");
  }
  return;
}


int **imatrix(int nr1,int nr2,int nl1,int nl2){
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

void free_imatrix(int **a, int nr1,int nr2,int nl1,int nl2){
  if(nl2!=nl1 && nr2!=nr1){
    for (int i = nr1 ; i <= nr2 ; i++) free((void *)(a[i]+nl1));
    free((void *)(a+nr1));
  }else{
    printf("Non Matrix dimension!\n");
  }
  return;
}

double *matrix_vector_product(double **a, double *b ,int n){
  int sta=1;
  int end=sta+n-1;
  double *c=dvector(sta,end);

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

//Lpノルム
double vector_norm1( double *a, int m, int n,double p){
  int i; 
  double norm = 0.0;
  double num=0;
  for ( i = m; i <= n; i++ ){
    num+=1.0;
    norm += pow(fabs(a[i]),p);
  }
  return pow(norm/num,1/p); 
}

//２本のベクトルからなす角を計算する
double arg(double *x1,double *x2,int M){
  int sta=1;
  int end=sta+M;
  return acos(inner_product(0,M,x1,x2)/(vector_norm1(x1,sta,end,2)*vector_norm1(x2,sta,end,2)));
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

double *CG(double **A,double *b,int M){//M次元正定値対称行列
  int sta=1;
  int end=sta+M-1;

  double *x=dvector(sta,end);
  double eps=pow(10,-6);//初期解

  double *r=dvector(sta,end);//-∇f(=r0)
  double *p=dvector(sta,end);//directionベクトル

  double *AT=matrix_vector_product(A,x,M);
  for(int i=sta;i<=end;i++){
    r[i]=b[i]-AT[i];//r0=b-AT
    p[i]=r[i];//d0=r0
  }

  free_dvector(AT,sta,end);

  double gamma=inner_product(sta,end,r,r);//γ=(r,r)

  //k番目について
  int count=0;
  while(count<pow(10,10)){//反復回数

    double *q=matrix_vector_product(A,p,M);//q=Ap
    double alpha=gamma/inner_product(sta,end,p,q);//γ/(p,q)

    for(int i=sta;i<=end;i++){
      x[i]=x[i]+alpha*p[i];//x=x+αp
      r[i]=r[i]-alpha*q[i];//r=r-αq
    }

    if(vector_norm1(r,sta,end,2.0)<eps*vector_norm1(b,sta,end,2.0)){
      break;
    }

    double delta=inner_product(sta,end,r,r);//δ=(r,r)

    double beta=delta/gamma;//δ/γ

    for(int i=sta;i<=end;i++)p[i]=r[i]+beta*p[i];
    
    gamma=delta;

    count+=1;
    free_dvector(q,sta,end);
  }

  free_dvector(r,sta,end);
  free_dvector(p,sta,end);

  return x;
}

/*===============疎行列格納形式=====================*/
typedef struct CRS{
  double *A;//non-zeroの要素配列
  int *ia;//行方向の要素
  int *ja;//列方向の要素

  int M;//行列のサイズ
  int A_l;//non-zero要素の個数
  int ia_l;//iaの長さ
  int ja_l;//jaの長さ
}CRS_t;


void CRS(CRS_t *CRS,double **A,int M){
  // printf("Make CRS of Matrix A\n");

  CRS->M=M;
  int sta=1;

  int al=0;
  for(int i=sta;i<sta+M;i++){
    for(int j=sta;j<sta+M;j++){
      if(A[i][j]!=0.0){
        al+=1;
      }
    }
  }

  CRS->A_l=al;//ベクトルの長さ
  CRS->ja_l=al;
  CRS->ia_l=M+1;
  int A_l=CRS->A_l;

  CRS->A=dvector(sta,A_l);
  CRS->ia=ivector(sta,M+1);
  CRS->ja=ivector(sta,A_l);

  // double *a=CRS->A;
  // int *ia=CRS->ia;
  // int *ja=CRS->ja;

  // double *p=&a[sta];//pでaを参照
  // int *p_j=&ja[sta];//pjでjaを参照
  // int *p_i=&ia[sta];//piでiaを参照

  double *p=&(CRS->A[sta]);//pでaを参照
  int *p_j=&(CRS->ja[sta]);//pjでjaを参照
  int *p_i=&(CRS->ia[sta]);//piでiaを参照

  int count=0;
  for(int i=sta;i<sta+M;i++){
    for(int j=sta;j<sta+M;j++){
      double aij=A[i][j];
      if(aij!=0.0){
        *(p+count)=aij;
        *(p_j+count)=j;
        // (CRS->A[count+sta])=aij;
        // (CRS->ja[count+sta])=j;
        count++;
      }//０でない値を格納
    }
  }

  //iaの作成
  int locate_ja=sta;
  int loc_pi=sta;
  for(int col=sta;col<sta+M;col++){
    *p_i=locate_ja;
    p_i++;

    int a_l_j=0;
    for(int j=sta;j<sta+M;j++){
      if(A[col][j]!=0.0){
        a_l_j++;//col行目の0でない個数のトータル
      }
    }
    locate_ja+=a_l_j;//staから各行ごとに0でないものの個数分追加していく
  }
  *p_i=A_l+1;
}

void free_CRS(CRS_t *CRS,int M,int N){
  int sta=1;
  int A_l=CRS->A_l;
  // printf("\nCRS_free: start ....\n"); fflush(stdout);
  free_dvector(CRS->A,sta,A_l);
  free_ivector(CRS->ia,sta,M);
  free_ivector(CRS->ja,sta,A_l);
}

double *matrix_vector_product_CRS(CRS_t *CRS_A,double *x,int n){

  // CRS_t CRS_A;
  // CRS(&CRS_A,A,n);
  // double *a,int *ja,int *ia,
  
  double *a=CRS_A->A;
  int *ja=CRS_A->ja;
  int *ia=CRS_A->ia;

  int sta=1;
  double *ret=dvector(sta,sta+n);
  //Aのij成分にアクセスする
  //i行目の成分はia[i]<k<ia[i+1]の間にある

  for(int col=sta;col<sta+n;col++){
    double sum=0.0;
    for(int k=ia[col];k<ia[col+1];k++){
      sum+=a[k]*x[ja[k]];
    }
    ret[col]=sum;
  }
  // free_CRS(&CRS_A,n,n);

  return ret;
}

double *CG_CRS(double **A,double *b,int M){//M次元正定値対称行列
  CRS_t CRS_A;
  CRS(&CRS_A,A,M);

  int sta=1;
  int end=sta+M-1;

  double *x=dvector(sta,end);
  double eps=pow(10,-7);//初期解

  double *r=dvector(sta,end);//-∇f(=r0)
  double *p=dvector(sta,end);//directionベクトル

  double *AT=matrix_vector_product_CRS(&CRS_A,x,M);
  for(int i=sta;i<=end;i++){
    r[i]=b[i]-AT[i];//r0=b-AT
    p[i]=r[i];//d0=r0
  }
  double gamma=inner_product(sta,end,r,r);//γ=(r,r)

  //k番目について
  int count=0;
  while(count<pow(10,10)){//反復回数
    double *q=matrix_vector_product_CRS(&CRS_A,p,M);//q=Ap
    double alpha=gamma/inner_product(sta,end,p,q);//γ/(p,q)

    for(int i=sta;i<=end;i++){
      x[i]=x[i]+alpha*p[i];//x=x+αp
      r[i]=r[i]-alpha*q[i];//r=r-αq
    }

    if(vector_norm1(r,sta,end,2.0)<eps*vector_norm1(b,sta,end,2.0)){
      break;
    }

    double delta=inner_product(sta,end,r,r);//δ=(r,r)

    double beta=delta/gamma;//δ/γ

    for(int i=sta;i<=end;i++)p[i]=r[i]+beta*p[i];
    
    gamma=delta;

    count+=1;
    free_dvector(q,sta,end);
  }

  free_dvector(AT,sta,end);
  free_dvector(r,sta,end);
  free_dvector(p,sta,end);
  free_CRS(&CRS_A,M,M);


  return x;
}
