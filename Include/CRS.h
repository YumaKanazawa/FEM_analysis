#include "functions.h"
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
  double eps=pow(10,-6);//初期解

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
