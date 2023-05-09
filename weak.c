#include "main_mesh.h"//メッシュ作成の関数

#define D 1.25*pow(10,-4)//Diffusion Coefficient
#define N 80
#define Nt (5/2)*sqrt(sqrt(2)*N)//時間方向の分割数
#define T_max M_PI
#define delta_t T_max/Nt //time step;

/*=========================移流拡散=================================*/
/*ここは手動で入力する*/
double *u(double x,double y){
    int dim=2;
    double *ret=dvector(1,dim);
    //ノイマン条件が与えられている境界上でu・n=0を満たすように作る
    //数値積分が面倒なので．．．
    ret[1]=-y;
    ret[2]=x;

    return ret;
}

//φi(u・∇φj)の定義(関数)
double phi_ij(mesh_t *mesh,int Kl,int i,int j,double x,double y){
    int dim=mesh->dim;
    double *U=u(x,y);

    double *Coord_j=coef_plate_grad(mesh,Kl,j);

    double phi_j_int=inner_product(1,dim,U,Coord_j);//u・∇φj
    double phi_i=phi(mesh,i,x,y,Kl);//φi

    free_dvector(Coord_j,0,dim);
    free_dvector(U,1,dim);

    return phi_i*phi_j_int;
}

double advect(mesh_t *mesh,int Kl,int i,int j){
    pfunc func=&phi_ij;
    int dim=mesh->dim;
    int **elnp= mesh->elnp;

    double S_Kl=area(mesh,Kl);//Kl番目の面積

    int ver1=elnp[Kl][i],ver2=elnp[Kl][j];//lを構成するi番目の節点番号
    double *C_i=coef_plate_grad(mesh,Kl,ver1);//φiの勾配
    double *C_j=coef_plate_grad(mesh,Kl,ver2);//φjの勾配

    double Int_ij=Int(mesh,i,j,Kl);
    double Int_div_ij=Int_div_five(func,mesh,Kl,ver1,ver2);
    double Int_div_ji=Int_div_five(func,mesh,Kl,ver2,ver1);

    double ret=Int_ij+0.5*delta_t*(Int_div_ij-Int_div_ji)+D*delta_t*inner_product(1,dim,C_i,C_j)*S_Kl;

    free_dvector(C_i,0,dim);
    free_dvector(C_j,0,dim);

    return ret;
}
/*=========================移流拡散=================================*/


/*=========================拡散====================================*/
//熱方程式の弱形式
double heat(mesh_t *mesh,int Kl,int i,int j){
    int dim=mesh->dim;
    int **elnp= mesh->elnp;

    double S_Kl=area(mesh,Kl);//Kl番目の面積

    int ver1=elnp[Kl][i],ver2=elnp[Kl][j];//lを構成するi番目の節点番号
    double *C_i=coef_plate_grad(mesh,Kl,ver1);//φiの勾配
    double *C_j=coef_plate_grad(mesh,Kl,ver2);//φjの勾配

    double Int_ij=Int(mesh,i,j,Kl);

    double ret=Int_ij+D*delta_t*inner_product(1,dim,C_i,C_j)*S_Kl;

    free_dvector(C_i,0,dim);
    free_dvector(C_j,0,dim);

    return ret;
}


/*=========================波動方程式=================================*/
//波動方程式の弱形式の離散化
double wave(mesh_t *mesh,int Kl,int i,int j){
    int dim=mesh->dim;
    int **elnp= mesh->elnp;

    double S_Kl=area(mesh,Kl);//Kl番目の面積

    int ver1=elnp[Kl][i],ver2=elnp[Kl][j];//lを構成するi番目の節点番号
    double *C_i=coef_plate_grad(mesh,Kl,ver1);//φiの勾配
    double *C_j=coef_plate_grad(mesh,Kl,ver2);//φjの勾配

    double Int_ij=Int(mesh,i,j,Kl);

    double ret=Int_ij+D*pow(delta_t,2)*inner_product(1,dim,C_i,C_j)*S_Kl;

    free_dvector(C_i,0,dim);
    free_dvector(C_j,0,dim);

    return ret;
}

/*=========================右辺の離散化=================================*/
//右辺のベクトル離散化
double RHS_right(double *u_p,int i,double x,double y){
    return u_p[i]+delta_t*f(x,y);
}