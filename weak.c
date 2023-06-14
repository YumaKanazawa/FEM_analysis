#include "main_mesh.h"//メッシュ作成の関数

#define D 0.1//(1.25*pow(10,-4))//Diffusion Coefficient
#define N 32
#define Nt (5*pow(N,2)/2)//時間方向の分割数
#define T_max M_PI
#define delta_t (T_max/Nt) //time step;

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

//移流項の弱形式の離散化の関数
double div_ij(mesh_t *mesh,int Kl,int i,int j,double x,double y){
    //(∇・u)φiφjの関数

    int dim=mesh->dim;
    double **npxy=mesh->npxy;

    double *U_i=u(npxy[i][1],npxy[i][2]);
    double *U_j=u(npxy[j][1],npxy[j][2]);

    double *Coord_i=coef_plate_grad(mesh,Kl,i);//∇φi
    double div_i=inner_product(1,dim,U_i,Coord_i);//π_{Vh} ∇・u
    double phi_i=phi(mesh,i,x,y,Kl);//φi

    double *Coord_j=coef_plate_grad(mesh,Kl,j);//∇φj
    double div_j=inner_product(1,dim,U_j,Coord_j);//π_{Vh} ∇・u
    double phi_j=phi(mesh,j,x,y,Kl);//φj

    free_dvector(U_i,1,dim);
    free_dvector(U_j,1,dim);
    free_dvector(Coord_j,0,dim);
    free_dvector(Coord_i,0,dim);

    return 0.5*(div_i*phi_i)+0.5*(div_j*phi_j);
}


//移流項の弱形式の離散化の関数
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


//移流項の弱形式の離散化の関数
double div_ij_SUPG(mesh_t *mesh,int Kl,int i,int j,double x,double y){
    int dim=mesh->dim;
    double **npxy=mesh->npxy;

    double *U_i=u(npxy[i][1],npxy[i][2]);

    double *Coord_i=coef_plate_grad(mesh,Kl,i);//∇φi
    double div_i=inner_product(1,dim,U_i,Coord_i);//π_{Vh} ∇・u
    double phi_i=phi(mesh,i,x,y,Kl);//φi

    double ret=(div_i*phi_i);//(∇・u)φi

    double *U=u(x,y);
    double *Coord_j=coef_plate_grad(mesh,Kl,j);
    double phi_j_int=inner_product(1,dim,U,Coord_j);//u・∇φj


    free_dvector(U,1,dim);
    free_dvector(U_i,1,dim);
    free_dvector(Coord_i,0,dim);
    free_dvector(Coord_j,0,dim);

    return ret*phi_j_int;
}

//移流項の弱形式の離散化の関数(SUPG)
double phi_ij_SUPG(mesh_t *mesh,int Kl,int i,int j,double x,double y){
    int dim=mesh->dim;
    double *U=u(x,y);

    double *Coord_i=coef_plate_grad(mesh,Kl,i);
    double *Coord_j=coef_plate_grad(mesh,Kl,j);

    double phi_j=inner_product(1,dim,U,Coord_j);//u・∇φj
    double phi_i=inner_product(1,dim,U,Coord_i);//u・∇φi

    free_dvector(Coord_j,0,dim);
    free_dvector(Coord_i,0,dim);
    free_dvector(U,1,dim);

    return phi_i*phi_j;
}

double advect(mesh_t *mesh,int Kl,int i,int j){
    pfunc func=&phi_ij;//φi(u・∇φj）
    pfunc func_SUPG=&phi_ij_SUPG;//
    pfunc div_u=&div_ij;//(∇・u)φiφj
    pfunc div_u_SUPG=&div_ij_SUPG;//(∇・u)φi(u・∇φj)


    int dim=mesh->dim;
    int **elnp= mesh->elnp;
    // double **npxy=mesh->npxy;
    // int np=mesh->np;

    double S_Kl=area(mesh,Kl);//Kl番目の面積

    int ver1=elnp[Kl][i],ver2=elnp[Kl][j];//lを構成するi番目の節点番号
    double *C_i=coef_plate_grad(mesh,Kl,ver1);//φiの勾配
    double *C_j=coef_plate_grad(mesh,Kl,ver2);//φjの勾配


    //第一項
    double Int_ij=Int(mesh,i,j,Kl);//φiφjの積分
    double Int_ij_SUPG=Int_div_five(func,mesh,Kl,ver1,ver2);//∇φi(u・∇φj)の積分

    //第二項
    double d_id_j=inner_product(1,dim,C_i,C_j)*S_Kl;
    //P１近似ということを考慮し，二回微分の項は消す
    
    //第三項
    double Int_div_ij=Int_div_five(div_u,mesh,Kl,ver1,ver2);//(∇・u)φiφjの積分
    double Int_div_ij_SUPG=Int_div_five(div_u_SUPG,mesh,Kl,ver1,ver2);//(∇・u)φi(u・∇φj)の積分

    double Int_u_dot_c=Int_div_five(func,mesh,Kl,ver1,ver2);//∇φi(u・∇φj)の積分
    double Int_u_dot_c_SUPG=Int_div_five(func_SUPG,mesh,Kl,ver1,ver2);//(u・∇φi)(u・∇φj)の積分


    double ret=
    Int_ij+delta_t*D*d_id_j+delta_t*(Int_div_ij+Int_u_dot_c)+
    SUPG*Int_ij_SUPG+
    SUPG*delta_t*(Int_div_ij_SUPG+Int_u_dot_c_SUPG);

    free_dvector(C_i,0,dim);
    free_dvector(C_j,0,dim);

    return ret;
}

//右辺のベクトル離散化
double RHS_advect(double *u_p,int i,double x,double y){
    return u_p[i]+delta_t*f(x,y);
}
/*========================================================================*/


/*=========================poisson=================================*/
double poisson(mesh_t *mesh,int Kl,int i,int j){
    int dim=mesh->dim;
    int **elnp= mesh->elnp;

    double S_Kl=area(mesh,Kl);//Kl番目の面積

    int ver1=elnp[Kl][i],ver2=elnp[Kl][j];//lを構成するi番目の節点番号
    double *C_i=coef_plate_grad(mesh,Kl,ver1);//φiの勾配
    double *C_j=coef_plate_grad(mesh,Kl,ver2);//φjの勾配

    double ret=D*inner_product(1,dim,C_i,C_j)*S_Kl;

    free_dvector(C_i,0,dim);
    free_dvector(C_j,0,dim);

    return ret;
}

double RHS_poisson(double *u_p,int i,double x,double y){
    /*注意分消すためだけの動作*/
    for(int j=i;j<=1;j++)u_p[i]=0.0;

    return f(x,y);
}

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
//右辺のベクトル離散化
double RHS_heat(double *u_p,int i,double x,double y){
    return u_p[i]+delta_t*f(x,y);
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


