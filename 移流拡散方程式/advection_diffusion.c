//advection_diffucion
#include "../main_mesh.h"//メッシュ作成の関数

//nonlinear of advection diffusion equation
double f(double x,double y){
    return 1.0;
}

//initial condition
double init(double x,double y){
    // double r=0.0,theta=pi/4;
    return exp(-50*(pow(x,2)+pow(y+1,2)));
}

//ディリクレ境界条件
double g(double x,double y){
    return init(x,y);
}

//ノイマン境界条件
double g1(double x,double y){
    return 0.0;
}

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

//φi(u・∇φj)の定義
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
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    // int n=mesh->n;
    // int np=mesh->np;
    // int ne=mesh->ne;
    // int nb=mesh->nb;
    // double **npxy;
    int **elnp;//,**bound;
    // npxy = mesh->npxy;
    elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

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

double heat(mesh_t *mesh,int Kl,int i,int j){
    // pfunc func=&phi_ij;//数値積分の関数定義

    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    // int n=mesh->n;
    // int np=mesh->np;
    // int ne=mesh->ne;
    // int nb=mesh->nb;
    // double **npxy;
    int **elnp;//,**bound;
    // npxy = mesh->npxy;
    elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

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

double wave(mesh_t *mesh,int Kl,int i,int j){
    // pfunc func=&phi_ij;//数値積分の関数定義

    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    // int n=mesh->n;
    // int np=mesh->np;
    // int ne=mesh->ne;
    // int nb=mesh->nb;
    // double **npxy;
    int **elnp;//,**bound;
    // npxy = mesh->npxy;
    elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

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

//右辺のベクトル離散化
double RHS_right(double *u_p,int i,double x,double y){
    return u_p[i]+delta_t*f(x,y);
}


int main(int argc,char *argv[]){
    weak weak_form=&advect;//ここを変える
    out pre_sol=&RHS_right;//右辺のベクトルの離散化

    if(argc < 3){printf("Usage:./PDE_name dim ../Mesh/mesh01.msh\n"); exit(1);}//実行の仕方
    /*==================構造体のデータ読み込み==========================*/
    mesh_t mesh;
    //alloc(and scan)
    alloc_scan_mesh(&mesh,argv[1],argv[2]);//argv[1]=次元の数　argv[2]=meshファイルの名前
    // int dim=mesh.dim;
    // int n=mesh.n;
    int np=mesh.np;
    // int ne=mesh.ne;
    // int nb=mesh.nb;
    // double **npxy;
    // int **elnp;
    // int **bound;
    // npxy = mesh.npxy;
    // elnp = mesh.elnp;
    // bound = mesh.bound;
    /*==============================================================*/

    //初期条件の代入
    double *u_old=dvector(1,mesh.np);
    for(int i=1;i<=mesh.np;i++){
        double x=mesh.npxy[i][1],y=mesh.npxy[i][2];
        u_old[i]=init(x,y);
    }

    double **A=Al(weak_form,&mesh);//剛性行列

    double *err=dvector(1,np);//空のベクトル
    Diriclet(&mesh,A,err);//ディリクレ境界条件の挿入(行列のみ)
    free_dvector(err,1,np);

    double **L=dmatrix(1,np,1,np);
    double **U=dmatrix(1,np,1,np);
    LU(A,np,L,U);  

    for(int T=0;T<100;T++){//時刻Tにおいて解を求める

        double *RHS=out_force(pre_sol,&mesh,u_old);//要素質量ベクトル

        double **ERR=dmatrix(1,np,1,np);
        Diriclet(&mesh,ERR,RHS);//ディリクレ境界条件の挿入(ベクトルのみ)
        free_dmatrix(ERR,1,np,1,np);

        // matrix_vector_print(mesh,A,RHS);

        /*Data output*/
        char str[200];
        sprintf(str,"figure/mesh%d.dat",T);
        // make_mesh_data_for_GLSC(mesh,u_old,str);//gnuplot用のファイル作成
        make_result_data_for_GLSC(&mesh,u_old,str);//GLSC用のデータ出力

        printf("t=%f,|u|=%f\n",delta_t*T,vector_norm1(u_old,1,mesh.np,1.0));
    
        double *u=LU_Decomp(L,U,RHS,mesh.np);//解の計算

        //解の更新
        printf("update solution \n");
        for(int i=1;i<=mesh.np;i++){
            u_old[i]=u[i];
        }
        
        free_dvector(u,1,mesh.np);
        free_dvector(RHS,1,mesh.np);
    }

    free_dmatrix(A,1,mesh.np,1,mesh.np);
    free_dmatrix(L,1,mesh.np,1,mesh.np);
    free_dmatrix(U,1,mesh.np,1,mesh.np);
    free_dvector(u_old,1,mesh.np);

    mesh_free(&mesh);//free

    return 0;
}

// typedef double (*test_func)(double x,double y);

// double a(double x,double y){
//     return x;
// }

// double da(test_func func,double x){
//     double h=0.001;
//     double y=0.0;
//     double df=func(x+h,y)-func(x-h,y);

//     return df/(2*h);
// }

// double Sa(test_func func,double a,double b){
//     double y=0.0;
//     double N=1000;
//     double sum=0.0;
//     for(int i=0;i<=N;i++){
//         double xi=a+(b-a)*(i/N);
//         sum+=func(xi,y);
//     }
//     return sum/N;
// }

// int main(void){
//     test_func func=&a;
//     printf("%f\n",da(a,1));
//     return 0;
// }