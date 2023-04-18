//Heat equation FEM
//メッシュによってスパイクが生じる
#include "../main_mesh.h"//メッシュ作成の関数


#define D 0.0//Diffusion param
#define delta_t 0.01//time step

//非線形項
double f(double x,double y,double t){
    return 0.0;//q2*t*(y-y*y)+2*D*(1+t*t);
}

// /*================厳密解==================*/
// double ue(double x,double y,double t){
//     return (1+t*t)*(y-y*y);
// }

//initial condition
double init(double x,double y){
    double x0=0.0,y0=0.0;
    
    double ret=1-pow(x-x0,2)-pow(y-y0,2);

    // // ret=exp(ret);
    // if(x0==x && y0==y){
    //     ret=1.0;
    // }else{
    //     ret=0.0;
    // }
    return ret;
}

//ディリクレ境界条件u=g
double g(double x,double y){
    return 0.0;
}

//ノイマン境界条件∂u/∂n=g1
double g1(double x,double y){
    return 0.0;
}

//要素剛性行列の作成
double **Al(mesh_t *mesh){
    printf("Al(要素剛性行列)\n");
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    int n=mesh->n;
    int np=mesh->np;
    int ne=mesh->ne;
    // int nb=mesh->nb;
    // double **npxy;
    int **elnp;//,**bound;
    // npxy = mesh->npxy;
    elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

    double **A=dmatrix(1,np,1,np);//要素剛性行列Aの初期化

    for(int l=1;l<=ne;l++){
        //要素KlがPi,Pj,Pkで構成されている
        double S_Kl=area(mesh,l);//Kl番目の要素の面積
        /*===================剛性行列の成分計算==========================*/
        for(int i=1;i<=n;i++){
            for(int j=1;j<=n;j++){
                int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目,j番目の節点番号

                double *C_i=coef_plate_grad(mesh,l,ver1);//φiの勾配
                double *C_j=coef_plate_grad(mesh,l,ver2);//φjの勾配

                /*======ここを問題の弱形式に応じて変更する========*/
                double Int=0.0;//φiφjの積分
                if(i==j){
                    Int=S_Kl/6.0;
                }else{
                    Int=S_Kl/12.0;
                }

                double Discreate_WeakForm=Int+D*delta_t*inner_product(1,dim,C_i,C_j)*S_Kl;//φiφj+(∇φi)・(∇φj)の積分
                /*==========================================*/
                A[ver1][ver2]+=Discreate_WeakForm;//弱形式の離散結果

                free_dvector(C_i,0,dim);
                free_dvector(C_j,0,dim);
            }
        }
    }
    return A;
}

double *out_force(mesh_t *mesh,double *u_old,double t){
    printf("out_force(外力項の離散化)\n");
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    int n=mesh->n;
    int np=mesh->np;
    int ne=mesh->ne;
    int nb=mesh->nb;
    double **npxy;
    int **elnp;
    int **bound;
    npxy = mesh->npxy;
    elnp = mesh->elnp;
    bound = mesh->bound;
    /*==============================================================*/

    double *return_vector=dvector(1,np);//右辺のベクトル

    /*========外力項の計算=========*/
    for(int l=1;l<=ne;l++){
        //要素KlがPi,Pj,Pkで構成されている
        double S_Kl=area(mesh,l);//l番目の面積
       
        
        double **M=dmatrix(1,n,1,n);//要素質量行列
        // double *f_vector=dvector(1,n);//ベクトルfの初期化

        double *u_old_vector=dvector(1,dim+1);

        /*====================剛性行列の作成==========================*/
        for(int i=1;i<=n;i++){  
            for(int j=1;j<=n;j++){
                // int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目の節点番号
                /*==========∫φiφjdxの計算結果代入==============*/
                double Int;
                if(i==j){
                    Int=S_Kl/6.0;
                }else{
                    Int=S_Kl/12.0;
                }
                /*==========================================*/
                M[i][j]=Int;
            }
        }
        /*========================================================*/

        /*==================要素質量外力ベクトル======================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            double x=npxy[ver1][1],y=npxy[ver1][2];//l番目の要素のx,y座標

            u_old_vector[i]=u_old[ver1]+delta_t*f(x,y,t);//∫φiφjdx u^{n-1}+Δt∫φiφjdx f(外力項と1時刻前の計算結果の積分)
        }
        // double *rhs=matrix_vector_product(M,f_vector,n);
        double *rhs=matrix_vector_product(M,u_old_vector,n);
        /*=========================================================*/

        /*================返すベクトルに値を代入していく================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            return_vector[ver1]+=rhs[i];//∫φiφjdx u^{n-1}+Δt∫φiφjdx f(外力項と1時刻前の計算結果の積分)
        }
        /*=========================================================*/

        // free_dvector(f_vector,1,n);
        free_dvector(u_old_vector,1,dim+1);
        free_dmatrix(M,1,n,1,n);
        // free_dvector(rhs,1,np);
        free_dvector(rhs,1,n);
    }
    

    /*=====================ノイマン条件の反映=====================*/
    for(int b=1;b<=nb;b++){
        int Gamma=2;//ノイマン条件の付与
        if(bound[b][dim+1]==Gamma){//ノイマン条件の与えられているb番目の接点
            int Neumman_dim=2;//２点間の次元を考える
            double **B=dmatrix(1,Neumman_dim,1,Neumman_dim);
            double *Neumman_vector=dvector(1,Neumman_dim);//g1を反映するベクトル

            int ver1=bound[b][1],ver2=bound[b][2];//ver1番目とver2の節点が境界上にある
           
            double Lij=L(mesh,ver1,ver2);///ver1番目とver2の節点上の距離

            for(int i=1;i<=Neumman_dim;i++){
                for(int j=1;j<=Neumman_dim;j++){
                    double value;
                    if(i==j){
                        value=2.0;
                    }else{
                        value=1.0;
                    }
                    B[i][j]=Lij*value/6.0;
                }
            }

            for(int i=1;i<=Neumman_dim;i++){
                int ver=bound[b][i];
                Neumman_vector[i]=D*delta_t*g1(npxy[ver][1],npxy[ver][2]);
            }

            double *rhs_Neumman=matrix_vector_product(B,Neumman_vector,Neumman_dim);

            for(int i=1;i<=Neumman_dim;i++){
                int ver1=bound[b][i];
                return_vector[ver1]+=rhs_Neumman[i];
            }

            free_dvector(rhs_Neumman,1,Neumman_dim);
            free_dmatrix(B,1,Neumman_dim,1,Neumman_dim);
            free_dvector(Neumman_vector,1,Neumman_dim);
        }
    }    
    /*=========================================================*/
    
    /*===========================*/
    return return_vector;
}

int main(int argc,char *argv[]){
    if(argc < 3){printf("Usage:./Heat dim mesh.msh\n"); exit(1);}//実行の仕方

    /*==================構造体のデータ読み込み==========================*/
    mesh_t mesh;
    //alloc(and scan)
    alloc_scan_mesh(&mesh,argv[1],argv[2]);//argv[1]=次元の数　argv[2]=meshファイルの名前
    int dim=mesh.dim;
    int n=mesh.n;
    int np=mesh.np;
    int ne=mesh.ne;
    int nb=mesh.nb;
    double **npxy;
    int **elnp;
    int **bound;
    npxy = mesh.npxy;
    elnp = mesh.elnp;
    bound = mesh.bound;
    /*==============================================================*/
    
    // /*初期条件の代入*/
    double *u_old=dvector(1,np);
    for(int i=1;i<=np;i++){
        double x=npxy[i][1],y=npxy[i][2];
        u_old[i]=init(x,y);
    }

    double **A=Al(&mesh);//剛性行列

    double *err=dvector(1,np);//空のベクトル
    Diriclet(&mesh,A,err);//ディリクレ境界条件の挿入(行列のみ)
    free_dvector(err,1,np);

    double **L=dmatrix(1,np,1,np),**U=dmatrix(1,np,1,np);
    LU(A,np,L,U);//係数行列の分解


    for(int T=0;T<100;T++){//時刻Tにおいて解を求める
        double *RHS=out_force(&mesh,u_old,T*delta_t);//要素質量ベクトル

        double **ERR=dmatrix(1,np,1,np);
        Diriclet(&mesh,ERR,RHS);//ディリクレ境界条件の挿入(ベクトルのみ)
        free_dmatrix(ERR,1,np,1,np);
        
        // make_coef_matrix(&mesh,A,u_old,T);//ファイル書き出し

        // matrix_vector_print(mesh,A,RHS);

        /*========================Data output============================*/
        char str[200];
        sprintf(str,"figure/mesh%d.dat",T);
        // make_mesh_data_for_gnuplot(mesh,u_old,str);//gnuplot用のファイル作成
        make_result_data_for_GLSC(&mesh,u_old,str);//GLSC用のデータ出力
        /*===============================================================*/

        printf("t=%f,|u|=%f\n",delta_t*T,vector_norm1(u_old,1,mesh.np));

        /*============================解の計算=========================*/
        double *u=LU_Decomp(L,U,RHS,np);//解の計算
        /*============================================================*/

        /*============================解の更新=========================*/
        for(int i=1;i<=np;i++)u_old[i]=u[i];
        /*============================================================*/
        free_dvector(u,1,np);
        free_dvector(RHS,1,np);
        
    }

    free_dvector(u_old,1,np);
    free_dmatrix(A,1,np,1,np);
    free_dmatrix(L,1,np,1,np);
    free_dmatrix(U,1,np,1,np);


    mesh_free(&mesh);// free

    return 0;
}
