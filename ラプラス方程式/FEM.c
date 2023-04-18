//poisson equation FEM

#include "../main_mesh.h"//メッシュ作成の関数

//ポアソン方程式の右辺
double f(double x,double y){
    return 1.0;
}

//ディリクレ境界条件
double g(double x,double y){
    return 0.0;
}

//ノイマン境界条件
double g1(double x,double y){
    return 0.0;
}

// //厳密解(球対称，Δu=1(domain),u=0(border))
// double u_exa(double x,double y){
//     double a=sqrt(pow(x,2)+pow(y,2));
//     return (a-1);
// }

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
        double S_Kl=area(mesh,l);
        /*===================剛性行列の成分計算==========================*/
        // double **A_l=dmatrix(1,n,1,n);
        for(int i=1;i<=n;i++){
            for(int j=1;j<=n;j++){
                int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目の節点番号

                double *C_i=coef_plate_grad(mesh,l,ver1);//φiの勾配
                double *C_j=coef_plate_grad(mesh,l,ver2);//φjの勾配

                /*======ここを問題の弱形式に応じて変更する========*/
                double Discreate_WeakForm=inner_product(1,dim,C_i,C_j)*S_Kl;
                /*==========================================*/
                A[ver1][ver2]+=Discreate_WeakForm;//弱形式の離散結果
                // A_l[i][j]=Discreate_WeakForm;
                free_dvector(C_i,0,dim);
                free_dvector(C_j,0,dim);
            }
        }

        /*=================行と列のそれぞれの和=========================*/
        // double sum_col[4]={0.0};
        // for(int i=1;i<=n;i++){
        //     double sum=0.0;
        //     for(int j=1;j<=n;j++){
        //         sum+=A_l[i][j];
        //         sum_col[j]+=A_l[i][j];
        //         printf("(%d,%d),%0.2f,",i,j,A_l[i][j]);
        //     }
        //     printf("%0.2f\n",sum);
        //     // printf("\n");
        // }
        // printf("%0.2f,%0.2f,%0.2f\n",sum_col[1],sum_col[2],sum_col[3]);
        // printf("\n");
        /*===========================================================*/

        // free_dmatrix(A_l,1,n,1,n);
        /*=========================================================*/
    }
    return A;
}

double *out_force(mesh_t *mesh){
    printf("out_force(外力項の離散化)\n");
    /*==================構造体のデータ読み込み==========================*/
    // int dim=mesh->dim;
    int n=mesh->n;
    int np=mesh->np;
    int ne=mesh->ne;
    int nb=mesh->nb;
    double **npxy;
    int **elnp,**bound;
    npxy = mesh->npxy;
    elnp = mesh->elnp;
    bound = mesh->bound;
    /*==============================================================*/

    double *return_vector=dvector(1,np);

    /*========外力項の計算=========*/
    for(int l=1;l<=ne;l++){
        //要素KlがPi,Pj,Pkで構成されている
        double S_Kl=area(mesh,l);
        
        double **M=dmatrix(1,n,1,n);//要素質量行列
        double *f_vector=dvector(1,n);
        /*===================剛性行列の作成==========================*/
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
                M[i][j]+=Int;
            }
        }
        /*=========================================================*/

        /*==================要素質量外力ベクトル======================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            f_vector[i]=f(npxy[ver1][1],npxy[ver1][2]);
        }
        double *rhs=matrix_vector_product(M,f_vector,n);
        /*=========================================================*/


        /*================返すベクトルに値を代入していく================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            return_vector[ver1]+=rhs[i];
        }
        /*=========================================================*/

        free_dvector(f_vector,1,n);
        free_dmatrix(M,1,n,1,n);
        free_dvector(rhs,1,n);
    }

    /*=====================ノイマン条件の反映=====================*/
    for(int b=1;b<=nb;b++){
        int Gamma=2;
        if(bound[b][n]==Gamma){
            int Neumman_dim=2;
            double **B=dmatrix(1,Neumman_dim,1,Neumman_dim);//二次元上の直線
            double *Neumman_vector=dvector(1,Neumman_dim);//g1を反映すベクトル

            int ver1=bound[b][1],ver2=bound[b][2];
            double Lij=L(mesh,ver1,ver2);

            /*辺質量行列Bの作成*/
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
                Neumman_vector[i]=g1(npxy[ver][1],npxy[ver][2]);
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
    if(argc < 3){printf("Usage:./FEM dim mesh.msh\n"); exit(1);}//実行の仕方

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

    /*============================LU分解===============================*/
    double **A=Al(&mesh);//剛性行列

    double *err=dvector(1,np);//空のベクトル
    Diriclet(&mesh,A,err);//ディリクレ境界条件の挿入(行列のみ)
    free_dvector(err,1,np);

    double **L=LU(A,np,"L");
    double **U=LU(A,np,"U");
    /*==============================================================*/

    /*==============================================================*/
    double *RHS=out_force(&mesh);//要素質量ベクトル

    double **ERR=dmatrix(1,np,1,np);
    Diriclet(&mesh,ERR,RHS);//ディリクレ境界条件の挿入(ベクトルのみ)
    free_dmatrix(ERR,1,np,1,np);
    /*==============================================================*/

    double *u=LU_Decomp(L,U,RHS,mesh.np);//解の評価


    /*========================Data output============================*/
    int T=0;
    char str[200];
    sprintf(str,"figure/mesh%d.dat",T);

    printf("|u|=%f\n",vector_norm1(u,1,mesh.np));

    // main
    // make_mesh_data_for_gnuplot(mesh,u);//gnuplot用のファイル作成
    make_result_data_for_GLSC(&mesh,u,str);
    /*================================================================*/

    free_dvector(u,1,mesh.np);
    free_dmatrix(A,1,mesh.np,1,mesh.np);
    free_dvector(RHS,1,mesh.np);
    free_dmatrix(L,1,mesh.np,1,mesh.np);
    free_dmatrix(U,1,mesh.np,1,mesh.np);
    
    mesh_free(&mesh);// free
    return 0;
}