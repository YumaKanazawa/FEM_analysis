//advection_diffucion
#include "../main_mesh.h"//メッシュ作成の関数

#define delta_t 0.01//time step
#define D 0.0 //Diffusion Coefficient

//nonlinear of advection diffusion equation
double f(double x,double y){
    return 0.0;
}

//initial condition
double init(double x,double y){
    // double r=0.0,theta=pi/4;
    return 1.0-(x*x+y*y);
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
    ret[1]=x;
    ret[2]=-y;

    return ret;
}

double div_u(mesh_t *mesh,double x,double y){
    double h=1e-3;
    double *U=u(x,y),*U_x=u(x+h,y),*U_y=u(x,y+h);

    double dx=U_x[1]-U[1],dy=U_y[2]-U[2];
    // printf("%f\n",(dx+dy)/h);

    free_dvector(U,1,mesh->dim);
    free_dvector(U_x,1,mesh->dim);
    free_dvector(U_y,1,mesh->dim);

    return (dx+dy)/h;
}

//φiφjの積分
double Int(mesh_t *mesh,int i,int j,int Kl){
    double I=0.0;
    double S_Kl=area(mesh,Kl);//Kl番目の面積
    if(i==j){
        I=S_Kl/6;
    }else{
        I=S_Kl/12;
    }
    return I;
}

// //関数ポインタ
// typedef double (* Func)(mesh_t,int,double,double,int);
// Func df=phi(mesh,int,double,double,int);

double Int_div_five(mesh_t *mesh,int i,int j,int Kl){
    int dim=mesh->dim;
    double S=area(mesh,Kl);
    double **Coord=NumInt_deg_five(mesh,Kl);//積分点の定義,５次
    int N=7;

    double *Coord_i=coef_plate_grad(mesh,Kl,i);
    double *Coord_j=coef_plate_grad(mesh,Kl,j);

    double ret=0.0;
    double w;
    for(int num=1;num<=N;num++){
        double x_num_I=Coord[num][1],y_num_I=Coord[num][2];
        if(i==1){
            w=9.0/40.0;
        }else if(2<=i && i<=4){
            w=(155.0-sqrt(15))/1200.0;
        }else{
            w=(155.0+sqrt(15))/1200.0;
        }

        
        double *U=u(x_num_I,y_num_I);
        double div=div_u(mesh,x_num_I,y_num_I);


        double div_int=div*phi(mesh,i,x_num_I,y_num_I,Kl)*phi(mesh,j,x_num_I,y_num_I,Kl);//(∇・u)φiφj
        // printf("∇・u=%f,phi_i=%f,phi_j%f\n",div,phi(mesh,i,x_num_I,y_num_I,Kl),phi(mesh,j,x_num_I,y_num_I,Kl));
        double phi_i_int=inner_product(1,dim,U,Coord_i);//u・∇φi
        // printf("u・∇φj=%f\n",phi_j_int);
        double phi_j_int=inner_product(1,dim,U,Coord_j);//u・∇φj
        // printf("u・∇φj=%f\n",phi_j_int);
        double phi_i=phi(mesh,i,x_num_I,y_num_I,Kl);//φi
        double phi_j=phi(mesh,j,x_num_I,y_num_I,Kl);//φj

        double cover_int=phi_i_int*phi_j-div_int-phi_j_int*phi_i;
        ret+=w*cover_int;

        free_dvector(U,1,dim);
    }
    // if(ret!=0.0){
    //     printf("K=%d,ret=%f\n",Kl,area(mesh,Kl)*ret);
    // }

    free_dmatrix(Coord,1,N,1,mesh->np);
    free_dvector(Coord_i,0,dim);
    free_dvector(Coord_j,0,dim);

    return S*ret;
}


// double Int_div_two(mesh_t *mesh,int i,int j,int Kl){
//     double **Coord=NumInt_deg_two(mesh,Kl);//2次の数値積分
//     int N=mesh->dim+1;
//     double ret=0.0;
//     double w=1/(N);
//     for(int num=1;num<=N;num++){
//         double x_num_I=Coord[num][1],y_num_I=Coord[num][2];
//         ret+=w*div_u(x_num_I,y_num_I)*phi(mesh,i,x_num_I,y_num_I,Kl)*phi(mesh,j,x_num_I,y_num_I,Kl);
//     }
//     free_dmatrix(Coord,1,N,1,mesh->np);
//     return area(mesh,Kl)*ret;
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
        for(int i=1;i<=n;i++){
            for(int j=1;j<=n;j++){
                int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目の節点番号

                double *C_i=coef_plate_grad(mesh,l,ver1);//φiの勾配
                double *C_j=coef_plate_grad(mesh,l,ver2);//φjの勾配

                /*======ここを問題の弱形式に応じて変更する========*/
                double Int_ij=Int(mesh,i,j,l);
                double Int_div_ij=Int_div_five(mesh,i,j,l);
                

                double Discreate_WeakForm=Int_ij+(0.5)*delta_t*Int_div_ij+D*delta_t*inner_product(1,dim,C_i,C_j)*S_Kl;
                /*==========================================*/
                A[ver1][ver2]+=Discreate_WeakForm;//弱形式の離散結果
                
                free_dvector(C_i,0,dim);
                free_dvector(C_j,0,dim);
            }
        }
    }   
    return A;
}

double *out_force(mesh_t *mesh,double *u_old){
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
        double **M=dmatrix(1,n,1,n);//要素質量行列
        double *f_vector=dvector(1,n);
        double *u_old_vector=dvector(1,n);
        /*===================剛性行列の作成==========================*/
        for(int i=1;i<=n;i++){  
            for(int j=1;j<=n;j++){
                // int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目の節点番号
                /*==========∫φiφjdxの計算結果代入==============*/
                double I=Int(mesh,i,j,l);
                /*==========================================*/
                M[i][j]+=I;
            }
        }
        /*=========================================================*/

        /*==================要素質量外力ベクトル======================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            double x=npxy[ver1][1],y=npxy[ver1][2];
            f_vector[i]=f(x,y);
            u_old_vector[i]=u_old[ver1];
        }
        double *rhs=matrix_vector_product(M,f_vector,n);
        double *rhs_1=matrix_vector_product(M,u_old_vector,n);
        /*=========================================================*/

        /*================返すベクトルに値を代入していく================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            return_vector[ver1]+=(delta_t*rhs[i]+rhs_1[i]);
        }
        /*=========================================================*/

        free_dvector(f_vector,1,n);
        free_dvector(u_old_vector,1,n);
        free_dmatrix(M,1,n,1,n);
        free_dvector(rhs,1,np);
        free_dvector(rhs_1,1,np);
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

            for(int i=1;i<=Neumman_dim;i++){
                for(int j=1;j<=Neumman_dim;j++){
                    int value;
                    if(i==j){
                        value=2;
                    }else{
                        value=1;
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

    return return_vector;
}

int main(int argc,char *argv[]){
    if(argc < 3){printf("Usage:./advection_diffusion dim ../Mesh/mesh01.msh\n"); exit(1);}//実行の仕方
    /*==================構造体のデータ読み込み==========================*/
    mesh_t mesh;
    //alloc(and scan)
    alloc_scan_mesh(&mesh,argv[1],argv[2]);//argv[1]=次元の数　argv[2]=meshファイルの名前
    // int dim=mesh.dim;
    // int n=mesh.n;
    int np=mesh.np;
    // int ne=mesh.ne;
    // int nb=mesh.nb;
    double **npxy;
    // int **elnp;
    // int **bound;
    npxy = mesh.npxy;
    // elnp = mesh.elnp;
    // bound = mesh.bound;
    /*==============================================================*/

    //初期条件の代入
    double *u_old=dvector(1,mesh.np);
    for(int i=1;i<=mesh.np;i++){
        double x=mesh.npxy[i][1],y=mesh.npxy[i][2];
        u_old[i]=init(x,y);
    }

    double **A=Al(&mesh);//剛性行列

    double *err=dvector(1,np);//空のベクトル
    Diriclet(&mesh,A,err);//ディリクレ境界条件の挿入(行列のみ)
    // matrix_vector_print(mesh,A,err);
    free_dvector(err,1,np);

    double **L=dmatrix(1,np,1,np);
    double **U=dmatrix(1,np,1,np);
    LU(A,np,L,U);  

    

    for(int T=0;T<100;T++){//時刻Tにおいて解を求める
        double *RHS=out_force(&mesh,u_old);//要素質量ベクトル

        double **ERR=dmatrix(1,np,1,np);
        Diriclet(&mesh,ERR,RHS);//ディリクレ境界条件の挿入(ベクトルのみ)
        free_dmatrix(ERR,1,np,1,np);

        // matrix_vector_print(mesh,A,RHS);

        /*Data output*/
        char str[200];
        sprintf(str,"figure/mesh%d.dat",T);
        // make_mesh_data_for_GLSC(mesh,u_old,str);//gnuplot用のファイル作成
        make_result_data_for_GLSC(&mesh,u_old,str);//GLSC用のデータ出力

        printf("t=%f,|u|=%f\n",delta_t*T,vector_norm1(u_old,1,mesh.np));
    
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