//特性曲線用のプログラム
#include "../main_mesh.h"

#define D 0.01//(1.25*pow(10,-4))//Diffusion Coefficient
#define N 32
#define Nt (5*pow(N,2)/2)//時間方向の分割数
#define T_max M_PI
#define delta_t (T_max/Nt) //time step;

double *u(double x,double y){
    int dim=2;
    double *ret=dvector(1,dim);
    //ノイマン条件が与えられている境界上でu・n=0を満たすように作る
    //数値積分が面倒なので．．．
    ret[1]=-y;
    ret[2]=x;

    return ret;
}

//移流拡散の初期座標
#define x0 0.6
#define y0 0


//方程式の解析解
double u_exa(double x,double y,double t){
    double x1_t=x*cos(t)+y*sin(t);
    double x2_t=-x*sin(t)+y*cos(t);
    double sigma=0.01;

    double coef=(sigma+4*D*t);
    double exp_in=pow(x1_t-x0,2)+pow(x2_t-y0,2);

    return sigma/coef*exp(-exp_in/coef);
}

double phi_ij(mesh_t *mesh,int Kl,int i,int j,double x,double y){
    int n=mesh->n;
    return 0.0;
}

//nonlinear of advection diffusion equation
double f(double x,double y){
    return 0.0;
}

//initial condition
double init(double x,double y){
    return u_exa(x,y,0.0);
}

//ディリクレ境界条件
double g(double x,double y){
    return 0.0*(x+y);
}

//ノイマン境界条件
double g1(double x,double y){
    return 0.0*(x+y);
}

//熱方程式の弱形式
double charactic_advection(mesh_t *mesh,int Kl,int i,int j){
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


//積分の近似
double upstream_Integration(mesh_t *mesh,int Kl,int i,int j){
    int dim=mesh->dim;
    // double **npxy;
    // int **elnp,**bound;
    // npxy = mesh->npxy;
    // elnp = mesh->elnp;

    int N_5=7;
    double w;
    double ret=0.0;
    double S=area(mesh,Kl);

    for(int num=1;num<=N_5;num++){
        double *x=Pi(mesh,Kl,num);//Klのnum番目の積分点の座標
        double *U=u(x[1],x[2]);//積分点上の流速

        int upstream=search_past_point(mesh,x[1],x[2],U,delta_t);//xの上流点が含まれる要素番号

        free_dvector(U,1,dim);

        if(num==1){
            w=9.0/40.0;
        }else if(2<=num && num<=4){
            w=(155.0-sqrt(15))/1200.0;
        }else{
            w=(155.0+sqrt(15))/1200.0;
        }

        double *x1=Pi(mesh,upstream,num);//Klのnum番目の上流点の座標
        
        double phi_i=phi(mesh,i,x1[1],x1[2],Kl);
        double phi_j=phi(mesh,j,x[1],x[2],Kl);
        ret+=w*phi_i*phi_j;
        
        free_dvector(x,1,dim);
        free_dvector(x1,1,dim);
    }

    return S*ret;
}

//右辺の行列
double *out_force_charcterize(mesh_t *mesh,double *u_old){
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

    for(int l=1;l<=ne;l++){
        //要素KlがPi,Pj,Pkで構成されている
        double **M=dmatrix(1,n,1,n);//要素質量行列
        double *f_vector=dvector(1,n);
        double *u_old_vector=dvector(1,n);

        /*========外力項の計算=========*/
        /*===================剛性行列の作成==========================*/
        for(int i=1;i<=n;i++){  
            for(int j=1;j<=n;j++){
                // int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目の節点番号
                /*==========∫φiφjdxの計算結果代入==============*/
                double I=Int(mesh,i,j,l);//i->Pi,j->Pjとして考えている．
                /*==========================================*/
                M[i][j]+=I;
            }
        }

        /*=========================================================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];//l番目の要素のi番目の接点
            double x=npxy[ver1][1],y=npxy[ver1][2];
            f_vector[i]=delta_t*f(x,y);
        }

        double *rhs=matrix_vector_product(M,f_vector,n);//外力項の計算結果
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            return_vector[ver1]+=rhs[i];
        }
        /*=============================*/


        /*===================上流点の要素内における積分の計算======================*/
        double **M_up=dmatrix(1,n,1,n);//要素質量行列
        for(int i=1;i<=n;i++){  
            for(int j=1;j<=n;j++){
                // int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目の節点番号
                /*==========∫φiφjdxの計算結果代入==============*/
                double I_up=upstream_Integration(mesh,l,i,j);//i->Pi,j->Pjとして考えている．
                /*==========================================*/
                M_up[i][j]+=I_up;
            }
        }

        /*==================要素質量外力ベクトル======================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];//l番目の要素のi番目の接点
            u_old_vector[i]=u_old[ver1];//上流点用の値格納，RHS(u_old,ver1,x,y);
        }
        double *rhs_1=matrix_vector_product(M_up,u_old_vector,n);
        /*=========================================================*/

        /*================返すベクトルに値を代入していく================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            return_vector[ver1]+=rhs_1[i];//上流点の位置はそれぞれで分けなくてはならない
        }
        /*=========================================================*/

        free_dvector(f_vector,1,n);
        free_dvector(u_old_vector,1,n);
        free_dmatrix(M,1,n,1,n);
        free_dmatrix(M_up,1,n,1,n);
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

    double **dummy=dmatrix(1,np,1,np);
    Diriclet(mesh,dummy,return_vector);
    free_dmatrix(dummy,1,np,1,np);

    return return_vector;
}


int main(int argc,char *argv[]){
    weak weak_form=&charactic_advection;//ここを変える
    printf("N=%d,Δt=%f\n",N,delta_t);
    if(argc < 3){printf("Usage:./charactrize dim ../Mesh/mesh01.msh\n"); exit(1);}//実行の仕方
    /*==================構造体のデータ読み込み==========================*/
    mesh_t mesh;
    //alloc(and scan)
    alloc_scan_mesh(&mesh,argv[1],argv[2]);//argv[1]=次元の数　argv[2]=meshファイルの名前
    int np=mesh.np;
    double **npxy=mesh.npxy;
    /*==============================================================*/


    // //初期条件の代入
    // double *u_old=dvector(1,np);
    // for(int i=1;i<=np;i++){
    //     double x=npxy[i][1],y=npxy[i][2];
    //     u_old[i]=init(x,y);
    // }
    // double **A=Al(weak_form,&mesh);//剛性行列(境界条件込み)

    // for(int T=0;T<=Nt;T+=1){//時刻Tにおいて解を求める
    //     double *RHS=out_force_charcterize(&mesh,u_old);
    //     printf("t=%d,|u|=%f\n",T,vector_norm1(u_old,1,np,1.0));

    //     char str[200];
    //     sprintf(str,"figure/mesh%d.dat",T);
    //     // make_mesh_data_for_gnuplot(mesh,u_old,str);//gnuplot用のファイル作成
    //     make_result_data_for_GLSC(&mesh,u_old,str);//GLSC用のデータ出力

    //     double *u=CG_CRS(A,RHS,np);//解の計算
    //     for(int i=1;i<=np;i++){
    //         u_old[i]=u[i];
    //     }

    //     free_dvector(u,1,np);
    //     free_dvector(RHS,1,np);
    // }

    // free_dmatrix(A,1,np,1,np);
    // free_dvector(u_old,1,np);


    int l=1;
    for(int i=1;i<=mesh.dim+1;i++){
        for(int j=1;j<=mesh.dim+1;j++){
            int ver1=mesh.elnp[l][i],ver2=mesh.elnp[l][j];
            printf("%f,%f",Int(&mesh,ver1,ver2,l),1./12);
        }
        printf("\n");
    }

    mesh_free(&mesh);//free;
    return 0;
}