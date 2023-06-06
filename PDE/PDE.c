//advection_diffucion
#include "../weak.c"

//移流拡散の初期座標
#define x0 0.25
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

//nonlinear of advection diffusion equation
double f(double x,double y){
    return 1.0;
}

//initial condition
double init(double x,double y){
    return u_exa(x,y,0.0);
}

//ディリクレ境界条件
double g(double x,double y){
    return 0.0;
}

//ノイマン境界条件
double g1(double x,double y){
    return 0.0;
}

//誤差のファイル書き出し
void error_write(double p,double err){
    char file_name[50];
    sprintf(file_name,"advection_diffusion_%d.txt",(int)p);

    FILE *file_open=fopen(file_name,"a");
    if(file_open==NULL){
        printf("file_open Error\n");
        exit(1);
    }
    double h=sqrt(2)/N;
    fprintf(file_open,"%f,%f\n",h,err);
    fclose(file_open);
}


//移流拡散(一回微分を含む方程式)
int main(int argc,char *argv[]){
    weak weak_form=&heat;//ここを変える
    out pre_sol=&RHS_heat;//右辺のベクトルの離散化

    // main_mesh.hのSUPGの値を0にする．
    printf("N=%d,Δt=%f\n",N,delta_t);
    if(argc < 3){printf("Usage:./PDE_name dim ../Mesh/mesh01.msh\n"); exit(1);}//実行の仕方
    /*==================構造体のデータ読み込み==========================*/
    mesh_t mesh;
    //alloc(and scan)
    alloc_scan_mesh(&mesh,argv[1],argv[2]);//argv[1]=次元の数　argv[2]=meshファイルの名前
    int np=mesh.np;
    /*==============================================================*/

    //初期条件の代入
    double *u_old=dvector(1,np);
    for(int i=1;i<=np;i++){
        double x=mesh.npxy[i][1],y=mesh.npxy[i][2];
        u_old[i]=init(x,y);
    }

    double **A=Al(weak_form,&mesh);//剛性行列(境界条件込み)

    // double **L=dmatrix(1,np,1,np);
    // double **U=dmatrix(1,np,1,np);
    // LU(A,np,L,U); 

    double L2n=0.0;//時刻nにおける解uのL2ノルム
    int Length=0;//点列の長さを入れる変数


    for(int T=0;T<=Nt;T+=1){//時刻Tにおいて解を求める

        double err_t=err_Lp(&mesh,u_old,2.0,T*delta_t);//誤差のL2ノルム
        // if(L2n<=err_t){
        //     L2n=err_t;
        // }
        L2n+=err_t;//ある時刻での解析解と数値解のLp誤差
        Length+=1;

        double *RHS=out_force(pre_sol,&mesh,u_old);

        /*============Data output=============*/
        char str[200];
        sprintf(str,"figure/mesh%d.dat",T);
        // make_mesh_data_for_gnuplot(mesh,u_old,str);//gnuplot用のファイル作成
        make_result_data_for_GLSC(&mesh,u_old,str);//GLSC用のデータ出力

        printf("t=%d,|u|=%f\n",T,vector_norm1(u_old,1,np,1.0));
        /*====================================*/

        /*============解の更新===========*/
        double *u=CG_CRS(A,RHS,np);//解の計算
        // double *u=LU_Decomp(L,U,RHS,np);
        printf("update solution \n");
        for(int i=1;i<=np;i++){
            u_old[i]=u[i];
            // printf("%f\n",u[i]);    
        }
        /*=============================*/
        
        free_dvector(u,1,np);
        free_dvector(RHS,1,np);

    }

    // printf("l^∞(l2)=%f\n",L2n);
    error_write(2,pow(L2n/Length,0.5));
    printf("l^2(l2)=%f\n",pow(L2n/Length,0.5));

    free_dmatrix(A,1,np,1,np);
    // free_dmatrix(L,1,np,1,np);
    // free_dmatrix(U,1,np,1,np);
    free_dvector(u_old,1,np);

    mesh_free(&mesh);//free

    return 0;
}


// //ポアソン方程式(時間微分を含まない)
// int main(int argc,char *argv[]){
//     weak weak_form=&poisson;//ここを変える
//     out pre_sol=&RHS_poisson;//右辺のベクトルの離散化

//     //main_mesh.hのSUPGの値を0にする．

//     printf("N=%d,Δt=%f\n",N,delta_t);
//     if(argc < 3){printf("Usage:./PDE_name dim ../Mesh/mesh01.msh\n"); exit(1);}//実行の仕方
//     /*==================構造体のデータ読み込み==========================*/
//     mesh_t mesh;
//     //alloc(and scan)
//     alloc_scan_mesh(&mesh,argv[1],argv[2]);//argv[1]=次元の数　argv[2]=meshファイルの名前
//     int np=mesh.np;
//     /*==============================================================*/

//     //初期条件の代入
//     double *u_old=dvector(1,np);
//     for(int i=1;i<=np;i++){
//         u_old[i]=0.0;
//     }

//     double **A=Al(weak_form,&mesh);//剛性行列(境界条件込み)
//     double *RHS=out_force(pre_sol,&mesh,u_old);

//     /*============解の更新===========*/
//     double *u=CG_CRS(A,RHS,np);//解の計算
//     // double *u=LU_Decomp(L,U,RHS,np);
//     for(int i=1;i<=np;i++){
//         u_old[i]=u[i];
//         printf("%f\n",u[i]);    
//     }
//     /*=============================*/
    
//     free_dvector(u,1,np);
//     free_dvector(RHS,1,np);
//     free_dmatrix(A,1,np,1,np);
//     free_dvector(u_old,1,np);

//     mesh_free(&mesh);//free

//     return 0;
// }
