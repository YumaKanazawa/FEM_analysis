//特性曲線用のプログラム
#include "../main_mesh.h"

#define D 0.01//(1.25*pow(10,-4))//Diffusion Coefficient
#define N 32
#define Nt (5*pow(N,2)/2)//時間方向の分割数
#define T_max M_PI
#define delta_t 0.07//(T_max/Nt) //time step;

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

//phi_i*phi_jの積
double phi_ij(mesh_t *mesh,int Kl,int i,int j,double x,double y){
    // int n=mesh->n;

    double phi_ij=phi(mesh,i,x,y,Kl)*phi(mesh,j,x,y,Kl);

    return phi_ij;
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
    return 0.0;
}

//ノイマン境界条件
double g1(double x,double y){
    return 0.0;
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

//上流点の座標x_{up}
double *up_point(mesh_t *mesh,double x,double y,double *u){
    double *ret=dvector(1,mesh->dim);
    ret[1]=x-delta_t*u[1],ret[2]=y-delta_t*u[2];
    return ret;
}

//上流点を含むベクトルの数値積分
double *upstream_Integration(mesh_t *mesh,double *u_old){//iの方を上流点の要素に変更
    printf("上流点の反映\n");
    int dim=mesh->dim;
    int **elnp=mesh->elnp;
    int ne=mesh->ne;
    int np=mesh->np;

    double *return_vector=dvector(1,np);//返す用のベクトル

    for(int Kl=1;Kl<=ne;Kl++){
        double S_Kl=area(mesh,Kl);//Kl番目の面積

        for(int j=1;j<=dim+1;j++){//Klを構成する３点のベクトルを返す
            double ret=0.0;//各行の線形和を入れる

            int ver2=elnp[Kl][j];//要素Klのj番目の節点番号
            
            int N_5=7;//積分点の個数
            for(int num=1;num<=N_5;num++){
                // printf("積分点%d\n",num);
                double w;//積分点の重み
                if(num==1){
                    w=9.0/40.0;
                }else if(2<=num && num<=4){
                    w=(155.0-sqrt(15))/1200.0;
                }else{
                    w=(155.0+sqrt(15))/1200.0;
                }

                double *x=Pi(mesh,Kl,num);//Klのnum番目の積分点の座標

                double *U=u(x[1],x[2]);//積分点における流速
                double *x_up=up_point(mesh,x[1],x[2],U);//上流点の座標
                int up_ele=search_past_point(mesh,x[1],x[2],U,delta_t);//上流点が含まれる要素番号

                double phi_i=0.0;
                for(int i=1;i<=dim+1;i++){//{u(Pi)φ(x)+u(Pj)φ(x)+u(Pk)φ(x)}v(x)
                    int ver1=elnp[up_ele][i];//要素Klのi番目の節点番号
                    phi_i+=u_old[ver1]*phi(mesh,ver1,x_up[1],x_up[2],up_ele);//Klの積分点における関数の値
                }

                double phi_j=phi(mesh,ver2,x[1],x[2],Kl);//Klの積分点における関数の値

                ret+=w*phi_i*phi_j;

                free_dvector(x,1,dim);
                free_dvector(U,1,dim);
                free_dvector(x_up,1,dim);
            }
            return_vector[ver2]+=ret*S_Kl;//各要素ごとに積分の値を加えていく
        }
    }
    return return_vector;
}

//外力項の離散化ベクトル
double *out_force_charcterize(mesh_t *mesh,double *u_old){
    printf("out_force(外力項の離散化)\n");
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
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
        // //要素KlがPi,Pj,Pkで構成されている
        // double **M=dmatrix(1,n,1,n);//要素質量行列
        // double *f_vector=dvector(1,n);

        // /*========外力項の計算=========*/
        // /*===================剛性行列の作成==========================*/
        // for(int i=1;i<=n;i++){  
        //     for(int j=1;j<=n;j++){
        //         // int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目の節点番号

        //         /*==========∫φiφjdxの計算結果代入==============*/
        //         double I=Int(mesh,i,j,l);//i->Pi,j->Pjとして考えている．
        //         /*==========================================*/

        //         M[i][j]+=I;
        //     }
        // }

        // /*=========================================================*/
        // for(int i=1;i<=n;i++){
        //     int ver1=elnp[l][i];//l番目の要素のi番目の接点
        //     double x=npxy[ver1][1],y=npxy[ver1][2];
        //     f_vector[i]=delta_t*f(x,y);
        // }

        // double *rhs=matrix_vector_product_CRS(M,f_vector,n);//外力項の計算結果


        // for(int i=1;i<=n;i++){
        //     int ver1=elnp[l][i];
        //     return_vector[ver1]+=rhs[i];
        // }
        // /*=============================*/

        // free_dvector(f_vector,1,n);
        // free_dmatrix(M,1,n,1,n);
        // free_dvector(rhs,1,n);

        double S_Kl=area(mesh,l);

        for(int j=1;j<=dim+1;j++){

            double ret=0.0;//各行の線形和を入れる

            int ver=elnp[l][j];//要素Klのj番目の節点番号
            
            int N_5=7;//積分点の個数
            for(int num=1;num<=N_5;num++){
                // printf("積分点%d\n",num);
                double w;//積分点の重み
                if(num==1){
                    w=9.0/40.0;
                }else if(2<=num && num<=4){
                    w=(155.0-sqrt(15))/1200.0;
                }else{
                    w=(155.0+sqrt(15))/1200.0;
                }

                double *x=Pi(mesh,l,num);//Klのnum番目の積分点の座標

                double *U=u(x[1],x[2]);//積分点における流速
                double *x_up=up_point(mesh,x[1],x[2],U);//上流点の座標
                int up_ele=search_past_point(mesh,x[1],x[2],U,delta_t);//上流点が含まれる要素番号

                double phi_i=0.0;
                for(int i=1;i<=dim+1;i++){//{u(Pi)φ(x)+u(Pj)φ(x)+u(Pk)φ(x)}v(x)
                    int ver1=elnp[up_ele][i];//要素Klのi番目の節点番号
                    phi_i+=u_old[ver1]*phi(mesh,ver1,x_up[1],x_up[2],up_ele);//Klの積分点における関数の値
                }

                double phi_j=phi(mesh,ver,x[1],x[2],l);//Klの積分点における関数の値

                ret+=w*phi_i*phi_j;

                free_dvector(x,1,dim);
                free_dvector(U,1,dim);
                free_dvector(x_up,1,dim);
            }

            //ret*S_Klが上流点の積分結果を格納したもの．各要素ごとに積分の値を加えていく

            //外力項の離散化について
            double mij=0.0;
            for(int i=1;i<=dim+1;i++){
                double I=Int(mesh,i,j,l);
                double x=npxy[ver][1],y=npxy[ver][2];
                mij+=delta_t*f(x,y)*I;
            }

            return_vector[ver]+=(mij+ret*S_Kl);//上流点の計算結果を反映させる
        }

    }

    // double *RHS=upstream_Integration(mesh,u_old);
    // for(int i=1;i<=np;i++)return_vector[i]+=RHS[i];
    // free_dvector(RHS,1,np);
    
    /*=====================ノイマン条件の反映=====================*/
    for(int b=1;b<=nb;b++){
        int Gamma=2;//ノイマン条件の判定
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

            double *rhs_Neumman=matrix_vector_product_CRS(B,Neumman_vector,Neumman_dim);

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

    //ディリクレ条件の反映
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

    //初期条件の代入
    double *u_old=dvector(1,np);
    for(int i=1;i<=np;i++){
        double x=npxy[i][1],y=npxy[i][2];
        u_old[i]=init(x,y);
    }

    double **A=Al(weak_form,&mesh);//剛性行列(境界条件込み)

    for(int T=0;T<=Nt;T+=1){//時刻Tにおいて解を求める
        double *RHS=out_force_charcterize(&mesh,u_old);
        printf("t=%d,|u|=%f\n",T,vector_norm1(u_old,1,np,1.0));

        char str[200];
        sprintf(str,"figure/mesh%d.dat",T);
        // make_mesh_data_for_gnuplot(mesh,u_old,str);//gnuplot用のファイル作成
        make_result_data_for_GLSC(&mesh,u_old,str);//GLSC用のデータ出力

        double *u=CG_CRS(A,RHS,np);//解の計算

        for(int i=1;i<=np;i++)u_old[i]=u[i];//解の更新

        free_dvector(u,1,np);
        free_dvector(RHS,1,np);
    }

    free_dmatrix(A,1,np,1,np);
    free_dvector(u_old,1,np);


    mesh_free(&mesh);//free;
    return 0;
}