#include "Include/functions.h"

#define D 0.02 //Diffusion Coefficient
#define delta_t 0.01//time step

typedef struct mesh{
    int dim, np, ne, nb;//データの次元，節点数，要素数，境界数
    double **npxy;//節点座標対応表
    int **elnp, **bound;//要素節点対応表
    int n;//n角形分割
    double Volume;//全領域での積分値
}mesh_t;

void print_mesh_data(mesh_t mesh);//メッシュデータの出力
void make_mesh_data_for_gnuplot(mesh_t mesh,double *u,char *str);//gnuplot用のデータの作成
void alloc_scan_mesh(mesh_t *mesh,char *s1,char *s2);//メッシュデータの読み込み，データ保存
void mesh_free(mesh_t *mesh);//free関数
double Volume(mesh_t *mesh);//全領域での積分値の計算
void drawney(mesh_t *mesh);//ドロネーのアルゴリズム
double Integer(int *alpfa);//数値積分その１
double *coord(mesh_t *mesh,int K,int x);//Kを構成するxベクトル
double **NumInt_deg_two(mesh_t *mesh,int K);//積分点の計算
double **NumInt_deg_five(mesh_t *mesh,int K);//積分点の計算
double *coef_plate_grad(mesh_t *mesh ,int K,int Pi);//要素KにおけるPiで１となる平面の勾配
double L(mesh_t *mesh,int i,int j);//２点間の距離
double **permu(double *a,double b,int N);//順列
void make_result_data_for_GLSC(mesh_t *mesh,double *u,char *str);//GLSC用のデータ
void make_coef_matrix(mesh_t *mesh,double **A,double *b,int t);//係数行列の作成
double area(mesh_t *mesh,int Kl);//Klの面積
double f(double x,double y);//Nonlinear function
double g(double x,double y);//Diriclet boundary condition's function
double g1(double x,double y);//Neumman boundary condition's function
double init(double x,double y);//initial condition
double *u(double x,double y);//Already given vector


typedef double (* pfunc)(mesh_t *mesh,int Kl,int i,int j,double x,double y);//被積分関数の定義
double Int_div_five(pfunc func,mesh_t *mesh,int Kl,int i,int j);//pfunc型の関数の数値積分

typedef double (*weak)(mesh_t *mesh,int Kl,int i,int j);//弱形式の入る関数ポインタ
double **Al(weak weak_form,mesh_t *mesh);

typedef double (*out)(double *u_p,int i,double x,double y);//右辺の離散化における関数
double *out_force(out RHS,mesh_t *mesh,double *u_old);


void alloc_scan_mesh(mesh_t *mesh,char *s1,char *s2){
    FILE *fp;
    int dim,np,ne,nb,i,j,k,dummy;//,n;(利用されていない)
    printf("alloc_scan_mesh: start .... "); fflush(stdout);

    dim = atoi(s1); mesh->dim = dim;//次元数を数字に変更してdimに格納
    mesh->n=3;//n角形で分割するとき

    
    if((fp=fopen(s2,"r"))==NULL){printf("Can’t open file: %s.\n",s2); exit(1);}//ファイルの読み込み
    

    fscanf(fp,"%d %d %d",&np,&ne,&nb); 
    mesh->np = np;
    mesh->ne = ne;
    mesh->nb = nb;//各種パラメータを変数に格納

    
    /*メモリの動的な確保*/
    mesh->npxy =dmatrix(1,np,1,dim); printf("npxy,"); fflush(stdout);
    mesh->elnp =imatrix(1,ne,1,dim+1); printf("elnp,"); fflush(stdout);
    mesh->bound=imatrix(1,nb,1,dim); printf("bound,"); fflush(stdout);


    /*行列の値格納*/
    for(i=1;i<=np;i++){
        for(j=1;j<=dim;j++){
            fscanf(fp,"%lf",&(mesh->npxy[i][j]));//座標の値
        }
        fscanf(fp,"%d",&dummy);
    }
    for(k=1;k<=ne;k++){
        for(i=1;i<=dim+1;i++){
            fscanf(fp,"%d",&(mesh->elnp[k][i]));//要素を構成する節点の集合
        }
        fscanf(fp,"%d",&dummy);
    }
    for(k=1;k<=nb;k++){
        for(i=1;i<=dim+1;i++){
            fscanf(fp,"%d",&(mesh->bound[k][i]));
        }
    }
    fclose(fp);
    printf(".... end\n"); fflush(stdout);
    return;
}

void mesh_free(mesh_t *mesh){
    int dim = mesh->dim, np = mesh->np, ne = mesh->ne, nb = mesh->nb;
    printf("mesh_free: start .... "); fflush(stdout);
    free_dmatrix(mesh->npxy,1,np,1,dim); printf("npxy,"); fflush(stdout);
    free_imatrix(mesh->elnp,1,ne,1,dim+1); printf("elnp,"); fflush(stdout);
    free_imatrix(mesh->bound,1,nb,1,dim+1); printf("bound,"); fflush(stdout);
    printf(".... end\n"); fflush(stdout);
    return;
}

void print_mesh_data(mesh_t mesh){
    int np=mesh.np,ne=mesh.ne,nb=mesh.nb,dim=mesh.dim,i,k;
    double **npxy;
    int **elnp, **bound;
    // FILE *fp;
    if(dim != 2){ printf("dim should be 2!\n"); exit(1);}
    npxy = mesh.npxy;
    elnp = mesh.elnp;
    bound = mesh.bound;
    // output the mesh data
    printf("np=%d, ne=%d, nb=%d\n",np,ne,nb);
    for(i=1;i<=np;i++){
        printf("%d: (%e, %e)\n",i,npxy[i][1],npxy[i][2]);
    }
    for(k=1;k<=ne;k++){
        printf("%d: (%d, %d, %d)\n",k,elnp[k][1],elnp[k][2],elnp[k][3]);
    }
    for(k=1;k<=nb;k++){
        printf("%d: (%d, %d); %d\n",k,bound[k][1],bound[k][2],bound[k][3]);
    }
    return;
}

//gnuplotのためのファイル作成
void make_mesh_data_for_gnuplot(mesh_t mesh,double *u,char *str){
    printf("make_mesh_data_for_gnuplot");
    int np=mesh.np,dim=mesh.dim,k;
    // int nb=mesh.nb;
    //int ne=mesh.ne;
    double **npxy;
    // int **bound;
    // int **elnp; 
    FILE *fp;
    if(dim != 2){ printf("dim should be 2!\n"); exit(1);}//二次元に対応してるのか？
    npxy = mesh.npxy;
    // elnp = mesh.elnp;
    // bound = mesh.bound;
    // make a file: mesh.dat

    if((fp=fopen(str,"w"))==NULL){//書き込み方式でファイルを開く
        printf("Can’t open file: %s.\n","mesh.dat");
        exit(1);
    }

    //近似解の出力
    for(k=1;k<=np;k++){
        double x=npxy[k][1];
        double y=npxy[k][2];
        fprintf(fp,"%f %f %f\n",x,y,u[k]);
    }
    
    fclose(fp);
    return;
}

//三角形の面積計算
double Volume(mesh_t *mesh){
    printf("Volume");
    int n=mesh->n;//要素の分割形状
    int dim=mesh->dim;//問題を考える次元

    // int k=0;
    int **elnp=mesh->elnp;
    double **npxy=mesh->npxy;

    double **A=dmatrix(1,dim+1,1,n);//行列のメモリ確保．要素番号の格納

    double area=0.0;
    for(int i=0;i<n;i++){
        A[0][i]=1.0;
        for(int j=1;j<dim+1;j++){
            int ver_n=elnp[i][j-1];//j番目の接点番号
            A[i][j]=npxy[ver_n][j-1];
            printf("%f\n",A[i][j]);
        }
        printf("\n");
        area+=1.0;
    }
    free_dmatrix(A,1,dim+1,1,n);
    return area;
}

//接点が領域に入ってるかどうかの確認
int count(mesh_t *mesh,double x,double y){
    // int np=mesh->np;
    int ne=mesh->ne;
    // int nb=mesh->nb;
    // int dim=mesh->dim;
    int n=mesh->n;
    // double **npxy;
    // int **elnp; 
    // int **bound;
    // npxy = mesh->npxy;
    // elnp = mesh->elnp;
    // bound = mesh->bound;

    int Kl;

    for(int K=1;K<=ne;K++){
        double *X1=dvector(0,n-1);
        double *X2=dvector(0,n-1);
        for(int i=0;i<n;i++){
            X1[i]=coord(mesh,K,1)[i]-x;
            X2[i]=coord(mesh,K,2)[i]-y;
        }

        double Arg=0.0;
        for(int i=0;i<n;i++){
            double ip=X1[i]*X1[(i+1)%n]+X2[i]*X2[(i+1)%n];
            double x1_l=X1[i]*X1[i]+X2[i]*X2[i];
            double x2_l=X1[(i+1)%n]*X1[(i+1)%n]+X2[(i+1)%n]*X2[(i+1)%n];
            if(x1_l!=0.0 && x2_l!=0.0){
                Arg+=acos(ip/sqrt(x1_l*x2_l));
                // printf("i=%d,%f\n",i,sqrt(x2_l));
            }else{
                Arg=0.0;
            }
        }
        // printf("Arg=%f\n",Arg);
    
        free_dvector(X1,0,n-1);
        free_dvector(X2,0,n-1);

        if(Arg==0.0 || Arg==2*pi){
            Kl=K;
            // printf("Arg=%f,K=%d\n",Arg,Kl);
        }
    }
    return Kl;
}


// void drawney(mesh_t *mesh){
//     int np=mesh->np,ne=mesh->ne,nb=mesh->nb,dim=mesh->dim,n=mesh->n;
//     double **npxy;
//     int **elnp, **bound;
//     npxy = mesh->npxy;
//     elnp = mesh->elnp;
//     bound = mesh->bound;

//     /*============1st========super Triangle============*/
//     //最大最小値の探索
//     double max=npxy[1][1];
//     double min=npxy[1][1];
//     for(int i=1;i<=np;i++){
//         for(int j=0;j<dim;j++){
//             if(max<=npxy[i][j]){
//                 max=npxy[i][j];
//             }
//             if(min>=npxy[i][j]){
//                 min=npxy[i][j];
//             }
//         }   
//     }
//     //上記４点によってなされる正方形に外接する円の情報
//     double L=0.5*sqrt(2)*fabs(max-min);
//     //三角形の内心=正方形の中心(x0,x0)
//     double x0=min+0.5*fabs(max-min);
//     double y0=x0;
//     //SuperTriangleを作る座標
//     double tx=x0,ty=y0+2.0*L;
//     double tx1=x0-sqrt(3)*L;
//     double ty1=y0-L;
//     double tx2=x0+sqrt(3)*L;
//     double ty2=y0-L;

//     Super Triangleの座標行列
//     double **Coord=dmatrix(1,n,1,dim);
//     Coord[1][1]=tx;Coord[1][2]=ty;
//     Coord[2][1]=tx1;Coord[2][2]=ty1;
//     Coord[3][1]=tx2;Coord[3][2]=ty2;
    
//     /*=================================================*/

//     /*==============2nd==============線を引く============*/

//     int data_number=1;
//     double *Xi=dvector(1,dim);
//     for(int i=1;i<=dim;i++)Xi[i]=npxy[data_number][i];
    
//     int C=count(Coord,Xi);

//     free_dmatrix(Coord,1,n,1,dim);
//     free_dvector(Xi,dim);

//     printf("%d\n",C);
    
//     /*==================================================*/
// }

double *coord(mesh_t *mesh,int K,int x){
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    int n=mesh->n;
    // int np=mesh->np;
    // int ne=mesh->ne;
    // int nb=mesh->nb;
    double **npxy;
    int **elnp;//**bound;
    npxy = mesh->npxy;
    elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/
    double *ret=dvector(1,n);

    for(int i=1;i<=dim+1;i++){
        int Pi=elnp[K][i];
        ret[i]=npxy[Pi][x];//Piを構成するx座標
    }
    return ret;
}

//順列の表現
double **permu(double *a,double b,int N){
    int sta=1;
    int end=sta+N;
    double **Ret=dmatrix(sta,end,sta,end);

    for(int j=sta;j<=end;j++){
        //Nはベクトルaの長さ
        //順列の計算
        double *ret=dvector(sta,end);
        for(int i=sta;i<=end;i++){
            if(i<j){
                ret[i]=a[i];
            }else if(i==j){
                ret[i]=b;
            }else{
                ret[i]=a[i-1];
            }
        }
        //順列の計算結果をコピー
        for(int k=sta;k<=end;k++){
            Ret[j][k]=ret[k];
        }

        free_dvector(ret,sta,end);
    }

    return Ret;
}

//φiφjの積分
double Int(mesh_t *mesh,int i,int j,int Kl){
    double I=0.0;
    double S_Kl=area(mesh,Kl);//Kl番目の面積
    if(i==j){
        I=S_Kl/6.0;
    }else{
        I=S_Kl/12.0;
    }
    return I;
}

//２次元の積分点
double **NumInt_deg_two(mesh_t *mesh,int K){
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    // int n=mesh->n;
    // int np=mesh->np;
    // int ne=mesh->ne;
    // int nb=mesh->nb;
    // double **npxy;
    // int **elnp,**bound;
    // npxy = mesh->npxy;
    // elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

    int N_2=dim+1;
    double p=(dim+2-sqrt(dim+2))/((dim+1)*(dim+2)),q=(dim+2+dim*sqrt(dim+2))/((dim+1)*(dim+2));

    double **Coord=dmatrix(1,N_2,1,dim);//戻り値用の行列
    /*=====================積分点の定義==========================*/
    //順列を二次元配列(行列)として表現する

    int Int_num=dim;//次元？順列の長さ

    double *p_perm=dvector(1,Int_num);
    for(int i=1;i<=Int_num;i++)p_perm[i]=p;

    double **coef=permu(p_perm,q,Int_num);//順列の入った行列

    int x=1,y=2;
    double *coor_x=coord(mesh,K,x);//x方向の座標格納,Piを構成するx座標
    double *coor_y=coord(mesh,K,y);//y方向の座標格納,Piを構成するy座標

    double *integer_point_x=matrix_vector_product(coef,coor_x,dim+1);//重心座標との線形和によって座標を表現する
    double *integer_point_y=matrix_vector_product(coef,coor_y,dim+1);//重心座標との線形和によって座標を表現する

    free_dmatrix(coef,1,dim+1,1,dim+1);
    free_dvector(coor_x,1,dim+1);free_dvector(coor_y,1,dim+1);
    free_dvector(p_perm,1,Int_num);
    /*==========================================================*/
    //返す用の変数格納
    for(int i=1;i<=N_2;i++){
        Coord[i][x]=integer_point_x[i];
        Coord[i][y]=integer_point_y[i];
    }
    return Coord;//積分点aのreturn
}

//２次元の積分点
double **NumInt_deg_five(mesh_t *mesh,int K){
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    // int n=mesh->n;
    // int np=mesh->np;
    // int ne=mesh->ne;
    // int nb=mesh->nb;
    // double **npxy;
    // int **elnp,**bound;
    // npxy = mesh->npxy;
    // elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

    int N_5=7;
    double p=(6-sqrt(15))/21,q=(9+2*sqrt(15))/21,r=(6+sqrt(15))/21,s=(9-2*sqrt(15))/21;
    double t=1./3.;

    double **Coord=dmatrix(1,N_5,1,dim);//戻り値用の行列
    /*=====================積分点の定義==========================*/
    //順列を二次元配列(行列)として表現する

    int Int_num=dim;

    int x=1,y=2;
    double *coor_x=coord(mesh,K,x);//x方向の座標格納,Piを構成するx座標
    double *coor_y=coord(mesh,K,y);//y方向の座標格納,Piを構成するy座標


    /*=====================i=1(t,t;t)====================*/
    //※全て同じ要素なので，順列は(t,t,t)のみ
    double *t_perm=dvector(1,Int_num+1);
    for(int i=1;i<=Int_num+1;i++)t_perm[i]=t;

    double a_1_x=inner_product(1,dim+1,t_perm,coor_x);
    double a_1_y=inner_product(1,dim+1,t_perm,coor_y);

    free_dvector(t_perm,1,Int_num);
    /*==================================================*/

    /*=====================i=2,3,4(p,p;q)===============*/
    double *p_perm=dvector(1,Int_num);
    for(int i=1;i<=Int_num;i++)p_perm[i]=p;
    double **coef=permu(p_perm,q,Int_num);


    double *bynaric_coord_p_x=matrix_vector_product(coef,coor_x,Int_num+1);
    double *bynaric_coord_p_y=matrix_vector_product(coef,coor_y,Int_num+1);
    free_dvector(p_perm,1,Int_num);
    free_dmatrix(coef,1,Int_num,1,Int_num);
    /*==================================================*/

    /*=====================i=5,6,7(r,r;s)===============*/
    double *r_perm=dvector(1,Int_num);
    for(int i=1;i<=Int_num;i++)r_perm[i]=r;
    coef=permu(r_perm,s,Int_num);

    double *bynaric_coord_r_x=matrix_vector_product(coef,coor_x,Int_num+1);
    double *bynaric_coord_r_y=matrix_vector_product(coef,coor_y,Int_num+1);
    free_dvector(r_perm,1,Int_num);
    free_dmatrix(coef,1,Int_num,1,Int_num);
    /*==================================================*/

    free_dvector(coor_x,1,dim+1);free_dvector(coor_y,1,dim+1);

    for(int i=1;i<=N_5;i++){
        double rx,ry;
        int number_p=1,number_r=1;
        if(i==1){
            rx=a_1_x;
            ry=a_1_y;

        }else if(2<=i && i<=4){
            rx=bynaric_coord_p_x[number_p];
            ry=bynaric_coord_p_y[number_p];
            number_p+=1;
        }else{
            rx=bynaric_coord_r_x[number_r];
            ry=bynaric_coord_r_y[number_r];
            number_r+=1;
        }
        Coord[i][x]=rx;
        Coord[i][y]=ry;
    }
    return Coord;//積分点aのreturn
}

//３点の面積
double S(double **Coord){
    //基準となる座標
    int x=1,y=2;
    double det=(Coord[2][x]*Coord[3][y]+Coord[3][x]*Coord[1][y]+Coord[1][x]*Coord[2][y])-(Coord[2][x]*Coord[1][y]+Coord[2][y]*Coord[3][x]+Coord[1][x]*Coord[3][y]);

    // double det=(Coord[2][x]-Coord[1][x])*(Coord[3][y]-Coord[1][y])+(Coord[3][x]-Coord[1][x])*(Coord[2][y]-Coord[1][y]);
    return 0.5*det;//行列式が正になるようにする
} 

//ver_number番目の基底関数のele_number番目の要素上での平面関数の勾配
double *coef_plate_grad(mesh_t *mesh ,int ele_number,int ver_number){
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    int n=mesh->n;
    // int np=mesh->np;
    // int ne=mesh->ne;
    // int nb=mesh->nb;
    double **npxy;
    int **elnp;
    // int **bound;
    npxy = mesh->npxy;
    elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

    /*============ele_numberを構成する座標行列の作成===============*/
    double **Pxy=dmatrix(1,n,1,dim);//要素を構成する節点座標
    for(int i=1;i<=n;i++){
        int ver=elnp[ele_number][i];//ele_numberを構成するi番目の節点番号
        for(int j=1;j<=dim;j++){
            Pxy[i][j]=npxy[ver][j];//要素を構成する節点座標
            // printf("%f,",Pxy[i][j]);
        }
        // printf("\n");
    }
    /*=========================================================*/
    
    /*===================一次平面の係数決定=====================*/
    //３点から面積計算
    double Sl=2*S(Pxy);//行列式
    // printf("detA=%f\n",2.0*Sl);
    double *C=dvector(0,dim);//return用の配列
    for(int i=1;i<=n;i++){
        int x=1,y=2;
        int ver=elnp[ele_number][i];//ele_number番目の要素を構成するi番目の節点番号
        if(ver==ver_number){
            int n1=i%3+1,n2=(i+1)%3+1;
            //クラメルの公式による連立方程式の解
            double xi=Pxy[n1][x],xi1=Pxy[n2][x];
            double yi=Pxy[n1][y],yi1=Pxy[n2][y];
            C[0]=(yi1*xi-yi*xi1)/Sl;
            C[1]=(yi-yi1)/Sl;
            C[2]=(xi1-xi)/Sl;

            // printf("i=%d,%f,%f,",i,Pxy[i][1],Pxy[i][2]);
            // printf("i=%d,(%d,%d)",i,n1,n2);
        }
        // printf("%f\n",C[0]+Pxy[i][1]*C[1]+Pxy[i][2]*C[2]);
    }
    // printf("\n");
    free_dmatrix(Pxy,1,n,1,dim);
    /*=======================================================*/
    return C;//係数ベクトルを返す
}

double L(mesh_t *mesh,int i,int j){
    /*==================構造体のデータ読み込み==========================*/
    // int dim=mesh->dim;
    // int n=mesh->n;
    // int np=mesh->np;
    // int ne=mesh->ne;
    // int nb=mesh->nb;
    double **npxy;
    // int **elnp;
    // int **bound;
    npxy = mesh->npxy;
    // elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

    double xi=npxy[i][1],yi=npxy[i][2];
    double xj=npxy[j][1],yj=npxy[j][2];

    double l=pow(xi-xj,2)+pow(yi-yj,2);

    return pow(l,0.5);
}

double area(mesh_t *mesh,int l){
    //l番目の要素の面積
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    int n=mesh->n;
    // int np=mesh->np;//,ne=mesh->ne,nb=mesh->nb;
    double **npxy;
    int **elnp;
    // int **bound;
    npxy = mesh->npxy;
    elnp = mesh->elnp;
    // bound = mesh->bound;
    /*==============================================================*/

    /*============Klを構成する座標行列の作成===============*/
    double **Pxy=dmatrix(1,n,1,dim);//要素を構成する節点座標
    for(int i=1;i<=n;i++){
        int ver=elnp[l][i];//lを構成するi番目の節点番号
        for(int j=1;j<=dim;j++){
            Pxy[i][j]=npxy[ver][j];//要素を構成する節点座標
        }
    }
    double S_Kl=S(Pxy);//要素Klの面積
    free_dmatrix(Pxy,1,n,1,dim);
    /*==================================================*/
    return S_Kl;
}



void Diriclet(mesh_t *mesh,double **A,double *b){
    printf("Dirichlet(ディリクレ境界条件の反映)\n");
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    // int n=mesh->n;
    int np=mesh->np;
    // int ne=mesh->ne;
    int nb=mesh->nb;
    double **npxy;
    // int **elnp;
    int **bound;
    npxy = mesh->npxy;
    // elnp = mesh->elnp;
    bound = mesh->bound;
    /*==============================================================*/

    int Gamma=2;//ノイマン境界の番号指定
    int Dir=1;


    /*=========================step1 fi=fi-ΣAijg(Pj)==================================*/
    for(int i=1;i<=nb;i++){

        int ver_Neum=-1;
        if(bound[i][dim+1]==Gamma){//ノイマン条件が与えられる番号の探索
            ver_Neum=bound[i][1];
        }

        if(ver_Neum!=-1){//節点番号がないときの処理
            double sum=0.0;//ΣAijg(Pj)の計算用の変数
            for(int j=1;j<=nb;j++){
                if(bound[j][dim+1]==Dir){//ディリクレ境界条件が与えられる節点
                    int b1=bound[j][1];
                    double x=npxy[b1][1],y=npxy[b1][2];

                    sum+=A[ver_Neum][b1]*g(x,y);
                }
            }
            b[ver_Neum]-=sum;//fi=fi-ΣAijg(Pj)
        }
    }
    /*================================================================================*/
    
    /*=================ディリクレ境界条件の挿入=================*/
    for(int i=1;i<=nb;i++){
        if(bound[i][dim+1]==Dir){//Γ1の判定
            int b1=bound[i][1];//Γ1上の節点番号//,b2;//Γ1上の節点番号用の変数
            //b1,b2の行と列をそれぞれ0にする
            b[b1]=A[b1][b1]*g(npxy[b1][1],npxy[b1][2]);

            for(int low_col=1;low_col<=np;low_col++){
                if(b1!=low_col){
                    A[b1][low_col]=0.0;
                    A[low_col][b1]=0.0;
                }
            }
           
        }
    }
    /*======================================================*/

    
    // for(int i=1;i<=nb;i++){
    //     if(bound[i][dim+1]==Dir){
    //         int b1=bound[i][1];
    //         for(int low_col=1;low_col<=np;low_col++){
    //             if(b1!=low_col){
    //                 A[b1][low_col]=0.0;
    //                 // A[low_col][b1]=0.0;
    //             }
    //         }
    //         A[b1][b1]=1.0;
    //         b[b1]=g(npxy[b1][1],npxy[b1][2]);
    //     }
    // }
}


//φi(u・∇φj)の積分(pfuncと同様の型の関数は積分可能)
double Int_div_five(pfunc func,mesh_t *mesh,int Kl,int i,int j){

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

        // double *U=u(x_num_I,y_num_I);
        // // double div=div_u(mesh,x_num_I,y_num_I);

        // // double div_int=div*phi(mesh,i,x_num_I,y_num_I,Kl)*phi(mesh,j,x_num_I,y_num_I,Kl);//(∇・u)φiφj
        // // printf("∇・u=%f,phi_i=%f,phi_j%f\n",div,phi(mesh,i,x_num_I,y_num_I,Kl),phi(mesh,j,x_num_I,y_num_I,Kl));

        // double phi_i_int=inner_product(1,dim,U,Coord_i);//u・∇φi
        // // printf("u・∇φj=%f,",phi_i_int);
        // double phi_j_int=inner_product(1,dim,U,Coord_j);//u・∇φj
        // // printf("u・∇φj=%f,",phi_j_int);
        // double phi_i=phi(mesh,i,x_num_I,y_num_I,Kl);//φi
        // double phi_j=phi(mesh,j,x_num_I,y_num_I,Kl);//φj

        double cover_int=func(mesh,Kl,i,j,x_num_I,y_num_I);
        // double cover_int=phi_i*phi_j_int-phi_j*phi_i_int;
        ret+=w*cover_int;

        // free_dvector(U,1,dim);
    }

    // free_dmatrix(Coord,1,N,1,mesh->np);
    // free_dvector(Coord_i,0,dim);
    // free_dvector(Coord_j,0,dim);

    return S*ret;
}

double phi(mesh_t *mesh,int i,double x,double y,int Kl){
    // int n=mesh->n;
    int dim=mesh->dim;
    //φiのKl番目上の平面
    double *coef=coef_plate_grad(mesh,Kl,i);
    double ret=coef[0]+coef[1]*x+coef[2]*y;
    free_dvector(coef,0,dim);

    if(ret<=0.0 && ret>1.0){
        ret=0.0;
    }

    return ret;
}

void matrix_vector_print(mesh_t mesh,double **A,double *RHS){
    for(int i=1;i<=mesh.np;i++){
        for(int j=1;j<=mesh.np;j++){
            if(A[i][j]>0){
                printf("+%0.2f,",A[i][j]);
            }else if(A[i][j]==0){
                printf("% 0.2f,",A[i][j]);
            }else{
                printf("%0.2f,",A[i][j]);
            }
        }
        printf(" %0.2f\n",RHS[i]);
    }
}

void make_result_data_for_GLSC(mesh_t *mesh,double *u,char *str){
    // int np=mesh.np;
    int dim=mesh->dim;
    int n=mesh->n;
    // int nb=mesh.nb;
    int ne=mesh->ne;
    double **npxy;
    // int **bound;
    int **elnp; 

    FILE *fp;
    if(dim != 2){ printf("dim should be 2!\n"); exit(1);}//二次元に対応してるのか？
    npxy = mesh->npxy;
    elnp = mesh->elnp;
    // bound = mesh.bound;
    // make a file: mesh.dat

    if((fp=fopen(str,"w"))==NULL){//書き込み方式でファイルを開く
        printf("Can’t open file: %s.\n","mesh.dat");
        exit(1);
    }

    // int count=0;
    //データの出力
    for(int l=1;l<=ne;l++){
        for(int i=1;i<=n;i++){
            // count+=1;
            int ele_k=elnp[l][i];
            //l番目の要素のi番目の接点番号
            double x=npxy[ele_k][1];
            double y=npxy[ele_k][2];
            fprintf(fp,"%f %f %f %d\n",x,y,u[ele_k],ele_k);
        }
    }
    // printf("total=%d\n",count);
    fclose(fp);
    return;
}

void search_ele_count(mesh_t *mesh,double *u){
    // int np=mesh.np;
    // int dim=mesh->dim;
    // int n=mesh->n;
    // int nb=mesh.nb;
    int ne=mesh->ne;
    // double **npxy;
    // int **bound;
    // int **elnp; 

    //ノイマン条件の場合，1番目の節点がその境界上にあれば発散する.その対処
    for(int ver_num=1;ver_num<=ne;ver_num++){
        if(u[ver_num]>=1.0){
            u[ver_num]=0.0;
            // double u_1=0.0;//u[1]の初期化
            // int count_ele=0;
            // for(int i=1;i<=mesh.ne;i++){
            //     for(int v=1;v<=mesh.n;v++){
            //         int l1=mesh.elnp[i][v];
            //         if(l1==1){//1番目の節点番号を含む要素がiに入る
            //             count_ele++;//1を含む要素の個数
            //             int l2=i%mesh.n+1,l3=(i+1)%mesh.n+1;//l1以外の要素確認
            //             int e_l2=mesh.elnp[i][l2],e_l3=mesh.elnp[i][l3];
            //             u[l1]=0.0,u[e_l2]=0.0,u[e_l3]=0.0;
            //             // u_1+=(u[e_l2]+u[e_l3]);//e_l2,e_l3番目の関数の値との平均を取る
            //         }
            //     }
            // }
            // u[1]=0.0;//u_1/(mesh.n*count_ele);//節点数=n×要素数
        }
    }
}

void make_coef_matrix(mesh_t *mesh,double **A,double *b,int t){
    int np=mesh->np;
    // int dim=mesh->dim;
    // int n=mesh->n;
    // int nb=mesh->nb;
    // int ne=mesh->ne;
    // double **npxy=mesh->npxy;
    // int **bound=mesh->bound;
    // int **elnp=mesh->elnp; 

    char str[200]="file/heat_matrix.csv";
    printf("matrix file name?");
    // scanf("%s",str);
    printf("Enterd\n");

    FILE *write_matrix;
    write_matrix=fopen(str,"w");

    if(write_matrix==NULL){
        printf("Can't open Matrix file!");
        exit(1);
    }

    for(int i=1;i<=np;i++){
        for(int j=1;j<=np;j++){
            fprintf(write_matrix,"%lf",A[i][j]);
            if(j!=np){
                fprintf(write_matrix,",");
            }
        }
        fprintf(write_matrix,"\n");
    }
    fclose(write_matrix);

    char str_v[200];
    printf("vector file name?");
    sprintf(str_v,"file/init_matrix%d.csv",t);
    printf("Enterd\n");

    FILE *write_vector;
    write_vector=fopen(str_v,"w");
    if(write_vector==NULL){
        printf("Can't open Vector file!");
        exit(1);
    }

    for(int j=1;j<=np;j++){
        fprintf(write_matrix,"%lf\n",b[j]);
    }
    fclose(write_vector);
}

//要素剛性行列の作成
double **Al(weak weak_form,mesh_t *mesh){
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
        // double S_Kl=area(mesh,l);
        /*===================剛性行列の成分計算==========================*/
        for(int i=1;i<=n;i++){
            for(int j=1;j<=n;j++){
                int ver1=elnp[l][i],ver2=elnp[l][j];//lを構成するi番目の節点番号

                // double *C_i=coef_plate_grad(mesh,l,ver1);//φiの勾配
                // double *C_j=coef_plate_grad(mesh,l,ver2);//φjの勾配

                // /*======ここを問題の弱形式に応じて変更する========*/
                // double Int_ij=Int(mesh,i,j,l);
                // double Int_div_ij=Int_div_five(func,mesh,l,ver1,ver2);
                // double Int_div_ji=Int_div_five(func,mesh,l,ver2,ver1);

                // //printf("(%d,%d),%.2f,%.2f,",i,j,Int_div_ij,Int_div_ji);

                // double Discreate_WeakForm=Int_ij+0.5*delta_t*(Int_div_ij-Int_div_ji)+D*delta_t*inner_product(1,dim,C_i,C_j)*S_Kl;
                double Discreate_WeakForm=weak_form(mesh,l,i,j);
                // /*==========================================*/
                A[ver1][ver2]+=Discreate_WeakForm;//弱形式の離散結果
                
                // free_dvector(C_i,0,dim);
                // free_dvector(C_j,0,dim);
            }
        }
    }   
    return A;
}


double *out_force(out RHS,mesh_t *mesh,double *u_old){
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
        // double *f_vector=dvector(1,n);
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
            // f_vector[i]=f(x,y);
            u_old_vector[i]=RHS(u_old,ver1,x,y);
        }
        // double *rhs=matrix_vector_product(M,f_vector,n);
        double *rhs_1=matrix_vector_product(M,u_old_vector,n);
        /*=========================================================*/

        /*================返すベクトルに値を代入していく================*/
        for(int i=1;i<=n;i++){
            int ver1=elnp[l][i];
            return_vector[ver1]+=rhs_1[i];
        }
        /*=========================================================*/

        // free_dvector(f_vector,1,n);
        free_dvector(u_old_vector,1,n);
        free_dmatrix(M,1,n,1,n);
        // free_dvector(rhs,1,np);
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

