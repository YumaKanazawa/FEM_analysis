#include "Include/functions.h"

#define SUPG 0.00//SUPG法の定数

typedef struct mesh{
    int dim, np, ne, nb;//データの次元，節点数，要素数，境界数
    double **npxy;//節点座標対応表
    int **elnp, **bound;//要素節点対応表
    int **elel;//要素要素対応表
    int n;//n角形分割
}mesh_t;

void print_mesh_data(mesh_t mesh);//メッシュデータの出力
void make_mesh_data_for_gnuplot(mesh_t mesh,double *u,char *str);//gnuplot用のデータの作成
void alloc_scan_mesh(mesh_t *mesh,char *s1,char *s2);//メッシュデータの読み込み，データ保存
void mesh_free(mesh_t *mesh);//free関数
int comp(int *a,int *b,int m,int n);//一致している個数
int comp_place(int *a,int *b,int m,int n);//一致していない個数
int ele_inside(mesh_t *mesh,int Kl,double x,double y);//要素Kl内に(x,y)が入ってるかどうか
double S(double **Pxy);//1要素の面積
void drawney(mesh_t *mesh);//ドロネーのアルゴリズム
double **NumInt_deg_five(mesh_t *mesh,int K);//積分点の計算
double **Centroid_coord(mesh_t *mesh);//積分点(重心座標系)
double *Pi(mesh_t *mesh,int Kl,int K_n);
double *coef_plate_grad(mesh_t *mesh ,int K,int Pi);//要素KにおけるPiで１となる平面の勾配
double L(mesh_t *mesh,int i,int j);//２点間の距離
void make_result_data_for_GLSC(mesh_t *mesh,double *u,char *str);//GLSC用のデータ
double area(mesh_t *mesh,int Kl);//Klの面積
double err_Lp(mesh_t *mesh,double *u,double p,double t);//解析解との誤差
double phi(mesh_t *mesh,int i,double x,double y,int Kl);//基底関数の値
double f(double x,double y);//Nonlinear function
double g(double x,double y);//Diriclet boundary condition's function
double g1(double x,double y);//Neumman boundary condition's function
double init(double x,double y);//initial condition
double *u(double x,double y);//Already given vector
double phi_ij(mesh_t *mesh,int Kl,int i,int j,double x,double y);//必要な数値積分用の関数
double u_exa(double x,double y,double t);
  
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
    int n=mesh->n;

    
    if((fp=fopen(s2,"r"))==NULL){printf("Can’t open file: %s.\n",s2); exit(1);}//ファイルの読み込み
    

    fscanf(fp,"%d %d %d",&np,&ne,&nb); 
    mesh->np = np;
    mesh->ne = ne;
    mesh->nb = nb;//各種パラメータを変数に格納

    /*メモリの動的な確保*/
    mesh->npxy =dmatrix(1,np,1,dim); printf("npxy,"); fflush(stdout);
    mesh->elnp =imatrix(1,ne,1,dim+1); printf("elnp,"); fflush(stdout);
    mesh->bound=imatrix(1,nb,1,dim); printf("bound,"); fflush(stdout);
    mesh->elel=imatrix(1,ne,1,dim+1);printf("elel,");fflush(stdout);

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


    for(int Kl=1;Kl<=ne;Kl++){//Kl番目の要素を考える
        int s=1,t=n;//比較用の位置
        int *Kl_v=ivector(s,t);//各要素の節点番号
        for(int i=s;i<=t;i++)Kl_v[i]=(mesh->elnp[Kl][(i-s)%n+1]);

        for(int Km=1;Km<=ne;Km++){
            int *Km_v=ivector(s,t);//各要素の節点番号
            for(int i=s;i<=t;i++)Km_v[i]=(mesh->elnp[Km][(i-s)%n+1]);

            int check=comp(Km_v,Kl_v,s,t);
            if(check==n-1){
                int elel_row=comp_place(Kl_v,Km_v,s,t);//二箇所が一致する時に要素同士は隣あう．
                mesh->elel[Kl][elel_row]=Km;//位置はelel_rowではなく，要素の反対側に来るように考えなくてはならない．
            }

            free_ivector(Km_v,s,t);
        }
        free_ivector(Kl_v,s,t);
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
    free_imatrix(mesh->elel,1,ne,1,dim+1);  printf("elel,"); fflush(stdout);
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

int ele_inside(mesh_t *mesh,int Kl,double x,double y){//入ってるかどうか．入っていなければ次に探索する方向を返す
    int n=mesh->n;
    int dim=mesh->dim;
    double **npxy;
    npxy = mesh->npxy;
    int **elnp;
    elnp=mesh->elnp;

    // Kl番目節点を通る行列の作成.Klのみに依存
    double **Pxy=dmatrix(1,n,1,dim);//要素を構成する節点座標
    for(int i=1;i<=n;i++){
        int ver=elnp[Kl][i];//ele_numberを構成するi番目の節点番号
        for(int j=1;j<=dim;j++){
            Pxy[i][j]=npxy[ver][j];//要素を構成する節点座標
        }
    }

    //３点から面積計算
    double Sl=2*S(Pxy);//行列式
    // printf("K=%d:",Kl);
    double *C=dvector(0,dim);//return用の配列
    int counter=0,next_position=0;
    for(int i=1;i<=n;i++){
        int x1=1,y1=2;

        int n1=i%3+1,n2=(i+1)%3+1;
        //クラメルの公式による連立方程式の解
        double xi=Pxy[n1][x1],xi1=Pxy[n2][x1];
        double yi=Pxy[n1][y1],yi1=Pxy[n2][y1];
        C[0]=(yi1*xi-yi*xi1)/Sl;
        C[1]=(yi-yi1)/Sl;
        C[2]=(xi1-xi)/Sl;

        double ret=C[0]+C[1]*x+C[2]*y;
        if(ret>=0.0){
            counter+=1;
        }else{
            next_position=i;
        }
        // printf("%f,",ret);

    }
    
    free_dvector(C,0,dim);
    free_dmatrix(Pxy,1,n,1,dim);

    
    int Return_value;
    if(counter==n){//全重心座標が>0ならnを返す．counter番目の要素に含まれる
        Return_value=-10;
    }else{//そうでなければ次の値を
        Return_value=next_position;
    }
    // printf(":%d\n",Return_value);

    return Return_value;
}

//接点が領域に入ってるかどうかの確認
int count(mesh_t *mesh,double x,double y){
    int ne=mesh->ne;   
    int Kl=-1;
    for(int K=1;K<=ne;K++){//K番目の要素
        int check_inside=ele_inside(mesh,K,x,y);
        if(check_inside==-10){//nに等しい時Kl内部に接点がある
            Kl=K;
            // break;
        }
    }
    return Kl;
}

int comp(int *a,int *b,int m,int n){//配列が一致している個数．(同じ要素がない配列のみ)
    int ret=0;
    for(int i=m;i<=n;i++){
        for(int j=m;j<=n;j++){
            if(a[i]==b[j]){
                ret+=1;
            }
        }
        
    }
    return ret;
}

int comp_place(int *a,int *b,int m,int n){//配列が一致していない個数．(同じ要素がない配列のみ)
    // int count=0;
    int ret=-1;
    int non_same;
    for(int i=m;i<=n;i++){
        non_same=0;
        for(int j=m;j<=n;j++){//a[i]が全ての要素と一致しない場合non_same=n-m+1となる．
            if(a[i]!=b[j]){
                non_same+=1;
            }
        }

        if(non_same==(n-m+1)){//全ての要素と一致しないものの探索．全てと一致する場合は-1を返す
            ret=i;
        }

    }

    return ret;
}

void drawney(mesh_t *mesh){
    printf("%d\n",mesh->n);
    // int np=mesh->np,ne=mesh->ne,nb=mesh->nb,dim=mesh->dim,n=mesh->n;
    // double **npxy;
    // int **elnp, **bound;
    // npxy = mesh->npxy;
    // elnp = mesh->elnp;
    // bound = mesh->bound;

    // /*============1st========super Triangle============*/
    // //最大最小値の探索
    // double max=npxy[1][1];
    // double min=npxy[1][1];
    // for(int i=1;i<=np;i++){
    //     for(int j=0;j<dim;j++){
    //         if(max<=npxy[i][j]){
    //             max=npxy[i][j];
    //         }
    //         if(min>=npxy[i][j]){
    //             min=npxy[i][j];
    //         }
    //     }   
    // }
    // //上記４点によってなされる正方形に外接する円の情報
    // double L=0.5*sqrt(2)*fabs(max-min);
    // //三角形の内心=正方形の中心(x0,x0)
    // double x0=min+0.5*fabs(max-min);
    // double y0=x0;
    // //SuperTriangleを作る座標
    // double tx=x0,ty=y0+2.0*L;
    // double tx1=x0-sqrt(3)*L;
    // double ty1=y0-L;
    // double tx2=x0+sqrt(3)*L;
    // double ty2=y0-L;

    // // Super Triangleの座標行列
    // double **Coord=dmatrix(1,n,1,dim);
    // Coord[1][1]=tx;Coord[1][2]=ty;
    // Coord[2][1]=tx1;Coord[2][2]=ty1;
    // Coord[3][1]=tx2;Coord[3][2]=ty2;
    
    // /*=================================================*/

    // /*==============2nd==============線を引く============*/

    // int data_number=1;
    // double *Xi=dvector(1,dim);
    // for(int i=1;i<=dim;i++)Xi[i]=npxy[data_number][i];
    
    // // int C=count(Coord,Xi);

    // free_dmatrix(Coord,1,n,1,dim);
    // free_dvector(Xi,1,dim);
    
    // /*==================================================*/
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

    // int dim=mesh->dim;
    // int **elnp=mesh->elnp;
    // double ret=0.0;

    // int N_5=7;//積分点の個数
    // for(int num=1;num<=N_5;num++){
    //     double w;//積分点の重み
    //     if(num==1){
    //         w=9.0/40.0;
    //     }else if(2<=num && num<=4){
    //         w=(155.0-sqrt(15))/1200.0;
    //     }else{
    //         w=(155.0+sqrt(15))/1200.0;
    //     }

    //     double *x=Pi(mesh,Kl,num);//Klのnum番目の積分点の座標
    //     int ver1=elnp[Kl][i],ver2=elnp[Kl][j];//要素Klのi,j番目の節点番号

    //     double phi_i=phi(mesh,ver1,x[1],x[2],Kl);//Klの積分点における関数の値
    //     double phi_j=phi(mesh,ver2,x[1],x[2],Kl);//Klの積分点における関数の値

    //     ret+=w*phi_i*phi_j;

    //     free_dvector(x,1,dim);
    // }

    // if(fabs(S_Kl*ret-I)>0.001){
    //     printf("Integration Error!%f\n",fabs(ret-I));
    //     exit(1);
    // }
    
    return I;
}

double **Centroid_coord(mesh_t *mesh){//積分点の重心座標系
    int dim=mesh->dim;
    int N_5=7;
    double p=(6-sqrt(15))/21,q=(9+2*sqrt(15))/21,r=(6+sqrt(15))/21,s=(9-2*sqrt(15))/21;
    double t=1./3.;

    double **t_perm=dmatrix(1,N_5,1,dim+1);

    for(int n=1;n<=N_5;n++){
        if(n==1){
            /*=====================i=1(t,t;t)====================*/
            for(int i=1;i<=dim+1;i++)t_perm[n][i]=t;
            /*===================================================*/
        }else if(2<=n && n<=4){
            /*=====================i=2,3,4(p,p;q)================*/
            for(int i=1;i<=dim+1;i++){
                double dummy;
                if(n-1==i){//対角成分(qの入る位置を考える)2,3,4を1,2,3に揃える
                    dummy=q;
                }else{
                    dummy=p;
                }
                t_perm[n][i]=dummy;//1,2,3番目に重みを格納していく
            }
            /*===================================================*/
        }else{
            for(int i=1;i<=dim+1;i++){
            /*=====================i=5,6,7(r,r;s)===============*/
                double dummy;
                if(n-4==i){//対角成分(qの入る位置を考える)5,6,7を1,2,3に揃える
                    dummy=s;
                }else{
                    dummy=r;
                }
                t_perm[n][i]=dummy;//1,2,3番目に重みを格納していく
            }
            /*==================================================*/
        }
    }
    return t_perm;
}

double *Pi(mesh_t *mesh,int Kl,int K_n){//Kl上のK_n個目の積分点の座標
    int dim=mesh->dim;
    double **npxy;
    int **elnp;
    npxy = mesh->npxy;
    elnp = mesh->elnp;

    int N_5=7;//積分点

    double *ret=dvector(1,dim);
    double **Coef=Centroid_coord(mesh);
    for(int i=1;i<=dim+1;i++){
        double x_i=npxy[elnp[Kl][i]][1];
        double y_i=npxy[elnp[Kl][i]][2];

        ret[1]+=Coef[K_n][i]*x_i;
        ret[2]+=Coef[K_n][i]*y_i;
    }
    free_dmatrix(Coef,1,N_5,1,dim+1);
    return ret;
}

//２次元の積分点5
double **NumInt_deg_five(mesh_t *mesh,int K){
    /*==================構造体のデータ読み込み==========================*/
    int dim=mesh->dim;
    /*==============================================================*/

    int N_5=7;

    double **Coord=dmatrix(1,N_5,1,dim);//戻り値用の行列
    /*=====================積分点の定義==========================*/
    //順列を二次元配列(行列)として表現する
    for(int i=1;i<=N_5;i++){
        double *coord_i=Pi(mesh,K,i);
        for(int j=1;j<=dim;j++){
            Coord[i][j]=coord_i[j];
        }
        free_dvector(coord_i,1,N_5);
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

//i,jの座標間の距離
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

//l番目の要素の面積
double area(mesh_t *mesh,int l){
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

//Diriclet境界条件の反映
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

//pfuncと同様の型の関数の積分
double Int_div_five(pfunc func,mesh_t *mesh,int Kl,int i,int j){

    int dim=mesh->dim;
    // double **npxy=mesh->npxy;
    // int **elnp=mesh->elnp;

    int N_5=7;
    double w;
    double ret=0.0;
    double S=area(mesh,Kl);

    //関数の値
    double **Integer_point_bynaritic=Centroid_coord(mesh);

    for(int num=1;num<=N_5;num++){
        double *x=Pi(mesh,Kl,num);//Klのnum番目の積分点の座標
        if(num==1){
            w=9.0/40.0;
        }else if(2<=num && num<=4){
            w=(155.0-sqrt(15))/1200.0;
        }else{
            w=(155.0+sqrt(15))/1200.0;
        }

        ret+=w*func(mesh,Kl,i,j,x[1],x[2]);//積分点上での関数の値をもちいて重みをかけて和を取る
        
        free_dvector(x,1,dim);
    }

    free_dmatrix(Integer_point_bynaritic,1,N_5,1,dim+1);


    return S*ret;
}

//基底関数phi
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

//GLSC出力用の関数
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

//要素剛性行列の作成
double **Al(weak weak_form,mesh_t *mesh){
    printf("Al(要素剛性行列)\n");
    /*==================構造体のデータ読み込み==========================*/
    // int dim=mesh->dim;
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

    double *err=dvector(1,np);
    Diriclet(mesh,A,err);
    free_dvector(err,1,np);
    
    return A;
}

//外力項の離散化
double *out_force(out RHS,mesh_t *mesh,double *u_old){
    pfunc func_SUPG=&phi_ij;

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
                double I=Int(mesh,i,j,l)+SUPG*Int_div_five(func_SUPG,mesh,l,i,j);
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
        double *rhs_1=matrix_vector_product_CRS(M,u_old_vector,n);
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

    double **dummy=dmatrix(1,np,1,np);
    Diriclet(mesh,dummy,return_vector);
    free_dmatrix(dummy,1,np,1,np);

    return return_vector;
}

//Lp誤差
double err_Lp(mesh_t *mesh,double *u,double p,double t){//時刻tにおける解析解との誤差
    double ret=0.0;
    double Lp=0.0;

    if(p==INFINITY){
        for(int i=1;i<=mesh->np;i++){
            double x=mesh->npxy[i][1],y=mesh->npxy[i][2];
            double lp_err=pow(u[i]-u_exa(x,y,t),p);
            if(Lp<=lp_err){
                Lp=lp_err;
            }
        }
        Lp=ret;
    }else{
        for(int i=1;i<=mesh->np;i++){
            double x=mesh->npxy[i][1],y=mesh->npxy[i][2];
            double lp_err=pow(u[i]-u_exa(x,y,t),p);
            Lp+=lp_err;
        }
        ret=pow(Lp/(mesh->np),1/p);
    }
    return ret;  
}

//流速に依存した上流点の探索
int search_past_point(mesh_t *mesh,double x_n,double y_n,double *u,double dt){

    // int ne=mesh->ne;
    int **elel=mesh->elel;
    // int n=mesh->n;
    // int dim=mesh->dim;

    int inside_Kl=count(mesh,x_n,y_n);//x_n,y_nを含む要素番号l0
    double x_p=x_n-dt*u[1],y_p=y_n-dt*u[2];//一個前の時刻の座標x-u*x,１時刻前の方に行くにはもうちょい工夫が必要

    int past_kl;//=count(mesh,x_p,y_p);
    //xp,ypの値が負になる方向の要素に探索する．
    int check_Kl_pi=ele_inside(mesh,inside_Kl,x_p,y_p);//この値に応じて次に探索する方向が決まる
    for(;;){
        if(check_Kl_pi==-10){//Klに入ってるか探索
            past_kl=inside_Kl;//inside_Kl要素にx_p,y_pがある
            // printf("Coordinate is inside %d\n",past_kl);q
            break;
        }else{
            inside_Kl=elel[inside_Kl][check_Kl_pi];//次に考えなくてはならない要素番号
            
            if(inside_Kl==-1){//要素番号が-1．上流点が外側にある場合
                // inside_Kl=past_kl;
                break;
            }

            check_Kl_pi=ele_inside(mesh,inside_Kl,x_p,y_p);
        }
    }
    return past_kl;
}

