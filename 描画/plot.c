#include "../main_mesh.h"//メッシュ作成の関数


//ディリクレ境界条件用の関数(ないとエラー吐くので...)
double g(double x,double y){
    return x+y;
}
double g1(double x,double y){
    return 0.0;
}
double phi_ij(mesh_t *mesh,int Kl,int i,int j,double x,double y){
    return 0.0;
}

#include<glsc3d_3.h>//glsc3Dを用いるときに使用

//標準座標系のパラメータ
#define back_x 1500
#define back_y 1000

//仮想座標系の範囲指定
#define x_max 1//グラフの横幅の長さ，最大値と最小値の差
#define x_min -1
#define y_max 1
#define y_min -1
#define z_max 0.5
#define z_min 0


#define base_x 10
#define base_y 10
#define size_x 1000
#define size_y 1000

/*============================================描画の関数================================================================*/
void plot_set(){
    /*--------ここは脳死でコピペ----------------*/
    g_enable_highdpi();
    g_set_antialiasing(4);
    g_init("標準座標",back_x,back_y);//標準座標系の作成。上部の文字日本語対応
    g_scr_color(1,1,1);//標準座標系の背景カラーrgb
    g_cls();//標準座標系の初期化
    
    // int pointnumber=2;//(pointnumber×pointnumber)の形に配置する
    // double def_pointx=base_x+size_x*1.1*(j%pointnumber);//グラフ描画の基準点(x)
    // double def_pointy=base_y+size_y*1.1*(j/pointnumber);//グラフ描画の基準点(y)

    // /*------------------------------標準座標系の定義---------------------------*/
    // g_line_color(0,0,0,1);
    // g_line_width(2);
    // g_def_scale_2D(j,x0,x1+x0,y_min,y_max,def_pointx,def_pointy,size_x,size_y);//標準座標系の定義
    // g_sel_scale(j);
    // g_box_2D(x0,x1+x0,y_min,y_max,G_YES,G_NO);
    // g_def_scale_3D(
    //     j,//グラフの番号
    //     x_min,x_max,//よくわからん
    //     y_min,y_max,//よくわからん
    //     z_min,z_max,//よくわからん
    //     x_min,x_max,//x_0_f~x_1_fの値を描画
    //     y_min,y_max,//y_0_f~y_1_f
    //     z_min,z_max,//z_0_f~z_1_f
    //     base_x,base_y,//枠の基準点
    //     size_x,size_y//枠の長さ
    //     );//標準座標系の定義

    // g_box_3D(x_min,x_max,//x_0_f~x_1_fの値を描画
    //     y_min,y_max,//y_0_f~y_1_f
    //     z_min,z_max,//z_0_f~z_1_f,
    //     G_YES,//枠線の描画
    //     G_YES//塗りつぶし
    //     );
    // /*---------------------------------------------------------*/
    //画像出力
    // g_capture_set("Plot");//保存するフォルダ名の指定も可能
}


void graph_shape(int dim){
    g_cls();
   int j=0;
//    g_text_standard(back_x-200,base_y+50,"Re：Reynolds Number");//kの描画
   if(dim==2){
        g_def_scale_2D(
            j,
            x_min,x_max,//よくわからん
            y_min,y_max,//よくわからん
            base_x,base_y,//枠の基準点
            size_x,size_y//枠の長さ
            );//標準座標系の定義

        g_sel_scale(j);

        g_box_2D(
            x_min,x_max,//よくわからん
            y_min,y_max,//よくわからん
            G_YES,//枠線の描画
            G_NO//塗りつぶし
            );//標準座標系の定義
   }else if(dim==3){

    // g_def_scale_3D(
    //     j,//グラフの番号
    //     x_min,x_max,//よくわからん
    //     y_min,y_max,//よくわからん
    //     z_min,z_max,//よくわからん
    //     x_min,x_max,//x_0_f~x_1_fの値を描画
    //     y_min,y_max,//y_0_f~y_1_f
    //     z_min,z_max,//z_0_f~z_1_f
    //     base_x,base_y,//枠の基準点
    //     size_x,size_y//枠の長さ
    //     );//標準座標系の定義

    g_def_scale_3D_fix(//数値と描画の範囲が同じ時
        j,//グラフの番号
        x_min,x_max,//x_0_f~x_1_fの値を描画
        y_min,y_max,//y_0_f~y_1_f
        z_min,z_max,//z_0_f~z_1_f
        base_x,base_y,//枠の基準点
        size_x,size_y//枠の長さ
        );//標準座標系の定義

    // g_vision(
    //     j,
    //     1,1,0.6,//視点位置
    //     0,0,1,//画面上方向の指定ベクトル
    //     1//拡大率
    // );
    g_sel_scale(j);

    g_box_3D(
        x_min,x_max,//x_0_f~x_1_fの値を描画
        y_min,y_max,//y_0_f~y_1_f
        z_min,z_max,//z_0_f~z_1_f,
        G_NO,//枠線の描画
        G_NO//塗りつぶし
        );

        //x方向の軸
        g_arrow_3D(
            -1.5, 0.0, 0,//始点
            1.0, 0.0, 0.0,//方向ベクトル
            5.0, //長さ
            0.25, //矢印
            G_YES, G_NO);

        //y方向の軸
        g_arrow_3D(
            0, -1.5, 0,//始点
            0.0, 1.0, 0.0,//方向ベクトル
            5.0, //長さ
            0.25, //矢印
            G_YES, G_NO);

        //z方向の軸
        g_arrow_3D(
            0, 0, -0.5,//始点
            0.0, 0.0, 1.0,//方向ベクトル
            5.0, //長さ
            0.25, //矢印
            G_YES, G_NO);

        g_circle_2D(
            0,0,//中心座標
            1,//半径
            G_YES,G_NO);
   }else{
        printf("visualize dimension is high!\n");
        exit(1);
   }

}

//描画用のデータ配列
double **Result_data(mesh_t *mesh,char *str){
    FILE *fp;
    
    if((fp=fopen(str,"r"))==NULL){
        printf("Can't open File!\n");
        exit(1);
    }

    double **ret=dmatrix(1,3*(mesh->ne),1,(mesh->n));

    int i=1;
    double data_x,data_y,ans;
    int Kl;
    while((fscanf(fp,"%lf %lf %lf %d",&data_x,&data_y,&ans,&Kl))!=EOF){
        ret[i][1]=data_x;
        ret[i][2]=data_y;
        ret[i][3]=ans;

        // printf("i=%d x=%f y=%f u(x,y)=%f \n",i,ret[i][1],ret[i][2],ret[i][3]);
        i++;
    }
    printf("total=%d,i=%d\n",3*(mesh->ne),i);
    return ret;
}

double u_exa(double x,double y,double t){
    return (1+t*t)*(x/y);
}

void Result_plot(mesh_t *mesh,int T,int dim,int waitime){
    double t=T*0.01;
    char str[200];

    sprintf(str,"../PDE/figure/mesh%d.dat",T);
    printf("%s\n",str);

    double **A=Result_data(mesh,str);

    for(int i=1;i<=3*(mesh->ne);i+=3){
        double x0=A[i][1],y0=A[i][2],z0=A[i][3];
        int i1=i+1;
        double x1=A[i1][1],y1=A[i1][2],z1=A[i1][3];
        int i2=i1+1;
        double x2=A[i2][1],y2=A[i2][2],z2=A[i2][3];

        double color=fabs(z0)+fabs(z1)+fabs(z2);

        if(dim==2){
            g_area_color(1-color,1-color,1-color,1);
            g_triangle_2D(
                x0,y0,
                x1,y1,
                x2,y2,
                G_NO,
                G_YES
            );

            // double color_exa=fabs(u_exa(x0,y0,t))+fabs(u_exa(x1,y1,t))+fabs(u_exa(x2,y2,t));
            // g_area_color(color,0,0,1);

            // g_triangle_2D(
            //     x0,y0,
            //     x1,y1,
            //     x2,y2,
            //     G_NO,
            //     G_YES
            // );


        }else{
            g_area_color(1,0,0,1);
            g_text_color(1,0,0,1);
            g_text_standard(100,100,"Num_sol");
            g_triangle_3D(
                x0,y0,z0,
                x1,y1,z1,
                x2,y2,z2,
                G_NO,
                G_YES
            );

            // g_area_color(0,0,1,1);
            // g_text_color(0,0,1,1);
            // g_text_standard(100,150,"u_exa");
            // g_triangle_3D(
            //     x0,y0,u_exa(x0,y0,t),
            //     x1,y1,u_exa(x1,y1,t),
            //     x2,y2,u_exa(x2,y2,t),
            //     G_NO,
            //     G_YES
            // );
        }
    }
    g_finish();//ここまでのg_関数を実行
    g_sleep(waitime);//描画の時間，-1の時は無限に

    free_dmatrix(A,1,(mesh->n)*(mesh->ne),1,(mesh->n));
}

int main(int argc,char *argv[]){
    mesh_t mesh;
    if(argc < 4){printf("Usage:./plot dim ../Mesh/mesh0.msh dim_v\n"); exit(1);}//実行の仕方

    //alloc(and scan)
    alloc_scan_mesh(&mesh,argv[1],argv[2]);//argv[1]=次元の数　argv[2]=meshファイルの名前

    plot_set();
    int waitime;
    printf("please enter waitime:");
    scanf("%d",&waitime);
    printf("Enterd!\n");

    int dim_v=atoi(argv[3]);
    
    for(int T=0;T<3140;T++){
        graph_shape(dim_v);
        
        char str[200];
        sprintf(str,"t=%d",T);
        g_text_standard(back_x-200,base_y+50,str);//kの描画
        Result_plot(&mesh,T,dim_v,waitime);
    }
    
    mesh_free(&mesh);//free

    return 0;
}

