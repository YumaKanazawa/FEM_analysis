//datファイルの動画化
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h> 

#include<glsc3d_3.h>//glsc3Dを用いるときに使用

//標準座標系のパラメータ
#define back_x 1500
#define back_y 1000

//仮想座標系の範囲指定
#define x0 0
#define x1 100//グラフの横幅の長さ，最大値と最小値の差
#define base_x 10
#define base_y 10
#define base_z 10
#define size_x 1500
#define size_y 900
#define size_z 800

#define N 1000
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

    //画像出力
    // g_capture_set("Plot");//保存するフォルダ名の指定も可能
}

void Plot_shape(int j){
    g_line_color(0,0,0,1);//色を黒くする
    g_line_width(2);//線の太さ

    double origin_x=-1,origin_y=-1,origin_z=-1;
    double l_x=2,l_y=2,l_z=2;
    g_def_scale_3D(
        j,
        origin_x,origin_x+l_x,
        origin_y,origin_y+l_y,
        origin_z,origin_z+l_z,
        origin_x,origin_x+l_x,
        origin_y,origin_y+l_y,
        origin_z,origin_z+l_z,
        base_x,base_y,size_x,size_y);//標準座標系の定義
    g_sel_scale(j);
    g_box_3D(
        origin_x,origin_x+l_x,
        origin_y,origin_y+l_y,
        origin_z,origin_z+l_z,
        G_YES,G_NO);

    // char y_M[50];
    // sprintf(y_M,"%f",y_max);
    // g_text_standard(base_x+size_x,base_y+20,y_M);//立軸の描画
    // char y_m[50];
    // sprintf(y_m,"%f",y_min);
    // g_text_standard(base_x+size_x,base_y+size_y,y_m);//横軸の描画
}


double **plot_data(int data_num){

    double u[N][3]={{0.0}};
    char str[200];
    sprintf(str,"mesh%d.dat",data_num);
    FILE *fp;
    fp=fopen(str,"r");//中間層行列データ

    double data_x,data_y,ans;
    int i=0;
    while(fscanf(fp,"%lf %lf %lf\n",&data_x,&data_y,&ans)!=NULL){
        u[i][0]=data_x;
        u[i][1]=data_y;
        u[i][2]=ans;
        i++;
    }
    return u;
}


int main(void){
    /*datファイルの読み込み*/
    int data_num;
    printf("Data_number?\n");

    scanf("%d",&data_num);

    double **ans=plot_data(data_num);

    plot_set();
    Plot_shape(0);

    g_finish();//ここまでのg_関数を実行
    g_sleep(-1);//描画の時間，-1の時は無限に
}