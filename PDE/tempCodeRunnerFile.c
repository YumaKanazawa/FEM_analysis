
    double *u=CG_CRS(A,RHS,np);//解の計算
    for(int i=1;i<=np;i++)printf("%f\n",u[i]);