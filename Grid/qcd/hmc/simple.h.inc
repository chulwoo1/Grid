    LaplacianRatParams gpar(1),mpar(1);
    gpar.offset = 1.;
    gpar.a0[0] = 500.;
    gpar.a1[0] = 0.;
    gpar.b0[0] = 1.;
    gpar.b1[0] = 2.;
    gpar.b2=1.;
    mpar.offset = 1.;
    mpar.a0[0] = -1.62443017922;
    mpar.a1[0] = -1.54707654538;
    mpar.b0[0] = 5.97654563524;
    mpar.b1[0] = 6.74194794773;
    mpar.b2=1.;
    for(int i=0;i<2;i++){
       gpar.a1[i] *=16.;
       gpar.b1[i] *=16.;
       mpar.a1[i] *=16.;
       mpar.b1[i] *=16.;
    }
    gpar.b2 *= 16.*16.;
    mpar.b2 *= 16.*16.;
