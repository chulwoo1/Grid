//    LaplacianRatParams gpar(2),mpar(2);
      LaplacianRatParams gpar(1,2),mpar(1,2);
//    int gpar_order=1;

//    int gpar_order=2;
//    gpar.offset = 1.;
    gpar.poly[0] = 1.;
    gpar.a0[0] =  1.999999999999265e+04;
    gpar.a1[0] = 0.;
    gpar.b0[0] = 4.999999999999944e-01;
    gpar.b1[0] = 1.;
    gpar.b2[0] = 0;
    gpar.a0[1] = -1.299999999999689e+04;
    gpar.a1[1] = -1.999999999999510e+04;
    gpar.b0[1] = 0.330;
    gpar.b1[1] = 1.15;
    gpar.b2[1]=1;

    mpar.poly[0] = 1.;
//    mpar.offset = 1.;
#if 1
    mpar.a0[0] =  -1.547136432463019e+00;
    mpar.a1[0] = 0.;
    mpar. b0[0] = 5.191768369836859e+00;
    mpar. b1[0] = 1.;
    mpar.b2[0] = 0;
    mpar.a0[1] =  -1.351197287420350e+01;
    mpar.a1[1] = 1.547136432463019e+00;
    mpar.b0[1] = 1.929304099580766e+01;
    mpar.b1[1] = -3.541768369836859e+00;
    mpar.b2[1]= 1;
#endif

