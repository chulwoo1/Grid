  LaplacianRatParams gpar(1),mpar(1);
gpar.offset =	-5	;
gpar.a0[0] =	-450	;
gpar.a1[0] =	0	;
gpar.b0[0] =	-17.	;
gpar.b1[0] =	-16.	;
gpar.b2=1.;
mpar.offset =	-0.2	;
mpar.a0[0] =	18	;
mpar.a1[0] =	0	;
mpar.b0[0] =	73	;
mpar.b1[0] =	-16.	;
mpar.b2=1.;
for(int i=0;i<1;i++){
gpar.a1[i] *=16.;
gpar.b1[i] *=16.;
mpar.a1[i] *=16.;
mpar.b1[i] *=16.;
}
gpar.b2 *= 16.*16.;
mpar.b2 *= 16.*16.;
