#include<stdint.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "tables.h"
#include "tables_g.h"

typedef union {
    double x;
    uint64_t i;
}u ;  

/*void Split(double a, double *ah, double *al){
    // We use algorithm of Veltkampf split
    unsigned long C =  1 << 27; 
    double z = C * a;
    *ah = z - ( z - a );
    *al = a - *ah;
    /*printf("ah = %la et al = %la \n", *ah, *al);*/
/*}

void FastTwoSum(double h, double l, double *x, double *y){
    *x = h + l;
    *y =  l - (*x - h);
}

void TwoSum(double h, double l, double *x, double *y){
    *x = h + l;
    double aprime = *x - l;
    double bprime = *x - aprime;
    double gammaa = h - aprime;
    double gammab = l - bprime;
    *y =  gammaa+gammab;
}

void FastTwoSum_modified(double a, double bh, double bl, double *x, double *y){
    double y1;
    FastTwoSum(a, bh, x, &y1);
    *y = y1+bl;
}

void TwoSum_modified(double a, double bh, double bl, double *x, double *y){
    double y1;
    TwoSum(a, bh, x, &y1);
    *y = y1+bl;
}

void FastTwoSum_modified2(double ah, double al, double bh, double bl, double *x, double *y){
    double y1;
    FastTwoSum(ah, bh, x, &y1);
    *y = y1+bl+al;
}

void TwoSum_modified2(double ah, double al, double bh, double bl, double *x, double *y){
    double y1;
    TwoSum(ah, bh, x, &y1);
    *y = y1+bl+al;
}*/

void Add112(double a, double b, double *s, double *t){
    *s = a+b ;
    double z = *s-a;
    *t = b-z;
}

void Add122(double a, double bh, double bl, double *s, double *t){
    double l;
    Add112(a, bh,s,&l);
    *t = l+bl;
}

void Add222(double ah, double al, double bh, double bl, double *s, double *t){
    double l,m;
    Add112(ah,bh,s,&l);
    m = l+al;
    *t = m+bl;
}


/* a is a double numbers
   bh,bm,bl is a triple-double numbers*/
void Add133(double a, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4;
    Add112(a,bh,rh,&t1); 
    Add112(t1,bm,&t2,&t3); 
    t4 = t3+bl;
    Add112(t2,t4,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

/* ah,al  is a double-double numbers
   bh,bm,bl is a triple-double numbers*/
void Add233(double ah, double al, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7;
    Add112(ah,bh,rh,&t1);
    Add112(al,bm,&t2,&t3);
    Add112(t1,t2,&t4,&t5);
    t6 = t3+bl;
    t7 = t6+t5;
    Add112(t4,t7,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

/*ah,am,al and bh,bm,bl are triple-double numbers*/
void Add333(double ah, double am, double al, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8;
    Add112(ah,bh,rh,&t1);
    Add112(am,bm,&t2,&t3);
    Add112(t1,t2,&t7,&t4);
    t6 = al+bl;
    t5 = t3+t4;
    t8 = t5+t6;
    Add112(t7,t8,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

void Mul112(double a, double b, double *r1, double *r2){
    *r1 = a * b;
    *r2 = __builtin_fma (a, b, -*r1);
}

void Mul122(double a, double bh, double bl, double *r1, double *r2){
    double t1,t2;
    Mul112(a,bh,&t1,&t2);
    double t3,t4;
    t3 = a*bl;
    t4 = t2+t3;
    Add112(t1,t4,r1,r2);
}

void Mul222(double ah, double al, double bh, double bl, double *r1, double *r2){
    double t1,t2;
    Mul112(ah, bh, &t1, &t2);
    double t3,t4;
    t3 = ah*bl;
    t4 = bh*al;
    double t5,t6;
    Add112(t3,t4,&t5,&t6);
    Add222(t1,t2,t5,t6,r1,r2);
}
/* a is a double numbers
   bh,bl is a double-double numbers*/
void Mul123(double a, double bh, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6;
    Mul112(a,bh,rh,&t1);
    Mul112(a,bl,&t2,&t3);
    Add112(t1,t2,&t5,&t4);
    t6 = t3+t4;
    Add112(t5,t6,rm,rl);
}

/* a is a double numbers
   bh,bm,bl is a triple-double numbers*/
void Mul133(double a, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    Mul112(a,bh,rh,&t2);
    Mul112(a,bm,&t3,&t4);
    t5 = a*bl;
    Add112(t2,t3,&t9,&t7);
    t8 = t4+t5;
    t10 = t7+t8;
    Add112(t9,t10,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

/* ah,al and bh,bl are  double-double numbers*/
void Mul223(double ah, double al, double bh, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    Mul112(ah,bh,rh,&t1);
    Mul112(ah,bl,&t2,&t3);
    Mul112(al,bh,&t4,&t5);
    t6 = al*bl;
    Add222(t2,t3,t4,t5,&t7,&t8);
    Add112(t1,t6,&t9,&t10);
    Add222(t7,t8,t9,t10,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

/* ah,al  is a double-double numbers
   bh,bm,bl is a triple-double numbers*/
void Mul233(double ah, double al, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    double t11,t12,t13,t14,t15,t16,t17,t18;
    Mul112(ah,bh,rh,&t1);
    Mul112(ah,bm,&t2,&t3);
    Mul112(ah,bl,&t4,&t5);
    Mul112(al,bh,&t6,&t7);
    Mul112(al,bm,&t8,&t9);
    t10 = al*bl;
    Add222(t2,t3,t4,t5,&t11,&t12);
    Add222(t6,t7,t8,t9,&t13,&t14);
    Add222(t11,t12,t13,t14,&t15,&t16);
    Add112(t1,t10,&t17,&t18);
    Add222(t17,t18,t15,t16,rm,rl); 
    /*rh,rm,rl is a triple-double numbers*/
}

void Mul333(double ah, double am, double al, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    double t11,t12,t13,t14,t15,t16,t17,t18,t19,t20;
    double t21,t22;
    Mul112(ah,bh,rh,&t1);
    Mul112(ah,bm,&t2,&t3);
    Mul112(am,bh,&t4,&t5);
    Mul112(am,bm,&t6,&t7);
    t8 = ah*bl;
    t9 = al*bh;
    t10 = am*bl;
    t11 = al*bm; 
    t12 = t8+t9;
    t13 = t10+t11;
    Add112(t1,t6,&t14,&t15);
    t16 = t7+t15;
    t17 = t12+t13;
    t18 = t16+t17;
    Add112(t14,t18,&t19,&t20);
    Add222(t2,t3,t4,t5,&t21,&t22);
    Add222(t21,t22,t19,t20,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}



/*void TwoSum_DEKKER(double a,double bh,double bl, double ch, double cl ,double *r1, double *r2){
    double x,y1,ah,al;
    Split(a,&ah,&al);
    TwoSum(ah, bh, &x, &y1);
    double y = y1+bl+al;
    DEKKER_Prod_Modified1( x, y, ch, cl, r1, r2);
}*/
int main(){
    /*printf("Test for Split\n");
    double x1 = 0x1.906f0ca2f82fbp+2;
    double h,l;
    Split(x1, &h, &l);
    printf(" x = %la , h = %la , l = %la\n" , x1, h, l);

    printf("Test for DEKKER_Prod\n");
    double a1 = 0x1.421f273fc712bp+9;
    double b1 = 0x1.e707af8e689a7p+9;
    double h1,l1;
    DEKKER_Prod(a1,b1,&h1,&l1);
    printf("a = %la, b = %la, h = %la, l = %la\n ", a1, b1, h1, l1);

    printf("Test for Mul12\n");
    double a2 = 0x1.421f273fc712bp+9;
    double b2 = 0x1.e707af8e689a7p+9;
    double h2,l2;
    Mul12(a2,b2,&h2,&l2);
    printf("a = %la, b = %la, h = %la, l = %la\n ", a2, b2, h2, l2);

    printf("Test for DEKKER_Prod_Modified2\n");
    double a3 = 0x1.c2148ef302f3fp+9;
    double b3 = 0x1.b9a87a94da0c1p+9;
    double b3h,b3l;
    Split(b3, &b3h, &b3l);
    double h3,l3;
    DEKKER_Prod_Modified2(a3,b3h,b3l,&h3,&l3);
    printf("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " , a3, b3, b3h, b3l, h3, l3);

    printf("Test for Mul122\n");
    double a4 = 0x1.c2148ef302f3fp+9;
    double b4 = 0x1.b9a87a94da0c1p+9;
    double b4h,b4l;
    Split(b4, &b4h, &b4l);
    double h4,l4;
    Mul122(a4,b4h,b4l,&h4,&l4);
    printf("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " , a4, b4, b4h, b4l, h4, l4);

    printf("Test for DEKKER_Prod_Modified1\n");
    double a5 = 0x1.39ad06459895p+9;
    double a5h,a5l;
    Split(a5, &a5h, &a5l);
    double b5 = 0x1.5bc29786125eap+7;
    double b5h,b5l;
    Split(b5, &b5h, &b5l);
    double h5,l5;
    DEKKER_Prod_Modified1(a5h,a5l,b5h,b5l,&h5, &l5);
    printf("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n " ,a5, b5, a5h, a5l, b5h, b5l, h5, l5);

    printf("Test for Mul22\n");
    double a6 = 0x1.39ad06459895p+9;
    double a6h,a6l;
    Split(a6, &a6h, &a6l);
    double b6 = 0x1.5bc29786125eap+7;
    double b6h,b6l;
    Split(b6, &b6h, &b6l);
    double h6,l6;
    Mul22(a6h,a6l,b6h,b6l,&h6, &l6);
    printf("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n " ,a6, b6, a6h, a6l, b6h, b6l, h6, l6);

    printf("Test for Mul133\n");
    double a7 = 0x1.5b2057ad879acp+8;
    double b7h = 0x1.c0b5b34724d3p+17;
    double b7m =-0x1.18440552f325bp-12;
    double b7l =0x1.f109p-66;
    double h7,m7,l7;
    Mul133(a7,b7h,b7m,b7l,&h7,&m7,&l7);
    printf("a = %la, bh = %la, bl = %la, bm = %la, h = %la, m = %la , l = %la\n " ,a7, b7h, b7m, b7l, h7, m7, l7);

    printf("Test for Mul23\n");
    double a8 = 0x1.d9bb8bb2434a7p+8;
    double a8h,a8l; 
    Split(a8,&a8h,&a8l);
    double b8 = 0x1.969146a9b9d69p+8;
    double b8h,b8l;
    Split(b8,&b8h,&b8l);
    double h8,m8,l8; 
    Mul23(a8h,a8l,b8h,b8l,&h8,&m8,&l8);
    printf("a = %la, ah = %la, al = %la,  b = %la, bh = %la, bl = %la, h = %la, m = %la , l = %la\n " ,a8, a8h, a8l, b8h, b8l, h8, m8, l8);

    printf("Test for Mul233\n");
    double a9 = 0x1.e716676d322c6p+9;
    double a9h,a9l;
    Split(a9,&a9h,&a9l);
    double b9h = 0x1.81260f7a7c404p+18;
    double b9m = 0x1.db484eca53ea2p-9;
    double b9l = -0x1.73c54ep-63;
    double h9,m9,l9;
    Mul233(a9h,a9l,b9h,b9m,b9l,&h9,&m9,&l9);
    printf("a = %la, ah = %la, al = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n ",a9,a9h,a9l,b9h,b9m,b9l,h9,m9,l9);

    printf("Test for Mul33\n");
    double a10h = 0x1.4b7b1763d5c36p+19;
    double a10m = -0x1.2eac725f98608p-7;
    double a10l = 0x1.d58524p-61;
    double b10h = 0x1.6906123c043bp+19;
    double b10m = -0x1.1b638b9cda375p-8;
    double b10l =-0x1.eae1fp-62;
    double h10,m10,l10; 
    Mul33(a10h,a10m,a10l,b10h,b10m,b10l,&h10,&m10,&l10);
    printf("ah = %la, am = %la, al = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n ",a10h,a10m,a10l,b10h,b10m,b10l,h10,m10,l10);

    printf("Test for FastTwoSum\n");
    double a11 = 0x1.09e5334316b0ap+9;
    double b11 = 0x1.6afe66d00be08p+8;
    double h11,l11;
    FastTwoSum(a11,b11, &h11, &l11);
    printf("a = %la, b = %la, h = %la, l = %la\n ", a11, b11, h11, l11);

    printf("Test for TwoSum\n");
    double a12 = 0x1.09e5334316b0ap+9;
    double b12 = 0x1.6afe66d00be08p+8;
    double h12,l12;
    TwoSum(a12,b12, &h12, &l12);
    printf("a = %la, b = %la, h = %la, l = %la\n ", a12, b12, h12, l12);

    printf("Test for Add12\n");
    double a13 = 0x1.09e5334316b0ap+9;
    double b13 = 0x1.6afe66d00be08p+8;
    double h13,l13;
    Add12(a13,b13, &h13, &l13);
    printf("a = %la, b = %la, h = %la, l = %la\n ", a13, b13, h13, l13);

    printf("Test for FastTwoSum_modified\n");
    double a14 = 0x1.01f03e5b06beep+9;
    double b14 = 0x1.3f20ef1cf5991p+7;
    double b14h,b14l;
    Split(b14, &b14h, &b14l);
    double h14, l14;
    FastTwoSum_modified(a14, b14h, b14l, &h14, &l14);
    printf("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " , a14, b14, b14h, b14l, h14, l14);

    printf("Test for TwoSum_modified\n");
    double a15 = 0x1.01f03e5b06beep+9;
    double b15 = 0x1.3f20ef1cf5991p+7;
    double b15h,b15l;
    Split(b15, &b15h, &b15l);
    double h15, l15;
    TwoSum_modified(a15, b15h, b15l, &h15, &l15);
    printf("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " , a15, b15, b15h, b15l, h15, l15);

    printf("Test for Add122\n");
    double a16 = 0x1.01f03e5b06beep+9;
    double b16 = 0x1.3f20ef1cf5991p+7;
    double b16h,b16l;
    Split(b16, &b16h, &b16l);
    double h16, l16;
    Add122(a16, b16h, b16l, &h16, &l16);
    printf("a = %la, b = %la, bh = %la, bl = %la , h = %la, l = %la\n " , a16, b16, b16h, b16l, h16, l16);
   

    printf("Test for FastTwoSum_modified2\n");
    double a17 = 0x1.1b33e8b39d9f1p+4;
    double a17h,a17l;
    Split(a17, &a17h, &a17l);
    double b17 = 0x1.e40fc62966a64p+9;
    double b17h,b17l;
    Split(b17, &b17h, &b17l);
    double h17,l17;
    FastTwoSum_modified2(a17h,a17l, b17h, b17l, &h17, &l17);
    printf("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n ", a17, b17, a17h, a17l, b17h, b17l, h17, l17);
    
    printf("Test for TwoSum_modified2\n");
    double a18 = 0x1.1b33e8b39d9f1p+4;
    double a18h,a18l;
    Split(a18, &a18h, &a18l);
    double b18 = 0x1.e40fc62966a64p+9;
    double b18h,b18l;
    Split(b18, &b18h, &b18l);
    double h18,l18;
    TwoSum_modified2(a18h,a18l, b18h, b18l, &h18, &l18);
    printf("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n ", a18, b18, a18h, a18l, b18h, b18l, h18, l18);

    printf("Test for Add22\n");
    double a19 = 0x1.1b33e8b39d9f1p+4;
    double a19h,a19l;
    Split(a19, &a19h, &a19l);
    double b19 = 0x1.e40fc62966a64p+9;
    double b19h,b19l;
    Split(b19, &b19h, &b19l);
    double h19,l19;
    Add22(a19h,a19l, b19h, b19l, &h19, &l19);
    printf("a = %la, b = %la, ah = %la, al = %la, bh = %la, bl = %la , h = %la, l = %la\n ", a19, b19, a19h, a19l, b19h, b19l, h19, l19);

    printf("Test for Add133\n");
    double a20 = 0x1.ce7079d3ea68p+7;
    double b20h = 0x1.23b1a4cf5f82p+18;
    double b20m = 0x1.80d37956f706fp-9;
    double b20l = -0x1.2ccap-67;
    double h20,l20,m20; 
    Add133(a20,b20h,b20m,b20l,&h20,&m20,&l20);
    printf("a = %la, bh = %la, bl = %la, bm = %la, h = %la, m = %la , l = %la\n " ,a20, b20h, b20m, b20l, h20, m20, l20);

    printf("Test for Add233\n");
    double a21 = 0x1.20a247f9999edp+9;
    double a21h,a21l;
    Split(a21,&a21h,&a21l);
    double b21h = 0x1.0dd04f65b5954p+18;
    double b21m = -0x1.9b0c2cc926733p-9;
    double b21l = -0x1.46ccd4p-65;
    double h21,l21,m21;
    Add233(a21h,a21l,b21h,b21m,b21l,&h21,&m21,&l21);
    printf("a = %la, ah = %la, al = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n ",a21,a21h,a21l,b21h,b21m,b21l,h21,m21,l21);

    printf("Test for Add33\n");
    double a22h = 0x1.2937ce70794fp+17;
    double a22m = 0x1.2ffbef6672a88p-11;
    double a22l = 0x1.6bd29cp-65;
    double b22h = 0x1.be1e53af7546p+17;
    double b22m = 0x1.031eb811ecc19p-9;
    double b22l = -0x1.16769p-64;
    double h22,m22,l22; 
    Mul33(a22h,a22m,a22l,b22h,b22m,b22l,&h22,&m22,&l22);
    printf("ah = %la, am = %la, al = %la, bh = %la, bm = %la, bl = %la, , h = %la, m = %la, l = %la\n ",a22h,a22m,a22l,b22h,b22m,b22l,h22,m22,l22);
*/
    printf("Test for Mul123\n");
    double a23 = 0x1.158922b181abfp+11;
    double b23h = 0x1.65580ca8ecfe5p+8;
    double b23l = 0x1.65580ca8ecfe5p-46;
    double h23,m23,l23; 
    Mul123(a23,b23h,b23l,&h23,&m23,&l23);
    printf("a = %la, bh = %la, bl = %la , h = %la, m = %la, l = %la\n " ,a23,b23h,b23l,h23,m23,l23);

    return 0;
}