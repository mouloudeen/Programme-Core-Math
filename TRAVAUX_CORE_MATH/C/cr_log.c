#include<stdint.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <time.h>
#include <mpfr.h>
#include <fenv.h>
#include "tables.h"
#include "tables_g.h"



typedef union {
    double x;
    uint64_t i;
}u ;  


/* Compute ulp(x)*/
double ulp_0_5(double x){
    double psi,a,b;
    psi = 0x1p-53 + 0x1p-54;
    a = __builtin_fma(x,psi,x);
    b = a-x;
    return fabs(0.5*b);
}   

/*s represents the sign 
 e  represents the exposant 
 m represents fraction */
void extract (double x, int *s, int *e, uint64_t *m ){
    u u;
    u.x = x;
    // We shift by 63 bits to find the most significant bit .
    *s = u.i >> 63;
    // We shift by 52 bits to remove the fraction bits.
    // We uses a mask to remove the sign bit.
    *e = (u.i >> 52) & 0x7ff;
    // We uses a mask to remove  the first 12 bits. 
    *m = u.i & 0xFFFFFFFFFFFFF;
   
}

    



void Add112(double a, double b, double *s, double *t){
    *s = a+b ;
    double z = *s-a;
    *t = b-z;
}

void Add112Cond(double a, double b, double *s, double *t){
    if (fabs(b) > fabs(a)){
        Add112(b,a,s,t);
    }
    else{
       Add112(a,b,s,t); 
    }
}

void Add122(double a, double bh, double bl, double *s, double *t){
    double l;
    Add112(a, bh,s,&l);
    *t = l+bl;
}

void Add122Cond(double a, double bh, double bl, double *s, double *t){
    double l;
    Add112Cond(a, bh,s,&l);
    *t = l+bl;
}

void Add222(double ah, double al, double bh, double bl, double *s, double *t){
    double l,m;
    Add112(ah,bh,s,&l);
    m = l+al;
    *t = m+bl;
}

void Add222Cond(double ah, double al, double bh, double bl, double *s, double *t){
    double l,m;
    Add112Cond(ah,bh,s,&l);
    m = l+al;
    *t = m+bl;
}

/* a is a double numbers
   bh,bm,bl is a triple-double numbers*/
void Add133(double a, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4;
    Add112Cond(a,bh,rh,&t1); 
    Add112(t1,bm,&t2,&t3); 
    t4 = t3+bl;
    Add112(t2,t4,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

/* ah,al  is a double-double numbers
   bh,bm,bl is a triple-double numbers*/
void Add233(double ah, double al, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7;
    Add112Cond(ah,bh,rh,&t1);
    Add112Cond(bm,al,&t2,&t3);
    Add112Cond(t1,t2,&t4,&t5);
    t6 = t3+bl;
    t7 = t6+t5;
    Add112Cond(t4,t7,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

/*ah,am,al and bh,bm,bl are triple-double numbers*/
void Add333(double ah, double am, double al, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8;
    Add112Cond(ah,bh,rh,&t1);
    Add112Cond(am,bm,&t2,&t3);
    Add112Cond(t1,t2,&t7,&t4);
    t6 = al+bl;
    t5 = t3+t4;
    t8 = t5+t6;
    Add112Cond(t7,t8,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

void Mul112(double a, double b, double *r1, double *r2){
    *r1 = a * b;
    *r2 = __builtin_fma (a, b, -*r1);
}

void Mul112Cond(double a, double b, double *r1, double *r2){
   if (fabs(b) > fabs(a)){
        Mul112(b,a,r1,r2);
    }
    else{
       Mul112(a,b,r1,r2); 
    }
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
    t5 = t3+t4;
    t6 = t2+t5;
    Add112Cond(t1,t6,r1,r2);
}
/* a is a double numbers
   bh,bl is a double-double numbers*/
void Mul123(double a, double bh, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6;
    Mul112(a,bh,rh,&t1);
    Mul112(a,bl,&t2,&t3);
    Add112Cond(t1,t2,&t5,&t4);
    t6 = t3+t4;
    Add112Cond(t5,t6,rm,rl);
}

/* a is a double numbers
   bh,bm,bl is a triple-double numbers*/
void Mul133(double a, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    Mul112(a,bh,rh,&t2);
    Mul112(a,bm,&t3,&t4);
    t5 = a*bl;
    Add112Cond(t2,t3,&t9,&t7);
    t8 = t4+t5;
    t10 = t7+t8;
    Add112Cond(t9,t10,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

/* ah,al and bh,bl are  double-double numbers*/
void Mul223(double ah, double al, double bh, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    Mul112(ah,bh,rh,&t1);
    Mul112(ah,bl,&t2,&t3);
    Mul112(bh,al,&t4,&t5);
    t6 = al*bl;
    Add222Cond(t2,t3,t4,t5,&t7,&t8);
    Add112Cond(t1,t6,&t9,&t10);
    Add222Cond(t7,t8,t9,t10,rm,rl);
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
    Mul112(bh,al,&t6,&t7);
    Mul112(al,bm,&t8,&t9);
    t10 = al*bl;
    Add222Cond(t2,t3,t4,t5,&t11,&t12);
    Add222Cond(t6,t7,t8,t9,&t13,&t14);
    Add222Cond(t11,t12,t13,t14,&t15,&t16);
    Add112Cond(t1,t10,&t17,&t18);
    Add222Cond(t17,t18,t15,t16,rm,rl); 
    /*rh,rm,rl is a triple-double numbers*/
}

void Mul333(double ah, double am, double al, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    double t11,t12,t13,t14,t15,t16,t17,t18,t19,t20;
    double t21,t22;
    Mul112Cond(ah,bh,rh,&t1);
    Mul112Cond(ah,bm,&t2,&t3);
    Mul112Cond(bh,am,&t4,&t5);
    Mul112Cond(am,bm,&t6,&t7);
    t8 = ah*bl;
    t9 = al*bh;
    t10 = am*bl;
    t11 = al*bm; 
    t12 = t8+t9;
    t13 = t10+t11;
    Add112Cond(t1,t6,&t14,&t15);
    t16 = t7+t15;
    t17 = t12+t13;
    t18 = t16+t17;
    Add112Cond(t14,t18,&t19,&t20);
    Add222Cond(t2,t3,t4,t5,&t21,&t22);
    Add222Cond(t21,t22,t19,t20,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

void TwoSum_DEKKER(double a,double bh,double bl, double ch, double cl ,double *r1, double *r2){
    double x,y;
    
    Add122(a, bh, bl, &x, &y);
    Mul222( x, y, ch, cl, r1, r2);
}



static double cr_log_fast_path(double x){ 
    int s;
    int e;
    uint64_t m;
    extract (x,&s,&e, &m);
    
    /* If x is a negative number return NAN*/
    if ((s == 1) && (e != 0)){
        return NAN;
    }
    /*If x is NAN return NAN*/
    if ((s == 0)  && (e == 0x7ff) && (m != 0)){
        return NAN;
    }
    /*If x=+0 ou -0 return INFINITY*/
    if (((s==1) && (e == 0) && (m == 0)) || ((s == 0) && (e == 0) && (m == 0))){
        return -(0x1p1023 +0x1p1023);
    }
    /*If x = +INFINITY return +INFINITY */
    if ((s == 0) && (e == 0x7ff)){
        return 0x1p1023 +0x1p1023;
    }
    /*The coefficients of the approximation polynomial of degree 7.  */
    double f7 =  0x1.2152a2de69894p-3;
    double f6 =  -0x1.555147415c204p-3;
    double f5 =  0x1.999997342c184p-3;
    double f4 =  -0x1.ffffffff57268p-3;
    double f3 =  0x1.55555555554cep-2;
    double f2 =  -0x1p-1;
    double f1 =  0x1p+0;

    if ( (x> 0x1.fe7814e49392fp-1) && (x<1)) {
    double cr_log_fast_1 = (((f7*(x-1)+f6)*(x-1)+f5)*(x-1)+f4);

    double cr_log_fast_1h,cr_log_fast_1l;
    
    
    
    double hr,lr;
    Add112(-1,x,&hr,&lr);
    
    double h1,l1;
    /* multiplication calculations of the truncated polynomial by (r1,r2)*/
    Mul122(cr_log_fast_1,hr,lr,&h1, &l1);
    
    double h2,l2;
    /* add f3  and multiply by (h,l)*/
    TwoSum_DEKKER(f3, h1, l1, hr,lr, &h2, &l2);
   
    double h3,l3;
    /* add f2 and multiply by (h,l) */
     TwoSum_DEKKER(f2, h2, l2, hr,lr, &h3, &l3);
    
    double h4,l4;
    /* add f1  and multiply by (h,l)*/
     TwoSum_DEKKER(f1, h3, l3, hr,lr, &h4, &l4);
     
    return (h4+l4);

    }
    
    if ((x> 1) && (x < 0x1.00068db8bac71p+0)){
    double cr_log_fast_1 = (((f7*(x-1)+f6)*(x-1)+f5)*(x-1)+f4);

    
    
    double hr,lr;
    Add112(x, -1, &hr, &lr);
    
    double h1,l1;
    /* multiplication calculations of the truncated polynomial by (r1,r2)*/
    Mul122(cr_log_fast_1,hr,lr,&h1, &l1);
    
    double h2,l2;
    /* add f3  and multiply by (h,l)*/
    TwoSum_DEKKER(f3, h1, l1, hr,lr, &h2, &l2);
   
    double h3,l3;
    /* add f2 and multiply by (h,l) */
     TwoSum_DEKKER(f2, h2, l2, hr,lr, &h3, &l3);
    
    double h4,l4;
    /* add f1  and multiply by (h,l)*/
     TwoSum_DEKKER(f1, h3, l3, hr,lr, &h4, &l4);
     
    return (h4+l4);

    }
    
    double m1;
    /*If x is a subnormal */
    if ((s ==0) && (e == 0) && (m !=0)) {
        uint64_t v = m;
        e = e - 1023;
        v = v*2;
        while (v < 0x10000000000000) {
            v *= 2;
            e--;
        }
        
        m1 =v *0x1p-52;
       
        u u;
        u.x = m1;
        m = u.i & 0xFFFFFFFFFFFFF;
    }
    /*If x is normal */
    else {
        m1 = 1 + m*0x1p-52;
       
        e = e - 1023;
    }
    /*We shift by 44 bits for to get the first 8 bits.*/
    uint64_t i = m>>44;
    
    /*k=8
    alpha_i_m is the double close to 2^k/(2^k+i) such that its log is very close to a double.*/
    double alpha_i_m = table_alpha_i_modified[(int)i];
    
   
    double hr,lr;
    /*multiply (hm1,lm1) by (halpha_i_m,lalpha_i_m)*/
    Mul112(m1,alpha_i_m,&hr,&lr);
    
    double h,l;
    /*add (hr,lr) by -1*/
    Add122(-1,hr,lr,&h,&l);
    
    
    
    double cr_log_fast_1 = (((f7*h+f6)*h+f5)*h+f4);
    
    double h1,l1;
    /* multiplication calculations of the truncated polynomial by (r1,r2)*/
    Mul122(cr_log_fast_1,h,l,&h1, &l1);
    
    
    double h2,l2;
    /* add f3  and multiply by (h,l)*/
    TwoSum_DEKKER(f3, h1, l1, h,l, &h2, &l2);
    
    double h3,l3;
    /* add f2 and multiply by (h,l) */
     TwoSum_DEKKER(f2, h2, l2, h,l, &h3, &l3);
    
    double h4,l4;
    /* add f1  and multiply by (h,l)*/
     TwoSum_DEKKER(f1, h3, l3, h,l, &h4, &l4);
    
    double log_alpha_i_m = table_log_alpha_i_modified[i];
   
    
    double h5,l5;
    /* add (hlog_alpha_i_m, llog_alpha_i_m)*/
    Add122(log_alpha_i_m, h4, l4, &h5, &l5);
    
    
    
    
    double h_log2 = 0x1.62e42fefa38p-1;
    double l_log2 = 0x1.ef35793c7673p-45;
    double e_hlog2,e_llog2;
    /* e*log(2) in (double,double) precision.*/
    Mul122(e, h_log2, l_log2, &e_hlog2, &e_llog2);
    
    double h6,l6;

    /*adding with (h10,l10)*/
    Add222(e_hlog2, e_llog2,h5,l5,&h6,&l6);
    
    return (h6+l6);
    

}

static double cr_log_accurate_path(double x/*, FILE *fic*/){
    int s;
    int e;
    uint64_t m;
    extract (x,&s,&e, &m);

    /*The coefficients in double,double of the approximation polynomial of degree 11.  */
    static const double G[11][2]= {
        {0x1.6d24dd22a92bp-4, 0x1.3e16deaf82f98p-58},   /* g11 decomposed into double,double*/
        {-0x1.99889bcf944ecp-4, -0x1.14f5635667f7p-58}, /* g10 decomposed into double,double*/       
        {0x1.c71c5b2ab9b1bp-4, -0x1.bba8402b8dab8p-58}, /* g9 decomposed into double,double*/
        {-0x1.ffffffed56645p-4, 0x1.fa475383e2fep-60},  /* g8 decomposed into double,double*/
        {0x1.249249248d599p-3, 0x1.a5c92da9d350cp-57},  /* g7 decomposed into double,double*/
        {-0x1.555555555553bp-3, 0x1.887e7d473e27p-57},  /* g6 decomposed into double,double*/
        {0x1.999999999999ap-3, -0x1.affa1a4be0d04p-57}, /* g5 decomposed into double,double*/
        {-0x1p-2, 0x1.590555132p-72},                   /* g4 decomposed into double,double*/
        {0x1.5555555555555p-2, 0x1.55555540b5114p-56},  /* g3 decomposed into double,double*/
        {-0x1p-1, 0x1.b78p-98},                         /* g2 decomposed into double,double*/
        {0x1p+0, 0x0p+0}                                /* g1 decomposed into double,double*/
    };  

    static const double S[14][2]= {
{-0x1.259816ddc6b9ep-4,0x1.f229b0ebc95e8p-58}, /* s14 decomposed into double,double*/
{0x1.3b33cbc0e63c5p-4,-0x1.31408a44f4e3p-59}, /* s13 decomposed into double,double*/
{-0x1.5555cadde3778p-4,-0x1.bda92435f3d6p-60}, /* s12 decomposed into double,double*/
{0x1.745d1820b70e5p-4,-0x1.dbf2039485fp-63}, /* s11 decomposed into double,double*/
{-0x1.9999999a94ebp-4,-0x1.50d87266ccc78p-58}, /* s10 decomposed into double,double*/
{0x1.c71c71c71d304p-4,-0x1.c82589f54b808p-58}, /* s9 decomposed into double,double*/
{-0x1.0000000000003p-3,-0x1.8d4472e3c5b8p-60}, /* s8 decomposed into double,double*/
{0x1.2492492492492p-3,0x1.28c8644888e58p-57}, /* s7 decomposed into double,double*/
{-0x1.5555555555555p-3,-0x1.555649ca82244p-57}, /* s6 decomposed into double,double*/
{0x1.999999999999ap-3,-0x1.9999997747bfp-57}, /* s5 decomposed into double,double*/
{-0x1p-2,-0x1.52ep-96}, /* s4 decomposed into double,double*/
{0x1.5555555555555p-2,0x1.5555555555556p-56}, /* s3 decomposed into double,double*/
{-0x1p-1,0x0p+0}, /* s2 decomposed into double,double*/
{0x1p+0,0x0p+0}, /* s1 decomposed into double,double*/
};

    if ((x> 0x1.fba5e353f7ceep-1) && (x<1)) {
        double hr,lr ;
        /*Add x by -1 result r*/
        Add112(-1,x,&hr,&lr);
       
        double hf = S[0][0];
        double lf = S[0][1];
        
        double h,l;
        /*Multiply (hf,lf) by r result (h,l)*/
        Mul222(hr,lr,hf,lf,&h,&l);
        
        double hfi,lfi,h1,l1;

        for (int i = 1;i<14;i++){
            hfi = S[i][0];
            lfi = S[i][1];

            /*Add (hfi,lfi) by (h,l) result (h1,l1)*/
            Add222(hfi,lfi,h,l,&h1,&l1);

            /*multiply (hr,lr) by (h1,l1) result (h,l)*/
            Mul222(hr,lr,h1,l1,&h,&l);
            
        }
        /*fprintf(fic,",(%la,%la) \n",h,l);*/
        return (h+l);
    }

    if ((x> 1) && (x <= 0x1.0178e5916f543p+0)){
        double hr,lr ;
        /*Add x by -1 result r*/
        Add112(x,-1,&hr,&lr);
       
        double hf = S[0][0];
        double lf = S[0][1];
        
        double h,l;
        /*Multiply (hf,lf) by r result (h,l)*/
        Mul222(hr,lr,hf,lf,&h,&l);
        
        double hfi,lfi,h1,l1;

        for (int i = 1;i<14;i++){
            hfi = S[i][0];
            lfi = S[i][1];

            /*Add (hfi,lfi) by (h,l) result (h1,l1)*/
            Add222(hfi,lfi,h,l,&h1,&l1);

            /*multiply (hr,lr) by (h1,l1) result (h,l)*/
            Mul222(hr,lr,h1,l1,&h,&l);
            
        }
        /*fprintf(fic,",(%la,%la) \n",h,l);*/
        return (h+l);
    }
    double m1;
    /*If x is a subnormal */
    if ((s ==0) && (e == 0) && (m !=0)) {
        uint64_t v = m;
        e = e - 1023;
        v = v*2;
        while (v < 0x10000000000000) {
            v *= 2;
            e--;
        }
        
        m1 =v *0x1p-52;
       
        u u;
        u.x = m1;
        m = u.i & 0xFFFFFFFFFFFFF;
    }
    /*If x is normal */
    else {
        m1 = 1 + m*0x1p-52;
       
        e = e - 1023;
    }

    /*We shift by 44 bits for to get the first 8 bits.*/
    uint64_t i = m>>44;

    /*k=8
    halpha_i,lalpha_i is the double,double of 2^k/(2^k+i).*/
    double halpha_i = table_alpha_i[(int)i][0];
    double lalpha_i = table_alpha_i[(int)i][1];
    
    double hr,lr;
    /*Multiply (halpha_i,lalpha_i_m) by m1 result (hr,lr) */
    Mul122(m1,halpha_i, lalpha_i,&hr,&lr);
    
    double h,l;
    /*Add (hr,lr) by (-1) result (h,l)*/
    Add122(-1,hr,lr,&h,&l);
    
    double hf = S[0][0];
    double lf = S[0][1];

    double h1,l1;
    /*Multiply (hf,lf) by (h,l) result (h1,l1)*/
    Mul222(hf,lf,h,l,&h1,&l1);
    
    double hfi,lfi,h2,l2;

    for (int i = 1;i<14;i++){
        hfi = S[i][0];
        lfi = S[i][1];

        /*Add (hfi,lfi) by (h1,l1) result (h2,l2)*/
        Add222(hfi,lfi,h1,l1,&h2,&l2);

        /*multiply (h,l) by (h2,l2) result (h1,l1)*/
        Mul222(h2,l2,h,l,&h1,&l1);
       
    }
    /*k=8
    hlog_alpha_i,llog_alpha_i is the double,double of log(2^k/(2^k+i)).*/
    double hlog_alpha_i = table_log_alpha_i[(int)i][0];
    double llog_alpha_i = table_log_alpha_i[(int)i][1];


    double h5,l5;
    /*Add (h1,l1) by (hlog_alpha_i,llog_alpha_i) result (h5,l5)*/
    Add222(hlog_alpha_i,llog_alpha_i,h1,l1,&h5,&l5);
   
    double h_log2 = 0x1.62e42fefa38p-1;
    double l_log2 = 0x1.ef35793c7673p-45;
    double e_hlog2,e_llog2;
    /* e*log(2) in (double,double) precision.*/
    Mul122(e, h_log2, l_log2, &e_hlog2, &e_llog2);
    
    double h6,l6;

    /*adding (e_hlog2,e_llog2) by (h5,l5) result (h6,l6)*/
    Add222(e_hlog2, e_llog2, h5, l5, &h6,&l6);
    /*fprintf(fic, ",(%la,%la) \n",h6,l6);*/
    return (h6+l6);


}

static double cr_log_accurate_advanced_path(double x/*,FILE *fich,FILE *ficm,FILE *ficl,FILE *ficml*/){
    int s;
    int e;
    uint64_t m;
    extract (x,&s,&e, &m);
    

    /*The coefficients in triple-double numbers of the approximation polynomial of degree 14. */
    static const double U[14][3]= {
        {-0x1.1d35b84c5cc56p-4, 0x1.780337dcf5705p-58, -0x1.b91450d691e2p-114}, /* u14 decomposed into triple-double*/
        {0x1.3afccbd40d20fp-4, 0x1.9c4dd4dab116p-58, 0x1.3a1e112af7efp-112}, /* u13 decomposed into triple-double*/
        {-0x1.55552b4de4c76p-4, -0x1.50aa78415c4bfp-59, 0x1.2fafbf95a1a8p-115}, /* u12 decomposed into triple-double*/
        {0x1.745d171352e2dp-4, -0x1.b9a4b0517c2c6p-58, -0x1.6128b0932efc8p-112}, /* u11 decomposed into triple-double*/
        {-0x1.999999996fed4p-4, 0x1.8ec6a7a38be73p-58, -0x1.85507f01f5868p-112}, /* u10 decomposed into triple-double*/
        {0x1.c71c71c71c59ap-4, 0x1.1e1945381ebebp-58, -0x1.b30f3065fce38p-112}, /* u9 decomposed into triple-double*/
        {-0x1.fffffffffffffp-4, -0x1.89d824c2891b7p-58, 0x1.b4f656bbfa58p-115}, /* u8 decomposed into triple-double*/
        {0x1.2492492492492p-3, 0x1.2438b8708dbd8p-57, 0x1.e00436fa1c45p-113}, /* u7 decomposed into triple-double*/
        {-0x1.5555555555555p-3, -0x1.5555440ddf7efp-57, 0x1.0fce30dd1d9p-116}, /* u6 decomposed into triple-double*/
        {0x1.999999999999ap-3, -0x1.9999999bc2648p-57, 0x1.4e3f5d1a3a7p-117}, /* u5 decomposed into triple-double*/
        {-0x1p-2, 0x1.4b03dac3ad4d9p-100, 0x1.78p-156}, /* u4 decomposed into triple-double*/
        {0x1.5555555555555p-2, 0x1.5555555555555p-56, 0x1.219f0c0c44756p-110}, /* u3 decomposed into triple-double*/
        {-0x1p-1, 0x1.9745a799cp-127, 0x0p+0}, /* u2 decomposed into triple-double*/
        {0x1p+0, -0x1.09f14p-143, 0x0p+0}, /* u1 decomposed into triple-double*/
    };

   if ((x> 0x1.fc65aa1908a66p-1) && (x<1)) {
        double hr,lr ;
        /*Add x by -1 result hr,lr*/
        Add112(-1.0,x,&hr,&lr);
       /*printf("hr = %la lr = %la\n",hr,lr);*/

       double hf = U[0][0];
        double mf = U[0][1];
        double lf = U[0][2];
        
        double h,m,l;
        /*Multiply (hf,mf,lf) by r result (h,m,l)*/
        Mul233(hr,lr,hf,mf,lf,&h,&m,&l);
        /*printf("i = 0, h = %la, m = %la, l = %la \n",h,m,l);*/
        double hfi,mfi,lfi,h1,m1,l1;

        for (int i = 1;i<14;i++){
            hfi = U[i][0];
            mfi = U[i][1];
            lfi = U[i][2];
            /*Add (hfi,mfi,lfi) by (h,m,l) result (h1,m1,l1)*/
            Add333(hfi,mfi,lfi,h,m,l,&h1,&m1,&l1);

            /*multiply (hr,lr) by (h1,m1,l1) result (h,m,l)*/
            Mul233(hr,lr,h1,m1,l1,&h,&m,&l);
            /*printf("i = %d, h1 = %la, m1 = %la, l1 = %la, h = %la, m = %la, l = %la \n",i,h1,m1,l1,h,m,l);*/
       }
        /*printf("h = %la, m = %la, l = %la\n", h, m, l);*/
        return h+m;
    }
    
    if ((x> 1) && (x <= 0x1.0178e5916f543p+0)) {
        double hr,lr ;
        /*Add x by -1 result hr,lr*/
     Add112(x,-1.0,&hr,&lr);
        /*printf("hr = %la lr = %la\n",hr,lr);*/
        
        double hf = U[0][0];
        double mf = U[0][1];
        double lf = U[0][2];
        
        double h,m,l;
        /*Multiply (hf,mf,lf) by r result (h,m,l)*/
        Mul233(hr,lr,hf,mf,lf,&h,&m,&l);
        /*printf("i = 0, h = %la, m = %la, l = %la \n",h,m,l);*/
        double hfi,mfi,lfi,h1,m1,l1;

        for (int i = 1;i<14;i++){
            hfi = U[i][0];
            mfi = U[i][1];
            lfi = U[i][2];
            /*Add (hfi,mfi,lfi) by (h,m,l) result (h1,m1,l1)*/
            Add333(hfi,mfi,lfi,h,m,l,&h1,&m1,&l1);

            /*multiply (hr,lr) by (h1,m1,l1) result (h,m,l)*/
            Mul233(hr,lr,h1,m1,l1,&h,&m,&l);
            /*printf("i = %d, h1 = %la, m1 = %la, l1 = %la, h = %la, m = %la, l = %la \n",i,h1,m1,l1,h,m,l);*/
        }
        /*printf("h = %la, m = %la, l = %la\n", h, m, l);*/
         return h+m;
    }
    double m1;
    /*If x is a subnormal */
    if ((s ==0) && (e == 0) && (m !=0)) {
        uint64_t v = m;
        e = e - 1023;
        v = v*2;
        while (v < 0x10000000000000) {
            v *= 2;
            e--;
        }
        
        m1 =v *0x1p-52;
       
        u u;
        u.x = m1;
        m = u.i & 0xFFFFFFFFFFFFF;
    }
    /*If x is normal */
    else {
        m1= 1 + m*0x1p-52;
       
        e = e - 1023;
    }

    /*We shift by 44 bits for to get the first 8 bits.*/
    uint64_t i = m>>44;
    
    
    /*k=8
    halpha_i,lalpha_i is the double,double of 2^k/(2^k+i).*/
    double halpha_i = table_alpha_i_triple[(int)i][0];
    double malpha_i = table_alpha_i_triple[(int)i][1];
    double lalpha_i = table_alpha_i_triple[(int)i][2];

      double hr,mr,lr;
    /*Multiply (halpha_i, malpha_i,lalpha_i_m) by m1 result (hr, mr, lr) */
    Mul133(m1,halpha_i, malpha_i, lalpha_i,&hr,&mr,&lr);
    /*printf(" m1 = %la  halpha_i = %la , malpha_i = %la, lalpha_i = %la ,hr = %la, mr = %la, lr = %la\n",m1,halpha_i, malpha_i, lalpha_i,hr,mr,lr);*/
    double h,mprime,l;
    /*Add (hr,lr) by (-1) result (h,m,l)*/
    Add133(-1,hr,mr,lr,&h,&mprime,&l);
    /*printf("h = %la, mprime = %la, l = %la\n", h,mprime,l);*/
    double hf = U[0][0];
    double mf = U[0][1];
    double lf = U[0][2];
    double h1, m1prime, l1;
    /*Multiply (hf,mf,lf) by (h,m,l) result (h1,m1,l1)*/
    Mul333(h,mprime,l,hf,mf,lf,&h1,&m1prime,&l1);
    /*printf("h1 = %la, m1prime = %la, l1 = %la\n", h1, m1prime,l1);*/
    double hfi,mfi,lfi,h2,m2,l2;

    for (int i = 1;i<14;i++){
        hfi = U[i][0];
        mfi = U[i][1];
        lfi = U[i][2];
        /*Add (hfi,mfi,lfi) by (h1,m1,l1) result (h2,m2,l2)*/
        Add333(hfi,mfi,lfi,h1,m1prime,l1,&h2,&m2,&l2);

        /*multiply (h,m,l) by (h2,m2,l2) result (h1,m1,l1)*/
        Mul333(h,mprime,l,h2,m2,l2,&h1,&m1prime,&l1);
       /*printf("i = %d, h1 = %la, m1prime = %la, l1 = %la, h2 = %la, m2 = %lal2 = %la\n", i,h1,m1prime,l1,h2,m2,l2);*/
    }
    /*k=8
    hlog_alpha_i,llog_alpha_i is the double,double of log(2^k/(2^k+i)).*/
    double hlog_alpha_i = table_log_alpha_i_triple[(int)i][0];
    double mlog_alpha_i = table_log_alpha_i_triple[(int)i][1];
    double llog_alpha_i = table_log_alpha_i_triple[(int)i][2];


    double h5,m5,l5;
    /*Add (h1,m1,l1) by (hlog_alpha_i,mlog_alpha_i,llog_alpha_i) result (h5,m5,l5)*/
    Add333(hlog_alpha_i,mlog_alpha_i,llog_alpha_i,h1,m1prime,l1,&h5,&m5,&l5);
    /*printf("h5 = %la ,m5 = %la, l5 = %la , hlog_alpha_i = %la, mlog_alpha_i = %la, llog_alpha_i = %la\n",h5,m5,l5,hlog_alpha_i,mlog_alpha_i, llog_alpha_i);*/
    double h_log2 =0x1.62e42fefa39efp-1;
    double m_log2 =0x1.abc9e3b39803fp-56;
    double l_log2 = 0x1.7b57a079a1934p-111;
    double e_hlog2,e_mlog2,e_llog2;
    /* e*log(2) in triple-double precision.*/
    Mul133(e, h_log2, m_log2, l_log2, &e_hlog2, &e_mlog2, &e_llog2);
    
    double h6,m6,l6;

    /*adding (e_hlog2,e_mlog2,e_llog2) by (h5,m5,l5) result (h6,m6,l6)*/
    Add333(e_hlog2, e_mlog2, e_llog2, h5, m5, l5, &h6, &m6, &l6);
    /*printf("h = %la, m = %la, l = %la\n", h6, m6, l6);*/
    /*printf("h6 = %la, m6 = %la, l6 =%la\n",h6,m6,l6);*/
    
    double ml;
   ml = m6+l6;
    
    
    
    double m62 = 2*m6;
   if ( (fabs(m6) == ulp_0_5(h6)) && m6>0){
        if (l6>0){
            ml = nextafter(m6,m62);
        }
        if (l6<0){
             ml= nextafter(m6,0);
        }
    }
    if ( (fabs(m6) == ulp_0_5(h6)) && (m6<0)) {
        if (l6<0){
            ml = nextafter(m6,m62);
        }
        if (l6>0){
            ml= nextafter(m6,0);
        }
    }
    
   
   /* fprintf(fich,"%la\n",h6);
    fprintf(ficm,"%la\n",m6);
    fprintf(ficl,"%la\n",l6);
    fprintf(ficml,"%la\n", ml);*/



    return h6+ml;

   /* fprintf(fich,"%la\n",h6);
    fprintf(ficm,"%la\n",m6);
    fprintf(ficl,"%la\n",l6);
    fprintf(ficml,"%la\n", ml);*/



    


}
int main(){

    
  FILE *fic1 = fopen("log.wc", "r"); 
    if (fic1 == NULL)
        exit(1);
   FILE *fic2 = fopen("../sage/logDN.wc", "w");
    if (fic2 == NULL)
        exit(1);
    FILE *fic3 = fopen("../sage/logDU.wc", "w");
    if (fic3 == NULL)
        exit(1);
   FILE *fic4 = fopen("../sage/logDD.wc", "w");
    if (fic4 == NULL)
        exit(1);
    FILE *fic5 = fopen("../sage/logDZ.wc", "w");
    if (fic5 == NULL)
        exit(1);
    
    
    printf("Worst cases test with MPFR with RNDN\n");
    
    int j = 0;
    int n = 0;
    mpfr_t x2,x3;
    mpfr_init2(x2,53);
    mpfr_init2(x3,53);
    double x,x1;
    
   while (!feof(fic1)){
        n+=1;
        fscanf(fic1,"%la",&x);
        
        x1 = cr_log_accurate_advanced_path(x);
        mpfr_set_d(x2,x,MPFR_RNDN);
        mpfr_log (x3,x2,MPFR_RNDN);
        /*printf("x = %la , cr_log(x) = %la, mpfr_log = %la \n",x,x1,mpfr_get_d(x3,MPFR_RNDN));*/

       if (x1!= mpfr_get_d(x3,MPFR_RNDN)){
            j=j+1;
            fprintf(fic2,"%la\n",  x);
            
        }
        
    }
    
    printf("%d errors on %d cases \n", j, n-1);
    
    
    printf("Worst cases test with MPFR with RNDU\n");
    fesetround (FE_UPWARD);
    j = 0;
    n = 0;
    rewind(fic1);
    while (!feof(fic1)){
        n+=1;
        fscanf(fic1,"%la",&x);

        x1 = cr_log_accurate_advanced_path(x);
        mpfr_set_d(x2,x,MPFR_RNDU);
        
        mpfr_log (x3,x2,MPFR_RNDU);
       
        /*printf("x = %la , cr_log(x) = %la, mpfr_log = %la \n",x,x1,mpfr_get_d(x3,MPFR_RNDU));*/
       if (x1!= mpfr_get_d(x3,MPFR_RNDU)){
            j=j+1;
            fprintf(fic3," %la \n", x);
        }
        
    }
    
   printf("%d errors on %d cases \n", j, n-1);
    
   printf("Worst cases test with MPFR with RNDD\n");
   fesetround (FE_DOWNWARD);
    j = 0;
    n = 0;
    rewind(fic1);
    while (!feof(fic1)){
        n+=1;
        fscanf(fic1,"%la",&x);
        
        x1 = cr_log_accurate_advanced_path(x);
        mpfr_set_d(x2,x,MPFR_RNDD);
        mpfr_log (x3,x2,MPFR_RNDD);
        /*printf("x = %la , cr_log(x) = %la, mpfr_log = %la \n",x,x1,mpfr_get_d(x3,MPFR_RNDN));*/
       if (x1!= mpfr_get_d(x3,MPFR_RNDD)){
            j=j+1;
            fprintf(fic4," %la \n", x);
        }
        
    }
    
   printf("%d errors on %d cases \n", j, n-1);
    
    printf("Worst cases test with MPFR with RNDZ\n");
    fesetround (FE_TOWARDZERO);
    j = 0;
    n = 0;
    rewind(fic1);
    while (!feof(fic1)){
        n+=1;
        fscanf(fic1,"%la",&x);
        
        x1 = cr_log_accurate_advanced_path(x);
        mpfr_set_d(x2,x,MPFR_RNDZ);
        mpfr_log (x3,x2,MPFR_RNDZ);
        /*printf("x = %la , cr_log(x) = %la, mpfr_log = %la \n",x,x1,mpfr_get_d(x3,MPFR_RNDN));*/
      if (x1!= mpfr_get_d(x3,MPFR_RNDZ)){
            j=j+1;
            fprintf(fic5," %la \n", x);
        }
        
    }
    printf("%d errors on %d cases \n", j, n-1);

    mpfr_clear(x2);
    mpfr_clear(x3);
    fclose(fic1);
    fclose(fic2);
    fclose(fic3);
    fclose(fic4);
    fclose(fic5);

/*FILE *fic1 = fopen("../sage/logDN.wc", "r");
    if (fic1 == NULL)
        exit(1);
FILE *fic1ml = fopen("../sage/mllogDN.wc", "w");
if (fic1ml == NULL)
        exit(1);
    FILE *fic1h = fopen("../sage/hlogDN.wc","w");
    if (fic1h == NULL)
         exit(1);
    FILE *fic1m = fopen("../sage/mlogDN.wc","w");
    if (fic1m == NULL)
         exit(1);
    FILE *fic1l = fopen("../sage/llogDN.wc","w");
    if (fic1l == NULL)
        exit(1);
   /* FILE *fic2 = fopen("../sage/logDU.wc", "r");
    if (fic2 == NULL)
        exit(1);
    FILE *fic2h = fopen("../sage/hlogDU.wc","w");
    if (fic2h == NULL)
        exit(1);
    FILE *fic2m = fopen("../sage/mlogDU.wc","w");
    if (fic2m == NULL)
        exit(1);
    FILE *fic2l = fopen("../sage/llogDU.wc","w");
    if (fic2l == NULL)
        exit(1);
    FILE *fic3 = fopen("../sage/logDD.wc", "r");
    if (fic3 == NULL)
        exit(1);
    FILE *fic3h = fopen("../sage/hlogDD.wc","w");
    if (fic3h == NULL)
        exit(1);
    FILE *fic3m = fopen("../sage/mlogDD.wc","w");
    if (fic3m == NULL)
        exit(1);
    FILE *fic3l = fopen("../sage/llogDD.wc","w");
    if (fic3l == NULL)
        exit(1);
    FILE *fic4 = fopen("../sage/logDZ.wc", "r");
    if (fic4 == NULL)
        exit(1);
    FILE *fic4h = fopen("../sage/hlogDZ.wc","w");
    if (fic4h == NULL)
        exit(1);
    FILE *fic4m = fopen("../sage/mlogDZ.wc","w");
    if (fic4m == NULL)
        exit(1);
    FILE *fic4l = fopen("../sage/llogDZ.wc","w");
    if (fic4l == NULL)
        exit(1);*/
    
    /*printf("Worst cases test with MPFR with RNDN\n");
    
    double x,x1;
    int n =0;
    while (!feof(fic1)){
        
        fscanf(fic1,"%la",&x);
        
        x1 = cr_log_accurate_advanced_path(x,fic1h,fic1m,fic1l,fic1ml);
        
        n+=1;
    }
    printf("n = %d\n",n-1);
    
    
   /* printf("Worst cases test with MPFR with RNDU\n");
    n =0;
    while (!feof(fic2)){
        
        fscanf(fic2,"%la",&x);
        
        x1 = cr_log_accurate_advanced_path(x,fic2h,fic2m,fic2l);
       
        n+=1;
    }
    printf("n = %d\n",n-1);
   
    
    printf("Worst cases test with MPFR with RNDD\n");
    n =0;
    while (!feof(fic3)){
        
        fscanf(fic3,"%la",&x);
        
        x1 = cr_log_accurate_advanced_path(x,fic3h,fic3m,fic3l);
        
        n+=1;
    }
    printf("n = %d\n",n-1);
    
    
    printf("Worst cases test with MPFR with RNDZ\n");
    n =0;
    while (!feof(fic4)){
        
        fscanf(fic4,"%la",&x);
        
        x1 = cr_log_accurate_advanced_path(x,fic4h,fic4m,fic4l);
        
        n+=1;
        
    }
    printf("n = %d\n",n-1);*/
   

   
    /*fclose(fic1);
    fclose(fic1h);
    fclose(fic1m);
    fclose(fic1l);
    fclose(fic1ml);*/

    /*fclose(fic2);
    fclose(fic2h);
    fclose(fic2m);
    fclose(fic2l);*/

    /*fclose(fic3);
    fclose(fic3h);
    fclose(fic3m);
    fclose(fic3l);*/

    /*fclose(fic4);
    fclose(fic4h);
    fclose(fic4m);
    fclose(fic4l);*/
    
    /*double x,x1;
    x= 0x1.ff78281f6daf5p-1;
    x1 = cr_log_accurate_advanced_path(x);
    printf("%la \n",x1);*/
    return 0;

   
}

/* gcc cr_log.c -lmpfr
   ./a.out*/
