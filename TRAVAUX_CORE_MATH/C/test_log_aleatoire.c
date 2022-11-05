#include<stdint.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <time.h>
#include <mpfr.h>
#include <fenv.h>
#include "tables.h"
#include "tables_g.h"
#include "wc.h"



typedef union {
    double x;
    uint64_t i;
}u ;  

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
    Add112Cond(t1,bm,&t2,&t3); 
    t4 = t3+bl;
    Add112Cond(t2,t4,rm,rl);
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
    Add112Cond(t3,t4,&t5,&t6);
    Add222Cond(t1,t2,t5,t6,r1,r2);
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

   if ((x> 0x1.fba5e353f7ceep-1) && (x<1)) {
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
    
    return h6+m6;


}
int main(){

   FILE *fic = fopen("resultat.txt", "w");
    if (fic == NULL)
        exit(1); 
  
    int n =1000000000;
    
    fprintf(fic,"Random test with MPFR with RNDN\n");
    int j = 0;
    srand(time(NULL));
    mpfr_t x2,x3;
    mpfr_init2(x2,53);
    mpfr_init2(x3,53);
    for (int i =0;i<n;i++){
        double x = 1 +(double)rand() / ((double)RAND_MAX/1000);
        double x1 = cr_log_accurate_advanced_path(x);
        mpfr_set_d(x2,x,MPFR_RNDN);
        mpfr_log (x3,x2,MPFR_RNDN);
        /*printf("x = %la , cr_log(x) = %la, mpfr_log = %la \n",x,x1,mpfr_get_d(x3,MPFR_RNDN));*/
        if (x1!= mpfr_get_d(x3,MPFR_RNDN)){
            j=j+1;
            printf("x = %la %d %e\n", x,j, (double) j / (double) (i+1));
        }
    }
    fprintf(fic,"%d errors on %d cases \n", j, n);

    
    
    fprintf(fic,"Random test with MPFR with RNDU\n");
    fesetround (FE_UPWARD);
    j=0;
    for (int i =0;i<n;i++){
        double x = 1 +(double)rand() / ((double)RAND_MAX/1000);
        double x1 = cr_log_accurate_advanced_path(x);
        mpfr_set_d(x2,x,MPFR_RNDU);
        mpfr_log (x3,x2,MPFR_RNDU);
        /*printf("x = %la , cr_log(x) = %la, mpfr_log = %la \n",x,x1,mpfr_get_d(x3,MPFR_RNDN));*/
        if (x1!= mpfr_get_d(x3,MPFR_RNDU)){
            j=j+1;
            printf("x = %la %d %e\n", x,j, (double) j / (double) (i+1));
        }
    }
    fprintf(fic,"%d errors on %d cases \n", j, n);


    
   fprintf(fic,"Random test with MPFR with RNDD\n");
   fesetround (FE_DOWNWARD);
    j = 0;
    for (int i =0;i<n;i++){
        double x = 1 +(double)rand() / ((double)RAND_MAX/1000);
        double x1 = cr_log_accurate_advanced_path(x);
        mpfr_set_d(x2,x,MPFR_RNDD);
        mpfr_log (x3,x2,MPFR_RNDD);
        /*printf("x = %la , cr_log(x) = %la, mpfr_log = %la \n",x,x1,mpfr_get_d(x3,MPFR_RNDN));*/
        if (x1!= mpfr_get_d(x3,MPFR_RNDD)){
            j=j+1;
            printf("x = %la %d %e\n", x,j, (double) j / (double) (i+1));
        }
    }
    fprintf(fic,"%d errors on %d cases \n", j, n);
  
    
    fprintf(fic,"Random test with MPFR with RNDZ\n");
    fesetround (FE_TOWARDZERO);
    j = 0;
   
   for (int i =0;i<n;i++){
        double x = 1 +(double)rand() / ((double)RAND_MAX/1000);
        double x1 = cr_log_accurate_advanced_path(x);
        mpfr_set_d(x2,x,MPFR_RNDZ);
        mpfr_log (x3,x2,MPFR_RNDZ);
        /*printf("x = %la , cr_log(x) = %la, mpfr_log = %la \n",x,x1,mpfr_get_d(x3,MPFR_RNDN));*/
        if (x1!= mpfr_get_d(x3,MPFR_RNDZ)){
            j=j+1;
            printf("x = %la %d %e\n", x,j, (double) j / (double) (i+1));
        }
    }
     fprintf(fic,"%d errors on %d cases \n", j, n);

    mpfr_clear(x2);
    mpfr_clear(x3);
   

    fclose(fic);
    return 0;

   
}

/* gcc test_log_aleatoire.c -lmpfr
   ./a.out*/
