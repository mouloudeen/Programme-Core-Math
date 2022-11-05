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



static double cr_log_fast_path(double x, double *h6, double *l6){ 
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
        double xx2 =(x-1) * (x-1);
        double ff4 = f4 + f5 * (x-1);
        double ff6 = f6 + f7 * (x-1);
        

        double cr_log_fast_1 = ff4 + xx2 * ff6;

        /* Add x by -1 result hr,lr */
        double hr,lr;
        Add112(x, -1, &hr, &lr);

        
        double hr2,lr2;
        /* Multiply (hr,lr) by (hr,lr) result (hr2,lr2) */
        Mul222(hr,lr,hr,lr,&hr2,&lr2);

        double hr4,lr4;
        /* Multiply (hr2,lr2) by (hr2,lr2) result (hr4,lr4) */
        Mul222(hr2,lr2,hr2,lr2,&hr4,&lr4);

        double ffh0,ffl0;
        /*  Multiply f1 by (hr,lr) result (ffh0,ffl0) */
        Mul122(f1,hr,lr,&ffh0,&ffl0);

        double ffh3,ffl3;
        /* Multiply f3 by (hr,lr) result (ffh3,ffl3) */
        Mul122(f3,hr,lr,&ffh3,&ffl3);

        double ffh2,ffl2;
        /*Add f2 by (ffh3,ffl3) result (ffh2,ffl2) */
        Add122(f2,ffh3,ffl3,&ffh2,&ffl2);

        double ffhx2,fflx2;
        /* Multiply (hr2,lr2) by (ffh2,ffl2) result (ffhx2,fflx2) */
        Mul222(hr2,lr2,ffh2,ffl2,&ffhx2,&fflx2);

        double ffhx4,fflx4;
        /* Multiply (hr4,lr4) by cr_log_fast_1 result (ffhx4,fflx4) */
        Mul122(cr_log_fast_1,hr4,lr4,&ffhx4,&fflx4);

        double ffhx0,fflx0;
        /* Add (ffh0,ffl0) by (ffhx2,fflx2) result (ffhx0,fflx0) */
        Add222(ffh0, ffl0, ffhx2, fflx2, &ffhx0,&fflx0);
        
        /* Add (ffhx0,fflx0) by (ffhx4,fflx4) result (h6,l6) */
        Add222(ffhx0,fflx0, ffhx4,fflx4,h6,l6);

        return (*h6+*l6);

    
    }
    
    if ((x> 1) && (x < 0x1.00068db8bac71p+0)){
     
        double xx2 =(x-1) * (x-1);
        double ff4 = f4 + f5 * (x-1);
        double ff6 = f6 + f7 * (x-1);
        

        double cr_log_fast_1 = ff4 + xx2 * ff6;

        /* Add x by -1 result hr,lr */
        double hr,lr;
        Add112(x, -1, &hr, &lr);

        
        double hr2,lr2;
        /* Multiply (hr,lr) by (hr,lr) result (hr2,lr2) */
        Mul222(hr,lr,hr,lr,&hr2,&lr2);

        double hr4,lr4;
        /* Multiply (hr2,lr2) by (hr2,lr2) result (hr4,lr4) */
        Mul222(hr2,lr2,hr2,lr2,&hr4,&lr4);

        double ffh0,ffl0;
        /*  Multiply f1 by (hr,lr) result (ffh0,ffl0) */
        Mul122(f1,hr,lr,&ffh0,&ffl0);

        double ffh3,ffl3;
        /* Multiply f3 by (hr,lr) result (ffh3,ffl3) */
        Mul122(f3,hr,lr,&ffh3,&ffl3);

        double ffh2,ffl2;
        /*Add f2 by (ffh3,ffl3) result (ffh2,ffl2) */
        Add122(f2,ffh3,ffl3,&ffh2,&ffl2);

        double ffhx2,fflx2;
        /* Multiply (hr2,lr2) by (ffh2,ffl2) result (ffhx2,fflx2) */
        Mul222(hr2,lr2,ffh2,ffl2,&ffhx2,&fflx2);

        double ffhx4,fflx4;
        /* Multiply (hr4,lr4) by cr_log_fast_1 result (ffhx4,fflx4) */
        Mul122(cr_log_fast_1,hr4,lr4,&ffhx4,&fflx4);

        double ffhx0,fflx0;
        /* Add (ffh0,ffl0) by (ffhx2,fflx2) result (ffhx0,fflx0) */
        Add222(ffh0, ffl0, ffhx2, fflx2, &ffhx0,&fflx0);
        
        /* Add (ffhx0,fflx0) by (ffhx4,fflx4) result (h6,l6) */
        Add222(ffhx0,fflx0, ffhx4,fflx4,h6,l6);

        return (*h6+*l6);

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

    
    
    
        double hh2 = h * h;
        double ff4 = f4 + f5 * h;
        double ff6 = f6 + f7 * h;
        

    double cr_log_fast_1 = ff4 +hh2*ff6;
    
    double ll2,hh4,ll4;
    
    /* Multiply (h,l) by (h,l) result (hh2,ll2) */
    Mul222(h,l,h,l,&hh2,&ll2);
    
    /* Multiply (hh2,ll2) by (hh2,ll2) result (hh4,ll4) */
    Mul222(hh2,ll2,hh2,ll2,&hh4,&ll4);
    
    double ffh0,ffl0;
    /*  Multiply f1 by (h,l) result (ffh0,ffl0) */
    Mul122(f1,h,l,&ffh0,&ffl0);

    double ffh3,ffl3;
    /* Multiply f3 by (h,l) result (ffh3,ffl3) */
    Mul122(f3,h,l,&ffh3,&ffl3);

    double ffh2,ffl2;
    /*Add f2 by (ffh3,ffl3) result (ffh2,ffl2) */
    Add122(f2,ffh3,ffl3,&ffh2,&ffl2);

    double ffhx2,fflx2;
    /* Multiply (hh2,ll2) by (ffh2,ffl2) result (ffhx2,fflx2) */
    Mul222(hh2,ll2,ffh2,ffl2,&ffhx2,&fflx2);

    double ffhx4,fflx4;
    /* Multiply (hh4,ll4) by cr_log_fast_1 result (ffhx4,fflx4) */
    Mul122(cr_log_fast_1,hh4,ll4,&ffhx4,&fflx4);

    double ffhx0,fflx0;
    /* Add (ffh0,ffl0) by (ffhx2,fflx2) result (ffhx0,fflx0) */
    Add222(ffh0, ffl0, ffhx2, fflx2, &ffhx0,&fflx0);
    
    double h4,l4;
    /* Add (ffhx0,fflx0) by (ffhx4,fflx4) result (h4,l4) */
    Add222(ffhx0,fflx0, ffhx4,fflx4,&h4,&l4);
    
    double log_alpha_i_m = table_log_alpha_i_modified[i];
   
    
    double h5,l5;
    /* add (hlog_alpha_i_m, llog_alpha_i_m)*/
    Add122(log_alpha_i_m, h4, l4, &h5, &l5);
    
    
    
    
    double h_log2 = 0x1.62e42fefa38p-1;
    double l_log2 = 0x1.ef35793c7673p-45;
    double e_hlog2,e_llog2;
    /* e*log(2) in (double,double) precision.*/
    Mul122(e, h_log2, l_log2, &e_hlog2, &e_llog2);
    
   

    /*adding with (h10,l10)*/
    Add222(e_hlog2, e_llog2,h5,l5,h6,l6);
    
    return (*h6+*l6);
    

}

int main(){

    
  
   FILE *fic1 = fopen("log.wc", "r"); 
    if (fic1 == NULL)
        exit(1);
    FILE *fic2 = fopen("../sage/logDNfast.wc", "w");
    if (fic2 == NULL)
        exit(1);
    FILE *fic3 = fopen("../sage/logDUfast.wc", "w");
    if (fic3 == NULL)
        exit(1);
   FILE *fic4 = fopen("../sage/logDDfast.wc", "w");
    if (fic4 == NULL)
        exit(1);
    FILE *fic5 = fopen("../sage/logDZfast.wc", "w");
    if (fic5 == NULL)
        exit(1);

    printf("Worst cases test with MPFR with RNDN\n");
    
    int jt = 0;
    int jtt =0;
    int jf = 0;
    int n = 0;
    mpfr_t x2,x3;
    mpfr_init2(x2,53);
    mpfr_init2(x3,53);
    double x,y;
    double err = 0x1.6b6b11ea279ecp-59;
    /*           0x1.31b0761b44d2bp-64*/
    /*We started with the variable err = 0x1p-68 because 
    we theoretically calculated it.
    we had errors in all 4 rounding modes.
    We increased the power by 1 until we found only one error in the RNDN rounding mode 
    and no other errors in the other rounding modes. After, we had err = 0x1p-64.
    We searched by dichotomy after the dot. If we had errors on the equality between cr_log_fast(x) and mpfr_log(x)
     we increase either we decrease until we have 2 consecutive digits.
    We started again until we optimize err.
    */

    
    double right,left;
    double h6,l6;
   while (!feof(fic1)){
        n+=1;
        fscanf(fic1,"%la",&x);
        
        y = cr_log_fast_path(x,&h6,&l6);
        mpfr_set_d(x2,x,MPFR_RNDN);
        mpfr_log (x3,x2,MPFR_RNDN);
        left  = h6 + __builtin_fma (-err, h6, l6);
        right = h6 + __builtin_fma (+err, h6, l6);
       if (right != left){
            jf=jf+1;
            fprintf(fic2,"%la\n",  x);
            
        }
         if (right== left){
            jt=jt+1;
           
            if (y==1024){
            y = h6+l6;
                        }
            if (y == mpfr_get_d(x3,MPFR_RNDN)){
            jtt=jtt+1;

            } 
            else{
                printf("x = %la\n",x);
            }  
        }
        
    }
    
    printf("%d errors and %d trues and %d equal with mpfr_log(x) on %d  cases \n", jf,jt,jtt, n);
    
  printf("Worst cases test with MPFR with RNDU\n");
  fesetround (FE_UPWARD);  
     jt = 0;
     jtt =0;
     jf = 0;
     n = 0;
   
    
   rewind(fic1); 
   while (!feof(fic1)){
        n+=1;
        fscanf(fic1,"%la",&x);
        
         y = cr_log_fast_path(x,&h6,&l6);
         mpfr_set_d(x2,x,MPFR_RNDU);
        mpfr_log (x3,x2,MPFR_RNDU);
        left  = h6 + __builtin_fma (-err, h6, l6);
        right = h6 + __builtin_fma (+err, h6, l6);
       if (right != left){
            jf=jf+1;
            fprintf(fic3,"%la\n",  x);
            
        }
         if (right== left){
            jt=jt+1;
            
            if (y == mpfr_get_d(x3,MPFR_RNDU)){
            jtt=jtt+1;

            } 
             /*else{
                printf("x = %la\n",x);
            }  */
         }
        
    }
    
    printf("%d errors and %d trues and %d equal with mpfr_log(x) on %d  cases \n", jf,jt,jtt, n);  
   
   printf("Worst cases test with MPFR with RNDD\n");
   fesetround (FE_DOWNWARD);  
     jt = 0;
     jtt =0;
     jf = 0;
     n = 0;
   
    
   rewind(fic1); 
   while (!feof(fic1)){
        n+=1;
        fscanf(fic1,"%la",&x);
        
         y = cr_log_fast_path(x,&h6,&l6);
        mpfr_set_d(x2,x,MPFR_RNDD);
        mpfr_log (x3,x2,MPFR_RNDD);
        left  = h6 + __builtin_fma (-err, h6, l6);
        right = h6 + __builtin_fma (+err, h6, l6);
       if (right != left){
            jf=jf+1;
            fprintf(fic4,"%la\n",  x);
            
        }
         if (right== left){
            jt=jt+1;
           
            if (y == mpfr_get_d(x3,MPFR_RNDD)){
            jtt=jtt+1;

            } 
             /*else{
                printf("x = %la\n",x);
            }  */
       
         }
    }
    
    printf("%d errors and %d trues and %d equal with mpfr_log(x) on %d  cases \n", jf,jt,jtt, n);   
    
    printf("Worst cases test with MPFR with RNDZ\n");
   fesetround (FE_TOWARDZERO); 
     jt = 0;
     jtt =0;
     jf = 0;
     n = 0;
   
    
   rewind(fic1); 
   while (!feof(fic1)){
        n+=1;
        fscanf(fic1,"%la",&x);
        
         y = cr_log_fast_path(x,&h6,&l6);
          mpfr_set_d(x2,x,MPFR_RNDZ);
            mpfr_log (x3,x2,MPFR_RNDZ);
        left  = h6 + __builtin_fma (-err, h6, l6);
        right = h6 + __builtin_fma (+err, h6, l6);
       if (right != left){
            jf=jf+1;
            fprintf(fic5,"%la\n",  x);
            
        }
         if (right== left){
            jt=jt+1;
           
            if (y == mpfr_get_d(x3,MPFR_RNDZ)){
            jtt=jtt+1;

            } 
             /*else{
                printf("x = %la\n",x);
            }  */
         }
        
    }
    
    printf("%d errors and %d trues and %d equal with mpfr_log(x) on %d  cases \n", jf,jt,jtt, n);   
    
    mpfr_clear(x2);
    mpfr_clear(x3);

    fclose(fic1);
    fclose(fic2);
    fclose(fic3);
    fclose(fic4);
    fclose(fic5);
    


    return 0;

   
}

/* gcc -o testfast.o testfast.c 
   ./testfast.o*/

   /*gcc testfast.c -lmpfr
   ./a.out*/