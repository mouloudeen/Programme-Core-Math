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
static double cr_log_accurate_path(double x, double *h6, double *l6/*, FILE *fic*/){
    int s;
    int e;
    uint64_t m;
    extract (x,&s,&e, &m);
/* Special case*/
    if (x == 0x1.fb85251a3f26fp-1){
        *h6 = -0x1.1ff9b8e8b38bep-7;
        *l6 = -0x1.0b393919c1fa3p-109;
        return (*h6+*l6);
    }

    if (x == 0x1.fc65aa1908a66p-1){
        *h6 = -0x1.cecc4ad8d358bp-8;
        *l6 = -0x1.65a43e3cf2b61p-107;
        return (*h6+*l6);
    }

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

    
    if ((x> 0x1.fe73451b9c74fp-1) && (x<1)) {
        double hr,lr ;
        /*Add x by -1 result r*/
        Add112(-1,x,&hr,&lr);
       
        double hf = G[0][0];
        double lf = G[0][1];
        
        double h,l;
        /*Multiply (hf,lf) by r result (h,l)*/
        Mul222(hr,lr,hf,lf,&h,&l);
        
        double hfi,lfi,h1,l1;

        for (int i = 1;i<11;i++){
            hfi = G[i][0];
            lfi = G[i][1];

            /*Add (hfi,lfi) by (h,l) result (h1,l1)*/
            Add222(hfi,lfi,h,l,&h1,&l1);

            /*multiply (hr,lr) by (h1,l1) result (h,l)*/
            Mul222(hr,lr,h1,l1,&h,&l);
            
        }
        *h6 = h;
        *l6 = l;
        /*fprintf(fic,",(%la,%la) \n",h,l);*/
        return (h+l);
    }

    if ((x> 1) && (x <= 0x1.0178e5916f543p+0)){
        double hr,lr ;
        /*Add x by -1 result r*/
        Add112(x,-1,&hr,&lr);
       
        double hf = G[0][0];
        double lf = G[0][1];
        
        double h,l;
        /*Multiply (hf,lf) by r result (h,l)*/
        Mul222(hr,lr,hf,lf,&h,&l);
        
        double hfi,lfi,h1,l1;

        for (int i = 1;i<11;i++){
            hfi = G[i][0];
            lfi = G[i][1];

            /*Add (hfi,lfi) by (h,l) result (h1,l1)*/
            Add222(hfi,lfi,h,l,&h1,&l1);

            /*multiply (hr,lr) by (h1,l1) result (h,l)*/
            Mul222(hr,lr,h1,l1,&h,&l);
            
        }
        *h6 = h;
        *l6 = l;
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
    
    double hf = G[0][0];
    double lf = G[0][1];

    double h1,l1;
    /*Multiply (hf,lf) by (h,l) result (h1,l1)*/
    Mul222(hf,lf,h,l,&h1,&l1);
    
    double hfi,lfi,h2,l2;

    for (int i = 1;i<11;i++){
        hfi = G[i][0];
        lfi = G[i][1];

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
   
    double h_log2 = 0x1.62e42fefa39efp-1;
    double l_log2 = 0x1.abc9e3b39803fp-56;
    double e_hlog2,e_llog2;
    /* e*log(2) in (double,double) precision.*/
    Mul122(e, h_log2, l_log2, &e_hlog2, &e_llog2);
    
    

    /*adding (e_hlog2,e_llog2) by (h5,l5) result (h6,l6)*/
    Add222(e_hlog2, e_llog2, h5, l5, h6,l6);
    /*fprintf(fic, ",(%la,%la) \n",h6,l6);*/
    return (*h6+*l6);


}


int main(){

    
  
   
    FILE *fic2 = fopen("../sage/logDNfast.wc", "r");
    if (fic2 == NULL)
        exit(1);
    FILE *fic3 = fopen("../sage/logDUfast.wc", "r");
    if (fic3 == NULL)
        exit(1);
   FILE *fic4 = fopen("../sage/logDDfast.wc", "r");
    if (fic4 == NULL)
        exit(1);
    FILE *fic5 = fopen("../sage/logDZfast.wc", "r");
    if (fic5 == NULL)
        exit(1);

    FILE *fic6 = fopen("../sage/logDNaccurate.wc", "w");
    if (fic6 == NULL)
        exit(1);
    FILE *fic7 = fopen("../sage/logDUaccurate.wc", "w");
    if (fic7 == NULL)
        exit(1);
   FILE *fic8 = fopen("../sage/logDDaccurate.wc", "w");
    if (fic8 == NULL)
        exit(1);
    FILE *fic9 = fopen("../sage/logDZaccurate.wc", "w");
    if (fic9 == NULL)
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
    double err = 0x1.810c9a86fc45ep-99;
    
    double right,left;
    double h6,l6;
   while (!feof(fic2)){
        n+=1;
        fscanf(fic2,"%la",&x);
        
         y = cr_log_accurate_path(x,&h6,&l6);
         mpfr_set_d(x2,x,MPFR_RNDN);
            mpfr_log (x3,x2,MPFR_RNDN);
        left  = h6 + __builtin_fma (-err, h6, l6);
        right = h6 + __builtin_fma (+err, h6, l6);
       if (right != left){
            jf=jf+1;
            fprintf(fic6,"%la\n",  x);
            
        }
         if (right== left){
            jt=jt+1;
            
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
   
    
   
   while (!feof(fic3)){
        n+=1;
        fscanf(fic3,"%la",&x);
        
        y = cr_log_accurate_path(x,&h6,&l6);
         mpfr_set_d(x2,x,MPFR_RNDU);
            mpfr_log (x3,x2,MPFR_RNDU);
        left  = h6 + __builtin_fma (-err, h6, l6);
        right = h6 + __builtin_fma (+err, h6, l6);
       if (right != left){
            jf=jf+1;
            fprintf(fic7,"%la\n",  x);
            
        }
         if (right== left){
            jt=jt+1;
           
            if (y == mpfr_get_d(x3,MPFR_RNDU)){
            jtt=jtt+1;

            } 
            else{
                printf("x = %la\n",x);
            }  
        }
        
    }
    
    printf("%d errors and %d trues and %d equal with mpfr_log(x) on %d  cases \n", jf,jt,jtt, n);  
   
   printf("Worst cases test with MPFR with RNDD\n");
   fesetround (FE_DOWNWARD);  
     jt = 0;
     jtt =0;
     jf = 0;
     n = 0;
   
    
    
   while (!feof(fic4)){
        n+=1;
        fscanf(fic4,"%la",&x);
        
         y = cr_log_accurate_path(x,&h6,&l6);
          mpfr_set_d(x2,x,MPFR_RNDD);
            mpfr_log (x3,x2,MPFR_RNDD);
        left  = h6 + __builtin_fma (-err, h6, l6);
        right = h6 + __builtin_fma (+err, h6, l6);
       if (right != left){
            jf=jf+1;
            fprintf(fic8,"%la\n",  x);
            
        }
         if (right== left){
            jt=jt+1;
           
            if (y == mpfr_get_d(x3,MPFR_RNDD)){
            jtt=jtt+1;

            } 
            else{
                printf("x = %la\n",x);
            }  
        }
        
    }
    
    printf("%d errors and %d trues and %d equal with mpfr_log(x) on %d  cases \n", jf,jt,jtt, n);   
    
    printf("Worst cases test with MPFR with RNDZ\n");
   fesetround (FE_TOWARDZERO); 
     jt = 0;
     jtt =0;
     jf = 0;
     n = 0;
   
    
   
   while (!feof(fic5)){
        n+=1;
        fscanf(fic5,"%la",&x);
        
         y = cr_log_accurate_path(x,&h6,&l6);
          mpfr_set_d(x2,x,MPFR_RNDZ);
            mpfr_log (x3,x2,MPFR_RNDZ);
        left  = h6 + __builtin_fma (-err, h6, l6);
        right = h6 + __builtin_fma (+err, h6, l6);
       if (right != left){
            jf=jf+1;
            fprintf(fic9,"%la\n",  x);
            
        }
         if (right== left){
            jt=jt+1;
           
            if (y == mpfr_get_d(x3,MPFR_RNDZ)){
            jtt=jtt+1;

            } 
            else{
                printf("x = %la\n",x);
            }  
        }
        
    }
    
    printf("%d errors and %d trues and %d equal with mpfr_log(x) on %d  cases \n", jf,jt,jtt, n);   
    
    mpfr_clear(x2);
    mpfr_clear(x3);

    
    fclose(fic2);
    fclose(fic3);
    fclose(fic4);
    fclose(fic5);
    fclose(fic6);
    fclose(fic7);
    fclose(fic8);
    fclose(fic9);


    return 0;

   
}



   /*gcc testaccurate.c -lmpfr
   ./a.out*/