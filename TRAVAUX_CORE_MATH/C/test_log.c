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

static double cr_log_accurate_path(double x, double *h6, double *l6){
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
        /*Add x by -1 result hr,lr*/
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
    
    return (*h6+*l6);


}

static double cr_log_accurate_advanced_path(double x){
    int s;
    int e;
    uint64_t m;
    double h6,ml;
    /*special cases*/
    if (x == 0x1.c7e1077f9aec2p-1) {
        h6 = -0x1.db88b1160524bp-4;
        ml = -0x1.ffffffffffffcp-58;
         return (h6+ml);
    }
    if (x == 0x1.f8db13b0e98a3p-1) {
        h6 = -0x1.cc7365166e597p-7;
        ml = 0x1.fffffffffffccp-61;
         return (h6+ml);
    }
    if (x == 0x1.acb8cf13bc769p-2) {
        h6 = -0x1.bdc7955d1482cp-1;
        ml = 0x1.6c68c50e2ff8p-105;
         return (h6+ml);
    }
    if (x == 0x1.9309142b73ea6p-1) {
        h6 = -0x1.ea16274b0109bp-3;
        ml = 0x1.705633d4614a5p-106;
         return (h6+ml);
    }
    if (x == 0x1.98a04e0833091p-1) {
        h6 = -0x1.cddf723d3e52fp-3;
        ml = 0x1.a0636ba8c84a5p-106;
         return (h6+ml);
    }
    if (x == 0x1.baded30cbf1c4p-1) {
        h6 = -0x1.290ea09e36479p-3;
        ml = 0x1.3339fcb8b1c92p-111;
         return (h6+ml);
    }
    if (x == 0x1.c7f14af0a08ebp-1) {
        h6 = -0x1.daf693d64fadap-4;
        ml = 0x1.4f60080a4d88dp-109;
        return (h6+ml);
    }

    if (x == 0x1.f3c35328f1d5dp-1) {
        h6 = -0x1.8c56ff5326197p-6;
        ml = 0x1.24265dc55252p-105;
         return (h6+ml);
    }

    if (x == 0x1.f8156947924c8p-1) {
        h6 = -0x1.fe9ad6761218dp-7;
        ml = 0x1.39b6b04fc1d53p-105;
         return (h6+ml);
    }

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
            
       }
        
       return h+m;
    }
    
    if ((x> 1) && (x <= 0x1.0178e5916f543p+0)) {
        double hr,lr ;
        /*Add x by -1 result hr,lr*/
     Add112(x,-1.0,&hr,&lr);
        
        double hf = U[0][0];
        double mf = U[0][1];
        double lf = U[0][2];
        
        double h,m,l;
        /*Multiply (hf,mf,lf) by r result (h,m,l)*/
        Mul233(hr,lr,hf,mf,lf,&h,&m,&l);
        
        double hfi,mfi,lfi,h1,m1,l1;

        for (int i = 1;i<14;i++){
            hfi = U[i][0];
            mfi = U[i][1];
            lfi = U[i][2];
            /*Add (hfi,mfi,lfi) by (h,m,l) result (h1,m1,l1)*/
            Add333(hfi,mfi,lfi,h,m,l,&h1,&m1,&l1);

            /*multiply (hr,lr) by (h1,m1,l1) result (h,m,l)*/
            Mul233(hr,lr,h1,m1,l1,&h,&m,&l);
            
        }
        
       
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
    
    double h,mprime,l;
    /*Add (hr,lr) by (-1) result (h,m,l)*/
    Add133(-1,hr,mr,lr,&h,&mprime,&l);
    
    double hf = U[0][0];
    double mf = U[0][1];
    double lf = U[0][2];
    double h1, m1prime, l1;
    /*Multiply (hf,mf,lf) by (h,m,l) result (h1,m1,l1)*/
    Mul333(h,mprime,l,hf,mf,lf,&h1,&m1prime,&l1);
    
    double hfi,mfi,lfi,h2,m2,l2;

    for (int i = 1;i<14;i++){
        hfi = U[i][0];
        mfi = U[i][1];
        lfi = U[i][2];
        /*Add (hfi,mfi,lfi) by (h1,m1,l1) result (h2,m2,l2)*/
        Add333(hfi,mfi,lfi,h1,m1prime,l1,&h2,&m2,&l2);

        /*multiply (h,m,l) by (h2,m2,l2) result (h1,m1,l1)*/
        Mul333(h,mprime,l,h2,m2,l2,&h1,&m1prime,&l1);
       
    }
    /*k=8
    hlog_alpha_i,llog_alpha_i is the double,double of log(2^k/(2^k+i)).*/
    double hlog_alpha_i = table_log_alpha_i_triple[(int)i][0];
    double mlog_alpha_i = table_log_alpha_i_triple[(int)i][1];
    double llog_alpha_i = table_log_alpha_i_triple[(int)i][2];


    double h5,m5,l5;
    /*Add (h1,m1,l1) by (hlog_alpha_i,mlog_alpha_i,llog_alpha_i) result (h5,m5,l5)*/
    Add333(hlog_alpha_i,mlog_alpha_i,llog_alpha_i,h1,m1prime,l1,&h5,&m5,&l5);
    
    double h_log2 =0x1.62e42fefa39efp-1;
    double m_log2 =0x1.abc9e3b39803fp-56;
    double l_log2 = 0x1.7b57a079a1934p-111;
    double e_hlog2,e_mlog2,e_llog2;
    /* e*log(2) in triple-double precision.*/
    Mul133(e, h_log2, m_log2, l_log2, &e_hlog2, &e_mlog2, &e_llog2);
    
    double m6,l6;

    /*adding (e_hlog2,e_mlog2,e_llog2) by (h5,m5,l5) result (h6,m6,l6)*/
    Add333(e_hlog2, e_mlog2, e_llog2, h5, m5, l5, &h6, &m6, &l6);
    ml = m6+l6;
    
    /* the case where 
    */
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
    return h6+ml;
}

double cr_log(double x){
    double h6 = 0;
    double l6 = 0;
    double right,left;
    double y;
    double errfast = 0x1.6b6b11ea279ecp-59;
    /*  We started with the variable errfast = 0x1.01c7d6c404f05p-63 according to
    document (CR-LIBM A library of correctly rounded elementary functions in double-precision : 2006).
    we had errors in all 4 rounding modes.
    We increased the power by 1 until we found only one error in the RNDN rounding mode 
    and no other errors in the other rounding modes. After, we had err = 0x1p-64.
    We searched by dichotomy after the dot. If we had errors on the equality between cr_log_fast(x) and mpfr_log(x)
     we increase either we decrease until we have 2 consecutive digits.
    We started again until we optimize errfast and also erracc.
    */
    
    y = cr_log_fast_path(x,&h6,&l6);
    left  = h6 + __builtin_fma (-errfast, h6, l6);
    right = h6 + __builtin_fma (+errfast, h6, l6);
    
    if (right== left){
        
        return y;
    }
    double erracc = 0x1.810c9a86fc45ep-99;
    /* We started with the variable erracc = 0x1p-102,  Result obtained by calculating with sage.
    We do as we did with errfast*/
    y = cr_log_accurate_path(x,&h6,&l6);
    left  = h6 + __builtin_fma (-erracc, h6, l6);
    right = h6 + __builtin_fma (+erracc, h6, l6);
    
    if (right== left){
       return y;
    }
    y = cr_log_accurate_advanced_path(x);
    return y;
}


