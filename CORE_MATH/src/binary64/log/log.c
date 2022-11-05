/* Correctly-rounded logarithm of binary64 value.
Copyright (c) 2022 Sid Ali Zitouni Terki.
This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include<stdint.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include <mpfr.h>
#include <fenv.h>

typedef union {
    double x;
    uint64_t i;
}u ;  

/* Compute 1/2*ulp(x)*/
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



/*ah,am,al and bh,bm,bl are triple-double numbers*/
void Add333(double ah, double am, double al, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t1,t2,t3,t4,t5,t6,t7,t8;
    Add112(ah,bh,rh,&t1);
    Add112Cond(am,bm,&t2,&t3);
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
    t5 = t3+t4;
    t6 = t2+t5;
    Add112(t1,t6,r1,r2);
}


/* a is a double numbers
   bh,bm,bl is a triple-double numbers*/
void Mul133(double a, double bh, double bm, double bl, double *rh, double *rm, double *rl){
    double t2,t3,t4,t5,t7,t8,t9,t10;
    Mul112(a,bh,rh,&t2);
    Mul112(a,bm,&t3,&t4);
    t5 = a*bl;
    Add112(t2,t3,&t9,&t7);
    t8 = t4+t5;
    t10 = t7+t8;
    Add112(t9,t10,rm,rl);
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
    Add222(t4,t5,t2,t3,&t21,&t22);
    Add222(t21,t22,t19,t20,rm,rl);
    /*rh,rm,rl is a triple-double numbers*/
}

/* table_alpha_i 
 alpha_i = {2^k}/{2^k+i}  with 0<=i<2*k with k=8.
 alpha_i = {h,l} to calculate in {double, double} precision
 h is the main value and l is the error value.
*/

static const double table_alpha_i[256][2]={
    //{h, l} value for each alpha_i
    {0x1p+0, 0x0p+0}, /* i= 0 */
    {0x1.fe01fe01fe02p-1, -0x1.fe01fe01fe02p-57}, /* i= 1 */
    {0x1.fc07f01fc07fp-1, 0x1.fc07f01fc07fp-57}, /* i= 2 */
    {0x1.fa11caa01fa12p-1, -0x1.aaff02f71aaffp-56}, /* i= 3 */
    {0x1.f81f81f81f82p-1, -0x1.f81f81f81f82p-55}, /* i= 4 */
    {0x1.f6310aca0dbb5p-1, 0x1.d2e19807d8c43p-55}, /* i= 5 */
    {0x1.f44659e4a4271p-1, 0x1.5fc17734c36b8p-55}, /* i= 6 */
    {0x1.f25f644230ab5p-1, 0x1.94ed8175c78b3p-58}, /* i= 7 */
    {0x1.f07c1f07c1f08p-1, -0x1.f07c1f07c1f08p-56}, /* i= 8 */
    {0x1.ee9c7f8458e02p-1, -0x1.163807ba71fe1p-57}, /* i= 9 */
    {0x1.ecc07b301eccp-1, 0x1.ecc07b301eccp-55}, /* i= 10 */
    {0x1.eae807aba01ebp-1, -0x1.7f8545fe1518p-57}, /* i= 11 */
    {0x1.e9131abf0b767p-1, 0x1.503d226357e17p-56}, /* i= 12 */
    {0x1.e741aa59750e4p-1, 0x1.9b1f67bb7ac41p-55}, /* i= 13 */
    {0x1.e573ac901e574p-1, -0x1.4dbf86a314dcp-55}, /* i= 14 */
    {0x1.e3a9179dc1a73p-1, 0x1.fa5504b926bb1p-56}, /* i= 15 */
    {0x1.e1e1e1e1e1e1ep-1, 0x1.e1e1e1e1e1e1ep-57}, /* i= 16 */
    {0x1.e01e01e01e01ep-1, 0x1.e01e01e01e01ep-61}, /* i= 17 */
    {0x1.de5d6e3f8868ap-1, 0x1.1c077975b8fe2p-55}, /* i= 18 */
    {0x1.dca01dca01dcap-1, 0x1.dca01dca01dcap-61}, /* i= 19 */
    {0x1.dae6076b981dbp-1, -0x1.9f89467e251ap-57}, /* i= 20 */
    {0x1.d92f2231e7f8ap-1, -0x1.2f2231e7f89b4p-55}, /* i= 21 */
    {0x1.d77b654b82c34p-1, -0x1.ba03aef6ca97p-55}, /* i= 22 */
    {0x1.d5cac807572b2p-1, 0x1.d5cac807572b2p-61}, /* i= 23 */
    {0x1.d41d41d41d41dp-1, 0x1.075075075075p-55}, /* i= 24 */
    {0x1.d272ca3fc5b1ap-1, 0x1.ae01d272ca3fcp-55}, /* i= 25 */
    {0x1.d0cb58f6ec074p-1, 0x1.96b1edd80e866p-56}, /* i= 26 */
    {0x1.cf26e5c44bfc6p-1, 0x1.b2347768073cap-57}, /* i= 27 */
    {0x1.cd85689039b0bp-1, -0x1.76fc64f52edf9p-56}, /* i= 28 */
    {0x1.cbe6d9601cbe7p-1, -0x1.34ff1a0c934ffp-56}, /* i= 29 */
    {0x1.ca4b3055ee191p-1, 0x1.ca4b3055ee191p-61}, /* i= 30 */
    {0x1.c8b265afb8a42p-1, 0x1.c8b265afb8a42p-61}, /* i= 31 */
    {0x1.c71c71c71c71cp-1, 0x1.c71c71c71c71cp-55}, /* i= 32 */
    {0x1.c5894d10d4986p-1, -0x1.f00e2c4a6886ap-56}, /* i= 33 */
    {0x1.c3f8f01c3f8fp-1, 0x1.c3f8f01c3f8fp-57}, /* i= 34 */
    {0x1.c26b5392ea01cp-1, 0x1.35a9c97500e13p-56}, /* i= 35 */
    {0x1.c0e070381c0ep-1, 0x1.c0e070381c0ep-55}, /* i= 36 */
    {0x1.bf583ee868d8bp-1, -0x1.41876d370b5bcp-57}, /* i= 37 */
    {0x1.bdd2b899406f7p-1, 0x1.2b899406f74aep-55}, /* i= 38 */
    {0x1.bc4fd65883e7bp-1, 0x1.d1239464aa169p-56}, /* i= 39 */
    {0x1.bacf914c1badp-1, -0x1.bacf914c1badp-55}, /* i= 40 */
    {0x1.b951e2b18ff23p-1, 0x1.5c3a9ce01b952p-55}, /* i= 41 */
    {0x1.b7d6c3dda338bp-1, 0x1.579fc90527845p-56}, /* i= 42 */
    {0x1.b65e2e3beee05p-1, 0x1.18d4559e6507bp-56}, /* i= 43 */
    {0x1.b4e81b4e81b4fp-1, -0x1.f92c5f92c5f93p-55}, /* i= 44 */
    {0x1.b37484ad806cep-1, -0x1.6f6a4ff2645bep-56}, /* i= 45 */
    {0x1.b2036406c80d9p-1, 0x1.b2036406c80d9p-61}, /* i= 46 */
    {0x1.b094b31d922a4p-1, -0x1.7a821cb9dfe4fp-57}, /* i= 47 */
    {0x1.af286bca1af28p-1, 0x1.af286bca1af28p-55}, /* i= 48 */
    {0x1.adbe87f94905ep-1, 0x1.adbe87f94905ep-61}, /* i= 49 */
    {0x1.ac5701ac5701bp-1, -0x1.d47f29d47f29dp-56}, /* i= 50 */
    {0x1.aaf1d2f87ebfdp-1, -0x1.578e97c3f5fe5p-55}, /* i= 51 */
    {0x1.a98ef606a63bep-1, -0x1.f959c427e5671p-55}, /* i= 52 */
    {0x1.a82e65130e159p-1, -0x1.6937821239fe5p-55}, /* i= 53 */
    {0x1.a6d01a6d01a6dp-1, 0x1.a6d01a6d01a6dp-61}, /* i= 54 */
    {0x1.a574107688a4ap-1, 0x1.566e4d604f05cp-57}, /* i= 55 */
    {0x1.a41a41a41a41ap-1, 0x1.069069069069p-55}, /* i= 56 */
    {0x1.a2c2a87c51cap-1, 0x1.3a11fe5d3d578p-55}, /* i= 57 */
    {0x1.a16d3f97a4b02p-1, -0x1.7a4b01a16d3f9p-55}, /* i= 58 */
    {0x1.a01a01a01a01ap-1, 0x1.a01a01a01a01ap-61}, /* i= 59 */
    {0x1.9ec8e951033d9p-1, 0x1.d2a2067b23a54p-57}, /* i= 60 */
    {0x1.9d79f176b682dp-1, 0x1.cab347dfb2792p-56}, /* i= 61 */
    {0x1.9c2d14ee4a102p-1, -0x1.8f4bac46d7bfap-55}, /* i= 62 */
    {0x1.9ae24ea5510dap-1, 0x1.20e71f4c3cfd9p-55}, /* i= 63 */
    {0x1.999999999999ap-1, -0x1.999999999999ap-55}, /* i= 64 */
    {0x1.9852f0d8ec0ffp-1, 0x1.9eb43c9c4fc03p-56}, /* i= 65 */
    {0x1.970e4f80cb872p-1, 0x1.f01970e4f80ccp-55}, /* i= 66 */
    {0x1.95cbb0be377aep-1, -0x1.b57f9a8d13d07p-55}, /* i= 67 */
    {0x1.948b0fcd6e9ep-1, 0x1.948b0fcd6e9ep-55}, /* i= 68 */
    {0x1.934c67f9b2ce6p-1, 0x1.934c67f9b2ce6p-61}, /* i= 69 */
    {0x1.920fb49d0e229p-1, -0x1.533d406483ed2p-56}, /* i= 70 */
    {0x1.90d4f120190d5p-1, -0x1.dbfcde561dbfdp-58}, /* i= 71 */
    {0x1.8f9c18f9c18fap-1, -0x1.f3831f3831f38p-56}, /* i= 72 */
    {0x1.8e6527af1373fp-1, 0x1.c031cca4f5e27p-59}, /* i= 73 */
    {0x1.8d3018d3018d3p-1, 0x1.8d3018d3018d3p-61}, /* i= 74 */
    {0x1.8bfce8062ff3ap-1, 0x1.8bfce8062ff3ap-61}, /* i= 75 */
    {0x1.8acb90f6bf3aap-1, -0x1.721ed7e75346fp-55}, /* i= 76 */
    {0x1.899c0f601899cp-1, 0x1.ec0313381ec03p-58}, /* i= 77 */
    {0x1.886e5f0abb04ap-1, -0x1.ad38b7f3bc8dp-55}, /* i= 78 */
    {0x1.87427bcc092b9p-1, -0x1.1937c8faa6975p-57}, /* i= 79 */
    {0x1.8618618618618p-1, 0x1.8618618618618p-55}, /* i= 80 */
    {0x1.84f00c2780614p-1, -0x1.fe7b0ff3d87fap-56}, /* i= 81 */
    {0x1.83c977ab2beddp-1, 0x1.4731fcf86d10bp-56}, /* i= 82 */
    {0x1.82a4a0182a4ap-1, 0x1.82a4a0182a4ap-57}, /* i= 83 */
    {0x1.8181818181818p-1, 0x1.8181818181818p-57}, /* i= 84 */
    {0x1.8060180601806p-1, 0x1.8060180601806p-61}, /* i= 85 */
    {0x1.7f405fd017f4p-1, 0x1.7f405fd017f4p-55}, /* i= 86 */
    {0x1.7e225515a4f1dp-1, 0x1.b9d7b26106b7ap-57}, /* i= 87 */
    {0x1.7d05f417d05f4p-1, 0x1.7d05f417d05f4p-57}, /* i= 88 */
    {0x1.7beb3922e017cp-1, -0x1.4c6dd1fe8414cp-57}, /* i= 89 */
    {0x1.7ad2208e0ecc3p-1, 0x1.516324fe852dep-55}, /* i= 90 */
    {0x1.79baa6bb6398bp-1, 0x1.bd9a30b10f7e2p-55}, /* i= 91 */
    {0x1.78a4c8178a4c8p-1, 0x1.78a4c8178a4c8p-57}, /* i= 92 */
    {0x1.77908119ac60dp-1, 0x1.a0a44f387b3b7p-56}, /* i= 93 */
    {0x1.767dce434a9b1p-1, 0x1.767dce434a9b1p-61}, /* i= 94 */
    {0x1.756cac201756dp-1, -0x1.4f7fa2a4d4f8p-55}, /* i= 95 */
    {0x1.745d1745d1746p-1, -0x1.745d1745d1746p-56}, /* i= 96 */
    {0x1.734f0c541fe8dp-1, -0x1.3c31507fa32c4p-55}, /* i= 97 */
    {0x1.724287f46debcp-1, 0x1.724287f46debcp-59}, /* i= 98 */
    {0x1.713786d9c7c09p-1, -0x1.62cb5b9545f3p-55}, /* i= 99 */
    {0x1.702e05c0b817p-1, 0x1.702e05c0b817p-56}, /* i= 100 */
    {0x1.6f26016f26017p-1, -0x1.b3fd21b3fd21bp-58}, /* i= 101 */
    {0x1.6e1f76b4337c7p-1, -0x1.a75461405b87ep-56}, /* i= 102 */
    {0x1.6d1a62681c861p-1, -0x1.3f77161b18f55p-59}, /* i= 103 */
    {0x1.6c16c16c16c17p-1, -0x1.f49f49f49f49fp-56}, /* i= 104 */
    {0x1.6b1490aa31a3dp-1, -0x1.c5d9b4d4be0ccp-60}, /* i= 105 */
    {0x1.6a13cd153729p-1, 0x1.0f8ed9cfe95ecp-55}, /* i= 106 */
    {0x1.691473a88d0cp-1, -0x1.691473a88d0cp-56}, /* i= 107 */
    {0x1.6816816816817p-1, -0x1.fa5fa5fa5fa6p-55}, /* i= 108 */
    {0x1.6719f3601671ap-1, -0x1.93fd31cc193fdp-58}, /* i= 109 */
    {0x1.661ec6a5122f9p-1, 0x1.661ec6a5122f9p-61}, /* i= 110 */
    {0x1.6524f853b4aa3p-1, 0x1.cf2bf20c8e4ccp-56}, /* i= 111 */
    {0x1.642c8590b2164p-1, 0x1.642c8590b2164p-56}, /* i= 112 */
    {0x1.63356b88ac0dep-1, 0x1.63356b88ac0dep-61}, /* i= 113 */
    {0x1.623fa7701624p-1, -0x1.623fa7701624p-55}, /* i= 114 */
    {0x1.614b36831ae94p-1, -0x1.5640dccf0211fp-55}, /* i= 115 */
    {0x1.6058160581606p-1, -0x1.fa7e9fa7e9fa8p-55}, /* i= 116 */
    {0x1.5f66434292dfcp-1, -0x1.e32c9c7b89f3ap-57}, /* i= 117 */
    {0x1.5e75bb8d015e7p-1, 0x1.6ee340579d6eep-55}, /* i= 118 */
    {0x1.5d867c3ece2a5p-1, 0x1.a485cd7b900afp-56}, /* i= 119 */
    {0x1.5c9882b931057p-1, 0x1.310572620ae4cp-56}, /* i= 120 */
    {0x1.5babcc647fa91p-1, 0x1.4339b8056eaf3p-55}, /* i= 121 */
    {0x1.5ac056b015acp-1, 0x1.5ac056b015acp-55}, /* i= 122 */
    {0x1.59d61f123ccaap-1, 0x1.bb1a57cf5de3ap-56}, /* i= 123 */
    {0x1.58ed2308158edp-1, 0x1.1840ac7691841p-56}, /* i= 124 */
    {0x1.580560158056p-1, 0x1.580560158056p-57}, /* i= 125 */
    {0x1.571ed3c506b3ap-1, -0x1.7749b79f7f547p-55}, /* i= 126 */
    {0x1.56397ba7c52e2p-1, -0x1.40d5e3ed48db4p-57}, /* i= 127 */
    {0x1.5555555555555p-1, 0x1.5555555555555p-55}, /* i= 128 */
    {0x1.54725e6bb82fep-1, 0x1.54725e6bb82fep-61}, /* i= 129 */
    {0x1.5390948f40febp-1, -0x1.c84a47a07f563p-56}, /* i= 130 */
    {0x1.52aff56a8054bp-1, -0x1.00a957fab5403p-55}, /* i= 131 */
    {0x1.51d07eae2f815p-1, 0x1.d07eae2f8151dp-57}, /* i= 132 */
    {0x1.50f22e111c4c5p-1, 0x1.b79bf81a52ebap-55}, /* i= 133 */
    {0x1.5015015015015p-1, 0x1.5015015015015p-61}, /* i= 134 */
    {0x1.4f38f62dd4c9bp-1, -0x1.eefa1b7fac31cp-55}, /* i= 135 */
    {0x1.4e5e0a72f0539p-1, 0x1.e0a72f0539783p-55}, /* i= 136 */
    {0x1.4d843bedc2c4cp-1, -0x1.c029b0877db86p-55}, /* i= 137 */
    {0x1.4cab88725af6ep-1, 0x1.d3d137e0cfeb3p-55}, /* i= 138 */
    {0x1.4bd3edda68fe1p-1, -0x1.bde4c79d7d156p-57}, /* i= 139 */
    {0x1.4afd6a052bf5bp-1, -0x1.fad40a57eb503p-55}, /* i= 140 */
    {0x1.4a27fad76014ap-1, 0x1.3fd6bb00a514p-56}, /* i= 141 */
    {0x1.49539e3b2d067p-1, -0x1.5de8d81edfd6dp-57}, /* i= 142 */
    {0x1.488052201488p-1, 0x1.488052201488p-55}, /* i= 143 */
    {0x1.47ae147ae147bp-1, -0x1.eb851eb851eb8p-57}, /* i= 144 */
    {0x1.46dce34596066p-1, 0x1.28382df70ff5dp-56}, /* i= 145 */
    {0x1.460cbc7f5cf9ap-1, 0x1.c051832f1fd74p-57}, /* i= 146 */
    {0x1.453d9e2c776cap-1, 0x1.453d9e2c776cap-61}, /* i= 147 */
    {0x1.446f86562d9fbp-1, -0x1.1be1958b67ebcp-57}, /* i= 148 */
    {0x1.43a2730abee4dp-1, 0x1.db5698f7c8601p-57}, /* i= 149 */
    {0x1.42d6625d51f87p-1, -0x1.064e2febd299ep-57}, /* i= 150 */
    {0x1.420b5265e5951p-1, 0x1.1ed21562c078cp-56}, /* i= 151 */
    {0x1.4141414141414p-1, 0x1.4141414141414p-57}, /* i= 152 */
    {0x1.40782d10e6566p-1, 0x1.909638551fecp-59}, /* i= 153 */
    {0x1.3fb013fb013fbp-1, 0x1.3fb013fb013fbp-61}, /* i= 154 */
    {0x1.3ee8f42a5af07p-1, -0x1.2ff608b85ead3p-56}, /* i= 155 */
    {0x1.3e22cbce4a902p-1, 0x1.f1165e7254814p-55}, /* i= 156 */
    {0x1.3d5d991aa75c6p-1, -0x1.10bc6f92e7d36p-55}, /* i= 157 */
    {0x1.3c995a47babe7p-1, 0x1.1013c995a47bbp-55}, /* i= 158 */
    {0x1.3bd60d9232955p-1, -0x1.f4e57985dc38cp-55}, /* i= 159 */
    {0x1.3b13b13b13b14p-1, -0x1.3b13b13b13b14p-55}, /* i= 160 */
    {0x1.3a524387ac822p-1, 0x1.83fd8b5b78f0ap-55}, /* i= 161 */
    {0x1.3991c2c187f63p-1, 0x1.b8f4f9e027324p-56}, /* i= 162 */
    {0x1.38d22d366088ep-1, -0x1.030e0d7107f15p-55}, /* i= 163 */
    {0x1.3813813813814p-1, -0x1.fb1fb1fb1fb2p-55}, /* i= 164 */
    {0x1.3755bd1c945eep-1, -0x1.f030a5658c773p-56}, /* i= 165 */
    {0x1.3698df3de0748p-1, -0x1.ab1232f514a02p-55}, /* i= 166 */
    {0x1.35dce5f9f2af8p-1, 0x1.0f21493ab4599p-56}, /* i= 167 */
    {0x1.3521cfb2b78c1p-1, 0x1.a90e7d95bc60ap-56}, /* i= 168 */
    {0x1.34679ace01346p-1, 0x1.e6b3804d19e6bp-55}, /* i= 169 */
    {0x1.33ae45b57bcb2p-1, -0x1.f3fb3146e92a1p-57}, /* i= 170 */
    {0x1.32f5ced6a1dfap-1, 0x1.32f5ced6a1dfap-61}, /* i= 171 */
    {0x1.323e34a2b10bfp-1, 0x1.9b8396ba9de81p-55}, /* i= 172 */
    {0x1.3187758e9ebb6p-1, 0x1.3187758e9ebb6p-61}, /* i= 173 */
    {0x1.30d190130d19p-1, 0x1.30d190130d19p-57}, /* i= 174 */
    {0x1.301c82ac4026p-1, 0x1.c82ac4026039p-56}, /* i= 175 */
    {0x1.2f684bda12f68p-1, 0x1.2f684bda12f68p-55}, /* i= 176 */
    {0x1.2eb4ea1fed14bp-1, 0x1.5e012eb4ea1ffp-57}, /* i= 177 */
    {0x1.2e025c04b8097p-1, 0x1.2e025c04b8097p-61}, /* i= 178 */
    {0x1.2d50a012d50ap-1, 0x1.2d50a012d50ap-57}, /* i= 179 */
    {0x1.2c9fb4d812cap-1, -0x1.2c9fb4d812cap-55}, /* i= 180 */
    {0x1.2bef98e5a3711p-1, -0x1.76eb7f1f0c4d5p-60}, /* i= 181 */
    {0x1.2b404ad012b4p-1, 0x1.2b404ad012b4p-55}, /* i= 182 */
    {0x1.2a91c92f3c105p-1, 0x1.fc804aa4724bdp-56}, /* i= 183 */
    {0x1.29e4129e4129ep-1, 0x1.04a7904a7904ap-55}, /* i= 184 */
    {0x1.293725bb804a5p-1, -0x1.1b488ff6b646dp-56}, /* i= 185 */
    {0x1.288b01288b013p-1, -0x1.dd3fb5dd3fb5ep-55}, /* i= 186 */
    {0x1.27dfa38a1ce4dp-1, 0x1.be1f34963f911p-55}, /* i= 187 */
    {0x1.27350b8812735p-1, 0x1.71024e6a17102p-58}, /* i= 188 */
    {0x1.268b37cd60127p-1, -0x1.d320ca7fb65d3p-55}, /* i= 189 */
    {0x1.25e22708092f1p-1, 0x1.3840497889c2p-57}, /* i= 190 */
    {0x1.2539d7e9177b2p-1, 0x1.ca2a615c34b06p-57}, /* i= 191 */
    {0x1.2492492492492p-1, 0x1.2492492492492p-55}, /* i= 192 */
    {0x1.23eb79717605bp-1, 0x1.ccaf9ba70e41p-56}, /* i= 193 */
    {0x1.23456789abcdfp-1, 0x1.23456789abcdfp-61}, /* i= 194 */
    {0x1.22a0122a0122ap-1, 0x1.22a0122a0122ap-61}, /* i= 195 */
    {0x1.21fb78121fb78p-1, 0x1.21fb78121fb78p-57}, /* i= 196 */
    {0x1.21579804855e6p-1, 0x1.21579804855e6p-61}, /* i= 197 */
    {0x1.20b470c67c0d9p-1, -0x1.e2adac8bd766ap-55}, /* i= 198 */
    {0x1.2012012012012p-1, 0x1.2012012012012p-61}, /* i= 199 */
    {0x1.1f7047dc11f7p-1, 0x1.1f7047dc11f7p-55}, /* i= 200 */
    {0x1.1ecf43c7fb84cp-1, 0x1.787008f67a1e4p-56}, /* i= 201 */
    {0x1.1e2ef3b3fb874p-1, 0x1.0c4c0478bbcedp-55}, /* i= 202 */
    {0x1.1d8f5672e4abdp-1, -0x1.f17fb89c2a634p-55}, /* i= 203 */
    {0x1.1cf06ada2811dp-1, -0x1.f2a4bafdc61f3p-58}, /* i= 204 */
    {0x1.1c522fc1ce059p-1, -0x1.32889b7cf21ep-56}, /* i= 205 */
    {0x1.1bb4a4046ed29p-1, 0x1.1bb4a4046ed29p-61}, /* i= 206 */
    {0x1.1b17c67f2bae3p-1, -0x1.37d830a8161dep-55}, /* i= 207 */
    {0x1.1a7b9611a7b96p-1, 0x1.1a7b9611a7b96p-57}, /* i= 208 */
    {0x1.19e0119e0119ep-1, 0x1.19e0119e0119ep-61}, /* i= 209 */
    {0x1.19453808ca29cp-1, 0x1.19453808ca29cp-59}, /* i= 210 */
    {0x1.18ab083902bdbp-1, -0x1.1adc5e4974c32p-55}, /* i= 211 */
    {0x1.1811811811812p-1, -0x1.fb9fb9fb9fbap-55}, /* i= 212 */
    {0x1.1778a191bd684p-1, 0x1.8045de28646f6p-57}, /* i= 213 */
    {0x1.16e0689427379p-1, -0x1.4b2a7c2fee92p-57}, /* i= 214 */
    {0x1.1648d50fc3201p-1, 0x1.648d50fc32011p-57}, /* i= 215 */
    {0x1.15b1e5f75270dp-1, 0x1.15b1e5f75270dp-59}, /* i= 216 */
    {0x1.151b9a3fdd5c9p-1, -0x1.a3fdd5c8cb804p-56}, /* i= 217 */
    {0x1.1485f0e0acd3bp-1, 0x1.a31b011485f0ep-55}, /* i= 218 */
    {0x1.13f0e8d344724p-1, 0x1.c0677a574f39bp-57}, /* i= 219 */
    {0x1.135c81135c811p-1, 0x1.ae4089ae4089bp-56}, /* i= 220 */
    {0x1.12c8b89edc0acp-1, -0x1.0a3272d9e52a6p-55}, /* i= 221 */
    {0x1.12358e75d3033p-1, 0x1.a82ad85e4269p-55}, /* i= 222 */
    {0x1.11a3019a74826p-1, 0x1.ebb0e6e1895a5p-55}, /* i= 223 */
    {0x1.1111111111111p-1, 0x1.1111111111111p-57}, /* i= 224 */
    {0x1.107fbbe01108p-1, -0x1.107fbbe01108p-55}, /* i= 225 */
    {0x1.0fef010fef011p-1, -0x1.0fef010fef011p-61}, /* i= 226 */
    {0x1.0f5edfab325a2p-1, -0x1.5fef0a12054cep-55}, /* i= 227 */
    {0x1.0ecf56be69c9p-1, -0x1.0ecf56be69c9p-56}, /* i= 228 */
    {0x1.0e40655826011p-1, -0x1.bf9aa7d9fef1cp-57}, /* i= 229 */
    {0x1.0db20a88f4696p-1, -0x1.9cf8a021b6415p-55}, /* i= 230 */
    {0x1.0d24456359e3ap-1, -0x1.69a8bd3d80c9ep-56}, /* i= 231 */
    {0x1.0c9714fbcda3bp-1, -0x1.f79b47582192ep-56}, /* i= 232 */
    {0x1.0c0a7868b4171p-1, -0x1.c669c021814f1p-55}, /* i= 233 */
    {0x1.0b7e6ec259dc8p-1, -0x1.b2ad73fbd2064p-55}, /* i= 234 */
    {0x1.0af2f722eecb5p-1, 0x1.c48fe6f938d4cp-55}, /* i= 235 */
    {0x1.0a6810a6810a7p-1, -0x1.fbd65fbd65fbdp-55}, /* i= 236 */
    {0x1.09ddba6af836p-1, 0x1.09ddba6af836p-57}, /* i= 237 */
    {0x1.0953f39010954p-1, -0x1.8dfded5818dfep-58}, /* i= 238 */
    {0x1.08cabb37565e2p-1, 0x1.08cabb37565e2p-61}, /* i= 239 */
    {0x1.0842108421084p-1, 0x1.0842108421084p-56}, /* i= 240 */
    {0x1.07b9f29b8eae2p-1, -0x1.8fb5d3b3c43fep-55}, /* i= 241 */
    {0x1.073260a47f7c6p-1, 0x1.b3eb701073261p-55}, /* i= 242 */
    {0x1.06ab59c7912fbp-1, 0x1.87f3aff7caa53p-55}, /* i= 243 */
    {0x1.0624dd2f1a9fcp-1, -0x1.89374bc6a7efap-57}, /* i= 244 */
    {0x1.059eea0727586p-1, 0x1.8c84dab2d7a2p-55}, /* i= 245 */
    {0x1.05197f7d73404p-1, 0x1.465fdf5cd0105p-57}, /* i= 246 */
    {0x1.04949cc1664c5p-1, 0x1.e27b2a3e17696p-55}, /* i= 247 */
    {0x1.041041041041p-1, 0x1.041041041041p-55}, /* i= 248 */
    {0x1.038c6b78247fcp-1, -0x1.c635bc123fdf9p-58}, /* i= 249 */
    {0x1.03091b51f5e1ap-1, 0x1.3bb3194be3abp-55}, /* i= 250 */
    {0x1.02864fc7729e9p-1, -0x1.d089575a61f4ep-56}, /* i= 251 */
    {0x1.0204081020408p-1, 0x1.0204081020408p-57}, /* i= 252 */
    {0x1.0182436517a37p-1, 0x1.4bf1eae05078bp-55}, /* i= 253 */
    {0x1.010101010101p-1, 0x1.010101010101p-57}, /* i= 254 */
    {0x1.008040201008p-1, 0x1.008040201008p-55} /* i= 255 */
};

/* table_log_alpha_i 
 log_alpha_i = -log{2^k}/{2^k+i}  with 0<=i<2*k with k=8.
 log_alpha_i= {h,l} to calculate in {double, double} precision
 h is the main value and l is the error value.
*/

static const double table_log_alpha_i[256][2]={
    //{h, l} value for each log(alpha_i)
    {0x0p+3, 0x0p+3}, /* i= 0 */
    {0x1.ff00aa2b10bcp-9, 0x1.2821ad5a6d353p-63}, /* i= 1 */
    {0x1.fe02a6b106789p-8, -0x1.e44b7e3711ebfp-67}, /* i= 2 */
    {0x1.7dc475f810a77p-7, -0x1.16d7687d3df21p-62}, /* i= 3 */
    {0x1.fc0a8b0fc03e4p-7, -0x1.83092c59642a1p-62}, /* i= 4 */
    {0x1.3cea44346a575p-6, -0x1.0cb5a902b3a1cp-62}, /* i= 5 */
    {0x1.7b91b07d5b11bp-6, -0x1.5b602ace3a51p-60}, /* i= 6 */
    {0x1.b9fc027af9198p-6, -0x1.0ae69229dc868p-64}, /* i= 7 */
    {0x1.f829b0e7833p-6, 0x1.33e3f04f1ef23p-60}, /* i= 8 */
    {0x1.1b0d98923d98p-5, -0x1.e9ae889bac481p-60}, /* i= 9 */
    {0x1.39e87b9febd6p-5, -0x1.5bfa937f551bbp-59}, /* i= 10 */
    {0x1.58a5bafc8e4d5p-5, -0x1.ce55c2b4e2b72p-59}, /* i= 11 */
    {0x1.77458f632dcfcp-5, 0x1.18d3ca87b9296p-59}, /* i= 12 */
    {0x1.95c830ec8e3ebp-5, 0x1.f5a0e80520bf2p-59}, /* i= 13 */
    {0x1.b42dd711971bfp-5, -0x1.eb9759c130499p-60}, /* i= 14 */
    {0x1.d276b8adb0b52p-5, 0x1.1e3c53257fd47p-61}, /* i= 15 */
    {0x1.f0a30c01162a6p-5, 0x1.85f325c5bbacdp-59}, /* i= 16 */
    {0x1.075983598e471p-4, 0x1.80da5333c45b8p-59}, /* i= 17 */
    {0x1.16536eea37ae1p-4, -0x1.79da3e8c22cdap-60}, /* i= 18 */
    {0x1.253f62f0a1417p-4, -0x1.c125963fc4cfdp-62}, /* i= 19 */
    {0x1.341d7961bd1d1p-4, -0x1.b599f227becbbp-58}, /* i= 20 */
    {0x1.42edcbea646fp-4, 0x1.ddd4f935996c9p-59}, /* i= 21 */
    {0x1.51b073f06183fp-4, 0x1.a49e39a1a8be4p-58}, /* i= 22 */
    {0x1.60658a93750c4p-4, -0x1.388458ec21b6ap-58}, /* i= 23 */
    {0x1.6f0d28ae56b4cp-4, -0x1.906d99184b992p-58}, /* i= 24 */
    {0x1.7da766d7b12cdp-4, -0x1.eeedfcdd94131p-58}, /* i= 25 */
    {0x1.8c345d6319b21p-4, -0x1.4a697ab3424a9p-61}, /* i= 26 */
    {0x1.9ab42462033adp-4, -0x1.2099e1c184e8ep-59}, /* i= 27 */
    {0x1.a926d3a4ad563p-4, 0x1.942f48aa70ea9p-58}, /* i= 28 */
    {0x1.b78c82bb0eda1p-4, 0x1.0878cf0327e21p-61}, /* i= 29 */
    {0x1.c5e548f5bc743p-4, 0x1.5d617ef8161b1p-60}, /* i= 30 */
    {0x1.d4313d66cb35dp-4, 0x1.790dd951d90fap-58}, /* i= 31 */
    {0x1.e27076e2af2e6p-4, -0x1.61578001e0162p-60}, /* i= 32 */
    {0x1.f0a30c01162a6p-4, 0x1.85f325c5bbacdp-58}, /* i= 33 */
    {0x1.fec9131dbeabbp-4, -0x1.5746b9981b36cp-58}, /* i= 34 */
    {0x1.0671512ca596ep-3, 0x1.50c647eb86499p-58}, /* i= 35 */
    {0x1.0d77e7cd08e59p-3, 0x1.9a5dc5e9030acp-57}, /* i= 36 */
    {0x1.14785846742acp-3, 0x1.a28813e3a7f07p-57}, /* i= 37 */
    {0x1.1b72ad52f67ap-3, 0x1.483023472cd74p-58}, /* i= 38 */
    {0x1.2266f190a5acbp-3, 0x1.f547bf1809e88p-57}, /* i= 39 */
    {0x1.29552f81ff523p-3, 0x1.301771c407dbfp-57}, /* i= 40 */
    {0x1.303d718e47fd3p-3, -0x1.6b9c7d96091fap-63}, /* i= 41 */
    {0x1.371fc201e8f74p-3, 0x1.de6cb62af18ap-58}, /* i= 42 */
    {0x1.3dfc2b0ecc62ap-3, -0x1.ab3a8e7d81017p-58}, /* i= 43 */
    {0x1.44d2b6ccb7d1ep-3, 0x1.9f4f6543e1f88p-57}, /* i= 44 */
    {0x1.4ba36f39a55e5p-3, 0x1.68981bcc36756p-57}, /* i= 45 */
    {0x1.526e5e3a1b438p-3, -0x1.746ff8a470d3ap-57}, /* i= 46 */
    {0x1.59338d9982086p-3, -0x1.65d22aa8ad7cfp-58}, /* i= 47 */
    {0x1.5ff3070a793d4p-3, -0x1.bc60efafc6f6ep-58}, /* i= 48 */
    {0x1.66acd4272ad51p-3, -0x1.0900e4e1ea8b2p-58}, /* i= 49 */
    {0x1.6d60fe719d21dp-3, -0x1.caae268ecd179p-57}, /* i= 50 */
    {0x1.740f8f54037a5p-3, -0x1.b264062a84cdbp-58}, /* i= 51 */
    {0x1.7ab890210d909p-3, 0x1.be36b2d6a0608p-59}, /* i= 52 */
    {0x1.815c0a14357ebp-3, -0x1.4be48073a0564p-58}, /* i= 53 */
    {0x1.87fa06520c911p-3, -0x1.bf7fdbfa08d9ap-57}, /* i= 54 */
    {0x1.8e928de886d41p-3, -0x1.569d851a5677p-57}, /* i= 55 */
    {0x1.9525a9cf456b4p-3, 0x1.d904c1d4e2e26p-57}, /* i= 56 */
    {0x1.9bb362e7dfb83p-3, 0x1.575e31f003e0cp-57}, /* i= 57 */
    {0x1.a23bc1fe2b563p-3, 0x1.93711b07a998cp-59}, /* i= 58 */
    {0x1.a8becfc882f19p-3, -0x1.e8c37918c39ebp-58}, /* i= 59 */
    {0x1.af3c94e80bff3p-3, -0x1.398cff3641985p-58}, /* i= 60 */
    {0x1.b5b519e8fb5a4p-3, 0x1.ba27fdc19e1ap-57}, /* i= 61 */
    {0x1.bc286742d8cd6p-3, 0x1.4fce744870f55p-58}, /* i= 62 */
    {0x1.c2968558c18c1p-3, -0x1.73dee38a3fb6bp-57}, /* i= 63 */
    {0x1.c8ff7c79a9a22p-3, -0x1.4f689f8434012p-57}, /* i= 64 */
    {0x1.cf6354e09c5dcp-3, 0x1.239a07d55b695p-57}, /* i= 65 */
    {0x1.d5c216b4fbb91p-3, 0x1.6e443597e4d4p-57}, /* i= 66 */
    {0x1.dc1bca0abec7dp-3, 0x1.834c51998b6fcp-57}, /* i= 67 */
    {0x1.e27076e2af2e6p-3, -0x1.61578001e0162p-59}, /* i= 68 */
    {0x1.e8c0252aa5a6p-3, -0x1.6e03a39bfc89bp-59}, /* i= 69 */
    {0x1.ef0adcbdc5936p-3, 0x1.48637950dc20dp-57}, /* i= 70 */
    {0x1.f550a564b7b37p-3, 0x1.c5f6dfd018c37p-61}, /* i= 71 */
    {0x1.fb9186d5e3e2bp-3, -0x1.caaae64f21acbp-57}, /* i= 72 */
    {0x1.00e6c45ad501dp-2, -0x1.cb9568ff6feadp-57}, /* i= 73 */
    {0x1.0402594b4d041p-2, -0x1.28ec217a5022dp-57}, /* i= 74 */
    {0x1.071b85fcd590dp-2, 0x1.d1707f97bde8p-58}, /* i= 75 */
    {0x1.0a324e27390e3p-2, 0x1.7dcfde8061c03p-56}, /* i= 76 */
    {0x1.0d46b579ab74bp-2, 0x1.03ec81c3cbd92p-57}, /* i= 77 */
    {0x1.1058bf9ae4ad5p-2, 0x1.89fa0ab4cb31dp-58}, /* i= 78 */
    {0x1.136870293a8bp-2, 0x1.7b66298edd24ap-56}, /* i= 79 */
    {0x1.1675cababa60ep-2, 0x1.ce63eab883717p-61}, /* i= 80 */
    {0x1.1980d2dd4236fp-2, 0x1.9d3d1b0e4d147p-56}, /* i= 81 */
    {0x1.1c898c16999fbp-2, -0x1.0e5c62aff1c44p-60}, /* i= 82 */
    {0x1.1f8ff9e48a2f3p-2, -0x1.c9fdf9a0c4b07p-56}, /* i= 83 */
    {0x1.22941fbcf7966p-2, -0x1.76f5eb09628afp-56}, /* i= 84 */
    {0x1.2596010df763ap-2, -0x1.0f76c57075e9ep-58}, /* i= 85 */
    {0x1.2895a13de86a3p-2, 0x1.7ad24c13f040ep-56}, /* i= 86 */
    {0x1.2b9303ab89d25p-2, -0x1.896b5fd852ad4p-56}, /* i= 87 */
    {0x1.2e8e2bae11d31p-2, -0x1.8f4cdb95ebdf9p-56}, /* i= 88 */
    {0x1.31871c9544185p-2, -0x1.51acc4c09b379p-60}, /* i= 89 */
    {0x1.347dd9a987d55p-2, -0x1.4dd4c580919f8p-57}, /* i= 90 */
    {0x1.3772662bfd85bp-2, -0x1.b5629d8117de7p-59}, /* i= 91 */
    {0x1.3a64c556945eap-2, -0x1.c68651945f97cp-57}, /* i= 92 */
    {0x1.3d54fa5c1f71p-2, -0x1.e3265c6a1c98dp-56}, /* i= 93 */
    {0x1.404308686a7e4p-2, -0x1.0bcfb6082ce6dp-56}, /* i= 94 */
    {0x1.432ef2a04e814p-2, -0x1.29931715ac903p-56}, /* i= 95 */
    {0x1.4618bc21c5ec2p-2, 0x1.f42decdeccf1dp-56}, /* i= 96 */
    {0x1.49006804009d1p-2, -0x1.9ffc341f177dcp-57}, /* i= 97 */
    {0x1.4be5f957778a1p-2, -0x1.259b35b04813dp-57}, /* i= 98 */
    {0x1.4ec973260026ap-2, -0x1.42a87d977dc5ep-56}, /* i= 99 */
    {0x1.51aad872df82dp-2, 0x1.3927ac19f55e3p-59}, /* i= 100 */
    {0x1.548a2c3add263p-2, -0x1.819cf7e308ddbp-57}, /* i= 101 */
    {0x1.5767717455a6cp-2, 0x1.526adb283660cp-56}, /* i= 102 */
    {0x1.5a42ab0f4cfe2p-2, -0x1.8ebcb7dee9a3dp-56}, /* i= 103 */
    {0x1.5d1bdbf5809cap-2, 0x1.4236383dc7fe1p-56}, /* i= 104 */
    {0x1.5ff3070a793d4p-2, -0x1.bc60efafc6f6ep-57}, /* i= 105 */
    {0x1.62c82f2b9c795p-2, 0x1.7b7af915300e5p-57}, /* i= 106 */
    {0x1.659b57303e1f3p-2, -0x1.f893d41c411f1p-56}, /* i= 107 */
    {0x1.686c81e9b14afp-2, -0x1.ddea0f7f58e3dp-57}, /* i= 108 */
    {0x1.6b3bb2235943ep-2, -0x1.da856ccd987b3p-56}, /* i= 109 */
    {0x1.6e08eaa2ba1e4p-2, -0x1.cfb1b39ca3a0fp-56}, /* i= 110 */
    {0x1.70d42e2789236p-2, -0x1.52cc811d78d59p-57}, /* i= 111 */
    {0x1.739d7f6bbd007p-2, -0x1.8c76ceb014b04p-56}, /* i= 112 */
    {0x1.7664e1239dbcfp-2, -0x1.f6d5d64f5daf8p-57}, /* i= 113 */
    {0x1.792a55fdd47a2p-2, 0x1.f057691fe9ed7p-56}, /* i= 114 */
    {0x1.7bede0a37afcp-2, -0x1.8783cb9801a5cp-56}, /* i= 115 */
    {0x1.7eaf83b82afc3p-2, 0x1.92ce979ed295p-56}, /* i= 116 */
    {0x1.816f41da0d496p-2, -0x1.2923ca04b701cp-56}, /* i= 117 */
    {0x1.842d1da1e8b17p-2, 0x1.24ec519784676p-56}, /* i= 118 */
    {0x1.86e919a330bap-2, 0x1.3f9b16feb7dd8p-59}, /* i= 119 */
    {0x1.89a3386c1425bp-2, -0x1.29639dfbbf0fbp-56}, /* i= 120 */
    {0x1.8c5b7c858b48bp-2, -0x1.e0ab4fdfa0595p-56}, /* i= 121 */
    {0x1.8f11e873662c7p-2, 0x1.f85da755a61a3p-56}, /* i= 122 */
    {0x1.91c67eb45a83ep-2, -0x1.e0e0ae234ae11p-56}, /* i= 123 */
    {0x1.947941c2116fbp-2, -0x1.16cc8bae0bbe4p-56}, /* i= 124 */
    {0x1.972a341135158p-2, 0x1.a5c09d24b70d9p-56}, /* i= 125 */
    {0x1.99d958117e08bp-2, -0x1.a2b6889dc3e72p-57}, /* i= 126 */
    {0x1.9c86b02dc0863p-2, -0x1.917eeb69dd421p-56}, /* i= 127 */
    {0x1.9f323ecbf984cp-2, -0x1.a92e513217f5cp-59}, /* i= 128 */
    {0x1.a1dc064d5b995p-2, 0x1.90128698ba0b8p-56}, /* i= 129 */
    {0x1.a484090e5bb0ap-2, 0x1.5fe535b875a75p-57}, /* i= 130 */
    {0x1.a72a4966bd9eap-2, 0x1.6a76b1a7d87c3p-58}, /* i= 131 */
    {0x1.a9cec9a9a084ap-2, -0x1.cadec02b436afp-56}, /* i= 132 */
    {0x1.ac718c258b0e4p-2, 0x1.8163d6f46f714p-59}, /* i= 133 */
    {0x1.af1293247786bp-2, 0x1.133844a15dc28p-58}, /* i= 134 */
    {0x1.b1b1e0ebdfc5bp-2, 0x1.a4479608a2c55p-56}, /* i= 135 */
    {0x1.b44f77bcc8f63p-2, -0x1.cd04495459c78p-56}, /* i= 136 */
    {0x1.b6eb59d3cf35ep-2, -0x1.8adbccd326a3cp-56}, /* i= 137 */
    {0x1.b9858969310fbp-2, 0x1.663ec53e23bc4p-56}, /* i= 138 */
    {0x1.bc1e08b0dad0ap-2, 0x1.09e8707055996p-56}, /* i= 139 */
    {0x1.beb4d9da71b7cp-2, -0x1.0f3c590a887cap-59}, /* i= 140 */
    {0x1.c149ff115f027p-2, -0x1.4cbcb90c06305p-56}, /* i= 141 */
    {0x1.c3dd7a7cdad4dp-2, 0x1.cecf052dea69bp-56}, /* i= 142 */
    {0x1.c66f4e3ff6ff8p-2, -0x1.82947258b688bp-58}, /* i= 143 */
    {0x1.c8ff7c79a9a22p-2, -0x1.4f689f8434012p-56}, /* i= 144 */
    {0x1.cb8e0744d7acap-2, -0x1.48879a214a2afp-61}, /* i= 145 */
    {0x1.ce1af0b85f3ebp-2, 0x1.edf4af2ab4267p-56}, /* i= 146 */
    {0x1.d0a63ae721e64p-2, 0x1.2acce112c40f2p-57}, /* i= 147 */
    {0x1.d32fe7e00ebd5p-2, 0x1.877b232fafa37p-56}, /* i= 148 */
    {0x1.d5b7f9ae2c684p-2, -0x1.a7be7f84ac06ap-57}, /* i= 149 */
    {0x1.d83e7258a2f3ep-2, 0x1.41456e8bb2511p-56}, /* i= 150 */
    {0x1.dac353e2c5954p-2, 0x1.18734b81a1bf8p-57}, /* i= 151 */
    {0x1.dd46a04c1c4a1p-2, -0x1.0467656d8b892p-56}, /* i= 152 */
    {0x1.dfc859906d5b5p-2, 0x1.01e1399f96398p-56}, /* i= 153 */
    {0x1.e24881a7c6c26p-2, 0x1.cbd8f45954a46p-58}, /* i= 154 */
    {0x1.e4c71a8687704p-2, 0x1.667923e1f5a8ep-57}, /* i= 155 */
    {0x1.e744261d68788p-2, -0x1.c825c90c344b9p-58}, /* i= 156 */
    {0x1.e9bfa659861f5p-2, 0x1.91bafc7dbe13p-56}, /* i= 157 */
    {0x1.ec399d2468ccp-2, 0x1.75cee53f35397p-58}, /* i= 158 */
    {0x1.eeb20c640ddf4p-2, 0x1.ac371d7c8f7f5p-57}, /* i= 159 */
    {0x1.f128f5faf06edp-2, -0x1.328df13bb38c3p-56}, /* i= 160 */
    {0x1.f39e5bc811e5cp-2, -0x1.97fc777bb19e5p-57}, /* i= 161 */
    {0x1.f6123fa7028acp-2, 0x1.8515b0f2db341p-56}, /* i= 162 */
    {0x1.f884a36fe9ec2p-2, 0x1.6315c9e0108p-57}, /* i= 163 */
    {0x1.faf588f78f31fp-2, -0x1.328260d8abcap-57}, /* i= 164 */
    {0x1.fd64f20f61572p-2, -0x1.adb0ac2cead1bp-57}, /* i= 165 */
    {0x1.ffd2e0857f498p-2, 0x1.565f40d9321afp-56}, /* i= 166 */
    {0x1.011fab125ff8ap-1, 0x1.810dd40845ddep-57}, /* i= 167 */
    {0x1.02552a5a5d0ffp-1, -0x1.cb1cb51408cp-56}, /* i= 168 */
    {0x1.0389eefce633bp-1, 0x1.e155c53483748p-56}, /* i= 169 */
    {0x1.04bdf9da926d2p-1, 0x1.97f304022c9dfp-55}, /* i= 170 */
    {0x1.05f14bd26459cp-1, 0x1.535b8ee4f9efep-58}, /* i= 171 */
    {0x1.0723e5c1cdf4p-1, 0x1.395e58e2445bbp-55}, /* i= 172 */
    {0x1.0855c884b450ep-1, 0x1.705826e49f318p-55}, /* i= 173 */
    {0x1.0986f4f573521p-1, -0x1.1b8095ac02f01p-55}, /* i= 174 */
    {0x1.0ab76bece14d2p-1, -0x1.fd6c935453f66p-56}, /* i= 175 */
    {0x1.0be72e4252a83p-1, -0x1.259da11330801p-55}, /* i= 176 */
    {0x1.0d163ccb9d6b8p-1, -0x1.f7b9a9a8bc30fp-57}, /* i= 177 */
    {0x1.0e44985d1cc8cp-1, -0x1.22a3442d2d384p-58}, /* i= 178 */
    {0x1.0f7241c9b497dp-1, 0x1.3a8443b9db19dp-55}, /* i= 179 */
    {0x1.109f39e2d4c97p-1, -0x1.0e09b27a4373ap-60}, /* i= 180 */
    {0x1.11cb81787ccf8p-1, 0x1.02387ab1fcc9p-55}, /* i= 181 */
    {0x1.12f719593efbcp-1, 0x1.4c048c671f435p-55}, /* i= 182 */
    {0x1.1422025243d45p-1, -0x1.ad0e24adb489ep-58}, /* i= 183 */
    {0x1.154c3d2f4d5eap-1, -0x1.59c33171a6876p-55}, /* i= 184 */
    {0x1.1675cababa60ep-1, 0x1.ce63eab883717p-60}, /* i= 185 */
    {0x1.179eabbd899a1p-1, -0x1.00e7c6417e0b4p-55}, /* i= 186 */
    {0x1.18c6e0ff5cf06p-1, 0x1.765142c2c671fp-58}, /* i= 187 */
    {0x1.19ee6b467c96fp-1, -0x1.9d1a11443f10cp-56}, /* i= 188 */
    {0x1.1b154b57da29fp-1, -0x1.011eb47db6a99p-57}, /* i= 189 */
    {0x1.1c3b81f713c25p-1, -0x1.0dac1c4c810e9p-55}, /* i= 190 */
    {0x1.1d610fe677003p-1, 0x1.09d58d91e58f2p-58}, /* i= 191 */
    {0x1.1e85f5e7040dp-1, 0x1.ef62cd2f9f1e3p-56}, /* i= 192 */
    {0x1.1faa34b87094cp-1, 0x1.817b8f7a193bp-58}, /* i= 193 */
    {0x1.20cdcd192ab6ep-1, -0x1.b2bf0bc229014p-55}, /* i= 194 */
    {0x1.21f0bfc65beecp-1, -0x1.e24f0c9187c92p-57}, /* i= 195 */
    {0x1.23130d7bebf43p-1, -0x1.f48725e374d6ep-55}, /* i= 196 */
    {0x1.2434b6f483934p-1, -0x1.debb8cf0f6d11p-57}, /* i= 197 */
    {0x1.2555bce98f7cbp-1, 0x1.e021d6d6881e7p-56}, /* i= 198 */
    {0x1.26762013430ep-1, -0x1.96a95781c6727p-56}, /* i= 199 */
    {0x1.2795e1289b11bp-1, -0x1.487c0c246978ep-57}, /* i= 200 */
    {0x1.28b500df60783p-1, -0x1.43f60605aaab3p-55}, /* i= 201 */
    {0x1.29d37fec2b08bp-1, -0x1.bd1949a2d1982p-56}, /* i= 202 */
    {0x1.2af15f02640adp-1, 0x1.cb064524acebp-57}, /* i= 203 */
    {0x1.2c0e9ed448e8cp-1, -0x1.1a158f3917586p-55}, /* i= 204 */
    {0x1.2d2b4012edc9ep-1, -0x1.51162c99b1cabp-55}, /* i= 205 */
    {0x1.2e47436e40268p-1, 0x1.0150861a4886bp-55}, /* i= 206 */
    {0x1.2f62a99509546p-1, 0x1.6c686739ffd99p-56}, /* i= 207 */
    {0x1.307d7334f10bep-1, 0x1.fb590a1f566dap-57}, /* i= 208 */
    {0x1.3197a0fa7fe6ap-1, 0x1.d6348fb97128fp-57}, /* i= 209 */
    {0x1.32b1339121d71p-1, 0x1.902ab5b3d916bp-56}, /* i= 210 */
    {0x1.33ca2ba328995p-1, -0x1.bf28b3205ede1p-56}, /* i= 211 */
    {0x1.34e289d9ce1d3p-1, 0x1.6eb92d885ce4fp-57}, /* i= 212 */
    {0x1.35fa4edd36eap-1, 0x1.27d4680964362p-60}, /* i= 213 */
    {0x1.37117b54747b6p-1, -0x1.d117edbdd9103p-56}, /* i= 214 */
    {0x1.38280fe58797fp-1, -0x1.015bd362a6e5dp-55}, /* i= 215 */
    {0x1.393e0d3562a1ap-1, -0x1.58eef67f2483ap-55}, /* i= 216 */
    {0x1.3a5373e7ebdfap-1, -0x1.cd8f775b8f76ep-55}, /* i= 217 */
    {0x1.3b68449fffc23p-1, -0x1.41c484f9e9b26p-55}, /* i= 218 */
    {0x1.3c7c7fff73206p-1, -0x1.be80db7025bedp-56}, /* i= 219 */
    {0x1.3d9026a7156fbp-1, -0x1.6fef670bd4b62p-55}, /* i= 220 */
    {0x1.3ea33936b2f5cp-1, -0x1.f099168a1360bp-55}, /* i= 221 */
    {0x1.3fb5b84d16f42p-1, 0x1.6d3a754172aefp-55}, /* i= 222 */
    {0x1.40c7a4880dce9p-1, 0x1.14f22de7fc9e1p-56}, /* i= 223 */
    {0x1.41d8fe84672aep-1, 0x1.9192f30bd1806p-55}, /* i= 224 */
    {0x1.42e9c6ddf80bfp-1, 0x1.657dc7a65061dp-56}, /* i= 225 */
    {0x1.43f9fe2f9ce67p-1, 0x1.e9c9ee6d83b86p-55}, /* i= 226 */
    {0x1.4509a5133bb0ap-1, 0x1.40fe2852d7b5ap-55}, /* i= 227 */
    {0x1.4618bc21c5ec2p-1, 0x1.f42decdeccf1dp-55}, /* i= 228 */
    {0x1.472743f33aaadp-1, 0x1.8d6cf012a2948p-56}, /* i= 229 */
    {0x1.48353d1ea88dfp-1, 0x1.cf57a2ecc07f4p-55}, /* i= 230 */
    {0x1.4942a83a2fc07p-1, 0x1.ed0c544652b5ap-55}, /* i= 231 */
    {0x1.4a4f85db03ebbp-1, 0x1.13dfa3d3761b6p-60}, /* i= 232 */
    {0x1.4b5bd6956e274p-1, -0x1.c87a06beea773p-55}, /* i= 233 */
    {0x1.4c679afccee3ap-1, -0x1.3a5c4c8b39e41p-55}, /* i= 234 */
    {0x1.4d72d3a39fdp-1, 0x1.1cd4d414e008dp-55}, /* i= 235 */
    {0x1.4e7d811b75bb1p-1, -0x1.8d3d9ea6e9ea9p-55}, /* i= 236 */
    {0x1.4f87a3f5026e9p-1, -0x1.e8ca8b1bcea9dp-55}, /* i= 237 */
    {0x1.50913cc01686bp-1, 0x1.2f2ce96c2d5b1p-55}, /* i= 238 */
    {0x1.519a4c0ba3446p-1, 0x1.9b32128e4a77fp-55}, /* i= 239 */
    {0x1.52a2d265bc5abp-1, -0x1.1883750ea4d0ap-57}, /* i= 240 */
    {0x1.53aad05b99b7dp-1, -0x1.55c8b052e2539p-55}, /* i= 241 */
    {0x1.54b2467999498p-1, -0x1.5baaf5d2f09f4p-55}, /* i= 242 */
    {0x1.55b9354b40bcdp-1, 0x1.e4197a357cb37p-56}, /* i= 243 */
    {0x1.56bf9d5b3f399p-1, 0x1.0471885cd8ff3p-55}, /* i= 244 */
    {0x1.57c57f336f191p-1, -0x1.e953a3bc88192p-55}, /* i= 245 */
    {0x1.58cadb5cd7989p-1, 0x1.849792ec98458p-56}, /* i= 246 */
    {0x1.59cfb25fae87ep-1, -0x1.172904559c6b6p-58}, /* i= 247 */
    {0x1.5ad404c359f2dp-1, -0x1.35955683f7196p-59}, /* i= 248 */
    {0x1.5bd7d30e71c73p-1, 0x1.bf8da6db2b45cp-57}, /* i= 249 */
    {0x1.5cdb1dc6c1765p-1, -0x1.cc2470e8a3df4p-55}, /* i= 250 */
    {0x1.5ddde57149923p-1, 0x1.dcfa37d75ef28p-55}, /* i= 251 */
    {0x1.5ee02a9241675p-1, 0x1.c358257f49082p-55}, /* i= 252 */
    {0x1.5fe1edad18919p-1, -0x1.ca8b610e18dbfp-55}, /* i= 253 */
    {0x1.60e32f44788d9p-1, -0x1.ac1bb52fa589bp-56}, /* i= 254 */
    {0x1.61e3efda46467p-1, -0x1.a1b727edefae3p-55} /* i= 255 */
};

/* table_alpha_i_modified:  
 alpha_i = {2^k}/{2^k+i}  with 0<=i<2*k with k=8.
 alpha_i_modified  is the nearest double precision alpha_i for
 log(alpha_i_modified) is close to a double precision.
*/

static const double table_alpha_i_modified[256]=
{
    0x1p+0,                /*i = 0 */  
    0xf.f00ff045af908p-4 , /* i = 1 */
    0xf.e03f8108eecf8p-4 , /* i = 2 */
    0xf.d08e552fabf2p-4 , /* i = 3 */
    0xf.c0fc0ff7c5c68p-4 , /* i = 4 */
    0xf.b188567de7a8p-4 , /* i = 5 */
    0xf.a232cf59f7108p-4 , /* i = 6 */
    0xf.92fb22280244p-4 , /* i = 7 */
    0xf.83e0f8d0da14p-4 , /* i = 8 */
    0xf.74e3fc5af46a8p-4 , /* i = 9 */
    0xf.6603d9cacd3fp-4 , /* i = 10 */
    0xf.57403d7a0b0d8p-4 , /* i = 11 */
    0xf.4898d66c8bdbp-4 , /* i = 12 */
    0xf.3a0d52e42cacp-4 , /* i = 13 */
    0xf.2b9d64d345b48p-4 , /* i = 14 */
    0xf.1d48bd26d36fp-4 , /* i = 15 */
    0xf.0f0f0fcca6db8p-4 , /* i = 16 */
    0xf.00f00f3c8cab8p-4 , /* i = 17 */
    0xe.f2eb72369a5b8p-4 , /* i = 18 */
    0xe.e500ee87dbap-4 , /* i = 19 */
    0xe.d7303ba994bc8p-4 , /* i = 20 */
    0xe.c97911ba76e4p-4 , /* i = 21 */
    0xe.bbdb2a9299fa8p-4 , /* i = 22 */
    0xe.ae56406e1f908p-4 , /* i = 23 */
    0xe.a0ea0f1d39b28p-4 , /* i = 24 */
    0xe.93965226004c8p-4 , /* i = 25 */
    0xe.865ac7f9613p-4 , /* i = 26 */
    0xe.79372e4bb12ap-4 , /* i = 27 */
    0xe.6c2b44ebb1dep-4 , /* i = 28 */
    0xe.5f36cb396cd5p-4 , /* i = 29 */
    0xe.525982f00277p-4 , /* i = 30 */
    0xe.45932d985cfap-4 , /* i = 31 */
    0xe.38e38e9768338p-4 , /* i = 32 */
    0xe.2c4a68b9cae68p-4 , /* i = 33 */
    0xe.1fc781122b158p-4 , /* i = 34 */
    0xe.135a9cbe24dc8p-4 , /* i = 35 */
    0xe.07038216ed6f8p-4 , /* i = 36 */
    0xd.fac1f77693bf8p-4 , /* i = 37 */
    0xd.ee95c5139bc1p-4 , /* i = 38 */
    0xd.e27eb2f1941p-4 , /* i = 39 */
    0xd.d67c8adc8fd78p-4 , /* i = 40 */
    0xd.ca8f15ad5d8p-4 , /* i = 41 */
    0xd.beb61f256a14p-4 , /* i = 42 */
    0xd.b2f171ffeb2c8p-4 , /* i = 43 */
    0xd.a740dac21d328p-4 , /* i = 44 */
    0xd.9ba425910a698p-4 , /* i = 45 */
    0xd.901b207390498p-4 , /* i = 46 */
    0xd.84a599039f4d8p-4 , /* i = 47 */
    0xd.79435f19ae22p-4 , /* i = 48 */
    0xd.6df43fe6c51cp-4 , /* i = 49 */
    0xd.62b80d8e75d1p-4 , /* i = 50 */
    0xd.578e97ec6a408p-4 , /* i = 51 */
    0xd.4c77b06b1fadp-4 , /* i = 52 */
    0xd.417328b8c24fp-4 , /* i = 53 */
    0xd.3680d3b234c48p-4 , /* i = 54 */
    0xd.2ba083d623d68p-4 , /* i = 55 */
    0xd.20d20d72c8dap-4 , /* i = 56 */
    0xd.16154408ffe3p-4 , /* i = 57 */
    0xd.0b69fd0119cap-4 , /* i = 58 */
    0xd.00d00d1c06a18p-4 , /* i = 59 */
    0xc.f6474acda179p-4 , /* i = 60 */
    0xc.ebcf8bd4ef26p-4 , /* i = 61 */
    0xc.e168a7b18a5dp-4 , /* i = 62 */
    0xc.d7127543c9cb8p-4 , /* i = 63 */
    0xc.cccccdb7387c8p-4 , /* i = 64 */
    0xc.c29786f44dc98p-4 , /* i = 65 */
    0xc.b8727c3becc6p-4 , /* i = 66 */
    0xc.ae5d861eecf98p-4 , /* i = 67 */
    0xc.a4587ea9d1878p-4 , /* i = 68 */
    0xc.9a633ff8547f8p-4 , /* i = 69 */
    0xc.907da5257dec8p-4 , /* i = 70 */
    0xc.86a78918a1f7p-4 , /* i = 71 */
    0xc.7ce0c839f47a8p-4 , /* i = 72 */
    0xc.73293da2f962p-4 , /* i = 73 */
    0xc.6980c6bd1db58p-4 , /* i = 74 */
    0xc.5fe74052e1268p-4 , /* i = 75 */
    0xc.565c87cfc5e2p-4 , /* i = 76 */
    0xc.4ce07b2b81aa8p-4 , /* i = 77 */
    0xc.4372f89593908p-4 , /* i = 78 */
    0xc.3a13de855e68p-4 , /* i = 79 */
    0xc.30c30ce878168p-4 , /* i = 80 */
    0xc.2780615890bap-4 , /* i = 81 */
    0xc.1e4bbd6997738p-4 , /* i = 82 */
    0xc.152500d96d5e8p-4 , /* i = 83 */
    0xc.0c0c0c5dae508p-4 , /* i = 84 */
    0xc.0300c051708cp-4 , /* i = 85 */
    0xb.fa02fea74904p-4 , /* i = 86 */
    0xb.f112a8c142e48p-4 , /* i = 87 */
    0xb.e82fa1005fc7p-4 , /* i = 88 */
    0xb.df59c93cee57p-4 , /* i = 89 */
    0xb.d691049c7ba7p-4 , /* i = 90 */
    0xb.cdd536063ac7p-4 , /* i = 91 */
    0xb.c52640da60b78p-4 , /* i = 92 */
    0xb.bc8408ead5cap-4 , /* i = 93 */
    0xb.b3ee725334e8p-4 , /* i = 94 */
    0xb.ab6561156d028p-4 , /* i = 95 */
    0xb.a2e8ba671e91p-4 , /* i = 96 */
    0xb.9a7862cf62638p-4 , /* i = 97 */
    0xb.92143fbed2738p-4 , /* i = 98 */
    0xb.89bc36e8d624p-4 , /* i = 99 */
    0xb.81702e56df19p-4 , /* i = 100 */
    0xb.79300b87be24p-4 , /* i = 101 */
    0xb.70fbb5dd3cc9p-4 , /* i = 102 */
    0xb.68d3136ac46c8p-4 , /* i = 103 */
    0xb.60b60bdf30afp-4 , /* i = 104 */
    0xb.58a48579770e8p-4 , /* i = 105 */
    0xb.509e68c1c4aep-4 , /* i = 106 */
    0xb.48a39d5083818p-4 , /* i = 107 */
    0xb.40b40b7e3ebd8p-4 , /* i = 108 */
    0xb.38cf9b19e41c8p-4 , /* i = 109 */
    0xb.30f63567439ep-4 , /* i = 110 */
    0xb.2927c2b04e228p-4 , /* i = 111 */
    0xb.21642cbd22118p-4 , /* i = 112 */
    0xb.19ab5c6263b88p-4 , /* i = 113 */
    0xb.11fd3bbea3bp-4 , /* i = 114 */
    0xb.0a59b443a338p-4 , /* i = 115 */
    0xb.02c0b06499288p-4 , /* i = 116 */
    0xa.fb321a2528078p-4 , /* i = 117 */
    0xa.f3addc9b96248p-4 , /* i = 118 */
    0xa.ec33e1f9cdb18p-4 , /* i = 119 */
    0xa.e4c4162ca482p-4 , /* i = 120 */
    0xa.dd5e63308a72p-4 , /* i = 121 */
    0xa.d602b5ab037p-4 , /* i = 122 */
    0xa.ceb0f8a03406p-4 , /* i = 123 */
    0xa.c769188974f4p-4 , /* i = 124 */
    0xa.c02b00c2cb4c8p-4 , /* i = 125 */
    0xa.b8f69e39dbc68p-4 , /* i = 126 */
    0xa.b1cbdd50a5d9p-4 , /* i = 127 */
    0xa.aaaaab18a3358p-4 , /* i = 128 */
    0xa.a392f3881d198p-4 , /* i = 129 */
    0xa.9c84a47d66458p-4 , /* i = 130 */
    0xa.957fab76bc42p-4 , /* i = 131 */
    0xa.8e83f5a2b3c9p-4 , /* i = 132 */
    0xa.879170b289cd8p-4 , /* i = 133 */
    0xa.80a80ab2592b8p-4 , /* i = 134 */
    0xa.79c7b18b37b18p-4 , /* i = 135 */
    0xa.72f053ccd1f48p-4 , /* i = 136 */
    0xa.6c21df86ffcb8p-4 , /* i = 137 */
    0xa.655c43c0c1b58p-4 , /* i = 138 */
    0xa.5e9f6ef778d3p-4 , /* i = 139 */
    0xa.57eb503d8cda8p-4 , /* i = 140 */
    0xa.513fd6dac0f98p-4 , /* i = 141 */
    0xa.4a9cf1e9479d8p-4 , /* i = 142 */
    0xa.4402910c48a6p-4 , /* i = 143 */
    0xa.3d70a45815f3p-4 , /* i = 144 */
    0xa.36e71a426d478p-4 , /* i = 145 */
    0xa.3065e40f7f738p-4 , /* i = 146 */
    0xa.29ecf17ae6808p-4 , /* i = 147 */
    0xa.237c32f86bb48p-4 , /* i = 148 */
    0xa.1d1398721a67p-4 , /* i = 149 */
    0xa.16b31321506f8p-4 , /* i = 150 */
    0xa.105a934ad1dfp-4 , /* i = 151 */
    0xa.0a0a0a27d3b2p-4 , /* i = 152 */
    0xa.03c168aa3b248p-4 , /* i = 153 */
    0x9.fd80a00b0cf4p-4 , /* i = 154 */
    0x9.f747a171989a8p-4 , /* i = 155 */
    0x9.f1165eae09dd8p-4 , /* i = 156 */
    0x9.eaecc8f858078p-4 , /* i = 157 */
    0x9.e4cad254296dp-4 , /* i = 158 */
    0x9.deb06ca4c5c28p-4 , /* i = 159 */
    0x9.d89d8a3b9e938p-4 , /* i = 160 */
    0x9.d2921c5924838p-4 , /* i = 161 */
    0x9.cc8e163a07848p-4 , /* i = 162 */
    0x9.c69169d6247e8p-4 , /* i = 163 */
    0x9.c09c0a05eabep-4 , /* i = 164 */
    0x9.baade908a137p-4 , /* i = 165 */
    0x9.b4c6fa199a29p-4 , /* i = 166 */
    0x9.aee72ff5029f8p-4 , /* i = 167 */
    0x9.a90e7dc344fa8p-4 , /* i = 168 */
    0x9.a33cd693373a8p-4 , /* i = 169 */
    0x9.9d722db01db38p-4 , /* i = 170 */
    0x9.97ae76cd29908p-4 , /* i = 171 */
    0x9.91f1a546d09ep-4 , /* i = 172 */
    0x9.8c3bac9976388p-4 , /* i = 173 */
    0x9.868c80bc0d758p-4 , /* i = 174 */
    0x9.80e415814bc78p-4 , /* i = 175 */
    0x9.7b425f44e482p-4 , /* i = 176 */
    0x9.75a75119b0188p-4 , /* i = 177 */
    0x9.7012e04d04eep-4 , /* i = 178 */
    0x9.6a8500ad947a8p-4 , /* i = 179 */
    0x9.64fda6f7298p-4 , /* i = 180 */
    0x9.5f7cc74619ee8p-4 , /* i = 181 */
    0x9.5a0256a5e60fp-4 , /* i = 182 */
    0x9.548e4991b3d28p-4 , /* i = 183 */
    0x9.4f209550885dp-4 , /* i = 184 */
    0x9.49b92df32012p-4 , /* i = 185 */
    0x9.44580955f5fp-4 , /* i = 186 */
    0x9.3efd1c653b32p-4 , /* i = 187 */
    0x9.39a85c86f3e6p-4 , /* i = 188 */
    0x9.3459be74b5d8p-4 , /* i = 189 */
    0x9.2f11386c06808p-4 , /* i = 190 */
    0x9.29cebf673a9b8p-4 , /* i = 191 */
    0x9.2492495944a78p-4 , /* i = 192 */
    0x9.1f5bcbaf50a68p-4 , /* i = 193 */
    0x9.1a2b3c6b86c5p-4 , /* i = 194 */
    0x9.1500915dd03bp-4 , /* i = 195 */
    0x9.0fdbc0b78c788p-4 , /* i = 196 */
    0x9.0abcc046ddc38p-4 , /* i = 197 */
    0x9.05a38666b20a8p-4 , /* i = 198 */
    0x9.009009216961p-4 , /* i = 199 */
    0x8.fb823f2d0837p-4 , /* i = 200 */
    0x8.f67a1e638718p-4 , /* i = 201 */
    0x8.f1779dcda1008p-4 , /* i = 202 */
    0x8.ec7ab3ab095b8p-4 , /* i = 203 */
    0x8.e7835709ccaa8p-4 , /* i = 204 */
    0x8.e2917e20f22ep-4 , /* i = 205 */
    0x8.dda520514ddb8p-4 , /* i = 206 */
    0x8.d8be341110b9p-4 , /* i = 207 */
    0x8.d3dcb0a043fdp-4 , /* i = 208 */
    0x8.cf008d0b975fp-4 , /* i = 209 */
    0x8.ca29c077c645p-4 , /* i = 210 */
    0x8.c55841e944398p-4 , /* i = 211 */
    0x8.c08c08e4ead7p-4 , /* i = 212 */
    0x8.bbc50caae73b8p-4 , /* i = 213 */
    0x8.b70344ae17348p-4 , /* i = 214 */
    0x8.b246a8925bap-4 , /* i = 215 */
    0x8.ad8f30113cb1p-4 , /* i = 216 */
    0x8.a8dcd21489ea8p-4 , /* i = 217 */
    0x8.a42f871c5ddcp-4 , /* i = 218 */
    0x8.9f8746b93b158p-4 , /* i = 219 */
    0x8.9ae408df9673p-4 , /* i = 220 */
    0x8.9645c50422aap-4 , /* i = 221 */
    0x8.91ac73cb1dd9p-4 , /* i = 222 */
    0x8.8d180ce8c31e8p-4 , /* i = 223 */
    0x8.8888893bdfeep-4 , /* i = 224 */
    0x8.83fddf1397e68p-4 , /* i = 225 */
    0x8.7f7808a9d20fp-4 , /* i = 226 */
    0x8.7af6fd76fa1d8p-4 , /* i = 227 */
    0x8.767ab62d7e5p-4 , /* i = 228 */
    0x8.72032ad456a5p-4 , /* i = 229 */
    0x8.6d90545549afp-4 , /* i = 230 */
    0x8.69222b37b392p-4 , /* i = 231 */
    0x8.64b8a82f972a8p-4 , /* i = 232 */
    0x8.6053c3512bafp-4 , /* i = 233 */
    0x8.5bf3762c7bea8p-4 , /* i = 234 */
    0x8.5797b935276cp-4 , /* i = 235 */
    0x8.534085709d9cp-4 , /* i = 236 */
    0x8.4eedd36684efp-4 , /* i = 237 */
    0x8.4a9f9c8f8f05p-4 , /* i = 238 */
    0x8.4655d9c986848p-4 , /* i = 239 */
    0x8.42108449e9908p-4 , /* i = 240 */
    0x8.3dcf94f276ad8p-4 , /* i = 241 */
    0x8.399305437ade8p-4 , /* i = 242 */
    0x8.355ace5374078p-4 , /* i = 243 */
    0x8.3126e9b086528p-4 , /* i = 244 */
    0x8.2cf7504c0cd7p-4 , /* i = 245 */
    0x8.28cbfc1291bd8p-4 , /* i = 246 */
    0x8.24a4e626e664p-4 , /* i = 247 */
    0x8.20820868006p-4 , /* i = 248 */
    0x8.1c635bc683fep-4 , /* i = 249 */
    0x8.1848dab0a0898p-4 , /* i = 250 */
    0x8.14327e57aa74p-4 , /* i = 251 */
    0x8.102040ba5f118p-4 , /* i = 252 */
    0x8.0c121b426713p-4 , /* i = 253 */
    0x8.0808082f8c4dp-4 , /* i = 254 */
    0x8.040201157cb4p-4 , /* i = 255 */
};

static const double table_log_alpha_i_modified[256] = 
{
    0x0p+3,    /* i = 0 */
    0x1.ff00a36c68bap-9 , /* i = 1 */
    0x1.fe02a600fba6ap-8 , /* i = 2 */
    0x1.7dc4747e38f8cp-7 , /* i = 3 */
    0x1.fc0a895299517p-7 , /* i = 4 */
    0x1.3cea437af59e3p-6 , /* i = 5 */
    0x1.7b91afa50fa62p-6 , /* i = 6 */
    0x1.b9fc021e8fca1p-6 , /* i = 7 */
    0x1.f829ae89ff9d9p-6 , /* i = 8 */
    0x1.1b0d981def938p-5 , /* i = 9 */
    0x1.39e87b06794f1p-5 , /* i = 10 */
    0x1.58a5babffb3f2p-5 , /* i = 11 */
    0x1.77458e6fe8fd4p-5 , /* i = 12 */
    0x1.95c830b9085d4p-5 , /* i = 13 */
    0x1.b42dd663f0133p-5 , /* i = 14 */
    0x1.d276b8357d102p-5 , /* i = 15 */
    0x1.f0a30a6e3398p-5 , /* i = 16 */
    0x1.07598319fc443p-4 , /* i = 17 */
    0x1.16536eabc6782p-4 , /* i = 18 */
    0x1.253f62b4b0549p-4 , /* i = 19 */
    0x1.341d790ee8c1dp-4 , /* i = 20 */
    0x1.42edcbbba1ca5p-4 , /* i = 21 */
    0x1.51b073b52e4dfp-4 , /* i = 22 */
    0x1.60658a5b70a5cp-4 , /* i = 23 */
    0x1.6f0d28265f998p-4 , /* i = 24 */
    0x1.7da766abfad8dp-4 , /* i = 25 */
    0x1.8c345d1a64a33p-4 , /* i = 26 */
    0x1.9ab4243456616p-4 , /* i = 27 */
    0x1.a926d32f3403fp-4 , /* i = 28 */
    0x1.b78c827c20b0fp-4 , /* i = 29 */
    0x1.c5e548ad99b34p-4 , /* i = 30 */
    0x1.d4313d48faf99p-4 , /* i = 31 */
    0x1.e270767859f47p-4 , /* i = 32 */
    0x1.f0a30bc7581c3p-4 , /* i = 33 */
    0x1.fec912e729dd5p-4 , /* i = 34 */
    0x1.0671511693a39p-3 , /* i = 35 */
    0x1.0d77e79bf57cp-3 , /* i = 36 */
    0x1.14785829189dbp-3 , /* i = 37 */
    0x1.1b72ad28b40a3p-3 , /* i = 38 */
    0x1.2266f17674dd8p-3 , /* i = 39 */
    0x1.29552f3a7c29ap-3 , /* i = 40 */
    0x1.303d717b37407p-3 , /* i = 41 */
    0x1.371fc1e12239ap-3 , /* i = 42 */
    0x1.3dfc2afbd88c2p-3 , /* i = 43 */
    0x1.44d2b69efab65p-3 , /* i = 44 */
    0x1.4ba36f23e0bf5p-3 , /* i = 45 */
    0x1.526e5e15f1282p-3 , /* i = 46 */
    0x1.59338d8bdd42p-3 , /* i = 47 */
    0x1.5ff3069339d91p-3 , /* i = 48 */
    0x1.66acd4163269ep-3 , /* i = 49 */
    0x1.6d60fe5778b7ep-3 , /* i = 50 */
    0x1.740f8f3bc1c4ap-3 , /* i = 51 */
    0x1.7ab890009c828p-3 , /* i = 52 */
    0x1.815c0a00b43ap-3 , /* i = 53 */
    0x1.87fa0625269e1p-3 , /* i = 54 */
    0x1.8e928dd3f40f4p-3 , /* i = 55 */
    0x1.9525a99d53067p-3 , /* i = 56 */
    0x1.9bb362d05f49fp-3 , /* i = 57 */
    0x1.a23bc1d47e855p-3 , /* i = 58 */
    0x1.a8becfb7c4dd3p-3 , /* i = 59 */
    0x1.af3c94bd22568p-3 , /* i = 60 */
    0x1.b5b519d5a5494p-3 , /* i = 61 */
    0x1.bc28671b93ddap-3 , /* i = 62 */
    0x1.c29685490551cp-3 , /* i = 63 */
    0x1.c8ff7be726545p-3 , /* i = 64 */
    0x1.cf6354c47198fp-3 , /* i = 65 */
    0x1.d5c216934bd09p-3 , /* i = 66 */
    0x1.dc1bc9ee3c487p-3 , /* i = 67 */
    0x1.e27076bb3896ap-3 , /* i = 68 */
    0x1.e8c0250f84031p-3 , /* i = 69 */
    0x1.ef0adc96e665dp-3 , /* i = 70 */
    0x1.f550a5557c423p-3 , /* i = 71 */
    0x1.fb918690c3443p-3 , /* i = 72 */
    0x1.00e6c44d38621p-2 , /* i = 73 */
    0x1.0402593f5af11p-2 , /* i = 74 */
    0x1.071b85f20b49ap-2 , /* i = 75 */
    0x1.0a324e1edbe5fp-2 , /* i = 76 */
    0x1.0d46b56bc58ap-2 , /* i = 77 */
    0x1.1058bf861b0bbp-2 , /* i = 78 */
    0x1.1368701d18e88p-2 , /* i = 79 */
    0x1.1675ca7e72f98p-2 , /* i = 80 */
    0x1.1980d2d3dc95cp-2 , /* i = 81 */
    0x1.1c898c113f221p-2 , /* i = 82 */
    0x1.1f8ff9dc8f3a2p-2 , /* i = 83 */
    0x1.22941fa1dcb5ap-2 , /* i = 84 */
    0x1.25960102d8a9p-2 , /* i = 85 */
    0x1.2895a13109876p-2 , /* i = 86 */
    0x1.2b9303a4cdaa4p-2 , /* i = 87 */
    0x1.2e8e2b97f1a63p-2 , /* i = 88 */
    0x1.31871c887ccb8p-2 , /* i = 89 */
    0x1.347dd99aa80d6p-2 , /* i = 90 */
    0x1.3772661d611acp-2 , /* i = 91 */
    0x1.3a64c54c5d804p-2 , /* i = 92 */
    0x1.3d54fa521613ep-2 , /* i = 93 */
    0x1.40430854f9e8fp-2 , /* i = 94 */
    0x1.432ef2993663ap-2 , /* i = 95 */
    0x1.4618bc0e536a5p-2 , /* i = 96 */
    0x1.490067f402f24p-2 , /* i = 97 */
    0x1.4be5f94dffc92p-2 , /* i = 98 */
    0x1.4ec9731cc7e9ap-2 , /* i = 99 */
    0x1.51aad856abf32p-2 , /* i = 100 */
    0x1.548a2c35ca1c2p-2 , /* i = 101 */
    0x1.5767715f7ce68p-2 , /* i = 102 */
    0x1.5a42ab009ea11p-2 , /* i = 103 */
    0x1.5d1bdbc9097f2p-2 , /* i = 104 */
    0x1.5ff306fc67046p-2 , /* i = 105 */
    0x1.62c82f231c71dp-2 , /* i = 106 */
    0x1.659b572bf381p-2 , /* i = 107 */
    0x1.686c81d3d0fd9p-2 , /* i = 108 */
    0x1.6b3bb21a5eb4bp-2 , /* i = 109 */
    0x1.6e08ea8c51734p-2 , /* i = 110 */
    0x1.70d42e20d9228p-2 , /* i = 111 */
    0x1.739d7f57c4c22p-2 , /* i = 112 */
    0x1.7664e119294dbp-2 , /* i = 113 */
    0x1.792a55e772551p-2 , /* i = 114 */
    0x1.7bede093f9997p-2 , /* i = 115 */
    0x1.7eaf83a39f588p-2 , /* i = 116 */
    0x1.816f41d40474ap-2 , /* i = 117 */
    0x1.842d1d8f155b2p-2 , /* i = 118 */
    0x1.86e919a1f5a64p-2 , /* i = 119 */
    0x1.89a33847afbdfp-2 , /* i = 120 */
    0x1.8c5b7c80ec503p-2 , /* i = 121 */
    0x1.8f11e863c567ap-2 , /* i = 122 */
    0x1.91c67eaf0f42ap-2 , /* i = 123 */
    0x1.947941a70f092p-2 , /* i = 124 */
    0x1.972a3408baf15p-2 , /* i = 125 */
    0x1.99d9580ae88c2p-2 , /* i = 126 */
    0x1.9c86b026d67e4p-2 , /* i = 127 */
    0x1.9f323ea2bc50bp-2 , /* i = 128 */
    0x1.a1dc063d6ea88p-2 , /* i = 129 */
    0x1.a484090d16a37p-2 , /* i = 130 */
    0x1.a72a49599df8bp-2 , /* i = 131 */
    0x1.a9cec996fa657p-2 , /* i = 132 */
    0x1.ac718c15b834fp-2 , /* i = 133 */
    0x1.af1293118a90ap-2 , /* i = 134 */
    0x1.b1b1e0e0f741bp-2 , /* i = 135 */
    0x1.b44f77a860969p-2 , /* i = 136 */
    0x1.b6eb59ca3f8ap-2 , /* i = 137 */
    0x1.b985895786876p-2 , /* i = 138 */
    0x1.bc1e08a2e4f54p-2 , /* i = 139 */
    0x1.beb4d9d2a43f4p-2 , /* i = 140 */
    0x1.c149ff050fb1bp-2 , /* i = 141 */
    0x1.c3dd7a76af7f2p-2 , /* i = 142 */
    0x1.c66f4e3b6db0dp-2 , /* i = 143 */
    0x1.c8ff7c47410f3p-2 , /* i = 144 */
    0x1.cb8e073c54616p-2 , /* i = 145 */
    0x1.ce1af0b049b5ep-2 , /* i = 146 */
    0x1.d0a63ade03af3p-2 , /* i = 147 */
    0x1.d32fe7c40c3f2p-2 , /* i = 148 */
    0x1.d5b7f9a30b77p-2 , /* i = 149 */
    0x1.d83e7242ed8fbp-2 , /* i = 150 */
    0x1.dac353d7c8aaep-2 , /* i = 151 */
    0x1.dd46a0403df12p-2 , /* i = 152 */
    0x1.dfc859826f3cp-2 , /* i = 153 */
    0x1.e24881935a12bp-2 , /* i = 154 */
    0x1.e4c71a7a2f704p-2 , /* i = 155 */
    0x1.e744260562806p-2 , /* i = 156 */
    0x1.e9bfa64b5c9e5p-2 , /* i = 157 */
    0x1.ec399d1b620c7p-2 , /* i = 158 */
    0x1.eeb20c5c46b9bp-2 , /* i = 159 */
    0x1.f128f5d2b802ep-2 , /* i = 160 */
    0x1.f39e5bbcc4c73p-2 , /* i = 161 */
    0x1.f6123f9452793p-2 , /* i = 162 */
    0x1.f884a3618a7d6p-2 , /* i = 163 */
    0x1.faf588db21eap-2 , /* i = 164 */
    0x1.fd64f200950f4p-2 , /* i = 165 */
    0x1.ffd2e073f241ap-2 , /* i = 166 */
    0x1.011fab0aa50dfp-1 , /* i = 167 */
    0x1.02552a50efc81p-1 , /* i = 168 */
    0x1.0389eef5995dbp-1 , /* i = 169 */
    0x1.04bdf9d9b03f9p-1 , /* i = 170 */
    0x1.05f14bcd5dcf4p-1 , /* i = 171 */
    0x1.0723e5b7815bp-1 , /* i = 172 */
    0x1.0855c87d0eed5p-1 , /* i = 173 */
    0x1.0986f4edf77eap-1 , /* i = 174 */
    0x1.0ab76be64b7a8p-1 , /* i = 175 */
    0x1.0be72e29ca74cp-1 , /* i = 176 */
    0x1.0d163cc60f105p-1 , /* i = 177 */
    0x1.0e449854ca7d5p-1 , /* i = 178 */
    0x1.0f7241c4d62e5p-1 , /* i = 179 */
    0x1.109f39d73673bp-1 , /* i = 180 */
    0x1.11cb81732787p-1 , /* i = 181 */
    0x1.12f7195144081p-1 , /* i = 182 */
    0x1.1422024d28689p-1 , /* i = 183 */
    0x1.154c3d1b0012ap-1 , /* i = 184 */
    0x1.1675cab5c019p-1 , /* i = 185 */
    0x1.179eabb9bc46p-1 , /* i = 186 */
    0x1.18c6e0faf7427p-1 , /* i = 187 */
    0x1.19ee6b373ab69p-1 , /* i = 188 */
    0x1.1b154b55be26ep-1 , /* i = 189 */
    0x1.1c3b81ed8d57dp-1 , /* i = 190 */
    0x1.1d610fdfcf158p-1 , /* i = 191 */
    0x1.1e85f5db7d086p-1 , /* i = 192 */
    0x1.1faa34b0a1064p-1 , /* i = 193 */
    0x1.20cdcd128a5a1p-1 , /* i = 194 */
    0x1.21f0bfc35333cp-1 , /* i = 195 */
    0x1.23130d7369739p-1 , /* i = 196 */
    0x1.2434b6ecd6c5ep-1 , /* i = 197 */
    0x1.2555bcde4b84ep-1 , /* i = 198 */
    0x1.2676200bf6c49p-1 , /* i = 199 */
    0x1.2795e117943d7p-1 , /* i = 200 */
    0x1.28b500d76af2fp-1 , /* i = 201 */
    0x1.29d37fe1eec75p-1 , /* i = 202 */
    0x1.2af15efdeed1fp-1 , /* i = 203 */
    0x1.2c0e9ec795707p-1 , /* i = 204 */
    0x1.2d2b400ec3464p-1 , /* i = 205 */
    0x1.2e474363e8d64p-1 , /* i = 206 */
    0x1.2f62a98fadac9p-1 , /* i = 207 */
    0x1.307d7330a1a49p-1 , /* i = 208 */
    0x1.3197a0f43e283p-1 , /* i = 209 */
    0x1.32b13385e0f9ep-1 , /* i = 210 */
    0x1.33ca2b9b97aabp-1 , /* i = 211 */
    0x1.34e289d17e731p-1 , /* i = 212 */
    0x1.35fa4ed693b6p-1 , /* i = 213 */
    0x1.37117b5180a81p-1 , /* i = 214 */
    0x1.38280fe0dec58p-1 , /* i = 215 */
    0x1.393e0d2169a4dp-1 , /* i = 216 */
    0x1.3a5373e2ed855p-1 , /* i = 217 */
    0x1.3b68449aaf09p-1 , /* i = 218 */
    0x1.3c7c7ff83d0c4p-1 , /* i = 219 */
    0x1.3d9026971df7fp-1 , /* i = 220 */
    0x1.3ea339339c62ep-1 , /* i = 221 */
    0x1.3fb5b8466ebccp-1 , /* i = 222 */
    0x1.40c7a4831d2bdp-1 , /* i = 223 */
    0x1.41d8fe5a5eaf2p-1 , /* i = 224 */
    0x1.42e9c6d97dfefp-1 , /* i = 225 */
    0x1.43f9fe25a5367p-1 , /* i = 226 */
    0x1.4509a50c4c73fp-1 , /* i = 227 */
    0x1.4618bc1405924p-1 , /* i = 228 */
    0x1.472743eeb1a6p-1 , /* i = 229 */
    0x1.48353d1b6b51bp-1 , /* i = 230 */
    0x1.4942a83350ed2p-1 , /* i = 231 */
    0x1.4a4f85c7ace69p-1 , /* i = 232 */
    0x1.4b5bd692ac991p-1 , /* i = 233 */
    0x1.4c679af6aa3dcp-1 , /* i = 234 */
    0x1.4d72d39c817c8p-1 , /* i = 235 */
    0x1.4e7d810ce7de1p-1 , /* i = 236 */
    0x1.4f87a3f174ae8p-1 , /* i = 237 */
    0x1.50913cbc75c74p-1 , /* i = 238 */
    0x1.519a4c080de12p-1 , /* i = 239 */
    0x1.52a2d25bd5c9bp-1 , /* i = 240 */
    0x1.53aad05642aafp-1 , /* i = 241 */
    0x1.54b24671f0a8ep-1 , /* i = 242 */
    0x1.55b93545ab578p-1 , /* i = 243 */
    0x1.56bf9d4da66e7p-1 , /* i = 244 */
    0x1.57c57f2ed474fp-1 , /* i = 245 */
    0x1.58cadb534a5f6p-1 , /* i = 246 */
    0x1.59cfb258e0b21p-1 , /* i = 247 */
    0x1.5ad404b1c1db3p-1 , /* i = 248 */
    0x1.5bd7d30d1e7b3p-1 , /* i = 249 */
    0x1.5cdb1dbe9dcc5p-1 , /* i = 250 */
    0x1.5ddde56a55bfcp-1 , /* i = 251 */
    0x1.5ee02a8406d28p-1 , /* i = 252 */
    0x1.5fe1eda6b7b37p-1 , /* i = 253 */
    0x1.60e32f3aa15d6p-1 , /* i = 254 */
    0x1.61e3efd509c9p-1 , /* i = 255 */
};

/* table_alpha_i
alpha_i = {2^k}/{2^k+i}  with 0<=i<2*k with k=8.
alpha_i = {h, m, l} to calculate in Triple-Double precision */
static const double table_alpha_i_triple[256][3] = {
//{h,m,l} for each alpha_i
{0x1p+0, 0x0p+0, 0x0p+0} /* i= 0 */
,{0x1.fe01fe01fe02p-1, -0x1.fe01fe01fe02p-57, 0x1.fe01fe01fe02p-113} /* i= 1 */
,{0x1.fc07f01fc07fp-1, 0x1.fc07f01fc07fp-57, 0x1.fc07f01fc07fp-113} /* i= 2 */
,{0x1.fa11caa01fa12p-1, -0x1.aaff02f71aaffp-56, -0x1.7b8d57f817b8dp-115} /* i= 3 */
,{0x1.f81f81f81f82p-1, -0x1.f81f81f81f82p-55, 0x1.f81f81f81f82p-109} /* i= 4 */
,{0x1.f6310aca0dbb5p-1, 0x1.d2e19807d8c43p-55, -0x1.35f244a8b479ap-109} /* i= 5 */
,{0x1.f44659e4a4271p-1, 0x1.5fc17734c36b8p-55, -0x1.38abf82ee6987p-109} /* i= 6 */
,{0x1.f25f644230ab5p-1, 0x1.94ed8175c78b3p-58, 0x1.a4807c97d9109p-114} /* i= 7 */
,{0x1.f07c1f07c1f08p-1, -0x1.f07c1f07c1f08p-56, 0x1.f07c1f07c1f08p-111} /* i= 8 */
,{0x1.ee9c7f8458e02p-1, -0x1.163807ba71fe1p-57, -0x1.63807ba71fe11p-113} /* i= 9 */
,{0x1.ecc07b301eccp-1, 0x1.ecc07b301eccp-55, 0x1.ecc07b301eccp-109} /* i= 10 */
,{0x1.eae807aba01ebp-1, -0x1.7f8545fe1518p-57, 0x1.eae807aba01ebp-111} /* i= 11 */
,{0x1.e9131abf0b767p-1, 0x1.503d226357e17p-56, -0x1.31abf0b7672ap-112} /* i= 12 */
,{0x1.e741aa59750e4p-1, 0x1.9b1f67bb7ac41p-55, -0x1.251d8079d06a9p-109} /* i= 13 */
,{0x1.e573ac901e574p-1, -0x1.4dbf86a314dcp-55, 0x1.e573ac901e574p-109} /* i= 14 */
,{0x1.e3a9179dc1a73p-1, 0x1.fa5504b926bb1p-56, -0x1.66f77f8715ba2p-110} /* i= 15 */
,{0x1.e1e1e1e1e1e1ep-1, 0x1.e1e1e1e1e1e1ep-57, 0x1.e1e1e1e1e1e1ep-113} /* i= 16 */
,{0x1.e01e01e01e01ep-1, 0x1.e01e01e01e01ep-61, 0x1.e01e01e01e01ep-121} /* i= 17 */
,{0x1.de5d6e3f8868ap-1, 0x1.1c077975b8fe2p-55, 0x1.a291c077975b9p-111} /* i= 18 */
,{0x1.dca01dca01dcap-1, 0x1.dca01dca01dcap-61, 0x1.dca01dca01dcap-121} /* i= 19 */
,{0x1.dae6076b981dbp-1, -0x1.9f89467e251ap-57, 0x1.dae6076b981dbp-111} /* i= 20 */
,{0x1.d92f2231e7f8ap-1, -0x1.2f2231e7f89b4p-55, -0x1.bb9c300ec9791p-110} /* i= 21 */
,{0x1.d77b654b82c34p-1, -0x1.ba03aef6ca97p-55, -0x1.619c8bf8a2127p-109} /* i= 22 */
,{0x1.d5cac807572b2p-1, 0x1.d5cac807572b2p-61, 0x1.d5cac807572b2p-121} /* i= 23 */
,{0x1.d41d41d41d41dp-1, 0x1.075075075075p-55, 0x1.d41d41d41d41dp-109} /* i= 24 */
,{0x1.d272ca3fc5b1ap-1, 0x1.ae01d272ca3fcp-55, 0x1.6c69ae01d272dp-109} /* i= 25 */
,{0x1.d0cb58f6ec074p-1, 0x1.96b1edd80e866p-56, -0x1.4e1227f179a54p-110} /* i= 26 */
,{0x1.cf26e5c44bfc6p-1, 0x1.b2347768073cap-57, -0x1.1a3bb4039e4ddp-111} /* i= 27 */
,{0x1.cd85689039b0bp-1, -0x1.76fc64f52edf9p-56, 0x1.b0ad12073615ap-111} /* i= 28 */
,{0x1.cbe6d9601cbe7p-1, -0x1.34ff1a0c934ffp-56, -0x1.a0c934ff1a0c9p-112} /* i= 29 */
,{0x1.ca4b3055ee191p-1, 0x1.ca4b3055ee191p-61, 0x1.ca4b3055ee191p-121} /* i= 30 */
,{0x1.c8b265afb8a42p-1, 0x1.c8b265afb8a42p-61, 0x1.c8b265afb8a42p-121} /* i= 31 */
,{0x1.c71c71c71c71cp-1, 0x1.c71c71c71c71cp-55, 0x1.c71c71c71c71cp-109} /* i= 32 */
,{0x1.c5894d10d4986p-1, -0x1.f00e2c4a6886ap-56, -0x1.30b83fc74ed66p-110} /* i= 33 */
,{0x1.c3f8f01c3f8fp-1, 0x1.c3f8f01c3f8fp-57, 0x1.c3f8f01c3f8fp-113} /* i= 34 */
,{0x1.c26b5392ea01cp-1, 0x1.35a9c97500e13p-56, 0x1.6a725d40384d7p-110} /* i= 35 */
,{0x1.c0e070381c0ep-1, 0x1.c0e070381c0ep-55, 0x1.c0e070381c0ep-109} /* i= 36 */
,{0x1.bf583ee868d8bp-1, -0x1.41876d370b5bcp-57, 0x1.338cab3fc815p-112} /* i= 37 */
,{0x1.bdd2b899406f7p-1, 0x1.2b899406f74aep-55, 0x1.3280dee95c4cap-110} /* i= 38 */
,{0x1.bc4fd65883e7bp-1, 0x1.d1239464aa169p-56, 0x1.bc4fd65883e7bp-117} /* i= 39 */
,{0x1.bacf914c1badp-1, -0x1.bacf914c1badp-55, 0x1.bacf914c1badp-109} /* i= 40 */
,{0x1.b951e2b18ff23p-1, 0x1.5c3a9ce01b952p-55, -0x1.d4e700dca8f16p-111} /* i= 41 */
,{0x1.b7d6c3dda338bp-1, 0x1.579fc90527845p-56, -0x1.19c59579fc905p-110} /* i= 42 */
,{0x1.b65e2e3beee05p-1, 0x1.18d4559e6507bp-56, 0x1.29f4036cbc5c7p-110} /* i= 43 */
,{0x1.b4e81b4e81b4fp-1, -0x1.f92c5f92c5f93p-55, 0x1.d0369d0369d03p-110} /* i= 44 */
,{0x1.b37484ad806cep-1, -0x1.6f6a4ff2645bep-56, 0x1.5b00d9ba4256cp-110} /* i= 45 */
,{0x1.b2036406c80d9p-1, 0x1.b2036406c80d9p-61, 0x1.b2036406c80d9p-121} /* i= 46 */
,{0x1.b094b31d922a4p-1, -0x1.7a821cb9dfe4fp-57, -0x1.ad3389b75705fp-111} /* i= 47 */
,{0x1.af286bca1af28p-1, 0x1.af286bca1af28p-55, 0x1.af286bca1af28p-109} /* i= 48 */
,{0x1.adbe87f94905ep-1, 0x1.adbe87f94905ep-61, 0x1.adbe87f94905ep-121} /* i= 49 */
,{0x1.ac5701ac5701bp-1, -0x1.d47f29d47f29dp-56, -0x1.1fca751fca752p-110} /* i= 50 */
,{0x1.aaf1d2f87ebfdp-1, -0x1.578e97c3f5fe5p-55, -0x1.438b41e0500d5p-109} /* i= 51 */
,{0x1.a98ef606a63bep-1, -0x1.f959c427e5671p-55, -0x1.3f2b3884fcacep-112} /* i= 52 */
,{0x1.a82e65130e159p-1, -0x1.6937821239fe5p-55, -0x1.f466bb3c7a9d7p-109} /* i= 53 */
,{0x1.a6d01a6d01a6dp-1, 0x1.a6d01a6d01a6dp-61, 0x1.a6d01a6d01a6dp-121} /* i= 54 */
,{0x1.a574107688a4ap-1, 0x1.566e4d604f05cp-57, 0x1.8b1ccf6f201a5p-112} /* i= 55 */
,{0x1.a41a41a41a41ap-1, 0x1.069069069069p-55, 0x1.a41a41a41a41ap-109} /* i= 56 */
,{0x1.a2c2a87c51cap-1, 0x1.3a11fe5d3d578p-55, 0x1.d71afd8bdc034p-110} /* i= 57 */
,{0x1.a16d3f97a4b02p-1, -0x1.7a4b01a16d3f9p-55, -0x1.e92c0685b4fe6p-109} /* i= 58 */
,{0x1.a01a01a01a01ap-1, 0x1.a01a01a01a01ap-61, 0x1.a01a01a01a01ap-121} /* i= 59 */
,{0x1.9ec8e951033d9p-1, 0x1.d2a2067b23a54p-57, 0x1.033d91d2a2068p-111} /* i= 60 */
,{0x1.9d79f176b682dp-1, 0x1.cab347dfb2792p-56, 0x1.5cdee3bc29fe6p-111} /* i= 61 */
,{0x1.9c2d14ee4a102p-1, -0x1.8f4bac46d7bfap-55, 0x1.c2d14ee4a101ap-109} /* i= 62 */
,{0x1.9ae24ea5510dap-1, 0x1.20e71f4c3cfd9p-55, 0x1.eb2282019ae25p-109} /* i= 63 */
,{0x1.999999999999ap-1, -0x1.999999999999ap-55, 0x1.999999999999ap-109} /* i= 64 */
,{0x1.9852f0d8ec0ffp-1, 0x1.9eb43c9c4fc03p-56, 0x1.852f0d8ec0ff3p-111} /* i= 65 */
,{0x1.970e4f80cb872p-1, 0x1.f01970e4f80ccp-55, -0x1.e360fe68f1b08p-109} /* i= 66 */
,{0x1.95cbb0be377aep-1, -0x1.b57f9a8d13d07p-55, -0x1.10a4dabfcd469p-110} /* i= 67 */
,{0x1.948b0fcd6e9ep-1, 0x1.948b0fcd6e9ep-55, 0x1.948b0fcd6e9ep-109} /* i= 68 */
,{0x1.934c67f9b2ce6p-1, 0x1.934c67f9b2ce6p-61, 0x1.934c67f9b2ce6p-121} /* i= 69 */
,{0x1.920fb49d0e229p-1, -0x1.533d406483ed2p-56, -0x1.d0e228d59857fp-110} /* i= 70 */
,{0x1.90d4f120190d5p-1, -0x1.dbfcde561dbfdp-58, 0x1.0d4f120190d4fp-113} /* i= 71 */
,{0x1.8f9c18f9c18fap-1, -0x1.f3831f3831f38p-56, -0x1.8f9c18f9c18fap-111} /* i= 72 */
,{0x1.8e6527af1373fp-1, 0x1.c031cca4f5e27p-59, -0x1.81f1fe719ad85p-115} /* i= 73 */
,{0x1.8d3018d3018d3p-1, 0x1.8d3018d3018d3p-61, 0x1.8d3018d3018d3p-121} /* i= 74 */
,{0x1.8bfce8062ff3ap-1, 0x1.8bfce8062ff3ap-61, 0x1.8bfce8062ff3ap-121} /* i= 75 */
,{0x1.8acb90f6bf3aap-1, -0x1.721ed7e75346fp-55, -0x1.2818acb90f6bfp-112} /* i= 76 */
,{0x1.899c0f601899cp-1, 0x1.ec0313381ec03p-58, 0x1.3381ec0313382p-114} /* i= 77 */
,{0x1.886e5f0abb04ap-1, -0x1.ad38b7f3bc8dp-55, -0x1.ea89f6cd69c5cp-109} /* i= 78 */
,{0x1.87427bcc092b9p-1, -0x1.1937c8faa6975p-57, 0x1.4a20187427bccp-113} /* i= 79 */
,{0x1.8618618618618p-1, 0x1.8618618618618p-55, 0x1.8618618618618p-109} /* i= 80 */
,{0x1.84f00c2780614p-1, -0x1.fe7b0ff3d87fap-56, 0x1.3c0309e0184fp-112} /* i= 81 */
,{0x1.83c977ab2beddp-1, 0x1.4731fcf86d10bp-56, -0x1.95f6e94731fdp-110} /* i= 82 */
,{0x1.82a4a0182a4ap-1, 0x1.82a4a0182a4ap-57, 0x1.82a4a0182a4ap-113} /* i= 83 */
,{0x1.8181818181818p-1, 0x1.8181818181818p-57, 0x1.8181818181818p-113} /* i= 84 */
,{0x1.8060180601806p-1, 0x1.8060180601806p-61, 0x1.8060180601806p-121} /* i= 85 */
,{0x1.7f405fd017f4p-1, 0x1.7f405fd017f4p-55, 0x1.7f405fd017f4p-109} /* i= 86 */
,{0x1.7e225515a4f1dp-1, 0x1.b9d7b26106b7ap-57, -0x1.6047a66ff40efp-111} /* i= 87 */
,{0x1.7d05f417d05f4p-1, 0x1.7d05f417d05f4p-57, 0x1.7d05f417d05f4p-113} /* i= 88 */
,{0x1.7beb3922e017cp-1, -0x1.4c6dd1fe8414cp-57, -0x1.b747fa10531b7p-111} /* i= 89 */
,{0x1.7ad2208e0ecc3p-1, 0x1.516324fe852dep-55, -0x1.1c1d986a8b192p-112} /* i= 90 */
,{0x1.79baa6bb6398bp-1, 0x1.bd9a30b10f7e2p-55, 0x1.f5abe570e046dp-109} /* i= 91 */
,{0x1.78a4c8178a4c8p-1, 0x1.78a4c8178a4c8p-57, 0x1.78a4c8178a4c8p-113} /* i= 92 */
,{0x1.77908119ac60dp-1, 0x1.a0a44f387b3b7p-56, -0x1.68e4dc0eaba51p-110} /* i= 93 */
,{0x1.767dce434a9b1p-1, 0x1.767dce434a9b1p-61, 0x1.767dce434a9b1p-121} /* i= 94 */
,{0x1.756cac201756dp-1, -0x1.4f7fa2a4d4f8p-55, 0x1.756cac201756dp-109} /* i= 95 */
,{0x1.745d1745d1746p-1, -0x1.745d1745d1746p-56, 0x1.745d1745d1746p-111} /* i= 96 */
,{0x1.734f0c541fe8dp-1, -0x1.3c31507fa32c4p-55, 0x1.8a83fd1961e75p-110} /* i= 97 */
,{0x1.724287f46debcp-1, 0x1.724287f46debcp-59, 0x1.724287f46debcp-117} /* i= 98 */
,{0x1.713786d9c7c09p-1, -0x1.62cb5b9545f3p-55, -0x1.431095fe8ec88p-109} /* i= 99 */
,{0x1.702e05c0b817p-1, 0x1.702e05c0b817p-56, 0x1.702e05c0b817p-111} /* i= 100 */
,{0x1.6f26016f26017p-1, -0x1.b3fd21b3fd21bp-58, -0x1.fe90d9fe90dap-113} /* i= 101 */
,{0x1.6e1f76b4337c7p-1, -0x1.a75461405b87ep-56, 0x1.2979907269d52p-111} /* i= 102 */
,{0x1.6d1a62681c861p-1, -0x1.3f77161b18f55p-59, 0x1.22f1066af6badp-114} /* i= 103 */
,{0x1.6c16c16c16c17p-1, -0x1.f49f49f49f49fp-56, -0x1.27d27d27d27d2p-110} /* i= 104 */
,{0x1.6b1490aa31a3dp-1, -0x1.c5d9b4d4be0ccp-60, -0x1.dc8afddf6127p-115} /* i= 105 */
,{0x1.6a13cd153729p-1, 0x1.0f8ed9cfe95ecp-55, 0x1.975646b7de0e2p-110} /* i= 106 */
,{0x1.691473a88d0cp-1, -0x1.691473a88d0cp-56, 0x1.691473a88d0cp-111} /* i= 107 */
,{0x1.6816816816817p-1, -0x1.fa5fa5fa5fa6p-55, 0x1.6816816816817p-109} /* i= 108 */
,{0x1.6719f3601671ap-1, -0x1.93fd31cc193fdp-58, -0x1.8e60c9fe98e61p-113} /* i= 109 */
,{0x1.661ec6a5122f9p-1, 0x1.661ec6a5122f9p-61, 0x1.661ec6a5122f9p-121} /* i= 110 */
,{0x1.6524f853b4aa3p-1, 0x1.cf2bf20c8e4ccp-56, -0x1.43a9810bdbba4p-110} /* i= 111 */
,{0x1.642c8590b2164p-1, 0x1.642c8590b2164p-56, 0x1.642c8590b2164p-111} /* i= 112 */
,{0x1.63356b88ac0dep-1, 0x1.63356b88ac0dep-61, 0x1.63356b88ac0dep-121} /* i= 113 */
,{0x1.623fa7701624p-1, -0x1.623fa7701624p-55, 0x1.623fa7701624p-109} /* i= 114 */
,{0x1.614b36831ae94p-1, -0x1.5640dccf0211fp-55, -0x1.a38950bbaff4fp-112} /* i= 115 */
,{0x1.6058160581606p-1, -0x1.fa7e9fa7e9fa8p-55, 0x1.6058160581606p-111} /* i= 116 */
,{0x1.5f66434292dfcp-1, -0x1.e32c9c7b89f3ap-57, -0x1.59e8aa3588944p-111} /* i= 117 */
,{0x1.5e75bb8d015e7p-1, 0x1.6ee340579d6eep-55, 0x1.a02bceb771a03p-110} /* i= 118 */
,{0x1.5d867c3ece2a5p-1, 0x1.a485cd7b900afp-56, -0x1.e60f04c756b2ep-111} /* i= 119 */
,{0x1.5c9882b931057p-1, 0x1.310572620ae4cp-56, 0x1.0572620ae4c41p-110} /* i= 120 */
,{0x1.5babcc647fa91p-1, 0x1.4339b8056eaf3p-55, 0x1.91fea454339b8p-111} /* i= 121 */
,{0x1.5ac056b015acp-1, 0x1.5ac056b015acp-55, 0x1.5ac056b015acp-109} /* i= 122 */
,{0x1.59d61f123ccaap-1, 0x1.bb1a57cf5de3ap-56, 0x1.6f73810360975p-112} /* i= 123 */
,{0x1.58ed2308158edp-1, 0x1.1840ac7691841p-56, -0x1.4e25b9efd4e26p-110} /* i= 124 */
,{0x1.580560158056p-1, 0x1.580560158056p-57, 0x1.580560158056p-113} /* i= 125 */
,{0x1.571ed3c506b3ap-1, -0x1.7749b79f7f547p-55, -0x1.2c3af94c65dd2p-112} /* i= 126 */
,{0x1.56397ba7c52e2p-1, -0x1.40d5e3ed48db4p-57, 0x1.966442d73a26cp-112} /* i= 127 */
,{0x1.5555555555555p-1, 0x1.5555555555555p-55, 0x1.5555555555555p-109} /* i= 128 */
,{0x1.54725e6bb82fep-1, 0x1.54725e6bb82fep-61, 0x1.54725e6bb82fep-121} /* i= 129 */
,{0x1.5390948f40febp-1, -0x1.c84a47a07f563p-56, -0x1.ed6e17e02a721p-110} /* i= 130 */
,{0x1.52aff56a8054bp-1, -0x1.00a957fab5403p-55, 0x1.6a8054abfd5aap-109} /* i= 131 */
,{0x1.51d07eae2f815p-1, 0x1.d07eae2f8151dp-57, 0x1.fab8be054742p-115} /* i= 132 */
,{0x1.50f22e111c4c5p-1, 0x1.b79bf81a52ebap-55, -0x1.aa72824da7d0ap-109} /* i= 133 */
,{0x1.5015015015015p-1, 0x1.5015015015015p-61, 0x1.5015015015015p-121} /* i= 134 */
,{0x1.4f38f62dd4c9bp-1, -0x1.eefa1b7fac31cp-55, -0x1.3a4566caf77d1p-110} /* i= 135 */
,{0x1.4e5e0a72f0539p-1, 0x1.e0a72f0539783p-55, -0x1.8d0fac687d634p-109} /* i= 136 */
,{0x1.4d843bedc2c4cp-1, -0x1.c029b0877db86p-55, 0x1.da38053610efbp-109} /* i= 137 */
,{0x1.4cab88725af6ep-1, 0x1.d3d137e0cfeb3p-55, 0x1.51de36942462cp-109} /* i= 138 */
,{0x1.4bd3edda68fe1p-1, -0x1.bde4c79d7d156p-57, -0x1.946a49e22ff5ap-112} /* i= 139 */
,{0x1.4afd6a052bf5bp-1, -0x1.fad40a57eb503p-55, 0x1.a814afd6a052cp-109} /* i= 140 */
,{0x1.4a27fad76014ap-1, 0x1.3fd6bb00a514p-56, -0x1.4a27fad76014ap-111} /* i= 141 */
,{0x1.49539e3b2d067p-1, -0x1.5de8d81edfd6dp-57, -0x1.630e2697cc8afp-111} /* i= 142 */
,{0x1.488052201488p-1, 0x1.488052201488p-55, 0x1.488052201488p-109} /* i= 143 */
,{0x1.47ae147ae147bp-1, -0x1.eb851eb851eb8p-57, -0x1.47ae147ae147bp-111} /* i= 144 */
,{0x1.46dce34596066p-1, 0x1.28382df70ff5dp-56, -0x1.b9c68b2c0cc4ap-110} /* i= 145 */
,{0x1.460cbc7f5cf9ap-1, 0x1.c051832f1fd74p-57, -0x1.978feb9f34381p-113} /* i= 146 */
,{0x1.453d9e2c776cap-1, 0x1.453d9e2c776cap-61, 0x1.453d9e2c776cap-121} /* i= 147 */
,{0x1.446f86562d9fbp-1, -0x1.1be1958b67ebcp-57, 0x1.be1958b67ebb9p-111} /* i= 148 */
,{0x1.43a2730abee4dp-1, 0x1.db5698f7c8601p-57, 0x1.0e89cc2afb934p-111} /* i= 149 */
,{0x1.42d6625d51f87p-1, -0x1.064e2febd299ep-57, 0x1.7547e1bbe6c74p-111} /* i= 150 */
,{0x1.420b5265e5951p-1, 0x1.1ed21562c078cp-56, 0x1.0fb98d85f9b5cp-110} /* i= 151 */
,{0x1.4141414141414p-1, 0x1.4141414141414p-57, 0x1.4141414141414p-113} /* i= 152 */
,{0x1.40782d10e6566p-1, 0x1.909638551fecp-59, -0x1.e0b4439959819p-113} /* i= 153 */
,{0x1.3fb013fb013fbp-1, 0x1.3fb013fb013fbp-61, 0x1.3fb013fb013fbp-121} /* i= 154 */
,{0x1.3ee8f42a5af07p-1, -0x1.2ff608b85ead3p-56, 0x1.e0db4027dd1e8p-110} /* i= 155 */
,{0x1.3e22cbce4a902p-1, 0x1.f1165e7254814p-55, -0x1.dd3431b56fd84p-111} /* i= 156 */
,{0x1.3d5d991aa75c6p-1, -0x1.10bc6f92e7d36p-55, 0x1.2987bf88fce69p-111} /* i= 157 */
,{0x1.3c995a47babe7p-1, 0x1.1013c995a47bbp-55, -0x1.062efec366a5cp-109} /* i= 158 */
,{0x1.3bd60d9232955p-1, -0x1.f4e57985dc38cp-55, -0x1.d9c1145b4bdffp-113} /* i= 159 */
,{0x1.3b13b13b13b14p-1, -0x1.3b13b13b13b14p-55, 0x1.3b13b13b13b14p-109} /* i= 160 */
,{0x1.3a524387ac822p-1, 0x1.83fd8b5b78f0ap-55, 0x1.beecf804e9491p-109} /* i= 161 */
,{0x1.3991c2c187f63p-1, 0x1.b8f4f9e027324p-56, -0x1.e9f3c04e6470bp-110} /* i= 162 */
,{0x1.38d22d366088ep-1, -0x1.030e0d7107f15p-55, -0x1.89785cde656c2p-109} /* i= 163 */
,{0x1.3813813813814p-1, -0x1.fb1fb1fb1fb2p-55, 0x1.3813813813814p-109} /* i= 164 */
,{0x1.3755bd1c945eep-1, -0x1.f030a5658c773p-56, 0x1.2d9b0f33afbbep-112} /* i= 165 */
,{0x1.3698df3de0748p-1, -0x1.ab1232f514a02p-55, -0x1.b4c6f9ef03a3dp-109} /* i= 166 */
,{0x1.35dce5f9f2af8p-1, 0x1.0f21493ab4599p-56, 0x1.da7a4026bb9ccp-112} /* i= 167 */
,{0x1.3521cfb2b78c1p-1, 0x1.a90e7d95bc60ap-56, -0x1.5bc609a90e7d9p-110} /* i= 168 */
,{0x1.34679ace01346p-1, 0x1.e6b3804d19e6bp-55, 0x1.c0268cf359c02p-110} /* i= 169 */
,{0x1.33ae45b57bcb2p-1, -0x1.f3fb3146e92a1p-57, -0x1.a70f9fd98a375p-114} /* i= 170 */
,{0x1.32f5ced6a1dfap-1, 0x1.32f5ced6a1dfap-61, 0x1.32f5ced6a1dfap-121} /* i= 171 */
,{0x1.323e34a2b10bfp-1, 0x1.9b8396ba9de81p-55, 0x1.91f1a515885fbp-110} /* i= 172 */
,{0x1.3187758e9ebb6p-1, 0x1.3187758e9ebb6p-61, 0x1.3187758e9ebb6p-121} /* i= 173 */
,{0x1.30d190130d19p-1, 0x1.30d190130d19p-57, 0x1.30d190130d19p-113} /* i= 174 */
,{0x1.301c82ac4026p-1, 0x1.c82ac4026039p-56, 0x1.56201301c82acp-110} /* i= 175 */
,{0x1.2f684bda12f68p-1, 0x1.2f684bda12f68p-55, 0x1.2f684bda12f68p-109} /* i= 176 */
,{0x1.2eb4ea1fed14bp-1, 0x1.5e012eb4ea1ffp-57, -0x1.75a750ff68a59p-112} /* i= 177 */
,{0x1.2e025c04b8097p-1, 0x1.2e025c04b8097p-61, 0x1.2e025c04b8097p-121} /* i= 178 */
,{0x1.2d50a012d50ap-1, 0x1.2d50a012d50ap-57, 0x1.2d50a012d50ap-113} /* i= 179 */
,{0x1.2c9fb4d812cap-1, -0x1.2c9fb4d812cap-55, 0x1.2c9fb4d812cap-109} /* i= 180 */
,{0x1.2bef98e5a3711p-1, -0x1.76eb7f1f0c4d5p-60, -0x1.e2b59a119309fp-115} /* i= 181 */
,{0x1.2b404ad012b4p-1, 0x1.2b404ad012b4p-55, 0x1.2b404ad012b4p-109} /* i= 182 */
,{0x1.2a91c92f3c105p-1, 0x1.fc804aa4724bdp-56, -0x1.f7d6037fb55b9p-113} /* i= 183 */
,{0x1.29e4129e4129ep-1, 0x1.04a7904a7904ap-55, 0x1.e4129e4129e41p-109} /* i= 184 */
,{0x1.293725bb804a5p-1, -0x1.1b488ff6b646dp-56, -0x1.11fed6c8da448p-111} /* i= 185 */
,{0x1.288b01288b013p-1, -0x1.dd3fb5dd3fb5ep-55, 0x1.6025116025116p-110} /* i= 186 */
,{0x1.27dfa38a1ce4dp-1, 0x1.be1f34963f911p-55, -0x1.eea9e56ae84e9p-110} /* i= 187 */
,{0x1.27350b8812735p-1, 0x1.71024e6a17102p-58, 0x1.39a85c40939a8p-112} /* i= 188 */
,{0x1.268b37cd60127p-1, -0x1.d320ca7fb65d3p-55, -0x1.0653fdb2e9906p-110} /* i= 189 */
,{0x1.25e22708092f1p-1, 0x1.3840497889c2p-57, 0x1.25e22708092f1p-112} /* i= 190 */
,{0x1.2539d7e9177b2p-1, 0x1.ca2a615c34b06p-57, 0x1.32f88e080494ep-111} /* i= 191 */
,{0x1.2492492492492p-1, 0x1.2492492492492p-55, 0x1.2492492492492p-109} /* i= 192 */
,{0x1.23eb79717605bp-1, 0x1.ccaf9ba70e41p-56, -0x1.23eb79717605bp-113} /* i= 193 */
,{0x1.23456789abcdfp-1, 0x1.23456789abcdfp-61, 0x1.23456789abcdfp-121} /* i= 194 */
,{0x1.22a0122a0122ap-1, 0x1.22a0122a0122ap-61, 0x1.22a0122a0122ap-121} /* i= 195 */
,{0x1.21fb78121fb78p-1, 0x1.21fb78121fb78p-57, 0x1.21fb78121fb78p-113} /* i= 196 */
,{0x1.21579804855e6p-1, 0x1.21579804855e6p-61, 0x1.21579804855e6p-121} /* i= 197 */
,{0x1.20b470c67c0d9p-1, -0x1.e2adac8bd766ap-55, -0x1.20b470c67c0d9p-114} /* i= 198 */
,{0x1.2012012012012p-1, 0x1.2012012012012p-61, 0x1.2012012012012p-121} /* i= 199 */
,{0x1.1f7047dc11f7p-1, 0x1.1f7047dc11f7p-55, 0x1.1f7047dc11f7p-109} /* i= 200 */
,{0x1.1ecf43c7fb84cp-1, 0x1.787008f67a1e4p-56, -0x1.1ecf43c7fb84cp-115} /* i= 201 */
,{0x1.1e2ef3b3fb874p-1, 0x1.0c4c0478bbcedp-55, -0x1.1e2ef3b3fb874p-115} /* i= 202 */
,{0x1.1d8f5672e4abdp-1, -0x1.f17fb89c2a634p-55, -0x1.b5437c5fee271p-109} /* i= 203 */
,{0x1.1cf06ada2811dp-1, -0x1.f2a4bafdc61f3p-58, 0x1.6d1408e78356dp-112} /* i= 204 */
,{0x1.1c522fc1ce059p-1, -0x1.32889b7cf21ep-56, 0x1.aa7b47a2b5085p-111} /* i= 205 */
,{0x1.1bb4a4046ed29p-1, 0x1.1bb4a4046ed29p-61, 0x1.1bb4a4046ed29p-121} /* i= 206 */
,{0x1.1b17c67f2bae3p-1, -0x1.37d830a8161dep-55, 0x1.1f842599285cep-109} /* i= 207 */
,{0x1.1a7b9611a7b96p-1, 0x1.1a7b9611a7b96p-57, 0x1.1a7b9611a7b96p-113} /* i= 208 */
,{0x1.19e0119e0119ep-1, 0x1.19e0119e0119ep-61, 0x1.19e0119e0119ep-121} /* i= 209 */
,{0x1.19453808ca29cp-1, 0x1.19453808ca29cp-59, 0x1.19453808ca29cp-117} /* i= 210 */
,{0x1.18ab083902bdbp-1, -0x1.1adc5e4974c32p-55, -0x1.baede8f9f8535p-109} /* i= 211 */
,{0x1.1811811811812p-1, -0x1.fb9fb9fb9fbap-55, 0x1.1811811811812p-109} /* i= 212 */
,{0x1.1778a191bd684p-1, 0x1.8045de28646f6p-57, -0x1.7be7fba21d79cp-111} /* i= 213 */
,{0x1.16e0689427379p-1, -0x1.4b2a7c2fee92p-57, 0x1.a2509cde3ad35p-111} /* i= 214 */
,{0x1.1648d50fc3201p-1, 0x1.648d50fc32011p-57, 0x1.923543f0c8046p-111} /* i= 215 */
,{0x1.15b1e5f75270dp-1, 0x1.15b1e5f75270dp-59, 0x1.15b1e5f75270dp-117} /* i= 216 */
,{0x1.151b9a3fdd5c9p-1, -0x1.a3fdd5c8cb804p-56, -0x1.51b9a3fdd5c8dp-110} /* i= 217 */
,{0x1.1485f0e0acd3bp-1, 0x1.a31b011485f0ep-55, 0x1.59a76d18d808ap-112} /* i= 218 */
,{0x1.13f0e8d344724p-1, 0x1.c0677a574f39bp-57, -0x1.49d5f64c87d09p-111} /* i= 219 */
,{0x1.135c81135c811p-1, 0x1.ae4089ae4089bp-56, -0x1.bf7651bf7651cp-112} /* i= 220 */
,{0x1.12c8b89edc0acp-1, -0x1.0a3272d9e52a6p-55, -0x1.7e1f20bce9fefp-109} /* i= 221 */
,{0x1.12358e75d3033p-1, 0x1.a82ad85e4269p-55, -0x1.dfddb94e3145ap-109} /* i= 222 */
,{0x1.11a3019a74826p-1, 0x1.ebb0e6e1895a5p-55, 0x1.2703bdba859c9p-110} /* i= 223 */
,{0x1.1111111111111p-1, 0x1.1111111111111p-57, 0x1.1111111111111p-113} /* i= 224 */
,{0x1.107fbbe01108p-1, -0x1.107fbbe01108p-55, 0x1.107fbbe01108p-109} /* i= 225 */
,{0x1.0fef010fef011p-1, -0x1.0fef010fef011p-61, 0x1.0fef010fef011p-121} /* i= 226 */
,{0x1.0f5edfab325a2p-1, -0x1.5fef0a12054cep-55, 0x1.686a010f5edfbp-109} /* i= 227 */
,{0x1.0ecf56be69c9p-1, -0x1.0ecf56be69c9p-56, 0x1.0ecf56be69c9p-111} /* i= 228 */
,{0x1.0e40655826011p-1, -0x1.bf9aa7d9fef1cp-57, 0x1.9560980439019p-115} /* i= 229 */
,{0x1.0db20a88f4696p-1, -0x1.9cf8a021b6415p-55, -0x1.1e8d2b3183affp-111} /* i= 230 */
,{0x1.0d24456359e3ap-1, -0x1.69a8bd3d80c9ep-56, 0x1.32fd5f255287ap-110} /* i= 231 */
,{0x1.0c9714fbcda3bp-1, -0x1.f79b47582192ep-56, -0x1.4fbcda3ac10c9p-111} /* i= 232 */
,{0x1.0c0a7868b4171p-1, -0x1.c669c021814f1p-55, 0x1.74be8f719a701p-110} /* i= 233 */
,{0x1.0b7e6ec259dc8p-1, -0x1.b2ad73fbd2064p-55, -0x1.3da62386cab5dp-109} /* i= 234 */
,{0x1.0af2f722eecb5p-1, 0x1.c48fe6f938d4cp-55, -0x1.98c40a6d7da76p-109} /* i= 235 */
,{0x1.0a6810a6810a7p-1, -0x1.fbd65fbd65fbdp-55, -0x1.97ef597ef597fp-109} /* i= 236 */
,{0x1.09ddba6af836p-1, 0x1.09ddba6af836p-57, 0x1.09ddba6af836p-113} /* i= 237 */
,{0x1.0953f39010954p-1, -0x1.8dfded5818dfep-58, 0x1.2a7e720212a7ep-114} /* i= 238 */
,{0x1.08cabb37565e2p-1, 0x1.08cabb37565e2p-61, 0x1.08cabb37565e2p-121} /* i= 239 */
,{0x1.0842108421084p-1, 0x1.0842108421084p-56, 0x1.0842108421084p-111} /* i= 240 */
,{0x1.07b9f29b8eae2p-1, -0x1.8fb5d3b3c43fep-55, 0x1.ee7ca6e3ab867p-112} /* i= 241 */
,{0x1.073260a47f7c6p-1, 0x1.b3eb701073261p-55, -0x1.6e020e64c149p-109} /* i= 242 */
,{0x1.06ab59c7912fbp-1, 0x1.87f3aff7caa53p-55, 0x1.c376824f018ap-111} /* i= 243 */
,{0x1.0624dd2f1a9fcp-1, -0x1.89374bc6a7efap-57, 0x1.26e978d4fdf3bp-112} /* i= 244 */
,{0x1.059eea0727586p-1, 0x1.8c84dab2d7a2p-55, -0x1.4706a488f12e8p-109} /* i= 245 */
,{0x1.05197f7d73404p-1, 0x1.465fdf5cd0105p-57, 0x1.97f7d73404146p-113} /* i= 246 */
,{0x1.04949cc1664c5p-1, 0x1.e27b2a3e17696p-55, -0x1.7aa7f3c908a6fp-109} /* i= 247 */
,{0x1.041041041041p-1, 0x1.041041041041p-55, 0x1.041041041041p-109} /* i= 248 */
,{0x1.038c6b78247fcp-1, -0x1.c635bc123fdf9p-58, 0x1.8d6f048ff7e3ap-114} /* i= 249 */
,{0x1.03091b51f5e1ap-1, 0x1.3bb3194be3abp-55, 0x1.03091b51f5e1ap-111} /* i= 250 */
,{0x1.02864fc7729e9p-1, -0x1.d089575a61f4ep-56, -0x1.0ea49b84cbfep-110} /* i= 251 */
,{0x1.0204081020408p-1, 0x1.0204081020408p-57, 0x1.0204081020408p-113} /* i= 252 */
,{0x1.0182436517a37p-1, 0x1.4bf1eae05078bp-55, 0x1.43e5d8c527bbap-109} /* i= 253 */
,{0x1.010101010101p-1, 0x1.010101010101p-57, 0x1.010101010101p-113} /* i= 254 */
,{0x1.008040201008p-1, 0x1.008040201008p-55, 0x1.008040201008p-109} /* i= 255 */
};


/* table_log_alpha_i
alpha_i = {2^k}/{2^k+i}  with 0<=i<2*k with k=8.
log_alpha_i = {h, m, l} to calculate in Triple-Double precision */
static const double table_log_alpha_i_triple[256][3] = {
//{h,m,l} for each alpha_i
{0x0p+3, 0x0p+3, 0x0p+3} /* i= 0 */
,{0x1.ff00aa2b10bcp-9, 0x1.2821ad5a6d353p-63, -0x1.12dcccb588a4ap-118} /* i= 1 */
,{0x1.fe02a6b106789p-8, -0x1.e44b7e3711ebfp-67, 0x1.a567b6587df34p-121} /* i= 2 */
,{0x1.7dc475f810a77p-7, -0x1.16d7687d3df21p-62, 0x1.a850a4a1800eap-117} /* i= 3 */
,{0x1.fc0a8b0fc03e4p-7, -0x1.83092c59642a1p-62, -0x1.52414fc416fc2p-116} /* i= 4 */
,{0x1.3cea44346a575p-6, -0x1.0cb5a902b3a1cp-62, 0x1.98d0797189a4dp-117} /* i= 5 */
,{0x1.7b91b07d5b11bp-6, -0x1.5b602ace3a51p-60, 0x1.dcd4f102a521dp-118} /* i= 6 */
,{0x1.b9fc027af9198p-6, -0x1.0ae69229dc868p-64, 0x1.9ffdb5331f453p-118} /* i= 7 */
,{0x1.f829b0e7833p-6, 0x1.33e3f04f1ef23p-60, -0x1.814544147acc9p-114} /* i= 8 */
,{0x1.1b0d98923d98p-5, -0x1.e9ae889bac481p-60, -0x1.f6acb8073198bp-114} /* i= 9 */
,{0x1.39e87b9febd6p-5, -0x1.5bfa937f551bbp-59, 0x1.c8d57ae1e11bdp-114} /* i= 10 */
,{0x1.58a5bafc8e4d5p-5, -0x1.ce55c2b4e2b72p-59, -0x1.33fb67ae4f6cep-114} /* i= 11 */
,{0x1.77458f632dcfcp-5, 0x1.18d3ca87b9296p-59, 0x1.63c9bf701b2a9p-116} /* i= 12 */
,{0x1.95c830ec8e3ebp-5, 0x1.f5a0e80520bf2p-59, -0x1.9e0ef8448a202p-113} /* i= 13 */
,{0x1.b42dd711971bfp-5, -0x1.eb9759c130499p-60, -0x1.6b5431d9cbf04p-116} /* i= 14 */
,{0x1.d276b8adb0b52p-5, 0x1.1e3c53257fd47p-61, 0x1.cecc7db99d86ap-117} /* i= 15 */
,{0x1.f0a30c01162a6p-5, 0x1.85f325c5bbacdp-59, -0x1.0ece597165991p-113} /* i= 16 */
,{0x1.075983598e471p-4, 0x1.80da5333c45b8p-59, -0x1.77ad5e5273f98p-116} /* i= 17 */
,{0x1.16536eea37ae1p-4, -0x1.79da3e8c22cdap-60, -0x1.b925bd6fa5998p-116} /* i= 18 */
,{0x1.253f62f0a1417p-4, -0x1.c125963fc4cfdp-62, -0x1.d2c3f5a497e44p-116} /* i= 19 */
,{0x1.341d7961bd1d1p-4, -0x1.b599f227becbbp-58, -0x1.15fbcbe26b491p-113} /* i= 20 */
,{0x1.42edcbea646fp-4, 0x1.ddd4f935996c9p-59, 0x1.7465d8f6866cfp-114} /* i= 21 */
,{0x1.51b073f06183fp-4, 0x1.a49e39a1a8be4p-58, 0x1.584bc9c7e09bcp-112} /* i= 22 */
,{0x1.60658a93750c4p-4, -0x1.388458ec21b6ap-58, 0x1.c66d48ed8883fp-112} /* i= 23 */
,{0x1.6f0d28ae56b4cp-4, -0x1.906d99184b992p-58, -0x1.bf31af3e109afp-112} /* i= 24 */
,{0x1.7da766d7b12cdp-4, -0x1.eeedfcdd94131p-58, 0x1.a115d17a663c2p-112} /* i= 25 */
,{0x1.8c345d6319b21p-4, -0x1.4a697ab3424a9p-61, -0x1.e547ecfe0df94p-115} /* i= 26 */
,{0x1.9ab42462033adp-4, -0x1.2099e1c184e8ep-59, -0x1.bb52cb975cbebp-115} /* i= 27 */
,{0x1.a926d3a4ad563p-4, 0x1.942f48aa70ea9p-58, 0x1.8f353ecfc45dap-113} /* i= 28 */
,{0x1.b78c82bb0eda1p-4, 0x1.0878cf0327e21p-61, -0x1.b0b1387f2d48fp-115} /* i= 29 */
,{0x1.c5e548f5bc743p-4, 0x1.5d617ef8161b1p-60, 0x1.da7659abe370ep-114} /* i= 30 */
,{0x1.d4313d66cb35dp-4, 0x1.790dd951d90fap-58, 0x1.20959368928d5p-113} /* i= 31 */
,{0x1.e27076e2af2e6p-4, -0x1.61578001e0162p-60, 0x1.55db94ebc4018p-116} /* i= 32 */
,{0x1.f0a30c01162a6p-4, 0x1.85f325c5bbacdp-58, -0x1.0ece597165991p-112} /* i= 33 */
,{0x1.fec9131dbeabbp-4, -0x1.5746b9981b36cp-58, -0x1.c4016e1d457eep-112} /* i= 34 */
,{0x1.0671512ca596ep-3, 0x1.50c647eb86499p-58, -0x1.e98f4812aa997p-113} /* i= 35 */
,{0x1.0d77e7cd08e59p-3, 0x1.9a5dc5e9030acp-57, -0x1.71dbd9a581398p-111} /* i= 36 */
,{0x1.14785846742acp-3, 0x1.a28813e3a7f07p-57, 0x1.bd933781e73cdp-112} /* i= 37 */
,{0x1.1b72ad52f67ap-3, 0x1.483023472cd74p-58, -0x1.81887026f66adp-112} /* i= 38 */
,{0x1.2266f190a5acbp-3, 0x1.f547bf1809e88p-57, 0x1.eea44ec5389a5p-111} /* i= 39 */
,{0x1.29552f81ff523p-3, 0x1.301771c407dbfp-57, -0x1.977b021b7c784p-111} /* i= 40 */
,{0x1.303d718e47fd3p-3, -0x1.6b9c7d96091fap-63, -0x1.5e72f6cc4e614p-117} /* i= 41 */
,{0x1.371fc201e8f74p-3, 0x1.de6cb62af18ap-58, -0x1.a2fc19b24ab16p-113} /* i= 42 */
,{0x1.3dfc2b0ecc62ap-3, -0x1.ab3a8e7d81017p-58, -0x1.b40efe811e153p-112} /* i= 43 */
,{0x1.44d2b6ccb7d1ep-3, 0x1.9f4f6543e1f88p-57, -0x1.f3be9a8337458p-111} /* i= 44 */
,{0x1.4ba36f39a55e5p-3, 0x1.68981bcc36756p-57, -0x1.04bfef68b5ce2p-116} /* i= 45 */
,{0x1.526e5e3a1b438p-3, -0x1.746ff8a470d3ap-57, 0x1.a6dbcc63b5444p-111} /* i= 46 */
,{0x1.59338d9982086p-3, -0x1.65d22aa8ad7cfp-58, 0x1.60e1f10db27cbp-112} /* i= 47 */
,{0x1.5ff3070a793d4p-3, -0x1.bc60efafc6f6ep-58, -0x1.140655471954p-113} /* i= 48 */
,{0x1.66acd4272ad51p-3, -0x1.0900e4e1ea8b2p-58, -0x1.80ab0a1bc6d9bp-112} /* i= 49 */
,{0x1.6d60fe719d21dp-3, -0x1.caae268ecd179p-57, -0x1.c825cda7da31dp-114} /* i= 50 */
,{0x1.740f8f54037a5p-3, -0x1.b264062a84cdbp-58, -0x1.0be957f10f5fbp-112} /* i= 51 */
,{0x1.7ab890210d909p-3, 0x1.be36b2d6a0608p-59, 0x1.91ff852536204p-117} /* i= 52 */
,{0x1.815c0a14357ebp-3, -0x1.4be48073a0564p-58, 0x1.435bddbbe732cp-112} /* i= 53 */
,{0x1.87fa06520c911p-3, -0x1.bf7fdbfa08d9ap-57, -0x1.0a5aa8fb49481p-112} /* i= 54 */
,{0x1.8e928de886d41p-3, -0x1.569d851a5677p-57, 0x1.c0d0e377c6294p-114} /* i= 55 */
,{0x1.9525a9cf456b4p-3, 0x1.d904c1d4e2e26p-57, -0x1.89d9afa096184p-111} /* i= 56 */
,{0x1.9bb362e7dfb83p-3, 0x1.575e31f003e0cp-57, 0x1.28792ae1aabc8p-112} /* i= 57 */
,{0x1.a23bc1fe2b563p-3, 0x1.93711b07a998cp-59, 0x1.3f1f8db36c599p-114} /* i= 58 */
,{0x1.a8becfc882f19p-3, -0x1.e8c37918c39ebp-58, 0x1.58b02842ae948p-114} /* i= 59 */
,{0x1.af3c94e80bff3p-3, -0x1.398cff3641985p-58, -0x1.a262591d1968bp-114} /* i= 60 */
,{0x1.b5b519e8fb5a4p-3, 0x1.ba27fdc19e1ap-57, 0x1.3dcf06e27bef1p-111} /* i= 61 */
,{0x1.bc286742d8cd6p-3, 0x1.4fce744870f55p-58, -0x1.e1d3c235b937cp-115} /* i= 62 */
,{0x1.c2968558c18c1p-3, -0x1.73dee38a3fb6bp-57, 0x1.f00f527d33467p-118} /* i= 63 */
,{0x1.c8ff7c79a9a22p-3, -0x1.4f689f8434012p-57, 0x1.a24ae3b2f53a1p-111} /* i= 64 */
,{0x1.cf6354e09c5dcp-3, 0x1.239a07d55b695p-57, 0x1.a1077102874fp-111} /* i= 65 */
,{0x1.d5c216b4fbb91p-3, 0x1.6e443597e4d4p-57, 0x1.c3c6ce7a257f4p-113} /* i= 66 */
,{0x1.dc1bca0abec7dp-3, 0x1.834c51998b6fcp-57, 0x1.dd2b51478112ep-113} /* i= 67 */
,{0x1.e27076e2af2e6p-3, -0x1.61578001e0162p-59, 0x1.55db94ebc4018p-115} /* i= 68 */
,{0x1.e8c0252aa5a6p-3, -0x1.6e03a39bfc89bp-59, 0x1.dee364d35208ap-113} /* i= 69 */
,{0x1.ef0adcbdc5936p-3, 0x1.48637950dc20dp-57, -0x1.eb052d7b3cbe3p-111} /* i= 70 */
,{0x1.f550a564b7b37p-3, 0x1.c5f6dfd018c37p-61, 0x1.98a014b61d51p-120} /* i= 71 */
,{0x1.fb9186d5e3e2bp-3, -0x1.caaae64f21acbp-57, -0x1.35f6dfd3ddd52p-111} /* i= 72 */
,{0x1.00e6c45ad501dp-2, -0x1.cb9568ff6feadp-57, 0x1.60709f1d0d49fp-113} /* i= 73 */
,{0x1.0402594b4d041p-2, -0x1.28ec217a5022dp-57, -0x1.0dddc4cf9a1f9p-111} /* i= 74 */
,{0x1.071b85fcd590dp-2, 0x1.d1707f97bde8p-58, 0x1.00ca1b7fa08dap-113} /* i= 75 */
,{0x1.0a324e27390e3p-2, 0x1.7dcfde8061c03p-56, 0x1.c51bc06b5f7c1p-113} /* i= 76 */
,{0x1.0d46b579ab74bp-2, 0x1.03ec81c3cbd92p-57, 0x1.7333da8be1a7dp-111} /* i= 77 */
,{0x1.1058bf9ae4ad5p-2, 0x1.89fa0ab4cb31dp-58, -0x1.eb31a74640ec7p-116} /* i= 78 */
,{0x1.136870293a8bp-2, 0x1.7b66298edd24ap-56, -0x1.4a5b394627b29p-113} /* i= 79 */
,{0x1.1675cababa60ep-2, 0x1.ce63eab883717p-61, 0x1.1f833e82521e1p-119} /* i= 80 */
,{0x1.1980d2dd4236fp-2, 0x1.9d3d1b0e4d147p-56, -0x1.8eb33aa901486p-110} /* i= 81 */
,{0x1.1c898c16999fbp-2, -0x1.0e5c62aff1c44p-60, -0x1.e623be88a509bp-115} /* i= 82 */
,{0x1.1f8ff9e48a2f3p-2, -0x1.c9fdf9a0c4b07p-56, 0x1.8cf23e43622b1p-110} /* i= 83 */
,{0x1.22941fbcf7966p-2, -0x1.76f5eb09628afp-56, -0x1.a168b2a9642c4p-111} /* i= 84 */
,{0x1.2596010df763ap-2, -0x1.0f76c57075e9ep-58, 0x1.82ce04d7e207dp-113} /* i= 85 */
,{0x1.2895a13de86a3p-2, 0x1.7ad24c13f040ep-56, 0x1.62d6a3aacbe58p-110} /* i= 86 */
,{0x1.2b9303ab89d25p-2, -0x1.896b5fd852ad4p-56, -0x1.0529c8be2b81bp-110} /* i= 87 */
,{0x1.2e8e2bae11d31p-2, -0x1.8f4cdb95ebdf9p-56, -0x1.864244294826fp-111} /* i= 88 */
,{0x1.31871c9544185p-2, -0x1.51acc4c09b379p-60, -0x1.19a07a2d2cc1ep-114} /* i= 89 */
,{0x1.347dd9a987d55p-2, -0x1.4dd4c580919f8p-57, 0x1.ee510a580b3b3p-111} /* i= 90 */
,{0x1.3772662bfd85bp-2, -0x1.b5629d8117de7p-59, 0x1.790d82b75e92p-113} /* i= 91 */
,{0x1.3a64c556945eap-2, -0x1.c68651945f97cp-57, 0x1.beb7a3cee7e03p-111} /* i= 92 */
,{0x1.3d54fa5c1f71p-2, -0x1.e3265c6a1c98dp-56, 0x1.229e62e452918p-111} /* i= 93 */
,{0x1.404308686a7e4p-2, -0x1.0bcfb6082ce6dp-56, -0x1.9ea6f9f60989cp-110} /* i= 94 */
,{0x1.432ef2a04e814p-2, -0x1.29931715ac903p-56, -0x1.3f95697c9bfc2p-110} /* i= 95 */
,{0x1.4618bc21c5ec2p-2, 0x1.f42decdeccf1dp-56, -0x1.77d446996dap-111} /* i= 96 */
,{0x1.49006804009d1p-2, -0x1.9ffc341f177dcp-57, 0x1.16c8675ad963dp-113} /* i= 97 */
,{0x1.4be5f957778a1p-2, -0x1.259b35b04813dp-57, 0x1.1eb953458673dp-112} /* i= 98 */
,{0x1.4ec973260026ap-2, -0x1.42a87d977dc5ep-56, -0x1.fcf3e64c8cd74p-110} /* i= 99 */
,{0x1.51aad872df82dp-2, 0x1.3927ac19f55e3p-59, 0x1.1d4f4f357cbfbp-115} /* i= 100 */
,{0x1.548a2c3add263p-2, -0x1.819cf7e308ddbp-57, -0x1.8294131dd7142p-111} /* i= 101 */
,{0x1.5767717455a6cp-2, 0x1.526adb283660cp-56, -0x1.7f83a3e5e6736p-111} /* i= 102 */
,{0x1.5a42ab0f4cfe2p-2, -0x1.8ebcb7dee9a3dp-56, 0x1.6f95d595cbf2ep-110} /* i= 103 */
,{0x1.5d1bdbf5809cap-2, 0x1.4236383dc7fe1p-56, 0x1.59f380b4a6b43p-112} /* i= 104 */
,{0x1.5ff3070a793d4p-2, -0x1.bc60efafc6f6ep-57, -0x1.140655471954p-112} /* i= 105 */
,{0x1.62c82f2b9c795p-2, 0x1.7b7af915300e5p-57, 0x1.7391362aee92cp-113} /* i= 106 */
,{0x1.659b57303e1f3p-2, -0x1.f893d41c411f1p-56, -0x1.3fe778dfe7cc6p-114} /* i= 107 */
,{0x1.686c81e9b14afp-2, -0x1.ddea0f7f58e3dp-57, 0x1.2c96f6f68e19dp-111} /* i= 108 */
,{0x1.6b3bb2235943ep-2, -0x1.da856ccd987b3p-56, 0x1.8378506ba0045p-114} /* i= 109 */
,{0x1.6e08eaa2ba1e4p-2, -0x1.cfb1b39ca3a0fp-56, -0x1.0fce95182c66ap-110} /* i= 110 */
,{0x1.70d42e2789236p-2, -0x1.52cc811d78d59p-57, -0x1.c657d4b4b3ef6p-114} /* i= 111 */
,{0x1.739d7f6bbd007p-2, -0x1.8c76ceb014b04p-56, -0x1.0d2a910f7918bp-111} /* i= 112 */
,{0x1.7664e1239dbcfp-2, -0x1.f6d5d64f5daf8p-57, 0x1.d4b7fcd3804aep-111} /* i= 113 */
,{0x1.792a55fdd47a2p-2, 0x1.f057691fe9ed7p-56, -0x1.fa980f34439f2p-110} /* i= 114 */
,{0x1.7bede0a37afcp-2, -0x1.8783cb9801a5cp-56, 0x1.30a6d4e7913d3p-112} /* i= 115 */
,{0x1.7eaf83b82afc3p-2, 0x1.92ce979ed295p-56, 0x1.0dc5832ff2fdcp-110} /* i= 116 */
,{0x1.816f41da0d496p-2, -0x1.2923ca04b701cp-56, -0x1.fd0a6e1849747p-112} /* i= 117 */
,{0x1.842d1da1e8b17p-2, 0x1.24ec519784676p-56, 0x1.a23c11851c7cep-110} /* i= 118 */
,{0x1.86e919a330bap-2, 0x1.3f9b16feb7dd8p-59, -0x1.45cedb41082dfp-113} /* i= 119 */
,{0x1.89a3386c1425bp-2, -0x1.29639dfbbf0fbp-56, 0x1.6cfff18ca06dp-110} /* i= 120 */
,{0x1.8c5b7c858b48bp-2, -0x1.e0ab4fdfa0595p-56, 0x1.2a2f1786f3a7dp-111} /* i= 121 */
,{0x1.8f11e873662c7p-2, 0x1.f85da755a61a3p-56, -0x1.9a18d00d0fc6fp-110} /* i= 122 */
,{0x1.91c67eb45a83ep-2, -0x1.e0e0ae234ae11p-56, 0x1.22fc55f6101c7p-110} /* i= 123 */
,{0x1.947941c2116fbp-2, -0x1.16cc8bae0bbe4p-56, -0x1.515b58cf688d8p-110} /* i= 124 */
,{0x1.972a341135158p-2, 0x1.a5c09d24b70d9p-56, 0x1.2107598781dc7p-110} /* i= 125 */
,{0x1.99d958117e08bp-2, -0x1.a2b6889dc3e72p-57, -0x1.16d1238da82edp-115} /* i= 126 */
,{0x1.9c86b02dc0863p-2, -0x1.917eeb69dd421p-56, 0x1.88a54f77fc355p-111} /* i= 127 */
,{0x1.9f323ecbf984cp-2, -0x1.a92e513217f5cp-59, 0x1.0c0cfa41ff669p-113} /* i= 128 */
,{0x1.a1dc064d5b995p-2, 0x1.90128698ba0b8p-56, 0x1.a892e1c78a129p-111} /* i= 129 */
,{0x1.a484090e5bb0ap-2, 0x1.5fe535b875a75p-57, -0x1.a6c6290af394ap-111} /* i= 130 */
,{0x1.a72a4966bd9eap-2, 0x1.6a76b1a7d87c3p-58, -0x1.ee4d9e07a81b8p-113} /* i= 131 */
,{0x1.a9cec9a9a084ap-2, -0x1.cadec02b436afp-56, -0x1.420f701b88eccp-111} /* i= 132 */
,{0x1.ac718c258b0e4p-2, 0x1.8163d6f46f714p-59, 0x1.f7e9fe1d457fbp-114} /* i= 133 */
,{0x1.af1293247786bp-2, 0x1.133844a15dc28p-58, 0x1.87134125f21c2p-115} /* i= 134 */
,{0x1.b1b1e0ebdfc5bp-2, 0x1.a4479608a2c55p-56, 0x1.d790ec4a16c08p-110} /* i= 135 */
,{0x1.b44f77bcc8f63p-2, -0x1.cd04495459c78p-56, -0x1.c437eb152cbdep-110} /* i= 136 */
,{0x1.b6eb59d3cf35ep-2, -0x1.8adbccd326a3cp-56, -0x1.dca18bc6bd6e1p-110} /* i= 137 */
,{0x1.b9858969310fbp-2, 0x1.663ec53e23bc4p-56, -0x1.8437e3152e77fp-110} /* i= 138 */
,{0x1.bc1e08b0dad0ap-2, 0x1.09e8707055996p-56, -0x1.2402cee15be62p-112} /* i= 139 */
,{0x1.beb4d9da71b7cp-2, -0x1.0f3c590a887cap-59, -0x1.b495a7c83dffcp-113} /* i= 140 */
,{0x1.c149ff115f027p-2, -0x1.4cbcb90c06305p-56, 0x1.34b43a830d5b7p-113} /* i= 141 */
,{0x1.c3dd7a7cdad4dp-2, 0x1.cecf052dea69bp-56, 0x1.82ed46395f605p-110} /* i= 142 */
,{0x1.c66f4e3ff6ff8p-2, -0x1.82947258b688bp-58, -0x1.07c424268805cp-112} /* i= 143 */
,{0x1.c8ff7c79a9a22p-2, -0x1.4f689f8434012p-56, 0x1.a24ae3b2f53a1p-110} /* i= 144 */
,{0x1.cb8e0744d7acap-2, -0x1.48879a214a2afp-61, 0x1.9e13827c5457cp-117} /* i= 145 */
,{0x1.ce1af0b85f3ebp-2, 0x1.edf4af2ab4267p-56, 0x1.2710c64600598p-110} /* i= 146 */
,{0x1.d0a63ae721e64p-2, 0x1.2acce112c40f2p-57, 0x1.87027a17f1c34p-111} /* i= 147 */
,{0x1.d32fe7e00ebd5p-2, 0x1.877b232fafa37p-56, -0x1.73aa590050815p-115} /* i= 148 */
,{0x1.d5b7f9ae2c684p-2, -0x1.a7be7f84ac06ap-57, -0x1.215d8bf93a178p-113} /* i= 149 */
,{0x1.d83e7258a2f3ep-2, 0x1.41456e8bb2511p-56, 0x1.d4a129983048fp-113} /* i= 150 */
,{0x1.dac353e2c5954p-2, 0x1.18734b81a1bf8p-57, 0x1.e1616e962bcf9p-112} /* i= 151 */
,{0x1.dd46a04c1c4a1p-2, -0x1.0467656d8b892p-56, 0x1.fe9f50684ce6cp-112} /* i= 152 */
,{0x1.dfc859906d5b5p-2, 0x1.01e1399f96398p-56, -0x1.14215547c2d4cp-110} /* i= 153 */
,{0x1.e24881a7c6c26p-2, 0x1.cbd8f45954a46p-58, 0x1.b1500f7c5d938p-113} /* i= 154 */
,{0x1.e4c71a8687704p-2, 0x1.667923e1f5a8ep-57, -0x1.958bdeb5faa65p-112} /* i= 155 */
,{0x1.e744261d68788p-2, -0x1.c825c90c344b9p-58, -0x1.2fed79c755684p-114} /* i= 156 */
,{0x1.e9bfa659861f5p-2, 0x1.91bafc7dbe13p-56, 0x1.a2e6a81cf3b6p-110} /* i= 157 */
,{0x1.ec399d2468ccp-2, 0x1.75cee53f35397p-58, -0x1.3dda340d7c50ap-118} /* i= 158 */
,{0x1.eeb20c640ddf4p-2, 0x1.ac371d7c8f7f5p-57, -0x1.ec6e2c3232e6fp-111} /* i= 159 */
,{0x1.f128f5faf06edp-2, -0x1.328df13bb38c3p-56, 0x1.d73d592445d0ap-110} /* i= 160 */
,{0x1.f39e5bc811e5cp-2, -0x1.97fc777bb19e5p-57, 0x1.de5246e8e04f1p-112} /* i= 161 */
,{0x1.f6123fa7028acp-2, 0x1.8515b0f2db341p-56, 0x1.2195120a66058p-110} /* i= 162 */
,{0x1.f884a36fe9ec2p-2, 0x1.6315c9e0108p-57, -0x1.956dcfe3d63fcp-112} /* i= 163 */
,{0x1.faf588f78f31fp-2, -0x1.328260d8abcap-57, -0x1.392b321d10e7bp-112} /* i= 164 */
,{0x1.fd64f20f61572p-2, -0x1.adb0ac2cead1bp-57, 0x1.56b699a6a9876p-113} /* i= 165 */
,{0x1.ffd2e0857f498p-2, 0x1.565f40d9321afp-56, 0x1.23719bce9f534p-111} /* i= 166 */
,{0x1.011fab125ff8ap-1, 0x1.810dd40845ddep-57, 0x1.e4aebfc09efa1p-111} /* i= 167 */
,{0x1.02552a5a5d0ffp-1, -0x1.cb1cb51408cp-56, -0x1.cb91b47473b3dp-112} /* i= 168 */
,{0x1.0389eefce633bp-1, 0x1.e155c53483748p-56, 0x1.c46213221b991p-120} /* i= 169 */
,{0x1.04bdf9da926d2p-1, 0x1.97f304022c9dfp-55, 0x1.a9b423911c3c4p-109} /* i= 170 */
,{0x1.05f14bd26459cp-1, 0x1.535b8ee4f9efep-58, 0x1.2da6f8cd96c9ap-112} /* i= 171 */
,{0x1.0723e5c1cdf4p-1, 0x1.395e58e2445bbp-55, -0x1.49a90b4515bdep-109} /* i= 172 */
,{0x1.0855c884b450ep-1, 0x1.705826e49f318p-55, -0x1.80d6fb1d01dc2p-110} /* i= 173 */
,{0x1.0986f4f573521p-1, -0x1.1b8095ac02f01p-55, 0x1.c089f89ad131cp-115} /* i= 174 */
,{0x1.0ab76bece14d2p-1, -0x1.fd6c935453f66p-56, -0x1.1795f418a9efep-112} /* i= 175 */
,{0x1.0be72e4252a83p-1, -0x1.259da11330801p-55, 0x1.a6d90d9beefcdp-110} /* i= 176 */
,{0x1.0d163ccb9d6b8p-1, -0x1.f7b9a9a8bc30fp-57, -0x1.9ef7ee909d097p-111} /* i= 177 */
,{0x1.0e44985d1cc8cp-1, -0x1.22a3442d2d384p-58, 0x1.a3f759ee145b4p-112} /* i= 178 */
,{0x1.0f7241c9b497dp-1, 0x1.3a8443b9db19dp-55, -0x1.27bf5e1f88b97p-109} /* i= 179 */
,{0x1.109f39e2d4c97p-1, -0x1.0e09b27a4373ap-60, -0x1.b7d38320cdf03p-117} /* i= 180 */
,{0x1.11cb81787ccf8p-1, 0x1.02387ab1fcc9p-55, 0x1.56e9eccf60493p-110} /* i= 181 */
,{0x1.12f719593efbcp-1, 0x1.4c048c671f435p-55, 0x1.6893b2757f501p-110} /* i= 182 */
,{0x1.1422025243d45p-1, -0x1.ad0e24adb489ep-58, -0x1.c16c0bc1a26a2p-113} /* i= 183 */
,{0x1.154c3d2f4d5eap-1, -0x1.59c33171a6876p-55, 0x1.53b4e8cc3cd07p-114} /* i= 184 */
,{0x1.1675cababa60ep-1, 0x1.ce63eab883717p-60, 0x1.1f833e82521e1p-118} /* i= 185 */
,{0x1.179eabbd899a1p-1, -0x1.00e7c6417e0b4p-55, -0x1.e54e3904f3714p-109} /* i= 186 */
,{0x1.18c6e0ff5cf06p-1, 0x1.765142c2c671fp-58, -0x1.355dad2cb4de5p-115} /* i= 187 */
,{0x1.19ee6b467c96fp-1, -0x1.9d1a11443f10cp-56, -0x1.5477c38afc9eap-111} /* i= 188 */
,{0x1.1b154b57da29fp-1, -0x1.011eb47db6a99p-57, 0x1.683fb14c9a0cp-112} /* i= 189 */
,{0x1.1c3b81f713c25p-1, -0x1.0dac1c4c810e9p-55, 0x1.f8efe9846f366p-109} /* i= 190 */
,{0x1.1d610fe677003p-1, 0x1.09d58d91e58f2p-58, -0x1.2dde3ee09695ap-114} /* i= 191 */
,{0x1.1e85f5e7040dp-1, 0x1.ef62cd2f9f1e3p-56, 0x1.7cb9f293d205ep-110} /* i= 192 */
,{0x1.1faa34b87094cp-1, 0x1.817b8f7a193bp-58, -0x1.c07ae7ea7aa51p-112} /* i= 193 */
,{0x1.20cdcd192ab6ep-1, -0x1.b2bf0bc229014p-55, 0x1.27a25206a44a1p-110} /* i= 194 */
,{0x1.21f0bfc65beecp-1, -0x1.e24f0c9187c92p-57, 0x1.5234d992b48aep-111} /* i= 195 */
,{0x1.23130d7bebf43p-1, -0x1.f48725e374d6ep-55, 0x1.48e379bf983ebp-113} /* i= 196 */
,{0x1.2434b6f483934p-1, -0x1.debb8cf0f6d11p-57, 0x1.e9df0af4351dep-111} /* i= 197 */
,{0x1.2555bce98f7cbp-1, 0x1.e021d6d6881e7p-56, 0x1.084750a06eb3p-112} /* i= 198 */
,{0x1.26762013430ep-1, -0x1.96a95781c6727p-56, -0x1.311e25567cac3p-111} /* i= 199 */
,{0x1.2795e1289b11bp-1, -0x1.487c0c246978ep-57, -0x1.fe56c1467b5e6p-119} /* i= 200 */
,{0x1.28b500df60783p-1, -0x1.43f60605aaab3p-55, -0x1.a97977c1a1956p-113} /* i= 201 */
,{0x1.29d37fec2b08bp-1, -0x1.bd1949a2d1982p-56, -0x1.28bcc0f82a9a6p-110} /* i= 202 */
,{0x1.2af15f02640adp-1, 0x1.cb064524acebp-57, 0x1.3fd09b70926d3p-116} /* i= 203 */
,{0x1.2c0e9ed448e8cp-1, -0x1.1a158f3917586p-55, -0x1.dab7eb5720f7bp-109} /* i= 204 */
,{0x1.2d2b4012edc9ep-1, -0x1.51162c99b1cabp-55, -0x1.e6ee8e84d6602p-109} /* i= 205 */
,{0x1.2e47436e40268p-1, 0x1.0150861a4886bp-55, -0x1.db5a61ad75a6fp-110} /* i= 206 */
,{0x1.2f62a99509546p-1, 0x1.6c686739ffd99p-56, -0x1.d81cad72edd27p-112} /* i= 207 */
,{0x1.307d7334f10bep-1, 0x1.fb590a1f566dap-57, -0x1.08f3fa47f6664p-111} /* i= 208 */
,{0x1.3197a0fa7fe6ap-1, 0x1.d6348fb97128fp-57, -0x1.08a857b49329ap-113} /* i= 209 */
,{0x1.32b1339121d71p-1, 0x1.902ab5b3d916bp-56, 0x1.9c56e84cd18b7p-114} /* i= 210 */
,{0x1.33ca2ba328995p-1, -0x1.bf28b3205ede1p-56, 0x1.e925f12836a8ep-111} /* i= 211 */
,{0x1.34e289d9ce1d3p-1, 0x1.6eb92d885ce4fp-57, -0x1.46d67110163eap-111} /* i= 212 */
,{0x1.35fa4edd36eap-1, 0x1.27d4680964362p-60, 0x1.248cc8c939424p-117} /* i= 213 */
,{0x1.37117b54747b6p-1, -0x1.d117edbdd9103p-56, -0x1.c1da9c99e4f6p-110} /* i= 214 */
,{0x1.38280fe58797fp-1, -0x1.015bd362a6e5dp-55, 0x1.ab9cc11bb5935p-113} /* i= 215 */
,{0x1.393e0d3562a1ap-1, -0x1.58eef67f2483ap-55, 0x1.c7b10b8be4f38p-111} /* i= 216 */
,{0x1.3a5373e7ebdfap-1, -0x1.cd8f775b8f76ep-55, -0x1.7f0b45615ae37p-110} /* i= 217 */
,{0x1.3b68449fffc23p-1, -0x1.41c484f9e9b26p-55, -0x1.7c524324c8d4ep-109} /* i= 218 */
,{0x1.3c7c7fff73206p-1, -0x1.be80db7025bedp-56, -0x1.006bcdebdbe0fp-111} /* i= 219 */
,{0x1.3d9026a7156fbp-1, -0x1.6fef670bd4b62p-55, 0x1.ed7013b2d2a96p-109} /* i= 220 */
,{0x1.3ea33936b2f5cp-1, -0x1.f099168a1360bp-55, -0x1.db19fd92dbb9ep-111} /* i= 221 */
,{0x1.3fb5b84d16f42p-1, 0x1.6d3a754172aefp-55, -0x1.937b130cc534bp-112} /* i= 222 */
,{0x1.40c7a4880dce9p-1, 0x1.14f22de7fc9e1p-56, -0x1.5da0ecfb398f5p-110} /* i= 223 */
,{0x1.41d8fe84672aep-1, 0x1.9192f30bd1806p-55, -0x1.0d58eede45763p-110} /* i= 224 */
,{0x1.42e9c6ddf80bfp-1, 0x1.657dc7a65061dp-56, -0x1.e9004fd2f0d6fp-111} /* i= 225 */
,{0x1.43f9fe2f9ce67p-1, 0x1.e9c9ee6d83b86p-55, 0x1.6d8376ee985fdp-109} /* i= 226 */
,{0x1.4509a5133bb0ap-1, 0x1.40fe2852d7b5ap-55, 0x1.2cfd3c8bc24e6p-109} /* i= 227 */
,{0x1.4618bc21c5ec2p-1, 0x1.f42decdeccf1dp-55, -0x1.77d446996dap-110} /* i= 228 */
,{0x1.472743f33aaadp-1, 0x1.8d6cf012a2948p-56, 0x1.80edce5db1351p-113} /* i= 229 */
,{0x1.48353d1ea88dfp-1, 0x1.cf57a2ecc07f4p-55, 0x1.2c307bef9e0cep-110} /* i= 230 */
,{0x1.4942a83a2fc07p-1, 0x1.ed0c544652b5ap-55, 0x1.182657bd7147p-109} /* i= 231 */
,{0x1.4a4f85db03ebbp-1, 0x1.13dfa3d3761b6p-60, 0x1.8b737b8c8ec58p-115} /* i= 232 */
,{0x1.4b5bd6956e274p-1, -0x1.c87a06beea773p-55, 0x1.95ff844550c6ep-109} /* i= 233 */
,{0x1.4c679afccee3ap-1, -0x1.3a5c4c8b39e41p-55, 0x1.8676c36226ef9p-109} /* i= 234 */
,{0x1.4d72d3a39fdp-1, 0x1.1cd4d414e008dp-55, -0x1.0112cb85b0ba3p-111} /* i= 235 */
,{0x1.4e7d811b75bb1p-1, -0x1.8d3d9ea6e9ea9p-55, 0x1.c34317af28812p-109} /* i= 236 */
,{0x1.4f87a3f5026e9p-1, -0x1.e8ca8b1bcea9dp-55, 0x1.d9ac37add819cp-110} /* i= 237 */
,{0x1.50913cc01686bp-1, 0x1.2f2ce96c2d5b1p-55, -0x1.2d0dc61275676p-112} /* i= 238 */
,{0x1.519a4c0ba3446p-1, 0x1.9b32128e4a77fp-55, -0x1.95b50c7e348bdp-112} /* i= 239 */
,{0x1.52a2d265bc5abp-1, -0x1.1883750ea4d0ap-57, -0x1.58412f6df095bp-112} /* i= 240 */
,{0x1.53aad05b99b7dp-1, -0x1.55c8b052e2539p-55, -0x1.515f6592356e7p-110} /* i= 241 */
,{0x1.54b2467999498p-1, -0x1.5baaf5d2f09f4p-55, -0x1.a5dae8aa5423bp-110} /* i= 242 */
,{0x1.55b9354b40bcdp-1, 0x1.e4197a357cb37p-56, 0x1.3188978609b1ep-110} /* i= 243 */
,{0x1.56bf9d5b3f399p-1, 0x1.0471885cd8ff3p-55, -0x1.8c8faa739028fp-110} /* i= 244 */
,{0x1.57c57f336f191p-1, -0x1.e953a3bc88192p-55, -0x1.cc564ea9b22dcp-111} /* i= 245 */
,{0x1.58cadb5cd7989p-1, 0x1.849792ec98458p-56, 0x1.544f1806acad7p-110} /* i= 246 */
,{0x1.59cfb25fae87ep-1, -0x1.172904559c6b6p-58, 0x1.8a981b81fd4e9p-112} /* i= 247 */
,{0x1.5ad404c359f2dp-1, -0x1.35955683f7196p-59, 0x1.08b073c08af03p-117} /* i= 248 */
,{0x1.5bd7d30e71c73p-1, 0x1.bf8da6db2b45cp-57, 0x1.8b103e22f031fp-111} /* i= 249 */
,{0x1.5cdb1dc6c1765p-1, -0x1.cc2470e8a3df4p-55, 0x1.5ec04a15c651dp-109} /* i= 250 */
,{0x1.5ddde57149923p-1, 0x1.dcfa37d75ef28p-55, 0x1.12508153fac48p-110} /* i= 251 */
,{0x1.5ee02a9241675p-1, 0x1.c358257f49082p-55, -0x1.0b39d60fb51b2p-112} /* i= 252 */
,{0x1.5fe1edad18919p-1, -0x1.ca8b610e18dbfp-55, 0x1.067cd2d04edf3p-111} /* i= 253 */
,{0x1.60e32f44788d9p-1, -0x1.ac1bb52fa589bp-56, 0x1.50cd45f38dd6bp-110} /* i= 254 */
,{0x1.61e3efda46467p-1, -0x1.a1b727edefae3p-55, -0x1.9e19fd1f774b7p-109} /* i= 255 */
};

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
/* x is double*/
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

/* log(x) in double*/
}
