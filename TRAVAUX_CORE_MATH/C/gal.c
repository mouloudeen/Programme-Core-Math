#include <stdio.h>
#include <mpfr.h>

/* 0 < i <256, find the smallest double-precision alpha >= 2^8/(2^8+i) such that
   log(alpha) in double precision is accurate to 71 bits */

int
main ()
{
  int K = 256; /* K = 2^k */
  mpfr_t alpha, t, u;
  mpfr_init2 (alpha, 53);
  mpfr_init2 (t, 53);
  mpfr_init2 (u, 71);
  for (int i = 1;i<256;i++){
    printf("i = %d\n",i);
  
    mpfr_set_ui (alpha, K, MPFR_RNDN);
    /* we round alpha upwards so that m*alpha >= 1 */
    mpfr_div_ui (alpha, alpha, K + i, MPFR_RNDU);
    /*mpfr_printf ("alpha0=%Ra\n", alpha);*/
    while (1)
    {
      mpfr_log (u, alpha, MPFR_RNDN);
      mpfr_set (t, u, MPFR_RNDN);
      if (mpfr_cmp (t, u) == 0)
        break;
      mpfr_nextabove (alpha);
    }
    mpfr_printf ("alpha=%Ra ,\n,", alpha);
  
  }
  mpfr_clear (alpha);
  mpfr_clear (t);
  mpfr_clear (u);
  return 0;
}
/*gcc gal.c -lmpfr
 time ./a.out */