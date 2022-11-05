static mpfr_rnd_t rnd2[] = { MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD };

static int rnd = 0; /* default is to nearest */

int ref_fesetround(int rounding_mode)
{
  rnd = rounding_mode;
  return 0;
}

void ref_init(void)
{
  mpfr_set_emin (-1073);
  mpfr_set_emax (1024);
}
