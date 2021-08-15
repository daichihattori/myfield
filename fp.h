#include "integ.h"
//init
void fp_init(fp_t *a) {
  integ_init(a->x0, FPLIMB);
}

//print
void fp_printf(char *str, fp_t *a) {
  integ_printf(str, a->x0, FPLIMB);
}

void fp_println(char *str, fp_t *a) {
  integ_println(str, a->x0, FPLIMB);
}

void fp_set(fp_t *ans, fp_t *a) {
  integ_set(ans->x0, a->x0, FPLIMB);
}

void fp_set_ui(fp_t *ans, unsigned long int UI) {
  integ_set_ui(ans->x0, FPLIMB, UI);
}

void fp_set_mpn(fp_t *ans, mp_limb_t *a) {
  integ_set(ans->x0, a, FPLIMB);
}

void fp_set_neg(fp_t *ans, fp_t *a) {
  integ_set_neg(ans->x0, a->x0, FPLIMB);
}

void fp_l1shift(fp_t *ans, fp_t *a) {
  integ_l1shift(ans->x0, a->x0, FPLIMB);
}

void fp_r1shift(fp_t *ans, fp_t *a) {
  integ_l1shift(ans->x0, a->x0, FPLIMB);
}

void fp_set_random(fp_t *ans, gmp_randstate_t state) {
  integ_set_random(ans->x0, prime, FPLIMB, state);
}

void pre_montgomery() {
  mp_limb_t tmp1[FPLIMB + 1], tmp2[FPLIMB2 + 2];
  mpz_t tmp_z;
  mpz_t R;
  mpz_t R3_z;
  mp_limb_t R2[FPLIMB2 + 2];

  mpz_init(tmp_z);
  mpz_init(R);
  mpz_init(R3_z);

  for (int i = 0; i < FPLIMB; i++) N[i] = prime[i];
  mpz_ui_pow_ui(R, 2, FPLIMB_BITS);
  mpz_invert(tmp_z, prime_z, R);
  mpz_sub(tmp_z, R, tmp_z);
  integ_set_mpz(tmp1, tmp_z, FPLIMB);
  Ni_neg = tmp1[0];

  integ_set_mpz(tmp1, R, FPLIMB);
  integ_mod(tmp1, tmp1, FPLIMB + 1);
  mpn_copyd(RmodP, tmp1, FPLIMB);

  mpz_pow_ui(R3_z, R, 3);
  mpz_mod(R3_z, R3_z, prime_z);
  integ_set_mpz(R3, R3_z, FPLIMB);

  mpz_clear(tmp_z);
  mpz_clear(R);
  mpz_clear(R3_z);
}

void fpr_to_fp_mod(fp_t *ans, fpr_t *a) {
  integ_mod_montgomery(ans->x0, a->x0);
}

void fp_to_fpr(fpr_t *ans, fp_t *a) {
  integ_to_montgomery(ans->x0, a->x0);
}

void fpd_add(fpd_t *ans, fpd_t *a, fpd_t *b) {
#ifdef DEBUG_COST_A
  cost_add_nonmod_double++;
#endif
  mpn_add_n(ans->x0, a->x0, b->x0, FPLIMB2);
}

void fpr_add_mpn(fpr_t *ans, fpr_t *a, mp_limb_t *b) {
#ifdef DEBUG_COST_A
  cost_add_nonmod++;
#endif
  mpn_add_n(ans->x0, a->x0, b, FPLIMB);
}

void fpr_add_ui(fpr_t *ans, fpr_t *a, unsigned long int UI) {
#ifdef DEBUG_COST_A
  cost_add_ui++;
#endif
  integ_add_ui(ans->x0, a->x0, FPLIMB, UI);
  if (mpn_cmp(ans->x0, prime, FPLIMB) > 0) mpn_sub_n(ans->x0, ans->x0, prime, FPLIMB);
}

// void fp_add_mpn(fp_t *ans, fp_t *a, mp_limb_t *b) {
// #ifdef DEBUG_COST_A
//   cost_add++;
// #endif
//   mpn_add_n(ans->x0, a->x0, b, FPLIMB);
//   if (mpn_cmp(ans->x0, prime, FPLIMB) > 0) mpn_sub_n(ans->x0, ans->x0, prime, FPLIMB);
// }

void fpr_sub_mod(fpr_t *ans, fpr_t *a, fpr_t *b) {
#ifdef DEBUG_COST_A
  cost_sub++;
#endif
#ifdef DEBUG_ASSERT
  assert((mpn_cmp(a->x0, prime, FPLIMB) > 0) || (mpn_cmp(b->x0, prime, FPLIMB) > 0));
#endif
  static mp_limb_t buf[FPLIMB];

  if (mpn_cmp(a->x0, b->x0, FPLIMB) < 0) {
    mpn_sub_n(buf, a->x0, b->x0, FPLIMB);
    mpn_add_n(ans->x0, prime, buf, FPLIMB);
  } else {
    mpn_sub_n(ans->x0, a->x0, b->x0, FPLIMB);
  }
}

void fpr_sub_while(fp_t *ans, fp_t *a, fp_t *b) {
#ifdef DEBUG_COST_A
  cost_sub_nonmod++;
#endif

  if (mpn_cmp(a->x0, b->x0, FPLIMB) < 0) {
    mpn_sub_n(ans->x0, b->x0, a->x0, FPLIMB);
    while (mpn_cmp(ans->x0, prime, FPLIMB) >= 0) {
      mpn_sub_n(ans->x0, ans->x0, prime, FPLIMB);
    }
    mpn_sub_n(ans->x0, prime, ans->x0, FPLIMB);
  } else {
    mpn_sub_n(ans->x0, a->x0, b->x0, FPLIMB);
  }
}

void fpd_sub_mod_prime2(fpd_t *ans, fpd_t *a, fpd_t *b) {
#ifdef DEBUG_COST_A
  cost_sub_nonmod_double++;
#endif
#ifdef DEBUG_ASSERT
  assert((mpn_cmp(a->x0 + FPLIMB, prime, FPLIMB) > 0) || (mpn_cmp(b->x0 + FPLIMB, prime, FPLIMB) > 0));
#endif
  static mp_limb_t buf[FPLIMB2];

  if (mpn_cmp(a->x0, b->x0, FPLIMB2) < 0) {
    mpn_sub_n(ans->x0, a->x0, b->x0, FPLIMB2);
    mpn_add_n(ans->x0 + FPLIMB, ans->x0 + FPLIMB, prime, FPLIMB);
  } else {
    mpn_sub_n(ans->x0, a->x0, b->x0, FPLIMB2);
  }
}

void fpr_sub_ui(fpr_t *ans, fpr_t *a, unsigned long int UI) {
#ifdef DEBUG_COST_A
  cost_sub_ui++;
#endif
  if (UI == 0)
    fpr_set(ans, a);
  else
    integ_sub_ui(ans->x0, a->x0, FPLIMB, UI);
}
void fp_sub_mpn(fp_t *ans, fp_t *a, mp_limb_t *b) {
#ifdef DEBUG_COST_A
  cost_sub++;
#endif
  static mp_limb_t buf[FPLIMB];

  if (mpn_cmp(a->x0, b, FPLIMB) < 0) {
    mpn_sub_n(buf, a->x0, b, FPLIMB);
    mpn_add_n(ans->x0, prime, buf, FPLIMB);
  } else {
    mpn_sub_n(ans->x0, a->x0, b, FPLIMB);
  }
}
//remove fp_mod
void fp_inv(fp_t *ans, fp_t *a) {
#ifdef DEBUG_COST_A
  cost_inv++;
#endif
  static mp_limb_t prime_tmp[FPLIMB], gp[FPLIMB], sp[FPLIMB], buf[FPLIMB];
  static mp_size_t buf_size;

  mpn_init(sp, FPLIMB);
  mpn_copyd(prime_tmp, prime, FPLIMB);

  mpn_add_n(buf, a->x0, prime, FPLIMB);
  mpn_gcdext(gp, sp, &buf_size, buf, FPLIMB, prime_tmp, FPLIMB);

  if (buf_size < 0) {
    mpn_sub_n(ans->x0, prime, sp, FPLIMB);
  } else {
    mpn_copyd(ans->x0, sp, FPLIMB);
  }
}

void fp_inv_montgomery(fp_t *ans, fp_t *a) {
#ifdef DEBUG_COST_A
  cost_mul--;
  cost_mod--;
#endif
  fp_inv(ans, a);
  mpn_mulmod_montgomery(ans->x0, FPLIMB, ans->x0, FPLIMB, R3, FPLIMB);
}

int fp_legendre(fp_t *a) {
  int i;
  mpz_t tmp1, tmp2;
  fp_t tmp1_fp;
  mpz_init(tmp1);
  mpz_init(tmp2);
  fp_init(&tmp1_fp);

  mpz_sub_ui(tmp1, prime_z, 1);
  mpz_tdiv_q_ui(tmp2, tmp1, 2);
  fp_pow(&tmp1_fp, a, tmp2);

  if (mpn_cmp_ui(tmp1_fp.x0, FPLIMB, 1) == 0)
    i = 1;
  else if (mpn_cmp_ui(tmp1_fp.x0, FPLIMB, 0) == 0)
    i = 0;
  else
    i = -1;

  mpz_clear(tmp1);
  mpz_clear(tmp2);

  return i;
}

int fp_isCNR(fp_t *a) {
  fp_t tmp;
  fp_init(&tmp);
  mpz_t exp;
  mpz_init(exp);

  mpz_sub_ui(exp, prime_z, 1);
  mpz_tdiv_q_ui(exp, exp, 3);
  fp_pow(&tmp, a, exp);

  if (fp_cmp_one(&tmp) == 0) {
    mpz_clear(exp);
    return 1;
  } else {
    mpz_clear(exp);
    return -1;
  }
}
void fp_sqrt(fp_t *ans, fp_t *a) {
  fp_t x, y, t, k, n, tmp;
  fp_init(&x);
  fp_init(&y);
  fp_init(&t);
  fp_init(&k);
  fp_init(&n);
  fp_init(&tmp);
  unsigned long int e, m;
  mpz_t exp, q, z, result;
  mpz_init(exp);
  mpz_init(q);
  mpz_init(z);
  mpz_init(result);
  gmp_randstate_t state1;
  gmp_randinit_default(state1);
  gmp_randseed_ui(state1, (unsigned long)time(NULL));
  fp_set_random(&n, state1);

  while (fp_legendre(&n) != -1) {
    fp_set_random(&n, state1);
  }
  mpz_sub_ui(q, prime_z, 1);
  mpz_mod_ui(result, q, 2);
  e = 0;
  while (mpz_cmp_ui(result, 0) == 0) {
    mpz_tdiv_q_ui(q, q, 2);
    mpz_mod_ui(result, q, 2);
    e++;
  }
  fp_pow(&y, &n, q);
  mpz_set_ui(z, e);
  mpz_sub_ui(exp, q, 1);
  mpz_tdiv_q_ui(exp, exp, 2);
  fp_pow(&x, a, exp);
  fp_mul(&tmp, &x, &x);
  fp_mul(&k, &tmp, a);
  fp_mul(&x, &x, a);
  while (mpn_cmp_ui(k.x0, FPLIMB, 1) != 0) {
    m = 1;
    mpz_ui_pow_ui(exp, 2, m);
    fp_pow(&tmp, &k, exp);
    while (mpn_cmp_ui(tmp.x0, FPLIMB, 1) != 0) {
      m++;
      mpz_ui_pow_ui(exp, 2, m);
      fp_pow(&tmp, &k, exp);
    }
    mpz_sub_ui(exp, z, m);
    mpz_sub_ui(exp, exp, 1);
    mpz_ui_pow_ui(result, 2, mpz_get_ui(exp));
    fp_pow(&t, &y, result);
    fp_mul(&y, &t, &t);
    mpz_set_ui(z, m);
    fp_mul(&x, &x, &t);
    fp_mul(&k, &k, &y);
  }
  fp_set(ans, &x);

  mpz_clear(exp);
  mpz_clear(q);
  mpz_clear(z);
  mpz_clear(result);
}

void fp_pow(fp_t *ans, fp_t *a, mpz_t scalar) {
  int i, length;
  length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);
  fp_t tmp;
  fp_init(&tmp);  //not need?

  fp_set(&tmp, a);

  for (i = 1; i < length; i++) {
    fp_mul(&tmp, &tmp, &tmp);
    if (binary[i] == '1') {
      fp_mul(&tmp, a, &tmp);
    }
  }
  fp_set(ans, &tmp);
}

void fp_pow_montgomery(fp_t *ans, fp_t *a, mpz_t scalar) {
  int length = (int)mpz_sizeinbase(scalar, 2);
  char binary[length + 1];
  mpz_get_str(binary, 2, scalar);
  fp_t tmp;
  fp_init(&tmp);  //not need?

  fp_set(&tmp, a);

  for (int i = 1; i < length; i++) {
    fp_mulmod_montgomery(&tmp, &tmp, &tmp);
    if (binary[i] == '1') {
      fp_mulmod_montgomery(&tmp, a, &tmp);
    }
  }
  fp_set(ans, &tmp);
}

//TODO: consider modification "return mpn_cmp()"
int fp_cmp(fp_t *a, fp_t *b) {
  if (!mpn_cmp(a->x0, b->x0, FPLIMB))
    return 0;
  else
    return 1;
}

int fp_cmp_ui(fp_t *a, unsigned long int UI) {
  if (!mpn_cmp_ui(a->x0, FPLIMB, UI))
    return 0;
  else
    return 1;
}

int fp_cmp_mpn(fp_t *a, mp_limb_t *b) {
  if (!mpn_cmp(a->x0, b, FPLIMB))
    return 0;
  else
    return 1;
}

int fp_cmp_zero(fp_t *a) {
  if (!mpn_cmp_ui(a->x0, FPLIMB, 0))
    return 0;
  else
    return 1;
}

int fp_cmp_one(fp_t *a) {
  if (!mpn_cmp_ui(a->x0, FPLIMB, 1))
    return 0;
  else
    return 1;
}

int fp_montgomery_trick_montgomery(fp_t *A_inv, fp_t *a, int n) {
  int i;
  fp_t ans[n], ALL_inv;

  fp_set(ans, a);

  for (i = 1; i < n; i++) {
    fp_mulmod_montgomery(&ans[i], &ans[i - 1], &a[i]);
  }
  fp_inv_montgomery(&ALL_inv, &ans[n - 1]);
  for (i = n - 1; i > 0; i--) {
    fp_mulmod_montgomery(&A_inv[i], &ALL_inv, &ans[i - 1]);
    fp_mulmod_montgomery(&ALL_inv, &ALL_inv, &a[i]);
  }

  fp_set(A_inv, &ALL_inv);
  return 0;
}

void fp_lshift_ui_nonmod_single(fp_t *ans, fp_t *a, int s) {
#ifdef DEBUG_COST_A
  cost_add_ui++;
#endif
  mpn_lshift(ans->x0, a->x0, FPLIMB, s);
}

void fp_lshift_ui_nonmod_double(fpd_t *ans, fpd_t *a, int s) {
#ifdef DEBUG_COST_A
  cost_sqr++;
#endif
  mpn_lshift(ans->x0, a->x0, FPLIMB2, s);
}

int fp_legendre_sqrt(fp_t *ans, fp_t *a) {
  //need to 4|(p+1)
  fp_t C, D, A_tmp;
  int i;

  //legendre
  fp_pow(&C, a, sqrt_power_z);
  fp_mul(&D, &C, &C);
  fp_mul(&D, &D, a);

  if (mpn_cmp_ui(D.x0, FPLIMB, 1) == 0)
    i = 1;
  else if (mpn_cmp_ui(D.x0, FPLIMB, 0) == 0)
    return 0;
  else
    return -1;

  //sqrt
  fp_mul(ans, &C, a);
  return 1;
}

int fp_legendre_sqrt_montgomery(fp_t *ans, fp_t *a) {
  //need to 4|(p+1)

  if (fp_cmp_zero(a) == 0) {
    fp_set(ans, a);
    return 0;
  }

  fp_t C, D, A_tmp;

  //legendre
  fp_pow_montgomery(&C, a, sqrt_power_z);
  fp_mulmod_montgomery(&D, &C, &C);
  fp_mulmod_montgomery(&D, &D, a);
  if (fp_cmp_mpn(&D, RmodP) != 0) {
    return -1;
  }

  //sqrt
  fp_mulmod_montgomery(ans, &C, a);
  return 1;
}