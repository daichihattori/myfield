#include "define.h"
void integ_init(mp_limb_t *a, mp_size_t size) {
  mpn_zero(a, size);
}

void integ_printf(char *str, mp_limb_t *a, mp_size_t size) {
  gmp_printf("%sNu", str, a, size);
}

void integ_println(char *str, mp_limb_t *a, mp_size_t size) {
  gmp_printf("%sNu\n", str, a, size);
}

void integ_set(mp_limb_t *ans, mp_limb_t *a, mp_size_t size) {
  mpn_copyd(ans, a, size);
}

void integ_set_ui(mp_limb_t *ans, mp_size_t size, unsigned long int ui) {
  unsigned long int i;

  ans[0] = ui;

  for (i = 1; i < size; i++) {
    ans[i] = 0;
  }
}

void integ_add_mod(mp_limb_t *ans, mp_limb_t *a, mp_limb_t *b) {
#ifdef DEBUG_COST_A
  cost_add++;
#endif
#ifdef DEBUG_ASSERT
  assert((mpn_cmp(a, prime, FPLIMB) > 0) || (mpn_cmp(b, prime, FPLIMB) > 0));
#endif
  mpn_add_n(ans, a, b, FPLIMB);
  if (mpn_cmp(ans, prime, FPLIMB) >= 0) mpn_sub_n(ans, ans, prime, FPLIMB);
}

void integ_add(mp_limb_t *ans, mp_limb_t *a, mp_limb_t *b, mp_size_t size) {
#ifdef DEBUG_COST_A
  cost_add++;
#endif
  mpn_add_n(ans, a, b, size);
}

void integ_add_mod_limb(mp_limb_t *ans, mp_limb_t *a, mp_limb_t *b) {
#ifdef DEBUG_COST_A
  cost_sub_nonmod_double++;
#endif
#ifdef DEBUG_ASSERT
  assert((mpn_cmp(a + FPLIMB, prime, FPLIMB) > 0) || (mpn_cmp(b + FPLIMB, prime, FPLIMB) > 0));
#endif
  static mp_limb_t buf[FPLIMB2];
  mpn_add_n(ans, a, b, FPLIMB);
  if (mpn_cmp(ans + FPLIMB, prime, FPLIMB) >= 0) mpn_sub_n(ans + FPLIMB, ans + FPLIMB, prime, FPLIMB);
}

void integ_sub_mod(mp_limb_t *ans, mp_limb_t *a, mp_limb_t *b) {
#ifdef DEBUG_COST_A
  cost_sub++;
#endif
#ifdef DEBUG_ASSERT
  assert((mpn_cmp(a, prime, FPLIMB) > 0) || (mpn_cmp(b, prime, FPLIMB) > 0));
#endif
  static mp_limb_t buf[FPLIMB];

  if (mpn_cmp(a, b, FPLIMB) < 0) {
    mpn_sub_n(buf, a, b, FPLIMB);
    mpn_add_n(a, prime, buf, FPLIMB);
  } else {
    mpn_sub_n(a, a, b, FPLIMB);
  }
}

void integ_sub_while(mp_limb_t *ans, mp_limb_t *a, mp_limb_t *b) {
#ifdef DEBUG_COST_A
  cost_sub_nonmod++;
#endif

  if (mpn_cmp(a, b, FPLIMB) < 0) {
    mpn_sub_n(ans, b, a, FPLIMB);
    while (mpn_cmp(ans, prime, FPLIMB) >= 0) {
      mpn_sub_n(ans, ans, prime, FPLIMB);
    }
    mpn_sub_n(ans, prime, ans, FPLIMB);
  } else {
    mpn_sub_n(ans, a, b, FPLIMB);
  }
}

void fp_sub_nonmod_double(fpd_t *ANS, fpd_t *A, fpd_t *B) {
#ifdef DEBUG_COST_A
  cost_sub_nonmod_double++;
#endif
#ifdef DEBUG_ASSERT
  assert((mpn_cmp(A->x0 + FPLIMB, prime, FPLIMB) > 0) || (mpn_cmp(B->x0 + FPLIMB, prime, FPLIMB) > 0));
#endif
  static mp_limb_t buf[FPLIMB2];

  if (mpn_cmp(A->x0, B->x0, FPLIMB2) < 0) {
    mpn_sub_n(ANS->x0, A->x0, B->x0, FPLIMB2);
    mpn_add_n(ANS->x0 + FPLIMB, ANS->x0 + FPLIMB, prime, FPLIMB);
  } else {
    mpn_sub_n(ANS->x0, A->x0, B->x0, FPLIMB2);
  }
}

void integ_l1shift(mp_limb_t *ans, mp_limb_t *a, mp_size_t size) {
#ifdef DEBUG_COST_A
  cost_add++;
#endif
  mpn_lshift(ans, a, size, 1);
  if (mpn_cmp(ans, prime, size) >= 0) mpn_sub_n(ans, ans, prime, size);
}

void integ_r1shift(mp_limb_t *ans, mp_limb_t *a) {
#ifdef DEBUG_COST_A
  cost_add++;
#endif
  if (a[0] & 1)
    mpn_add_n(ans, a, prime, FPLIMB);
  else
    mpn_copyd(ans, a, FPLIMB);
  mpn_rshift(ans, ans, FPLIMB, 1);
}

void integ_set_random(mp_limb_t *ans, mp_limb_t *max, mp_size_t size, gmp_randstate_t state) {
  mpz_t tmp;
  mpz_init(tmp);
  mpz_urandomm(tmp, state, prime_z);
  integ_set_mpz(ans, tmp, FPLIMB);
  mpz_clear(tmp);
}

void integ_set_mpz(mp_limb_t *ans, mpz_t a, mp_size_t size) {
  char *str;

  str = (char *)malloc(mpz_sizeinbase(a, 10) + 2);

  //gmp_printf("a=%Zd\n",a);
  str = mpz_get_str(str, 10, a);
  //printf("str1=%s",str);
  integ_set_char(ans, size, str);

  free(str);
}

void integ_set_neg(mp_limb_t *ans, mp_limb_t *a, mp_size_t size) {
  if (mpn_cmp_ui(a, 0, size) == 0)
    mpn_copyd(ans, a, size);
  else
    mpn_sub_n(ans, prime, a, size);
}

void integ_mod_normal(mp_limb_t *ans, mp_limb_t *a, mp_size_t size_a) {
  mp_limb_t dumy[size_a];
  mpn_tdiv_qr(dumy, ans, 0, a, size_a, prime, FPLIMB);
}

int mpn_cmp_ui(mp_limb_t *a, mp_size_t size, unsigned long int ui) {
  integ_set_ui(buf, size, ui);
  if (mpn_cmp(a, buf, size) == 0) {
    return 0;
  } else {
    return 1;
  }
}
void integ_lshift_ext(mp_limb_t *ans, mp_limb_t *a, mp_size_t size, long int L) {
  mp_limb_t tmp[size];

  mpn_copyd(tmp, a, size);
  while (L > 63) {
    mpn_lshift(tmp, tmp, size, 63);
    L = L - 63;
  }
  mpn_lshift(ans, tmp, size, L);
}
void integ_add_ui(mp_limb_t *ans, mp_limb_t *a, mp_size_t size, unsigned long int ui) {
  mp_limb_t buf[size];

  integ_set_ui(buf, size, ui);

  mpn_add_n(ans, a, buf, size);
}
void integ_sub_ui(mp_limb_t *ans, mp_limb_t *a, mp_size_t size, unsigned long int ui) {
  mp_limb_t buf[size];

  integ_set_ui(buf, size, ui);

  mpn_sub_n(ans, a, buf, size);
}
void integ_mul(mp_limb_t *ans, mp_limb_t *a, mp_size_t a_size, mp_limb_t *b, mp_size_t b_size) {
  mpn_mul(ans, a, a_size, b, b_size);
}
void integ_mul_ui(mp_limb_t *ans, mp_limb_t *a, mp_size_t size, unsigned long int ui) {
  mp_limb_t buf[size];

  integ_set_ui(buf, size, ui);

  mpn_mul_n(ans, a, buf, size);
}
void integ_sqr(mp_limb_t *ans, mp_limb_t *a, mp_size_t a_size) {
#ifdef DEBUG_COST_A
  cost_sqr++;
#endif
  mpn_sqr(ans, a, a_size);
}
void integ_invert(mp_limb_t *ans, mp_limb_t *a, mp_limb_t *p) {
  mp_limb_t prime_tmp[FPLIMB], gp[FPLIMB], sp[FPLIMB], tmp[FPLIMB];
  mp_size_t buf_size;

  integ_init(gp, FPLIMB);
  integ_init(sp, FPLIMB);
  integ_init(tmp, FPLIMB);
  integ_init(prime_tmp, FPLIMB);

  mpn_copyd(prime_tmp, p, FPLIMB);

  mpn_add_n(buf, a, p, FPLIMB);
  mpn_gcdext(gp, sp, &buf_size, buf, FPLIMB, prime_tmp, FPLIMB);

  if (buf_size < 0) {
    mpn_sub_n(tmp, p, sp, FPLIMB);
  } else {
    mpn_copyd(tmp, sp, FPLIMB);
  }

  integ_mod_normal(ans, tmp, FPLIMB);
}

void integ_mulmod_montgomery(mp_limb_t *ans, mp_limb_t *a, mp_limb_t *b) {
#ifdef DEBUG_COST_A
  cost_mul++;
  cost_mod++;
#endif

  static mp_limb_t T[FPLIMB2];
  mpn_zero(T, FPLIMB2);

  mpn_mul(T, a, FPLIMB, b, FPLIMB);
  for (int i = 0; i < FPLIMB; i++)
    T[i] = mpn_addmul_1(&T[i], prime, FPLIMB, T[i] * Ni_neg);

  mpn_add_n(ans, T + FPLIMB, T, FPLIMB);
  if (mpn_cmp(ans, prime, FPLIMB) != -1) mpn_sub_n(ans, ans, prime, FPLIMB);
}

void integ_sqrmod_montgomery(mp_limb_t *ans, mp_limb_t *a) {
  static mp_limb_t T[FPLIMB2];
  mpn_zero(T, FPLIMB2);

  mpn_sqr(T, a, FPLIMB);
  for (int i = 0; i < FPLIMB; i++)
    T[i] = mpn_addmul_1(&T[i], prime, FPLIMB, T[i] * Ni_neg);

  mpn_add_n(ans, T + FPLIMB, T, FPLIMB);
  if (mpn_cmp(ans, prime, FPLIMB) != -1) mpn_sub_n(ans, ans, prime, FPLIMB);
}

void integ_mod_montgomery(mp_limb_t *ans, mp_limb_t *a) {
  static mp_limb_t T[FPLIMB2];
  mpn_zero(T, FPLIMB2);

  mpn_copyd(T, a, FPLIMB);
  for (int i = 0; i < FPLIMB; i++)
    T[i] = mpn_addmul_1(&T[i], prime, FPLIMB, T[i] * Ni_neg);

  mpn_add_n(ans, T + FPLIMB, T, FPLIMB);
  while (mpn_cmp(ans, prime, FPLIMB) > 0) mpn_sub_n(ans, ans, prime, FPLIMB);
}

void integ_to_montgomery(mp_limb_t *ans, mp_limb_t *a) {
#ifdef DEBUG_COST_A
  //cost_mod++;
  cost_mod_nomal++;
#endif
  static int i;
  static mp_limb_t tmp[FPLIMB2];
  mpn_zero(tmp, FPLIMB2);
  for (i = FPLIMB; i < FPLIMB2; i++) tmp[i] = a[i - FPLIMB];
  integ_mod_normal(ans, tmp, FPLIMB2);
}

void integ_set_char(mp_limb_t *ans, mp_size_t mp_size, char *str) {
  unsigned long int i, sizeL;
  char *str_buf;
  mp_size_t size;

  sizeL = strlen(str);

  str_buf = (char *)malloc(sizeL * sizeof(char));

  for (i = 0; i < sizeL; i++) {
    str_buf[i] = str[i] - 48;
  }

  size = mpn_set_str(ans, (unsigned const char *)str_buf, sizeL, 10);
  for (i = size; i < mp_size; i++) {
    ans[i] = 0;
  }

  free(str_buf);
}
