#include "fp.h"

void fpr_init(fpr_t *a) {
  mpn_zero(a->x0, FPLIMB);
}

void fpr_set_neg(fp_t *ans, fp_t *a) {
  integ_set_neg(ans->x0, a->x0, FPLIMB);
}
void fpr_mod_printf(char *str, fpr_t *a) {
  static fp_t out;
  fpr_to_fp_mod(&out, a);
  gmp_printf("%s%Nu", str, out.x0, FPLIMB);
}

void fpr_mod_println(char *str, fpr_t *a) {
  static fp_t out;
  fpr_to_fp_mod(&out, a);
  gmp_printf("%s%Nu\n", str, out.x0, FPLIMB);
}

void fpr_set(fp_t *ans, fp_t *a) {
  mpn_copyd(ans->x0, a->x0, FPLIMB);
}

void fpr_set_mpn(fpr_t *ans, mp_limb_t *a) {
  fp_t tmp;
  mpn_copyd(ans->x0, a, FPLIMB);
  fp_to_fpr(ans, &tmp);
}

void fpr_l1shift(fpr_t *ans, fpr_t *a) {
  integ_l1shift(ans->x0, a->x0, FPLIMB);
}

void fpr_r1shift(fpr_t *ans, fpr_t *a) {
  integ_l1shift(ans->x0, a->x0, FPLIMB);
}

void fpr_mulmod(fpr_t *ans, fpr_t *a, fpr_t *b) {
  integ_mulmod_montgomery(ans->x0, a->x0, b->x0);
}

void fpr_sqrmod(fpr_t *ans, fpr_t *a) {
  integ_sqrmod_montgomery(ans->x0, a->x0);
}

void fpr_mul(fpd_t *ans, fpr_t *a, fpr_t *b) {
  integ_mul(ans->x0, a->x0, FPLIMB, b->x0, FPLIMB);
}

void fpr_mul_ui(fpr_t *ans, fpr_t *a, unsigned long int UI) {
  integ_mul_ui(ans->x0, a->x0, FPLIMB, UI);
}

void fpr_sqr(fpd_t *ans, fpr_t *a) {
  integ_sqr(ans->x0, a->x0, FPLIMB);
}

void fpr_add_mod(fpr_t *ans, fpr_t *a, fpr_t *b) {
#ifdef DEBUG_COST_A
  cost_add++;
#endif
#ifdef DEBUG_ASSERT
  assert((mpn_cmp(a->x0, prime, FPLIMB) > 0) || (mpn_cmp(b->x0, prime, FPLIMB) > 0));
#endif
  mpn_add_n(ans->x0, a->x0, b->x0, FPLIMB);
  if (mpn_cmp(ans->x0, prime, FPLIMB) >= 0) mpn_sub_n(ans->x0, ans->x0, prime, FPLIMB);
}

void fpr_add(fpr_t *ans, fpr_t *a, fpr_t *b) {
#ifdef DEBUG_COST_A
  cost_add_nonmod++;
#endif
  mpn_add_n(ans->x0, a->x0, b->x0, FPLIMB);
}