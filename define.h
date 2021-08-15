#ifndef DEFINE_H
#define DEFINE_H

#include <assert.h>
#include <gmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

/*************OPTION*************/

//bit
#define ARCBIT 64  //64bit processor
//#define ARCBIT 32 //32bit processor

//debug
//#define DEBUG_COST_A
#define DEBUG_ASSERT

#define PARAM_TAXONOMY_CHANGE_B

/* x-mood*/
//taxonomy cahnge b
#ifdef PARAM_TAXONOMY_CHANGE_B
#define X_MINUS
#define TWIST_PHI
#define CURVE_B_8_3_TYPE1
#define exp_type_2
#define PRIME_BIT 461
#define ORDER_BIT 308
#define X_BIT 77

#endif

/************************************/

#define scalar_t mpz_t

#define FPLIMB (PRIME_BIT / ARCBIT + 1)
#define FPLIMB2 FPLIMB * 2
#define FPLIMB_BITS FPLIMB* ARCBIT

#define FRLIMB (ORDER_BIT / ARCBIT + 1)
#define FRLIMB2 FRLIMB * 2
#define FRLIMB_BITS FRLIMB* ARCBIT

#define FXLIMB (X_BIT / ARCBIT + 1)
#define FXLIMB2 ((X_BIT * 2) / ARCBIT + 1)
#define FXLIMB_BITS FXLIMB* ARCBIT
#define bls12_X_length X_BIT
#define bls12_X2_length X_BIT - 1

int cost_add,
    cost_add_ui, cost_sub, cost_sub_ui, cost_mul, cost_set_neg, cost_r1shift, cost_sqr, cost_inv, cost_mod;
int cost_add_nonmod, cost_add_nonmod_double, cost_sub_nonmod, cost_sub_nonmod_double, cost_mod_nomal;

typedef struct {
  int add;
  int add_ui;
  int add_nonmod;
  int add_nonmod_double;
  int sub;
  int sub_ui;
  int sub_nonmod;
  int sub_nonmod_double;
  int mul;
  int set_neg;
  int r1shift;
  int sqr;
  int inv;
  int mod;
  int mod_nomal;
} cost;

/*============================================================================*/
/* Field                                                                      */
/*============================================================================*/
mp_limb_t buf[FPLIMB];

typedef struct {
  mp_limb_t x0[FPLIMB];
} fp_t;
typedef struct {
  mp_limb_t x0[FPLIMB];
} fpr_t;
typedef struct {
  mp_limb_t x0[FPLIMB2];
} fpd_t;
gmp_randstate_t state;
mpz_t X_z, prime_z, order_z, trace_z, X_abs_z;
mp_limb_t X_abs[FXLIMB], X2[FXLIMB2], prime[FPLIMB], order[FRLIMB], trace[FPLIMB];
mp_limb_t prime_carry[FPLIMB];
mp_limb_t prime2[FPLIMB2];
mp_limb_t epsilon1[FPLIMB], epsilon2[FPLIMB];
fp_t epsilon1_montgomery, epsilon2_montgomery;
mpz_t efp_total, efp12_total;
mp_limb_t curve_b[FPLIMB];
mpz_t sqrt_power_z;
mpz_t g1_power;
mpz_t X_mod_order_z;

mp_limb_t N[FPLIMB2], RmodP[FPLIMB], R3[FPLIMB];
mp_limb_t Ni_neg;  //Ni_neg=-N^(-1)
/*============================================================================*/
/* Test functions                                                             */
/*============================================================================*/
struct timeval tv_start, tv_end;
float MILLER_OPT_AFFINE, FINALEXP_OPT_AFFINE;
float MILLER_OPT_PROJECTIVE, FINALEXP_OPT_PROJECTIVE;

cost MILLER_OPT_AFFINE_COST, FINALEXP_OPT_AFFINE_COST;
cost MILLER_OPT_PROJECTIVE_COST, FINALEXP_OPT_PROJECTIVE_COST;

mpz_t to_g1_expo, to_g2_expo;
fp_t curve_b_montgomery;
fp_t inv2_montgomery;

#endif
