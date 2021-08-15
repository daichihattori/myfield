#include "fpr.h"

void fpd_init(fpd_t *a) {
  integ_init(a->x0, FPLIMB2);
}

void fpd_printf(char *str, fpd_t *a) {
  integ_printf(str, a->x0, FPLIMB2);
}

void fpd_println(char *str, fp_t *a) {
  integ_println(str, a->x0, FPLIMB2);
}