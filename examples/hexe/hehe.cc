// External function
extern "C" double v_mlr_he_he_(double *, double *, double *, double *);

// Public variable (as cache)
static double old_r(-1.0);
static double v(0.0);
static double dv(0.0);
static double d2v(0.0);

// Constant to tranform between units
const double cm2k = 1.438777356765961;

extern "C" double value(double r) {
  if (old_r != r){
    v_mlr_he_he_(&r, &v, &dv, &d2v);
    old_r = r;
  }
  return v * cm2k;
}

extern "C" double derivative(double r) {
  if (old_r != r){
    v_mlr_he_he_(&r, &v, &dv, &d2v);
    old_r = r;
  }
  return dv * cm2k;
}
