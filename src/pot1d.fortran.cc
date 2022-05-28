#include <tuple>
extern "C" {
void pot_(double *, double *, double *, double *);
void pot(double r, double *v, double *dv, double *d2v) {
  pot_(&r, v, dv, d2v);
  return;
}
}
