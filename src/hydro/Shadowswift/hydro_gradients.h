//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H

/**
 * @brief Initialize gradient variables.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part* p) {
  p->gradients.rho[0] = 0.0f;
  p->gradients.rho[1] = 0.0f;
  p->gradients.rho[2] = 0.0f;

  p->gradients.v[0][0] = 0.0f;
  p->gradients.v[0][1] = 0.0f;
  p->gradients.v[0][2] = 0.0f;

  p->gradients.v[1][0] = 0.0f;
  p->gradients.v[1][1] = 0.0f;
  p->gradients.v[1][2] = 0.0f;

  p->gradients.v[2][0] = 0.0f;
  p->gradients.v[2][1] = 0.0f;
  p->gradients.v[2][2] = 0.0f;

  p->gradients.P[0] = 0.0f;
  p->gradients.P[1] = 0.0f;
  p->gradients.P[2] = 0.0f;

  p->gradients.A[0] = 0.0f;
  p->gradients.A[1] = 0.0f;
  p->gradients.A[2] = 0.0f;

#ifdef SHADOWSWIFT_GRADIENTS_WLS
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->gradients.matrix_wls[i][j] = 0.f;
    }
  }
#endif
}

#ifdef SHADOWSWIFT_GRADIENTS

/**
 * @brief Gradients time extrapolation (makes scheme second order in time).
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_extrapolate_in_time(const struct part* p, const float* W,
                                    float dt, float* dW) {

  float drho[3], dvx[3], dvy[3], dvz[3], dP[3], dA[3];
  hydro_part_get_gradients(p, drho, dvx, dvy, dvz, dP, dA);
  const float div_v = dvx[0] + dvy[1] + dvz[2];

  dW[0] =
      -dt * (W[0] * div_v + W[1] * drho[0] + W[2] * drho[1] + W[3] * drho[2]);

  if (W[0] != 0.0f) {
    const float rho_inv = 1.f / W[0];
    dW[1] = -dt * (W[1] * div_v + rho_inv * dP[0]);
    dW[2] = -dt * (W[2] * div_v + rho_inv * dP[1]);
    dW[3] = -dt * (W[3] * div_v + rho_inv * dP[2]);
  } else {
    dW[1] = 0.0f;
    dW[2] = 0.0f;
    dW[3] = 0.0f;
  }
  dW[4] = -dt * (hydro_gamma * W[4] * div_v + W[1] * dP[0] + W[2] * dP[1] +
                 W[3] * dP[2]);
  dW[5] = -dt * (W[1] * dA[0] + W[2] * dA[1] + W[3] * dA[1]);
}

/**
 * @brief Extrapolate the given gradient over the given distance.
 *
 * @param gradient Gradient of a quantity.
 * @param dx Distance vector.
 * @return Change in the quantity after a displacement along the given distance
 * vector.
 */
__attribute__((always_inline)) INLINE static float
hydro_gradients_extrapolate_single_quantity(const float* gradient,
                                            const float* dx) {

  return gradient[0] * dx[0] + gradient[1] * dx[1] + gradient[2] * dx[2];
}

/**
 * @brief Extrapolate all the quantities over the given distance.
 *
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_extrapolate(
    const struct part* p, const float* dx, float* dW) {

  float drho[3], dvx[3], dvy[3], dvz[3], dP[3], dA[3];
  hydro_part_get_gradients(p, drho, dvx, dvy, dvz, dP, dA);

  dW[0] = hydro_gradients_extrapolate_single_quantity(drho, dx);
  dW[1] = hydro_gradients_extrapolate_single_quantity(dvx, dx);
  dW[2] = hydro_gradients_extrapolate_single_quantity(dvy, dx);
  dW[3] = hydro_gradients_extrapolate_single_quantity(dvz, dx);
  dW[4] = hydro_gradients_extrapolate_single_quantity(dP, dx);
  dW[5] = hydro_gradients_extrapolate_single_quantity(dA, dx);
}

__attribute__((always_inline)) INLINE static void
hydro_gradients_apply_extrapolation(float* W, const float* restrict dW) {
  if (-dW[0] > W[0]) {
    message(
        "Gradient extrapolation would lead to unphysical Density! "
        "Falling back to first order for this particle!");
  } else if (-dW[4] > W[4]) {
    message(
        "Gradient extrapolation would lead to unphysical Pressure! "
        "Falling back to first order for this particle!");
  } else if (-dW[5] > W[5]) {
    message(
        "Gradient extrapolation would lead to unphysical Entropy! "
        "Falling back to first order for this particle!");
  } else {
    W[0] += dW[0];
    W[1] += dW[1];
    W[2] += dW[2];
    W[3] += dW[3];
    W[4] += dW[4];
    W[5] += dW[5];
  }
}

#ifndef SHADOWSWIFT_GRADIENTS_WLS

/* Green-Gauss gradient estimates (default) */
#include "hydro_gradients_shadowswift.h"

#else

/* Weighted least square gradient estimates */
#include "hydro_gradients_wls.h"

#endif

#else

/* No gradients. Perfectly acceptable, but we have to provide empty functions */
#define HYDRO_GRADIENT_IMPLEMENTATION "No gradients (first order scheme)"

/**
 * @brief Update the gradient estimation for a particle using a given
 * neighbouring particle.
 *
 * @param pi Particle we are updating.
 * @param pj Particle we are using to update pi.
 */
__attribute__((always_inline)) INLINE void hydro_gradients_collect(
    struct part* restrict pi, const struct part* restrict pj,
    const double* restrict centroid, const double* restrict dx, double r,
    float surface_area) {}

/**
 * @brief Extrapolate all the quantities over the given distance.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_extrapolate(
    const struct part* p, const float* dx, float* dW) {
  dW[0] = 0.f;
  dW[1] = 0.f;
  dW[2] = 0.f;
  dW[3] = 0.f;
  dW[4] = 0.f;
}

/**
 * @brief Gradients reconstruction. Empty for no gradients.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    struct part* restrict pi, struct part* restrict pj, const float* dx,
    float r, const float* xij_i, float dt, float* Wi, float* Wj) {}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part* p) {}
#endif


#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_GRADIENTS_H
