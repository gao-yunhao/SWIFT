/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Loic Hausammann (loic.hausammann@epfl.ch)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_GEAR_FEEDBACK_IACT_H
#define SWIFT_GEAR_FEEDBACK_IACT_H

/* Local includes */
#include "feedback.h"
#include "hydro.h"
#include "random.h"
#include "timestep_sync_part.h"
#include <math.h>


/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * In GEAR, this function does nothing. What we need is the
 * star->density.wcount computed in runner_iact_nonsym_stars_density().
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sparticle.
 * @param pj Second particle (not updated).
 * @param xpj Extra particle data (not updated).
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time value
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_density(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct spart *si, const struct part *pj,
                                    const struct xpart *xpj,
                                    const struct cosmology *cosmo,
                                    const struct feedback_props *fb_props,
                                    const integertime_t ti_current) {}

/**
 * @brief Prepare the feedback by computing the required quantities (loop 1). 
 * Used for updating properties of star particles required for the feedback. 
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (updated).
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep1(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {

  /* Accumulate the sum in the numerator and denominator of f_plus and f_minus */
  double dx_ij_plus[3];
  double dx_ij_minus[3];
  double scalar_weight_j = feedback_compute_scalar_weight(r2, dx, hi, hj, si, pj,
							  dx_ij_plus, dx_ij_minus);

  si->feedback_data.f_plus_num[0] += scalar_weight_j*fabs(dx_ij_minus[0]);
  si->feedback_data.f_plus_num[1] += scalar_weight_j*fabs(dx_ij_minus[1]);
  si->feedback_data.f_plus_num[2] += scalar_weight_j*fabs(dx_ij_minus[2]);

  si->feedback_data.f_plus_denom[0] += scalar_weight_j*fabs(dx_ij_plus[0]);
  si->feedback_data.f_plus_denom[1] += scalar_weight_j*fabs(dx_ij_plus[1]);
  si->feedback_data.f_plus_denom[2] += scalar_weight_j*fabs(dx_ij_plus[2]);


  si->feedback_data.f_minus_num[0] += scalar_weight_j*fabs(dx_ij_plus[0]);
  si->feedback_data.f_minus_num[1] += scalar_weight_j*fabs(dx_ij_plus[1]);
  si->feedback_data.f_minus_num[2] += scalar_weight_j*fabs(dx_ij_plus[2]);

  si->feedback_data.f_minus_denom[0] += scalar_weight_j*fabs(dx_ij_minus[0]);
  si->feedback_data.f_minus_denom[1] += scalar_weight_j*fabs(dx_ij_minus[1]);
  si->feedback_data.f_minus_denom[2] += scalar_weight_j*fabs(dx_ij_minus[2]);
}


/**
 * @brief Prepare the feedback by computing the required quantities (loop 2). 
 * Used for updating properties of star particles required for the feedback.
 *
 * In GEAR, we update the enrichment weight. 
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (updated).
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep2(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {

  /* Now we can compute f_plus and f_minus for the star */
  double f_plus_i[3], f_minus_i[3], w_j[3];
  feedback_compute_vector_weight_non_normalized(r2, dx, hi, hj, si, pj, f_plus_i,
						f_minus_i, w_j);

  /* Accumulate w_j norm for later */
  const double w_j_norm_2 = w_j[0]*w_j[0] + w_j[1]*w_j[1] + w_j[2]*w_j[2];  
  si->feedback_data.enrichment_weight += sqrt(w_j_norm_2);
}


#ifdef EXTRA_STAR_LOOPS_2
/**
 * @brief Prepare the feedback by computing the required quantities (loop 3). 
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (updated).
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep3(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {
  message("prep3");
}

/**
 * @brief Prepare the feedback by computing the required quantities (loop 4). 
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (updated).
 * @param pj Second (gas) particle (not updated).
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_prep4(const float r2, const float dx[3],
                                  const float hi, const float hj,
                                  struct spart *si, const struct part *pj,
                                  const struct xpart *xpj,
                                  const struct cosmology *cosmo,
                                  const integertime_t ti_current) {
  message("prep4");
}
#endif /* EXTRA_STAR_LOOPS_2 */

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 * Used for updating properties of gas particles neighbouring a star particle
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (si - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First (star) particle (not updated).
 * @param pj Second (gas) particle.
 * @param xpj Extra particle data
 * @param cosmo The cosmological model.
 * @param fb_props Properties of the feedback scheme.
 * @param ti_current Current integer time used value for seeding random number
 * generator
 */
#ifndef EXTRA_STAR_LOOPS_2
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const struct phys_const* phys_const,
    const struct unit_system* us, const integertime_t ti_current) {

  const double E_ej = si->feedback_data.energy_ejected;

  /* Do we have supernovae? */
  if (E_ej == 0) {
    return;
  }

  if (si->feedback_data.mass_ejected < 200e-10) {
    return;
  }

  /* Finally, we can compute the w_j_bar. */
  double f_plus_i[3], f_minus_i[3], w_j[3];
  feedback_compute_vector_weight_non_normalized(r2, dx, hi, hj, si, pj, f_plus_i, f_minus_i, w_j);

  /* The normalized vector weight */
  double w_j_bar[3] =  {w_j[0]/si->feedback_data.enrichment_weight,
			w_j[1]/si->feedback_data.enrichment_weight,
			w_j[2]/si->feedback_data.enrichment_weight};
  const double w_j_bar_norm_2 = w_j_bar[0]*w_j_bar[0] + w_j_bar[1]*w_j_bar[1] + w_j_bar[2]*w_j_bar[2];
  const double w_j_bar_norm = sqrt(w_j_bar_norm_2);

  /* Here just get the feedback properties we want to distribute */
  const float mj = hydro_get_mass(pj);
  const double m_ej = si->feedback_data.mass_ejected;
  const double p_ej = sqrt(2*m_ej*E_ej) ;

  /* Distribute mass... (the max avoids to have dm=0 and NaN by dividing by dm) */
  double dm = max(w_j_bar_norm * m_ej, FLT_MIN);
  xpj->feedback_data.delta_mass += dm;

  /* ... metals */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    pj->chemistry_data.metal_mass[i] +=
	w_j_bar_norm * si->feedback_data.metal_mass_ejected[i];
  }

  /* ... momentum */
  const double dp[3] = {w_j_bar[0]*p_ej, w_j_bar[1]*p_ej, w_j_bar[2]*p_ej};
  const double dE = w_j_bar_norm * E_ej;

  /* Now boost to the 'laboratory' frame */
  double dp_prime[3] = {dp[0] + dm*si->v[0], dp[1] + dm*si->v[1], dp[2] + dm*si->v[2]};

  /* Note: The momentum is given to xpj later... We need to do some physics before. */

  /* ... Total energy */
  const double dp_norm_2 = dp[0]*dp[0] +  dp[1]*dp[1] +  dp[2]*dp[2];
  const double dp_prime_norm_2 = dp_prime[0]*dp_prime[0] +  dp_prime[1]*dp_prime[1] +  dp_prime[2]*dp_prime[2];
  const double dE_prime = dE + 1.0/(2.0*dm) * (dp_prime_norm_2 - dp_norm_2);

  /* ... internal energy */
  const double new_mass = mj + dm;

  /* Compute kinetic energy difference before and after SN */
  const double p_old_norm_2 = mj*mj*(xpj->v_full[0]*xpj->v_full[0] + xpj->v_full[1]*xpj->v_full[1] + xpj->v_full[2]*xpj->v_full[2]);
  const double p_new[3] = {mj*xpj->v_full[0] + dp_prime[0],
			   mj*xpj->v_full[1] + dp_prime[1],
			   mj*xpj->v_full[2] + dp_prime[2]};
  const double p_new_norm_2 = p_new[0]*p_new[0] + p_new[1]*p_new[1] + p_new[2]*p_new[2];

  const double E_kin_old = p_old_norm_2/(2.0*mj);
  const double E_kin_new = p_new_norm_2/(2.0*new_mass);
  const double dKE = E_kin_new - E_kin_old;

  const double U_old = xpj->u_full*mj;
  const double E_old = U_old + E_kin_old;
  const double E_new = E_old + dE_prime;
  const double U_new = E_new - E_kin_new;

  /* Compute the internal energy */
  double dU = U_new - U_old;
  /* message("dE_prime= %e, dU = %e,  dKE = %e, p_new_2 = %e, p_old_2 = %e", dE_prime, dU, dKE, p_new_norm_2, p_old_norm_2); */
  /* message("E_old = %e, E_new = %e, E_kin_old = %e, E_kin_new = %e, U_old = %e, U_new = %e", E_old, E_new, E_kin_old, E_kin_old, U_old, U_new); */

  /* --Now, we take into account for potentially unresolved energy-conserving
     phase of the SN explosion-- */

  /* Now, momentum given to the particle is not simply dp_prime. This is the
     momentum the location of the SN explosion. The ejecta must travel to the
     particle at distance r, i.e it must sweep up the mass mj of the particle
     as in a shock or a shell.
     This means that there is some work PdV that is done to reach the particle,
     work that converts thermal energy to kinetic energy. Therefore, the
     correct momentum to couple is proportional to the following factor: */
  const double PdV_work_fraction = sqrt(1 + mj/dm);

  /* At the end of the energy conserving phase, the blastwave reaches its
     terminal momentum. */
  const double p_terminal = feedback_get_SN_terminal_momentum(si, pj, xpj, phys_const, us);

  /* During the Taylor Sedov phase (energy conserving phase), dp_prime <=
     p_terminal. Thus the momentum coupling is irrelevant. We couple the
     thermal energy and leave the hydro solver correctly determines the PdV work. */

  /* If we can resolve the Taylor Sedov, then we give the right coupled
     momentum (which is by definition <= p_terminal). If we cannot resolve it,
     then we have reached the p_terminal. This is the upper limit of momentum,
     since afterwards the cooling is efficient and the thermal energy is radiated
     away.
     Thus, the factor to multiply dp_prime is: */
  const double p_factor = min(PdV_work_fraction, p_terminal/p_ej);

  /* if (p_factor == p_terminal/p_ej) { */
  /*   message("We do not resolve the Sedov-Taylor. Using p_terminal. p_t = %e, p_ej=%e, p_factor = %e, PdV_work_fraction=%e.",p_terminal, p_ej, p_factor, PdV_work_fraction); */
  /* } else { */
  /*   message("We do resolve the Sedov-Taylor. p_t = %e, p_ej = %e, p_factor=%e, p_t/p_ej = %e", p_terminal, p_ej, p_factor, p_terminal/p_ej); */
  /* } */

  dp_prime[0] *= p_factor;
  dp_prime[1] *= p_factor;
  dp_prime[2] *= p_factor;

  /* Compute the cooling radius */
  const double second_part = p_terminal*p_terminal/(p_ej*p_ej) - 1;
  const double r_cool = pow(3.0*m_ej*second_part/(4.0*M_PI*pj->rho), 1.0/3.0);

  /* If we do not resolve the Taylor-Sedov, we rescale the internal energy */
  if (r2 > r_cool*r_cool) {
    dU *= pow(sqrt(r2)/r_cool, -6.5);
    message("We do not resolve the Sedov-Taylor (r_cool = %e). Rescaling dU.", r_cool);
  } /* else we do not change dU */

  /* Now we can give momentum, thermal and kinetic energy to the xpart. */
  for (int i = 0; i < 3; i++) {
    xpj->feedback_data.delta_p[i] += dp_prime[i];
  }
  xpj->feedback_data.delta_u += dU/new_mass;
  xpj->feedback_data.delta_E_kin += dKE;
  /* xpj->feedback_data.number_SN += 1; */

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(pj);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);

  /*-------------------------------------------------------------------------*/
  /* Verify conservation things */
  si->feedback_data.delta_m_check += dm;
  si->feedback_data.delta_E_check += dE;
  si->feedback_data.delta_p_norm_check += sqrt(dp_norm_2);

  si->feedback_data.delta_p_check[0] += dp[0];
  si->feedback_data.delta_p_check[1] += dp[1];
  si->feedback_data.delta_p_check[2] += dp[2];

  /* message("Conservation check (star %lld): Sum dm_i = %e (m_ej), Sum dE_i = %e (e_ej), Sum |dp_i| = %e (p_ej), Sum dp_i = (%e, %e, %e) (0), m_ej = %e, E_ej = %e, p_ej = %e", si->id, si->feedback_data.delta_m_check, si->feedback_data.delta_E_check, si->feedback_data.delta_p_norm_check, si->feedback_data.delta_p_check[0], si->feedback_data.delta_p_check[1], si->feedback_data.delta_p_check[2], m_ej, e_sn, p_ej); */
}
#else /* With EXTRA_STARS_LOOPS_2 */

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_feedback_apply(
    const float r2, const float dx[3], const float hi, const float hj,
    struct spart *si, struct part *pj, struct xpart *xpj,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct feedback_props *fb_props, const struct phys_const* phys_const,
    const struct unit_system* us, const integertime_t ti_current) {

  const double E_ej = si->feedback_data.energy_ejected;

  /* Do we have supernovae? */
  if (E_ej == 0) {
    return;
  }

  /* Finally, we can compute the w_j_bar. */
  double f_plus_i[3], f_minus_i[3], w_j[3];
  feedback_compute_vector_weight_non_normalized(r2, dx, hi, hj, si, pj, f_plus_i, f_minus_i, w_j);

  /* The normalized vector weight */
  double w_j_bar[3] =  {w_j[0]/si->feedback_data.enrichment_weight,
			w_j[1]/si->feedback_data.enrichment_weight,
			w_j[2]/si->feedback_data.enrichment_weight};
  const double w_j_bar_norm_2 = w_j_bar[0]*w_j_bar[0] + w_j_bar[1]*w_j_bar[1] + w_j_bar[2]*w_j_bar[2];
  const double w_j_bar_norm = sqrt(w_j_bar_norm_2);

  /* Here just get the feedback properties we want to distribute */
  const float mj = hydro_get_mass(pj);
  const double m_ej = si->feedback_data.mass_ejected;
  const double p_ej = sqrt(2*m_ej*E_ej) ;

  /* Distribute mass... (the max avoids to have dm=0 and NaN by dividing by dm) */
  double dm = max(w_j_bar_norm * m_ej, FLT_MIN);
  xpj->feedback_data.delta_mass += dm;

  /* ... metals */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    pj->chemistry_data.metal_mass[i] +=
	w_j_bar_norm * si->feedback_data.metal_mass_ejected[i];
  }

  /* ... momentum */
  const double dp[3] = {w_j_bar[0]*p_ej, w_j_bar[1]*p_ej, w_j_bar[2]*p_ej};
  const double dE = w_j_bar_norm * E_ej;

  /* Now boost to the 'laboratory' frame */
  double dp_prime[3] = {dp[0] + dm*si->v[0], dp[1] + dm*si->v[1], dp[2] + dm*si->v[2]};

  /* Note: The momentum is given to xpj later... We need to do some physics before. */

  /* ... Total energy */
  const double dp_norm_2 = dp[0]*dp[0] +  dp[1]*dp[1] +  dp[2]*dp[2];
  const double dp_prime_norm_2 = dp_prime[0]*dp_prime[0] +  dp_prime[1]*dp_prime[1] +  dp_prime[2]*dp_prime[2];
  const double dE_prime = dE + 1.0/(2.0*dm) * (dp_prime_norm_2 - dp_norm_2);

  /* ... internal energy */
  const double new_mass = mj + dm;

  /* Compute kinetic energy difference before and after SN */
  const double p_old_norm_2 = mj*mj*(xpj->v_full[0]*xpj->v_full[0] + xpj->v_full[1]*xpj->v_full[1] + xpj->v_full[2]*xpj->v_full[2]);
  const double p_new[3] = {mj*xpj->v_full[0] + dp_prime[0],
			   mj*xpj->v_full[1] + dp_prime[1],
			   mj*xpj->v_full[2] + dp_prime[2]};
  const double p_new_norm_2 = p_new[0]*p_new[0] + p_new[1]*p_new[1] + p_new[2]*p_new[2];

  const double E_kin_old = p_old_norm_2/(2.0*mj);
  const double E_kin_new = p_new_norm_2/(2.0*new_mass);
  const double dKE = E_kin_new - E_kin_old;

  const double U_old = xpj->u_full*mj;
  const double E_old = U_old + E_kin_old;
  const double E_new = E_old + dE_prime;
  const double U_new = E_new - E_kin_new;

  /* Compute the internal energy */
  double dU = U_new - U_old;
  /* message("dE_prime= %e, dU = %e,  dKE = %e, p_new_2 = %e, p_old_2 = %e", dE_prime, dU, dKE, p_new_norm_2, p_old_norm_2); */
  /* message("E_old = %e, E_new = %e, E_kin_old = %e, E_kin_new = %e, U_old = %e, U_new = %e", E_old, E_new, E_kin_old, E_kin_old, U_old, U_new); */

  /* --Now, we take into account for potentially unresolved energy-conserving
     phase of the SN explosion-- */

  /* Now, momentum given to the particle is not simply dp_prime. This is the
     momentum the location of the SN explosion. The ejecta must travel to the
     particle at distance r, i.e it must sweep up the mass mj of the particle
     as in a shock or a shell.
     This means that there is some work PdV that is done to reach the particle,
     work that converts thermal energy to kinetic energy. Therefore, the
     correct momentum to couple is proportional to the following factor: */
  const double PdV_work_fraction = sqrt(1 + mj/dm);

  /* At the end of the energy conserving phase, the blastwave reaches its
     terminal momentum. */
  const double p_terminal = feedback_get_SN_terminal_momentum(si, pj, xpj, phys_const, us);

  /* During the Taylor Sedov phase (energy conserving phase), dp_prime <=
     p_terminal. Thus the momentum coupling is irrelevant. We couple the
     thermal energy and leave the hydro solver correctly determines the PdV work. */

  /* If we can resolve the Taylor Sedov, then we give the right coupled
     momentum (which is by definition <= p_terminal). If we cannot resolve it,
     then we have reached the p_terminal. This is the upper limit of momentum,
     since afterwards the cooling is efficient and the thermal energy is radiated
     away.
     Thus, the factor to multiply dp_prime is: */
  const double p_factor = min(PdV_work_fraction, p_terminal/p_ej);

  /* if (p_factor == p_terminal/p_ej) { */
  /*   message("We do not resolve the Sedov-Taylor. Using p_terminal. p_t = %e, p_ej=%e, p_factor = %e, PdV_work_fraction=%e.",p_terminal, p_ej, p_factor, PdV_work_fraction); */
  /* } else { */
  /*   message("We do resolve the Sedov-Taylor. p_t = %e, p_ej = %e, p_factor=%e, p_t/p_ej = %e", p_terminal, p_ej, p_factor, p_terminal/p_ej); */
  /* } */

  dp_prime[0] *= p_factor;
  dp_prime[1] *= p_factor;
  dp_prime[2] *= p_factor;

  /* Compute the cooling radius */
  const double second_part = p_terminal*p_terminal/(p_ej*p_ej) - 1;
  const double r_cool = pow(3.0*m_ej*second_part/(4.0*M_PI*pj->rho), 1.0/3.0);

  /* If we do not resolve the Taylor-Sedov, we rescale the internal energy */
  if (r2 > r_cool*r_cool) {
    dU *= pow(sqrt(r2)/r_cool, -6.5);
    message("We do not resolve the Sedov-Taylor (r_cool = %e). Rescaling dU.", r_cool);
  } /* else we do not change dU */

  /* Now we can give momentum, thermal and kinetic energy to the xpart. */
  for (int i = 0; i < 3; i++) {
    xpj->feedback_data.delta_p[i] += dp_prime[i];
  }
  xpj->feedback_data.delta_u += dU/new_mass;
  xpj->feedback_data.delta_E_kin += dKE;
  /* xpj->feedback_data.number_SN += 1; */

  /* Impose maximal viscosity */
  hydro_diffusive_feedback_reset(pj);

  /* Synchronize the particle on the timeline */
  timestep_sync_part(pj);

  /*-------------------------------------------------------------------------*/
  /* Verify conservation things */
  si->feedback_data.delta_m_check += dm;
  si->feedback_data.delta_E_check += dE;
  si->feedback_data.delta_p_norm_check += sqrt(dp_norm_2);

  si->feedback_data.delta_p_check[0] += dp[0];
  si->feedback_data.delta_p_check[1] += dp[1];
  si->feedback_data.delta_p_check[2] += dp[2];

  /* message("Conservation check (star %lld): Sum dm_i = %e (m_ej), Sum dE_i = %e (e_ej), Sum |dp_i| = %e (p_ej), Sum dp_i = (%e, %e, %e) (0), m_ej = %e, E_ej = %e, p_ej = %e", si->id, si->feedback_data.delta_m_check, si->feedback_data.delta_E_check, si->feedback_data.delta_p_norm_check, si->feedback_data.delta_p_check[0], si->feedback_data.delta_p_check[1], si->feedback_data.delta_p_check[2], m_ej, e_sn, p_ej); */
}


#endif /* EXTRA_STARS_LOOPS_2*/

#endif /* SWIFT_GEAR_FEEDBACK_IACT_H */
