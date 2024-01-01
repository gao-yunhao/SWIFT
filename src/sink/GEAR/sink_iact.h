/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_SINKS_IACT_H
#define SWIFT_GEAR_SINKS_IACT_H

/* Local includes */
#include "error.h"
#include "gravity.h"
#include "gravity_iact.h"
#include "hydro.h"
#include "sink/GEAR/sink.h"
#include "sink/GEAR/sink_part.h"
#include "timeline.h"
#include <math.h>
#include <stddef.h>


/**
 * @brief do sink computation after the runner_iact_density (symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_sink(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float cut_off_radius) {

  /* In order to prevent the formation of two sink particles at a distance
   * smaller than the sink cutoff radius, we keep only gas particles with
   * the smallest potential. */

  const float r = sqrtf(r2);

  /* if the distance is less than the cut off radius */
  if (r < cut_off_radius) {

    float potential_i = pi->gpart->potential;
    float potential_j = pj->gpart->potential;

    /* prevent the particle with the largest potential to form a sink */
    if (potential_i > potential_j) {
      pi->sink_data.can_form_sink = 0;
      return;
    }

    if (potential_j > potential_i) {
      pj->sink_data.can_form_sink = 0;
      return;
    }
  }
}

/**
 * @brief do sink computation after the runner_iact_density (non symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_sink(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H, const float cut_off_radius) {

  /* In order to prevent the formation of two sink particles at a distance
   * smaller than the sink cutoff radius, we keep only gas particles with
   * the smallest potential. */

  const float r = sqrtf(r2);

  if (r < cut_off_radius) {

    float potential_i = pi->gpart->potential;
    float potential_j = pj->gpart->potential;

    /* if the potential is larger
     * prevent the particle to form a sink */
    if (potential_i > potential_j) pi->sink_data.can_form_sink = 0;
  }
}

/**
 * @brief Compute sink-sink swallow interaction (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param rj Comoving cut off radius of particle j.
 * @param si First sink particle.
 * @param sj Second sink particle.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_sink_swallow(const float r2, const float dx[3],
                                      const float ri, const float rj,
                                      struct sink *restrict si,
                                      struct sink *restrict sj,
				      const int with_cosmology,
				      const struct cosmology *cosmo,
				      const struct gravity_props *grav_props) {
  /* Binding energy check */
  /* Compute the physical relative velocity between the particles */
  const float dv[3] = {(sj->v[0] - si->v[0]) * cosmo->a_inv,
		       (sj->v[1] - si->v[1]) * cosmo->a_inv,
		       (sj->v[2] - si->v[2]) * cosmo->a_inv};

  /* Kinetic energy of the gas */
  float E_kin_rel = 0.5f * (dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]) ;

  /* Compute the Newtonian or truncated potential the sink exherts onto the
     gas particle */
  const float eps = gravity_get_softening(si->gpart, grav_props);
  const float eps2 = eps * eps;
  const float eps_inv = 1.f / eps;
  const float eps_inv3 = eps_inv * eps_inv * eps_inv;
  const float si_mass = si->mass;
  const float sj_mass = sj->mass;

  float dummy, pot_ij, pot_ji;
  runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, si_mass, &dummy, &pot_ij);
  runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, sj_mass, &dummy, &pot_ji);

  /* Compute the potential energies */
  float E_pot_ij = grav_props->G_Newton * pot_ij * cosmo->a_inv;
  float E_pot_ji = grav_props->G_Newton * pot_ji * cosmo->a_inv;

  /* Mechanical energy of the pair i-j and j-i */
  float E_mec_si = E_kin_rel + E_pot_ij ;
  float E_mec_sj = E_kin_rel + E_pot_ji ;

  /* Now, check if one is bound to the other */
  if ((E_mec_si > 0) && (E_mec_sj > 0)) {
    return;
  }

  /* The sink with the smaller mass will be merged onto the one with the
   * larger mass.
   * To avoid rounding issues, we additionally check for IDs if the sink
   * have the exact same mass. */
  if ((sj->mass < si->mass) || (sj->mass == si->mass && sj->id < si->id)) {
      /* This particle is swallowed by the sink with the largest mass of all the
       * candidates wanting to swallow it (we use IDs to break ties)*/
      if ((sj->merger_data.swallow_mass < si->mass) ||
	  (sj->merger_data.swallow_mass == si->mass &&
	   sj->merger_data.swallow_id < si->id)) {
	// message("sink %lld wants to swallow sink particle %lld", si->id,
	// sj->id);
	sj->merger_data.swallow_id = si->id;
	sj->merger_data.swallow_mass = si->mass;
      }
    }

#ifdef DEBUG_INTERACTIONS_SINKS
  /* Update ngb counters */
  if (si->num_ngb_formation < MAX_NUM_OF_NEIGHBOURS_SINKS)
    si->ids_ngbs_formation[si->num_ngb_formation] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_formation;
#endif
}

/**
 * @brief Compute sink-gas swallow interaction (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param ri Comoving cut off radius of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param si First sink particle.
 * @param pj Second particle.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_sinks_gas_swallow(const float r2, const float dx[3],
                                     const float ri, const float hj,
                                     struct sink *restrict si,
                                     struct part *restrict pj,
				     const int with_cosmology,
				     const struct cosmology *cosmo,
				     const struct gravity_props *grav_props,
				     const struct sink_props* sink_properties) {

 
  // message("sink %lld wants to swallow gas particle %lld", si->id, pj->id);

  const float r = sqrtf(r2);
  const float f_acc_r_acc = sink_properties->f_acc * ri ;

  if (sink_properties->do_regulated_accretion) {
    /* Update the number of gas neighbours of this sink */
    si->N_neighbours += 1;

    /* Update the mass in the interaction zone */
    si->mass_interaction_current += hydro_get_mass(pj);

    /* Now, allocate the memory to add one element to the neighbour array and
       add this element to the array. */
    sink_add_part_neighbour_array(si->neighbour_array, pj);

  } else {

    /* If the gas falls within f_acc*r_acc, it is accreted without further check */
    if (r < f_acc_r_acc) {
      /* Check if a gas particle has not been already marked to be swallowed by
	 another sink particle. */
      if (pj->sink_data.swallow_id < si->id) {
	pj->sink_data.swallow_id = si->id;
      }
    } else if ((r >= f_acc_r_acc) && (r < ri)) /* f_acc*r_acc <= r <= r_acc, we perform other checks */ {

      /* Compute the physical relative velocity between the particles */
      const float dv[3] = {(pj->v[0] - si->v[0]) * cosmo->a_inv,
			   (pj->v[1] - si->v[1]) * cosmo->a_inv,
			   (pj->v[2] - si->v[2]) * cosmo->a_inv};

      /* Compute the physical distance between the particles */
      const float dx_physical[3] = {dx[0] * cosmo->a,
				    dx[1] * cosmo->a,
				    dx[2] * cosmo->a};
      const float r_physical = r * cosmo->a;


      /* Momentum check */
      /* Relative momentum of the gas */
      const float specific_angular_momentum_gas[3] = {dx_physical[1] * dv[2] - dx_physical[2] * dv[1],
						      dx_physical[2] * dv[0] - dx_physical[0] * dv[2],
						      dx_physical[0] * dv[1] - dx_physical[1] * dv[0]};
      const float L2_gas = specific_angular_momentum_gas[0]*specific_angular_momentum_gas[0] + specific_angular_momentum_gas[1]*specific_angular_momentum_gas[1] + specific_angular_momentum_gas[2]*specific_angular_momentum_gas[2];

      /* Keplerian angular speed squared */
      const float omega_acc_2 = grav_props->G_Newton*si->mass / (r_physical*r_physical*r_physical);

      /*Keplerian angular momentum squared */
      const float L2_acc = (si->r_cut*si->r_cut*si->r_cut*si->r_cut)*omega_acc_2;

      /* To be accreted, the gas momentum shoulb lower than the keplerian orbit momentum. */
      if (L2_gas > L2_acc) {
	return ;
      }

      /* Energy check */
      /* Kinetic energy of the gas */
      float E_kin_relative_gas = 0.5f * (dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]) ;

      /* Compute the Newtonian or truncated potential the sink exherts onto the
	 gas particle */
      const float eps = gravity_get_softening(si->gpart, grav_props);
      const float eps2 = eps * eps;
      const float eps_inv = 1.f / eps;
      const float eps_inv3 = eps_inv * eps_inv * eps_inv;
      const float sink_mass = si->mass;
      float dummy, pot_ij;
      runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, sink_mass, &dummy,
			       &pot_ij);

      /* Compute the potential energy that the sink exerts in the gas (do not
	 forget to convert to physical quantity)*/
      float E_pot_gas = grav_props->G_Newton * pot_ij * cosmo->a_inv;

      /* Mechanical energy of the pair sink-gas */
      float E_mec_sink_part = E_kin_relative_gas + E_pot_gas ;

      /* To be accreted, the gas must be gravitationally bound to the sink. */
      if (E_mec_sink_part >= 0) return;

      /* Most bound pair check */
      /* The pair gas-sink must be the most bound among all sinks */
      if (E_mec_sink_part >= pj->sink_data.E_mec_bound) {
	return ;
      }

      /* Since this pair gas-sink is the most bounf, keep track of the
	 E_mec_bound and set the swallow_id accordingly */
      pj->sink_data.E_mec_bound = E_mec_sink_part;
      pj->sink_data.swallow_id = si->id;

    }
  } /* End of if-else of do_regulated_accretion */

#ifdef DEBUG_INTERACTIONS_SINKS
  /* Update ngb counters */
  if (si->num_ngb_formation < MAX_NUM_OF_NEIGHBOURS_SINKS)
    si->ids_ngbs_formation[si->num_ngb_formation] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_formation;
#endif
}

/**
 * @brief Swallow the gas
 *
 * @param
 */
INLINE static void
runner_iact_nonsym_sinks_do_gas_swallow_regulated(struct engine* e, struct space *s,
						  struct cell* c) {

  const struct cosmology* cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const struct phys_const* phys_const = e->physical_constants;
  const struct sink_props* sink_props = e->sink_properties;

  /* Get a pointer to the sinks in the cell */
  struct sink *sinks = c->sinks.parts;
  const size_t nr_sink = c->sinks.count;

  /* Loop over all sinks in the cell */
  for (size_t i = 0; i < nr_sink ; ++i) {
    /* Get a handle on the sink. */
    struct sink *si = &sinks[i];
    sink_neighbour_array* neighbour_array = si->neighbour_array;

    /* Skip inactive particles */
    if (!sink_is_active(si, e)) continue;

    /* If sp does not have any neighbour, continue with the next sink */
    if (si->N_neighbours == 0) continue;

    /* Check if the neighbour array is sorted. If it is already sorted, the
       job has already been done. */
    if (neighbour_array->is_sorted) {
      /* The job has already been done, we can continue to the next sink. */
      message("Array already sorted");
      continue;
    }

    /* This should only happen if we read a sink from ICs. */
    if (si->mass_interaction_init == 0.0) {
      si->mass_interaction_init = si->mass_interaction_current;
    }

    /* Sort the neigbour array from the closest neigbour to the farthest from
       the sink. */
    sink_sort_neighbour_array(si);

    /* Compute the accretion timescale */
    const double t_rad = sink_compute_radial_accretion_timescale(si, cosmo, sink_props);
    const double t_disk = sink_compute_disk_accretion_timescale(si, cosmo, sink_props, phys_const);
    const double f = sink_compute_f_accretion_timescale(si, cosmo, sink_props, phys_const);
    double t_acc = pow(t_rad, 1-f) * pow(t_disk, f);
    message("t_rad = %e", t_rad);
    message("t_disk = %e", t_disk);
    message("f = %lf", f);
    message("t_acc = %e", t_acc);


    /* If the mass in the interaction zone exceeds the mass in the interaction
    zone at the creation of the sink, decrease t_acc to accrete the excess of
    mass more quickly. */
    if (si->mass_interaction_current > si->mass_interaction_init) {
      t_acc /= ( (si->mass_interaction_current/si->mass_interaction_init) *  (si->mass_interaction_current/si->mass_interaction_init));
    }

    /* Get particle time-step */
    double dt_sink = 0.0;
    const integertime_t ti_step_sink = get_integer_timestep(si->time_bin);
    const integertime_t ti_begin_sink = get_integer_time_begin(e->ti_current, si->time_bin);
    if (with_cosmology) {
      dt_sink = cosmology_get_delta_time(cosmo, ti_begin_sink,
					 ti_begin_sink + ti_step_sink);
    } else {
      dt_sink = get_timestep(si->time_bin, e->time_base);
    }

    /* Compute the accreted mass */
    double delta_M_remaining = si->mass_interaction_current * (1 - exp(- dt_sink/t_acc));

    /* Mass removed from a particle. Used to update the sink*/
    double delta_m_j = 0.0 ; 
    message("delta_M = %e", delta_M_remaining);
    message("N_neighbour = %i", si->N_neighbours);

    /* Compute the delta t for the special "Toomre criterion" */
    /* Give an explanation...*/
    const double r_cut_3 = si->r_cut*si->r_cut*si->r_cut*cosmo->a*cosmo->a*cosmo->a;
    const double dt_criterion = sink_props->tol_param*sqrt(r_cut_3/(phys_const->const_newton_G*si->mass_interaction_current));
    /* message("dt_criterion = %lf", dt_criterion); */

    /* Loop over the neighbour part */
    for (size_t j = 0; j < neighbour_array->size ; ++j) {
      /* Get a handle on the part. */
      struct part *pj = neighbour_array->part_neighbours[j];

      /* If there is no remaining mass, stop */
      if (delta_M_remaining <= 0.0) break;

      /* Ignore inhibited particles (they have already been removed!) */
      if (part_is_inhibited(pj, e)) {
	continue;
      }

      /* If this part has already been entirely swallowed, skip it */
      if (pj->sink_data.swallow_id >= 0){
	continue;
      }

      /* Get the pj timestep */
      double dt_pj = 0.0;
      const integertime_t ti_step_pj = get_integer_timestep(pj->time_bin);
      const integertime_t ti_begin_pj = get_integer_time_begin(e->ti_current, pj->time_bin);
      if (with_cosmology) {
	dt_pj = cosmology_get_delta_time(cosmo, ti_begin_pj,
					   ti_begin_pj + ti_step_pj);
      } else {
	dt_pj = get_timestep(pj->time_bin, e->time_base);
      }

      /* Lock the space as we are going to work directly on the neigbour list */
      lock_lock(&s->lock);

      /* Check if the timestep of the part is small enough to accrete it
	 entirely */
      if (dt_pj < dt_criterion) {
	/* Do not put the mass of the part to zero, it will be swallowed
	   properly in runner_do_sinks_gas_swallow(). If the mass is set to
	   zero, the hydro runs into problems. */

	/* Mark this part to be swallowed entirely */
	pj->sink_data.swallow_id = si->id;

	/* Update mass removed from the part to update the sink properties */
	delta_m_j = pj->mass ;

	/* Update the remaining mass to continue swallowing. */
	delta_M_remaining -= pj->mass;
	/* message("Criterion passed ! dt_part = %e", dt_pj); */
      } else {

	double mass_part = hydro_get_mass(pj);
	if (mass_part > delta_M_remaining) {
	  /* If the particle is *more massive* than delta_M_remaining, remove
	     delta_M_remaining from the part, The accretion is finished */
	  pj->mass -= delta_M_remaining;
	  pj->gpart->mass -= delta_M_remaining;

	  /* Nothing remains to be eaten. */
	  delta_M_remaining = 0.0; 

 	  /* Do not mark the part to be swallowed */

	} else { /* If the particle is *less massive* than delta_M, mark it to be swallowed. */
	  /* Do not put the mass of the part to zero, it will be swallowed
	   properly in runner_do_sinks_gas_swallow(). If the mass is set to
	   zero, the hydro runs into problems. */

	  /* Mark this part to be swallowed */
	  pj->sink_data.swallow_id = si->id;

	  /* Update mass removed from the part to update the sink properties */
	  delta_m_j = pj->mass ;

	  /* Update the remaining mass to continue swallowing. */
	  delta_M_remaining -= mass_part;
	} /* Enf of gas accretion treatment */
      } /* End of special criterion */

      /* Release the space as we are done updating the part */
      if (lock_unlock(&s->lock) != 0) error("Failed to unlock the space.");

      /* Update the sink properties */
      /* If the particle has been entirely swallowed, its sink_data.swallow_id
	 contains the id of sink that swallows it. In such case, here we do not
	 update the sink nor remove it; both things are done properly in
	 runner_do_sinks_gas_swallow(). */
      /* Finally, I changed my mind. We'll change depending on Yves' or
	 Matthieu's recommandation. I will update the sink here for all
	 cases and remove the part here, by copying the code and adapting
	 it. Then, to make sure the cell gets an updated value of
	 ti_beg_max, I won't end the runner_do_sinks_gas_swallow() after
	 the current fct. I will just mark the part as swalloed with
	 the fct sink_mark_part_as_swallowed(&p->sink_data). When leaving this
	 function, the part id will be -2, preventing it to be reswallowed
	 afterwards. But it will allow the cell to update the ti_beg_max. */
      /* Pay attention to update the gpart accordingly ! */

      /* Lock space again, we are updating the sink and the part. */
      lock_lock(&s->lock);

      sink_swallow_part_regulated_accretion(si, pj, NULL, cosmo, delta_m_j);

      /* Release the space as we are done updating the part. */
      if (lock_unlock(&s->lock) != 0)
	error("Failed to unlock the space.");

      /* Remove the particle flagged for swallowing */
      /* Get the ID of the sink that will swallow this part */
      const long long swallow_id = sink_get_part_swallow_id(&pj->sink_data);

      /* Has this particle been flagged for swallowing? */
      if (swallow_id >= 0) {

	/* If the gas particle is local, remove it */
	if (c->nodeID == e->nodeID) {

	  lock_lock(&e->s->lock);

	  /* Re-check that the particle has not been removed
	   * by another thread before we do the deed. */
	  if (!part_is_inhibited(pj, e)) {

	    /* Finally, remove the gas particle from the system
	     * Recall that the gpart associated with it is also removed
	     * at the same time. */
	    cell_remove_part(e, c, pj, NULL); /* For now xpart = NULL */
	  }

	  if (lock_unlock(&e->s->lock) != 0)
	    error("Failed to unlock the space!");
	}

	/* In any case, prevent the particle from being re-swallowed */
	sink_mark_part_as_swallowed(&pj->sink_data);
      } /* End of gas removal */

    }  /* End of neighbour loop */

    
    /* Quick angular momentum feedback */
    /* Do no do it in the neighbour loop because 1) it should have its own task
       and 2) Hubber et al 2013 does it a the end. Why ? Think about it. */
      /* Compute the amount of angular momentum transferred */

  } /* End of sink loop */

  return ;
}

#endif
