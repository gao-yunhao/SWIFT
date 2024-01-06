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
#ifndef SWIFT_GEAR_SINK_H
#define SWIFT_GEAR_SINK_H

#include <float.h>
#include <math.h>
#include <stdint.h>

/* Local includes */
#include "active.h"
#include "chemistry.h"
#include "cooling.h"
#include "error.h"
#include "feedback.h"
#include "gravity.h"
#include "gravity_iact.h"
#include "hydro.h"
#include "minmax.h"
#include "physical_constants.h"
#include "random.h"
#include "sink_part.h"
#include "sink_properties.h"

/**
 * @brief Computes the time-step of a given sink particle.
 *
 * @param sp Pointer to the sink-particle data.
 */
__attribute__((always_inline)) INLINE static float sink_compute_timestep(
    const struct sink* const sp) {

  return FLT_MAX;
}

/**
 * @brief Update the target mass of the sink particle.
 *
 * @param e The #engine
 * @param sink the sink particle.
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 */
INLINE static void sink_update_target_mass(struct sink* sink,
                                           const struct sink_props* sink_props,
                                           const struct engine* e, int rloop) {

  float random_number = random_unit_interval_part_ID_and_loop_idx(
      sink->id, rloop, e->ti_current, random_number_sink_formation);

  struct feedback_props* feedback_props = e->feedback_props;

  /* Pick the correct table. (if only one table, threshold is < 0) */

  const float metal =
      chemistry_get_sink_total_iron_mass_fraction_for_feedback(sink);
  const float threshold = feedback_props->metallicity_max_first_stars;

  struct stellar_model* model;
  double minimal_discrete_mass;

  if (metal >= threshold) {
    model = &feedback_props->stellar_model;
    minimal_discrete_mass = sink_props->minimal_discrete_mass_first_stars;
  } else {
    model = &feedback_props->stellar_model_first_stars;
    minimal_discrete_mass = sink_props->minimal_discrete_mass;
  }

  const struct initial_mass_function* imf = &model->imf;

  if (random_number < imf->sink_Pc) {
    // we are dealing with the continous part of the IMF
    sink->target_mass = imf->sink_stellar_particle_mass;
    sink->target_type = 0;
  } else {
    // we are dealing with the discrete part of the IMF
    random_number = random_unit_interval_part_ID_and_loop_idx(
        sink->id, rloop + 1, e->ti_current, random_number_sink_formation);
    double m =
        random_sample_power_law(minimal_discrete_mass, imf->mass_max,
                                imf->exp[imf->n_parts - 1], random_number);
    sink->target_mass = m;
    sink->target_type = 1;
  }
}


/*****************************************************************************/
/* Block size for the sink_neighbour_array */
#define NEIGHBOUR_ARRAY_PADDING 64

/**
 * @brief Create the sink_neighbour array.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * @param array The #sink_neighbour_array.
 */
INLINE static sink_neighbour_array* sink_construct_neighbour_array(sink_neighbour_array* array) {
  if (array != NULL) {
    /* Starts with size=0, allocated=0, is_sorted=0 and part_neighbour=NULL */
    sink_neighbour_array result = { 0, 0, 0, NULL };
    result.part_neighbours = calloc(NEIGHBOUR_ARRAY_PADDING, sizeof(struct part*));
    if (result.part_neighbours != NULL) {
      result.allocated = NEIGHBOUR_ARRAY_PADDING;
    } else {
      /* If the memory could not be allocated, return NULL */
      return NULL;
    }
    *array = result;
  }
  return array;
}

/**
 * @brief Delete the sink_neighbour array.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * @param array The #sink_neighbour_array.
 */
INLINE static void sink_delete_neighbour_array(sink_neighbour_array* array) {
  if ((array != NULL) && (array->part_neighbours != NULL)) {
    free(array->part_neighbours);
    array->size = 0;
    array->allocated = 0;
    array->is_sorted = 0;
  }
  return;
}

/**
 * @brief Increases the size of the the sink_neighbour array.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * @param array The #sink_neighbour_array.
 */
INLINE static sink_neighbour_array* sink_enlarge_neighbour_array(sink_neighbour_array* array) {
  if (array != NULL) {
    sink_neighbour_array result = *array;
    result.allocated += NEIGHBOUR_ARRAY_PADDING;
    if ((result.allocated > SIZE_MAX/sizeof(struct part*)) || ((result.part_neighbours = realloc(result.part_neighbours, result.allocated * sizeof(struct part*))) == NULL)) {
      /* If we could not reallocate, return NULL ; array has not been modified. */
      return NULL;
    }
    /* Final affectation */
    *array = result;
  }
  return array;
}

/**
 * @brief Add a  particle to the sink_neighbour array.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * @param array The #sink_neighbour_array.
 * @param p The  #part to add.
 */
INLINE static size_t sink_add_part_neighbour_array(sink_neighbour_array* array,
						   struct part* p) {
  if (array != NULL) {
    /* Check that the array has enough room for the element and allocate if
       necessary */
    while (array->size >= array->allocated) {
      if (sink_enlarge_neighbour_array(array) == NULL) {
	return 0;
      }
    } /* End of while */

    /* Now add the particle to the array */
    array->part_neighbours[array->size] = p;
    ++(array->size);
    return array->size;
  }

  /* If the array was NULL, its size is 0 */
  return 0;
}
/*****************************************************************************/


/**
 * @brief Initialises the sink-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon
 * @param sink_props The properties of the sink particles scheme.
 */
__attribute__((always_inline)) INLINE static void sink_first_init_sink(
								       struct sink* sp, const struct sink_props* sink_props, const struct engine* e) {

  if (sink_props->do_regulated_accretion) {
    sp->N_neighbours = 0;
    sp->mass_interaction_init = 0.f;
    sp->mass_interaction_current = 0.f;
    sp->neighbour_array = malloc(sizeof(sink_neighbour_array));
    message("Creating neighbour array after ICs...");

    if (sink_construct_neighbour_array(sp->neighbour_array) == NULL) {
      error("Could not construct sink neighbour_neighnour array. Aborting.");
    }

    /* The r_cut is defined by the SmoothingLength field in the ICs. This
    field is optional for sinks. If it has not been given, the value of
    sp->r_cut is set to 0.0. So, check that r_cut is not 0. */
    if (sp->r_cut == 0.) {
      sp->r_cut = sink_props->f_interaction * sink_props->cut_off_radius;
    } else {
      sp->r_cut *= sink_props->f_interaction;
    }
  } else {
    sp->r_cut = sink_props->cut_off_radius;
  }
  sp->time_bin = 0;

  sp->number_of_gas_swallows = 0;
  sp->number_of_direct_gas_swallows = 0;
  sp->number_of_sink_swallows = 0;
  sp->number_of_direct_sink_swallows = 0;
  sp->swallowed_angular_momentum[0] = 0.f;
  sp->swallowed_angular_momentum[1] = 0.f;
  sp->swallowed_angular_momentum[2] = 0.f;

  sink_mark_sink_as_not_swallowed(&sp->merger_data);

  /* Bug fix: Setup the target mass for sink formation after reading the
     ICs. Otherwise sink->target_mass = 0.0 and a sink present in the IC spawn
     a star of mass 0.0... */
  sink_update_target_mass(sp, sink_props, e, 0);
}

/**
 * @brief Initialisation of particle data before the hydro density loop.
 * Note: during initalisation (space_init)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sink_init_part(struct part* restrict p,
								 const struct sink_props* sink_props) {

  struct sink_part_data* cpd = &p->sink_data;

  if (sink_props->disable_sink_formation) {
    cpd->can_form_sink = 0;
  } else {
    cpd->can_form_sink = 1;
  }
  cpd->E_kin_neighbours = 0.f;
  cpd->E_int_neighbours = 0.f;
  cpd->E_rad_neighbours = 0.f;
  cpd->E_pot_self_neighbours = 0.f;
  cpd->E_pot_ext_neighbours = 0.f;
  cpd->E_mag_neighbours = 0.f;
  cpd->E_rot_neighbours = 0.f;
  cpd->potential = 0.f;
  cpd->E_mec_bound = 0.f; /* Gravitationally bound particles will have
			     E_mec_bound < 0. This is checked before comparing
			     any other value with this one. So no need to put
			     it to the max of float. */
  cpd->is_overlapping_sink = 0;
  cpd->mass_interaction_init = 0.f;
}

/**
 * @brief Initialisation of sink particle data before sink loops.
 * Note: during initalisation (space_init_sinks)
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sink_init_sink(
    struct sink* sp) {
#ifdef DEBUG_INTERACTIONS_SINKS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_accretion[i] = -1;
  sp->num_ngb_accretion = 0;

  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_merger[i] = -1;
  sp->num_ngb_merger = 0;

  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SINKS; ++i)
    sp->ids_ngbs_formation[i] = -1;
  sp->num_ngb_formation = 0;
#endif

  /* Reset the number of neighbour */
  sp->N_neighbours = 0;

  /* Free the old data about the neighbour array */
  sink_delete_neighbour_array(sp->neighbour_array);

  /* Allocates storage */
  sp->neighbour_array = malloc(sizeof(sink_neighbour_array));

  /* Create a new neigbour array */
  if (sink_construct_neighbour_array(sp->neighbour_array) == NULL) {
    error("Could not construct sink neighbour_neighnour array. Aborting.");
  }

  sp->mass_interaction_current = 0.f;

}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param sp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void sink_predict_extra(
    struct sink* restrict sp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The particle.
 */
__attribute__((always_inline)) INLINE static void sink_reset_predicted_values(
    struct sink* restrict sp) {}

/**
 * @brief Kick the additional variables
 *
 * @param sp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void sink_kick_extra(
    struct sink* sp, float dt) {}

/**
 * @brief Calculate if the gas has the potential of becoming
 * a sink.
 *
 * Return 0 if no sink formation should occur.
 * Note: called in runner_do_sink_formation
 *
 * @param sink_props the sink properties to use.
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cooling The cooling data struct.
 *
 */
INLINE static int sink_is_forming(
    const struct part* restrict p, const struct xpart* restrict xp,
    const struct sink_props* sink_props, const struct phys_const* phys_const,
    const struct cosmology* cosmo,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling,
    const struct entropy_floor_properties* restrict entropy_floor) {

  /* the particle is not elligible */
  if (!p->sink_data.can_form_sink) return 0;

  const struct sink_part_data* sink_data = &p->sink_data;

  const float temperature_max = sink_props->maximal_temperature;
  const float temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                    cosmo, cooling, p, xp);

  const float density_threshold = sink_props->density_threshold;
  const float density = hydro_get_physical_density(p, cosmo);

  const float div_v = hydro_get_div_v(p);

  const float h = p->h;
  const float sink_cut_off_radius = sink_props->cut_off_radius;

  double E_grav = sink_data->E_pot_self_neighbours;
  double E_tot = sink_data->E_kin_neighbours + sink_data->E_int_neighbours +
                 E_grav + sink_data->E_mag_neighbours;

  /* Density and temperature check */
  if (density <= density_threshold || temperature >= temperature_max) {
    return 0;
  }

  /* Contracting gas check */
  if ((sink_props->sink_formation_contracting_gas_check) && (div_v > 0)) {
    return 0;
  }

  /* Smoothing length check */
  if ((sink_props->sink_formation_smoothing_length_check) && (h >= 0.5 * sink_cut_off_radius)) {
    return 0;
  }

  /* Active neighbours check */
  /* This is checked on the fly in runner_do_sink_formation(). The part is
     flagged to not form sink through p->sink_data.can_form_sink */

  /* Jeans instability check */
  if ((sink_props->sink_formation_jeans_instability_check) && (sink_data->E_int_neighbours >= 0.5f * fabs(E_grav))) {
    return 0;
  }

  if ((sink_props->sink_formation_jeans_instability_check) && (sink_data->E_int_neighbours + sink_data->E_rot_neighbours >= fabs(E_grav))) {
    return 0;
  }

  /* Bound state check */
  if ((sink_props->sink_formation_bound_state_check) && (E_tot >= 0)) {
    return 0;
  }

  /* Minimum of the potential check */
  /* Done in density loop. The gas is then flagged through
     sink_data.can_form_sink to not form sink. The check is done at the
     beginning. */

  /* Overlapping existing sinks check */
  if (sink_data->is_overlapping_sink) {
    return 0;
  }

#ifdef SWIFT_DEBUG_CHECKS
  message("Gas particle %lld can form a sink !", p->id);
#endif
  return 1;
}

/**
 * @brief Decides whether a particle should be converted into a
 * sink or not.
 *
 * No SF should occur, so return 0.
 * Note: called in runner_do_sink_formation
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param sink_props The properties of the sink model.
 * @param e The #engine (for random numbers).
 * @param dt_sink The time-step of this particle
 * @return 1 if a conversion should be done, 0 otherwise.
 */
INLINE static int sink_should_convert_to_sink(
    const struct part* p, const struct xpart* xp,
    const struct sink_props* sink_props, const struct engine* e,
    const double dt_sink) {

  /* We do not use a stockastic approach.
   * Once elligible (sink_is_forming), the gas particle form a sink */

  return 1;
}

/**
 * @brief Copies the properties of the gas particle over to the
 * sink particle.
 *
 * Nothing to do here.
 *
 * @param e The #engine
 * @param p the gas particles.
 * @param xp the additional properties of the gas particles.
 * @param sink the new created sink  particle with its properties.
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 */
INLINE static void sink_copy_properties(
    struct part* p, const struct xpart* xp, struct sink* sink,
    const struct engine* e, const struct sink_props* sink_props,
    const struct cosmology* cosmo, const int with_cosmology,
    const struct phys_const* phys_const,
    const struct hydro_props* restrict hydro_props,
    const struct unit_system* restrict us,
    const struct cooling_function_data* restrict cooling) {

  /* First initialisation */
  sink_init_sink(sink);

  /* Set a smoothing length */
  if (sink_props->do_regulated_accretion) {
    sink->r_cut = sink_props->f_interaction*p->h;

    /* Create a new neigbour array */
    sink->neighbour_array = malloc(sizeof(sink_neighbour_array));
    if (sink_construct_neighbour_array(sink->neighbour_array) == NULL) {
      error("Could not construct sink neighbour_neighnour array. Aborting.");
    }

  } else {
    sink->r_cut = e->sink_properties->cut_off_radius;
  }

  /* Flag it as not swallowed */
  sink_mark_sink_as_not_swallowed(&sink->merger_data);

  /* Additional initialisation */

  /* setup the target mass for sink star formation */
  sink_update_target_mass(sink, sink_props, e, 0);

  /* Copy the chemistry properties */
  chemistry_copy_sink_properties(p, xp, sink);

  /* test the mass distribution */
  // for (int i=0;i<1000000;i++)
  //   {
  //     sink_update_target_mass(sink, sink_props, e, i);
  //     message("%g %d",sink->target_mass,sink->target_type);
  //   }
  // exit(-1);
}

/**
 * @brief Update the properties of a sink particles by swallowing
 * a gas particle.
 *
 * @param sp The #sink to update.
 * @param p The #part that is swallowed.
 * @param xp The #xpart that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_swallow_part(
    struct sink* sp, const struct part* p, const struct xpart* xp,
    const struct cosmology* cosmo) {

  /* Get the current dynamical masses */
  const float gas_mass = hydro_get_mass(p);
  const float sink_mass = sp->mass;

  /* store the mass of the sink part i */
  const float msp_old = sp->mass;

  /* Increase the dynamical mass of the sink. */
  sp->mass += gas_mass;
  sp->gpart->mass += gas_mass;

  /* Physical velocity difference between the particles */
  const float dv[3] = {(sp->v[0] - p->v[0]) * cosmo->a_inv,
                       (sp->v[1] - p->v[1]) * cosmo->a_inv,
                       (sp->v[2] - p->v[2]) * cosmo->a_inv};

  /* Physical distance between the particles */
  const float dx[3] = {(sp->x[0] - p->x[0]) * cosmo->a,
                       (sp->x[1] - p->x[1]) * cosmo->a,
                       (sp->x[2] - p->x[2]) * cosmo->a};

  /* Collect the swallowed angular momentum */
  sp->swallowed_angular_momentum[0] +=
      gas_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
  sp->swallowed_angular_momentum[1] +=
      gas_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
  sp->swallowed_angular_momentum[2] +=
      gas_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

  /* Update the sink momentum */
  const float sink_mom[3] = {sink_mass * sp->v[0] + gas_mass * p->v[0],
                             sink_mass * sp->v[1] + gas_mass * p->v[1],
                             sink_mass * sp->v[2] + gas_mass * p->v[2]};

  sp->v[0] = sink_mom[0] / sp->mass;
  sp->v[1] = sink_mom[1] / sp->mass;
  sp->v[2] = sink_mom[2] / sp->mass;
  sp->gpart->v_full[0] = sp->v[0];
  sp->gpart->v_full[1] = sp->v[1];
  sp->gpart->v_full[2] = sp->v[2];

  /* Update the sink metal masses fraction */
  chemistry_add_part_to_sink(sp, p, msp_old);

  /* This sink swallowed a gas particle */
  sp->number_of_gas_swallows++;
  sp->number_of_direct_gas_swallows++;

#ifdef SWIFT_DEBUG_CHECKS
  const float dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  message("sink %lld swallow gas particle %lld. "
	  "(Mass = %e, "
	  "Delta_v = [%f, %f, %f] U_V, "
	  "Delta_x = [%f, %f, %f] U_L, "
	  "Delta_v_rad = %f)", sp->id, p->id, sp->mass, -dv[0], -dv[1], -dv[2], -dx[0], -dx[1], -dx[2],
      (dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]) / dr);
#endif
}

/**
 * @brief Update the properties of a sink particles by swallowing
 * a sink particle.
 *
 * @param spi The #sink to update.
 * @param spj The #sink that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sink_swallow_sink(
    struct sink* spi, const struct sink* spj, const struct cosmology* cosmo) {

  /* Get the current dynamical masses */
  const float spi_dyn_mass = spi->mass;
  const float spj_dyn_mass = spj->mass;

  /* store the mass of the sink part i */
  const float mi_old = spi->mass;

  /* Increase the masses of the sink. */
  spi->mass += spj->mass;
  spi->gpart->mass += spj->mass;

  /* Collect the swallowed angular momentum */
  spi->swallowed_angular_momentum[0] += spj->swallowed_angular_momentum[0];
  spi->swallowed_angular_momentum[1] += spj->swallowed_angular_momentum[1];
  spi->swallowed_angular_momentum[2] += spj->swallowed_angular_momentum[2];

  /* Update the sink momentum */
  const float sink_mom[3] = {
      spi_dyn_mass * spi->v[0] + spj_dyn_mass * spj->v[0],
      spi_dyn_mass * spi->v[1] + spj_dyn_mass * spj->v[1],
      spi_dyn_mass * spi->v[2] + spj_dyn_mass * spj->v[2]};

  spi->v[0] = sink_mom[0] / spi->mass;
  spi->v[1] = sink_mom[1] / spi->mass;
  spi->v[2] = sink_mom[2] / spi->mass;
  spi->gpart->v_full[0] = spi->v[0];
  spi->gpart->v_full[1] = spi->v[1];
  spi->gpart->v_full[2] = spi->v[2];

  /* Update the sink metal masses fraction */
  chemistry_add_sink_to_sink(spi, spj, mi_old);

  /* This sink swallowed a sink particle */
  spi->number_of_sink_swallows++;
  spi->number_of_direct_sink_swallows++;

  /* Add all other swallowed particles swallowed by the swallowed sink */
  spi->number_of_sink_swallows += spj->number_of_sink_swallows;
  spi->number_of_gas_swallows += spj->number_of_gas_swallows;

#ifdef SWIFT_DEBUG_CHECKS
  message("sink %lld swallow sink particle %lld. New mass: %e.", spi->id, spj->id, spi->mass);
#endif
  
}

/**
 * @brief Should the sink spawn a star particle?
 *
 * @param e The #engine
 * @param sink the sink particle.
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 * @param us The internal unit system.
 */
INLINE static int sink_spawn_star(struct sink* sink, const struct engine* e,
                                  const struct sink_props* sink_props,
                                  const struct cosmology* cosmo,
                                  const int with_cosmology,
                                  const struct phys_const* phys_const,
                                  const struct unit_system* restrict us) {

  if (sink->mass > sink->target_mass * phys_const->const_solar_mass)
    return 1;
  else
    return 0;
}

/**
 * @brief Separate the #spart and #part by randomly moving both of them.
 *
 * @param e The #engine.
 * @param p The #part generating a star.
 * @param xp The #xpart generating a star.
 * @param sp The new #spart.
 */
INLINE static void sink_star_formation_separate_particles(
    const struct engine* e, struct sink* si, struct spart* sp) {
#ifdef SWIFT_DEBUG_CHECKS
  if (si->x[0] != sp->x[0] || si->x[1] != sp->x[1] || si->x[2] != sp->x[2]) {
    error(
        "Moving particles that are not at the same location."
        " (%g, %g, %g) - (%g, %g, %g)",
        si->x[0], si->x[1], si->x[2], sp->x[0], sp->x[1], sp->x[2]);
  }
#endif

  /* Move a bit the particle in order to avoid
     division by 0.
  */
  const float max_displacement = 0.1;
  const double delta_x =
      2.f * random_unit_interval(si->id, e->ti_current,
                                 (enum random_number_type)0) -
      1.f;
  const double delta_y =
      2.f * random_unit_interval(si->id, e->ti_current,
                                 (enum random_number_type)1) -
      1.f;
  const double delta_z =
      2.f * random_unit_interval(si->id, e->ti_current,
                                 (enum random_number_type)2) -
      1.f;

  sp->x[0] += delta_x * max_displacement * si->r_cut;
  sp->x[1] += delta_y * max_displacement * si->r_cut;
  sp->x[2] += delta_z * max_displacement * si->r_cut;

  /* Copy the position to the gpart */
  sp->gpart->x[0] = sp->x[0];
  sp->gpart->x[1] = sp->x[1];
  sp->gpart->x[2] = sp->x[2];

  /* Do the sink particle. */
  const double mass_ratio = sp->mass / si->mass;
  const double dx[3] = {mass_ratio * delta_x * max_displacement * si->r_cut,
                        mass_ratio * delta_y * max_displacement * si->r_cut,
                        mass_ratio * delta_z * max_displacement * si->r_cut};

  si->x[0] -= dx[0];
  si->x[1] -= dx[1];
  si->x[2] -= dx[2];

  /* Copy the position to the gpart */
  si->gpart->x[0] = si->x[0];
  si->gpart->x[1] = si->x[1];
  si->gpart->x[2] = si->x[2];
}

/**
 * @brief Copy the properties of the sink particle towards the new star.
 * This function also needs to update the sink particle.
 *
 * @param e The #engine
 * @param sink the sink particle.
 * @param sp The star particle.
 * @param sink_props the sink properties to use.
 * @param phys_const the physical constants in internal units.
 * @param cosmo the cosmological parameters and properties.
 * @param with_cosmology if we run with cosmology.
 * @param us The internal unit system.
 */
INLINE static void sink_copy_properties_to_star(
    struct sink* sink, struct spart* sp, const struct engine* e,
    const struct sink_props* sink_props, const struct cosmology* cosmo,
    const int with_cosmology, const struct phys_const* phys_const,
    const struct unit_system* restrict us) {

  sink_star_formation_separate_particles(e, sink, sp);

  /* set the mass */
  sp->mass = sink->target_mass * phys_const->const_solar_mass;

  /* set feedback type */
  sp->feedback_data.type = sink->target_type;

  /* Initialize the feedback */
  if (sp->feedback_data.type == 1)
    feedback_init_after_star_formation(sp, e->feedback_props);

  /* sph smoothing */
  sp->h = sink->r_cut;

  /* mass at birth */
  sp->sf_data.birth_mass = sp->mass;

  /* Store either the birth_scale_factor or birth_time depending  */
  if (with_cosmology) {
    sp->birth_scale_factor = cosmo->a;
  } else {
    sp->birth_time = e->time;
  }

  /* Store the tracers data */
  // sp->tracers_data = sink->tracers_data;

  /* Move over the splitting data */
  // sp->split_data = sink->split_data;

  ///* Store the birth density in the star particle */
  // sp->sf_data.birth_density = hydro_get_physical_density(p, cosmo);

  ///* Store the birth temperature*/
  // sp->sf_data.birth_temperature = cooling_get_temperature(
  //     phys_const, hydro_props, us, cosmo, cooling, p, xp);

  /* Copy the chemistry properties */
  chemistry_copy_sink_properties_to_star(sink, sp);

  /* Copy the progenitor id */
  sp->sf_data.progenitor_id = sink->id;
}

/**
 * @brief Store the gravitational potential of a particle by copying it from
 * its #gpart friend.
 *
 * @param p_data The sink data of a gas particle.
 * @param gp The part's #gpart.
 */
__attribute__((always_inline)) INLINE static void sink_store_potential_in_part(
    struct sink_part_data* p_data, const struct gpart* gp) {
  p_data->potential = gp->potential;
}

/**
 * @brief Compute the energies (kinetic, potential, etc ) of the gas particle
 * and all quantities required for the formation of a sink.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param p The #part.
 * @param xp The #xpart data of the particle #p.
 */
INLINE static void sink_prepare_part_sink_formation(struct engine* e, struct cell* c, struct part* restrict p, struct xpart* restrict xp) {
  const struct cosmology *cosmo = e->cosmology;
  const int count = c->hydro.count;
  const struct gravity_props *grav_props = e->gravity_properties;
  const struct sink_props *sink_props = e->sink_properties;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  const int with_self_grav = (e->policy & engine_policy_self_gravity);
  const float r_acc_p = sink_props->cut_off_radius*cosmo->a; /* Physical accretion radius of part p */

  /* No external potential for now */
  /* const int with_ext_grav = (e->policy & engine_policy_external_gravity); */
  /* const struct external_potential *potential = e->external_potential; */
  
  /* Loop over all particles to find the neighbours within r_acc. Then,
     compute all quantities you need.  */
  const float px[3] = {(float)(p->x[0] - c->loc[0]),
      (float)(p->x[1] - c->loc[1]),
      (float)(p->x[2] - c->loc[2])};

  /* Compute the physical velocity */
  const float v[3] = {(p->v[0]) * cosmo->a_inv,
		      (p->v[1]) * cosmo->a_inv,
		      (p->v[2]) * cosmo->a_inv};

  float E_rot_x = 0;
  float E_rot_y = 0;
  float E_rot_z = 0;

  /* Loop over the gas particles to find its neighbours */
  for (int i = 0; i < count; i++) {

    /*Get a handle on the part */
    struct part *restrict pi = &parts[i];
    struct xpart *restrict xpi = &xparts[i];

    /* If for some reason the particle has been flagged to not form sink,
       do not continue and save some computationnal ressources. */
    if (!p->sink_data.can_form_sink){
      break ;
    }

    /* Compute the pairwise physical distance */
    const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
			  (float)(pi->x[1] - c->loc[1]),
			  (float)(pi->x[2] - c->loc[2])};

    const float dx[3] = {(px[0] - pix[0]) * cosmo->a,
			 (px[1] - pix[1]) * cosmo->a,
			 (px[2] - pix[2]) * cosmo->a};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    /* Checks that this part is a neighbour */
    if ((r2 > r_acc_p*r_acc_p) || (r2 == 0.0)) {
      continue;
    }

    /* Do not form sinks if some neighbours are not active */
    if (!part_is_active(pi, e)) {
      p->sink_data.can_form_sink = 0 ; 
      continue ;
    }

    const float mi = hydro_get_mass(p);
    const float u_inter_i =
      hydro_get_drifted_physical_internal_energy(p, cosmo);

    /* Compute the relative physical velocity between p and pi */
    const float vi[3] = {(pi->v[0]) * cosmo->a_inv,
			 (pi->v[1]) * cosmo->a_inv,
			 (pi->v[2]) * cosmo->a_inv};
    const float dv[3] = {vi[0] - v[0], vi[1] - v[1], vi[2] - v[2]};

    /* Compute specific angular momentum between pk and pi */
    const float specific_angular_momentum[3] = {
      dx[1] * dv[2] - dx[2] * dv[1], dx[2] * dv[0] - dx[0] * dv[2],
      dx[0] * dv[1] - dx[1] * dv[0]};

    /* Updates the energies */
    p->sink_data.E_kin_neighbours += 0.5f * mi * ((vi[0] * vi[0] - v[0] * v[0]) +
						  (vi[1] * vi[1] - v[1] * v[1]) +
						  (vi[1] * vi[2] - v[2] * v[2]));
    p->sink_data.E_int_neighbours += mi * u_inter_i;
    p->sink_data.E_rad_neighbours += cooling_get_radiated_energy(xpi);

    /* Notice that we skip the potential of the current particle here
       instead of subtracting it later */
    if ((with_self_grav) && (pi != p))
      p->sink_data.E_pot_self_neighbours += 0.5 * mi * pi->sink_data.potential * cosmo->a_inv;

    /* No external potential for now */
    /* if (gpi != NULL && with_ext_grav)	 */
    /* p->sink_data.E_pot_ext_neighbours +=  mi *
     * external_gravity_get_potential_energy( */
    /* time, potential, phys_const, gpi); */

    /* Need to include mhd header */
    /* p->sink_data.E_mag_neighbours += mhd_get_magnetic_energy(p, xpi); */

    /* Compute rotation energies */
    E_rot_x += 0.5 * mi * specific_angular_momentum[0] *
      specific_angular_momentum[0] /
      sqrtf(dx[1] * dx[1] + dx[2] * dx[2]);
    E_rot_y += 0.5 * mi * specific_angular_momentum[1] *
      specific_angular_momentum[1] /
      sqrtf(dx[0] * dx[0] + dx[2] * dx[2]);
    E_rot_z += 0.5 * mi * specific_angular_momentum[2] *
      specific_angular_momentum[2] /
      sqrtf(dx[0] * dx[0] + dx[1] * dx[1]);

    /* Update the mass in the interaction zone of p */
    if (sink_props->do_regulated_accretion) {
      p->sink_data.mass_interaction_init += pi->mass;
    }
  } /* End of gas neighbour loop */

  p->sink_data.E_rot_neighbours +=
    sqrtf(E_rot_x * E_rot_x + E_rot_y * E_rot_y + E_rot_z * E_rot_z);

  /* Shall we reset the values of the energies for the next timestep? No, it is
     done in cell_drift.c and space_init.c, for active particles. The
     potential is set in runner_others.c->runner_do_end_grav_force() */


  /* Check that we are not forming a sink in the accretion radius of another
     one. The new sink may be swallowed by the older one.) */
  const int scount = c->sinks.count;
  struct sink *restrict sinks = c->sinks.parts;

  for (int i = 0; i < scount ; i++) {

    /* Do not continue if the gas cannot form sink for any reason */
    if (!p->sink_data.can_form_sink){
      break ;
    }
    
    /* Get a hold of the ith sinks in ci. */
    struct sink *restrict si = &sinks[i];
    float r_acc_si = si->r_cut*cosmo->a;  /* Physical accretion radius of sink si */

    /* Compute the pairwise physical distance */
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
			  (float)(si->x[1] - c->loc[1]),
			  (float)(si->x[2] - c->loc[2])};

    const float dx[3] = {(px[0] - six[0]) * cosmo->a,
			 (px[1] - six[1]) * cosmo->a,
			 (px[2] - six[2]) * cosmo->a};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    if (r2 < (r_acc_si + r_acc_p)*(r_acc_si + r_acc_p)) {
      p->sink_data.is_overlapping_sink = 1 ;
    }

    /* Hill Sphere criterion */
    /* If the check has not been disabled, do it */
    if (!sink_props->sink_formation_hill_density_check) {
      float rho_hill = 0.f;

      /* Acceleration difference between p and si */
      /* How do we convert from comoving to physical ??*/
      const float da[3] = {(p->a_hydro[0] - si->gpart->a_grav[0]),
			   (p->a_hydro[1] - si->gpart->a_grav[1]),
			   (p->a_hydro[2] - si->gpart->a_grav[2])};

      float dx_times_da = dx[0]*da[0] + dx[1]*da[1] + dx[2] *da[2];

      rho_hill = 3.f*sink_props->x_hill*(-dx_times_da)/(4.f*M_PI*grav_props->G_Newton*r2);

      if (p->rho < rho_hill) {
	p->sink_data.can_form_sink = 0;
      }
    } /* End of Hill sphere criterion */
  } /* End of sink neighbour loop */
}


/*****************************************************************************/

/**
 * @brief Compute the comoving distance squared between a #sink and a #part.
 *
 * @param si The #sink particle.
 * @param pj The gas #part.
 */
INLINE static double sink_compute_sink_gas_comoving_distance_squared(struct sink* si, struct part* pj) {
    /* Comoving distance between sink and gas i */
   double dx[3] = {(si->x[0] - pj->x[0]),
		   (si->x[1] - pj->x[1]),
		   (si->x[2] - pj->x[2])};
   double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
   return r2;
}

/**
 * @brief Sort the neighbour array of a #sink in ascending order for the
 * distance between the parts and the #sink.
 *
 * @param sink The #sink particle.
 */
INLINE static void sink_sort_neighbour_array(struct sink* sink) {

  sink_neighbour_array* neighbour_array = sink->neighbour_array;
  size_t size = neighbour_array->size;
  double r2_i = 0.0;
  int j = 0;

  /* Array of size <= 1 are already sorted. */
  if (size <= 1) {
    neighbour_array->is_sorted = 1;
    return;
  }

  for (size_t i = 1 ; i < size ; i++) {
    /* Get a handle on the part. */
    struct part *const pi = neighbour_array->part_neighbours[i];

    /* Comoving distance between sink and gas i */
    r2_i = sink_compute_sink_gas_comoving_distance_squared(sink, pi);

    /* Get a handle on the part. */
    j = i - 1;
    struct part* pj = neighbour_array->part_neighbours[j];
    while (j >= 0 && sink_compute_sink_gas_comoving_distance_squared(sink, pj) > r2_i) {
      /* message("j = %i", j); */
      neighbour_array->part_neighbours[j+1] = pj;
      j -= 1;

      /* Attention, when j=-1 at this stage, this points outside the array ! However,
      since the while loop first verify that j>=0, this should not result in
      segfault. */
      pj = neighbour_array->part_neighbours[j];
    }
    neighbour_array->part_neighbours[j+1] = pi;
  }

#ifdef SWIFT_DEBUG_CHECKS
  double r2 = 0.0;
  double r2_previous = 0.0;

  message("Verification that the array is well sorted...");
  for (size_t i = 1 ; i < size ; ++i) {

    /* Get a handle on the part. */
    struct part *const pi = neighbour_array->part_neighbours[i];
    struct part *const pi_minus_one = neighbour_array->part_neighbours[i-1];

    /* Comoving distance between sink and gas i */
    r2 = sink_compute_sink_gas_comoving_distance_squared(sink, pi);
    r2_previous = sink_compute_sink_gas_comoving_distance_squared(sink, pi_minus_one);

    /* If particle i-1 is farther than particle i, the sort failed */
    if (r2_previous > r2) {
      error("The sort was not performed properly. Particle i=%zu (%lld) is closer than particle (i-1)=%zu (%lld) (out of %i neighbours).", i+1, pi->id, i, pi_minus_one->id, sink->N_neighbours);
    }
  }
#endif
  return ;
}

/**
 * @brief Compute the $\mathcal{W}$ quantity that normalises the radial accretion
 * timescale and the disk accretion timescale.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * @param sp The #sink.
 * @param sink_props The sink properties structure.
 */
INLINE static double sink_compute_W_regulated_accretion(struct sink* restrict sp,
							const struct sink_props* sink_props) {
  const int N_neighbours = sp->N_neighbours;
  const struct part* parts_neighbours = *(sp->neighbour_array->part_neighbours);
  const double H_sink = sp->r_cut / sink_props->f_interaction;
  double W = 0.0;

  /* Loop over the neighbours */
  for (int i = 0 ; i < N_neighbours ; ++i) {
    /* Get a pointer to the particle. */
    const struct part* pi = &parts_neighbours[i];

    /* Comoving distance */
    const double dx[3] = {(sp->x[0] - pi->x[0]),
			  (sp->x[1] - pi->x[1]),
			  (sp->x[2] - pi->x[2])};
    const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    const double r = sqrt(r2);

    /* Kernel function. Notice that it uses h = sp->r_cut. */
    const float ui = r / H_sink; /* The ratio of comoving removes the
				       cosmo->a. That's why a was not used. */
    /* Kernel of the part i */
    double wi = 0.0;
    kernel_eval_double(ui, &wi);
    wi *= pow_dimension_plus_one(1.0/H_sink) ;

    /* Compute the \mathcal{W} */
    W += pi->mass * wi / pi->rho;
  } /* End of neighbour loop */
  return W;

}

/**
 * @brief Compute the radial accretion timescale.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * @param sp The #sink.
 * @param cosmo The #cosmology structure.
 * @param sink_props The sink properties structure.
 */
INLINE static double sink_compute_radial_accretion_timescale(struct sink* restrict sp,
							     const struct cosmology* cosmo,
							     const struct sink_props* sink_props) {
  const int N_neighbours = sp->N_neighbours;
  const struct part* parts_neighbours = *(sp->neighbour_array->part_neighbours);
  const double W = sink_compute_W_regulated_accretion(sp, sink_props);
  const double H_sink = sp->r_cut / sink_props->f_interaction ;
  double t_radial = 0.0;
  double t_rad_numerator = 0.0;
  double t_rad_denominator = 0.0;

  /* Compute \mathcal{W} */
  /* Loop over the neighbours */
  for (int i = 0 ; i < N_neighbours ; ++i) {
    /* Get a pointer to the particle. */
    const struct part* pi = &parts_neighbours[i];

    /* Comoving distance */
    const double dx[3] = {(pi->x[0] - sp->x[0]),
			  (pi->x[1] - sp->x[1]),
			  (pi->x[2] - sp->x[2])};
    const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    const double r = sqrt(r2);

    /* Relative velocity */
    const double dv[3] = {(pi->v[0] - sp->v[0]),
			  (pi->v[1] - sp->v[1]),
			  (pi->v[2] - sp->v[2])};

    /* Kernel function. Notice that it uses h = sp->r_cut. */
    const float ui = r / H_sink; /* The ratio of comoving removes the
				       cosmo->a. That's why a was not used. */
    /* Kernel of the part i */
    double wi = 0.0;
    kernel_eval_double(ui, &wi);
    wi *= pow_dimension_plus_one(1.0/H_sink) ;

    /* Compute the numerator of t_rad : Sum_j m_j.*/
    t_rad_numerator += pi->mass ;

    /* Compute the numerator of t_rad : */
    t_rad_denominator += r * fabs(dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2]) * pi->mass * wi ;
  } /* End of neighbour loop */

  t_rad_numerator *= W ;
  t_rad_denominator *= 4.0 * M_PI;
  
  /* Compute the timescale */
  t_radial = t_rad_numerator / (t_rad_denominator);

  return t_radial;
}

/**
 * @brief Compute the disk accretion timescale.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * @param sp The #sink.
 * @param cosmo The #cosmology structure.
 * @param sink_props The sink properties structure.
 * @param phys_const The structure containing the physical constants.
 */
INLINE static double sink_compute_disk_accretion_timescale(struct sink* restrict sp,
							   const struct cosmology* cosmo,
							   const struct sink_props* sink_props,
							   const struct phys_const* phys_const) {
  const int N_neighbours = sp->N_neighbours;
  const struct part *parts_neighbours = *(sp->neighbour_array->part_neighbours);
  const double W = sink_compute_W_regulated_accretion(sp, sink_props);
  const double H_sink = sp->r_cut / sink_props->f_interaction ;
  double t_disk = 0.0;

  /* Compute \mathcal{W} */
  /* Loop over the neighbours */
  for (int i = 0 ; i < N_neighbours ; ++i) {
    /* Get a pointer to the particle. */
    const struct part* pi = &parts_neighbours[i];

    /* Comoving distance */
    const double dx[3] = {(sp->x[0] - pi->x[0]),
			  (sp->x[1] - pi->x[1]),
			  (sp->x[2] - pi->x[2])};
    const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    const double r = sqrt(r2);

    /* Kernel function. Notice that it uses h = sp->r_cut. */
    const float ui = r / H_sink; /* The ratio of comoving removes the
				       cosmo->a. That's why a was not used. */
    /* Kernel of the part i */
    double wi = 0.0;
    kernel_eval_double(ui, &wi);
    wi *= pow_dimension_plus_one(1.0/H_sink) ;

    /* Compute disk timescale */
    double sound_speed = hydro_get_comoving_soundspeed(pi);
    t_disk += sqrt(r) * pi->mass * wi / (pi->rho * sound_speed*sound_speed);
  } /* End of neighbour loop */

  /* Multiply by the constant */
  t_disk *= sqrt(phys_const->const_newton_G * sp->mass) / (sink_props->alpha_AMT * W)  ;

  return t_disk;
}

/**
 * @brief Compute the f exponent to compute the accretion timescale. The
 * accretion timescale is given by: t_radial^(1-f) * t_disk^f.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * @param sp The #sink.
 * @param cosmo The #cosmology structure.
 * @param sink_props The sink properties structure.
 * @param phys_const The structure containing the physical constants.
 */
INLINE static double sink_compute_f_accretion_timescale(struct sink* restrict sp,
							const struct cosmology* cosmo,
							const struct sink_props* sink_props,
							const struct phys_const* phys_const) {
  const size_t N_neighbours = sp->neighbour_array->size;
  const struct part* parts_neighbours = *(sp->neighbour_array->part_neighbours);
  const double H_sink = sp->r_cut / sink_props->f_interaction * cosmo->a;
  double f = 0.0;

  double E_rot_denominator = 0.0 ;
  double angular_momentum[3] = { 0, 0, 0 };

  /* First compute the total angular momentum within the interaction zone */
  /* Loop over the neighbours */
  for (size_t i = 0 ; i < N_neighbours ; ++i) {
    /* Get a pointer to the particle. */
    const struct part* pi = &parts_neighbours[i];

    const double mass_part = hydro_get_mass(pi);

    /* Physical distance */
    const double dx[3] = {(pi->x[0] - sp->x[0])*cosmo->a,
			  (pi->x[1] - sp->x[1])*cosmo->a,
			  (pi->x[2] - sp->x[2])*cosmo->a};

  /* Compute the relative physical velocity between p and pi */
    const double dv[3] = {(pi->v[0] - sp->v[0])*cosmo->a,
			  (pi->v[1] - sp->v[1])*cosmo->a,
			  (pi->v[2] - sp->v[2])*cosmo->a};

    /* Compute specific angular momentum between sp and pi */
    angular_momentum[0] += mass_part*(dx[1] * dv[2] - dx[2] * dv[1]);
    angular_momentum[1] += mass_part*(dx[2] * dv[0] - dx[0] * dv[2]);
    angular_momentum[2] += mass_part*(dx[0] * dv[1] - dx[1] * dv[0]);

  } /* End of neighbour loop */

  /* Then compute the rotation energy */
  /* Loop over the neighbours */
  for (size_t i = 0 ; i < N_neighbours ; ++i) {
    /* Get a pointer to the particle. */
    const struct part* pi = &parts_neighbours[i];

    const double mass_part = hydro_get_mass(pi);

    /* Comoving distance */
    const double dx[3] = {(pi->x[0] - sp->x[0])*cosmo->a,
			  (pi->x[1] - sp->x[1])*cosmo->a,
			  (pi->x[2] - sp->x[2])*cosmo->a};

    /* Compute the denominator of E_rot */
    double scalar_product = dx[0]*angular_momentum[0] + dx[1]*angular_momentum[1] + dx[2]*angular_momentum[2];
    E_rot_denominator += mass_part * scalar_product*scalar_product;

  } /* End of neighbour loop */

  double angular_momentum_2 = angular_momentum[0]*angular_momentum[0] + angular_momentum[1]*angular_momentum[1] + angular_momentum[2]*angular_momentum[2];
  double E_rot = angular_momentum_2*angular_momentum_2 / (2.0* E_rot_denominator);

  /* Finally, compute the gravitational energy */
  double E_grav = 0.0 ;

  float dummy;
  float pot_sink_i_sink = 0.0;
  float pot_sink_i_i = 0.0;
  float pot_ij_i = 0.0 ;
  float pot_ij_j = 0.0;

  /* First neighbour loop */
  for (size_t i = 0; i < N_neighbours ; ++i) {
    /* Get a pointer to the particle. */
    const struct part* pi = &parts_neighbours[i];
    const float h_i = pi->h;

    /* Physical distance */
    const float dx_sink_i[3] = {(pi->x[0] - sp->x[0])*cosmo->a,
				(pi->x[1] - sp->x[1])*cosmo->a,
				(pi->x[2] - sp->x[2])*cosmo->a};
    const float r2_sink_i = dx_sink_i[0] * dx_sink_i[0] + dx_sink_i[1] * dx_sink_i[1] + dx_sink_i[2] * dx_sink_i[2];

    /* Compute the potentials */
    runner_iact_grav_pp_full(r2_sink_i, H_sink*H_sink, 1.f/H_sink, 1.f/(H_sink*H_sink*H_sink), sp->mass, &dummy, &pot_sink_i_sink);
    runner_iact_grav_pp_full(r2_sink_i, h_i*h_i, 1.f/h_i, 1.f/(h_i*h_i*h_i), pi->mass, &dummy, &pot_sink_i_i);

    /* Add the contribution sink-i to the gravitional energy */
    E_grav += 0.5*(pi->mass*pot_sink_i_sink + sp->mass*pot_sink_i_i);

    /* Second neighbour loop */
    for (size_t j = 0; j < N_neighbours ; ++j) {
      if (j != i){
	/* Get a pointer to the particle. */
	const struct part* pj = &parts_neighbours[j];
	const float h_j = pj->h;

	/* Physical distance */
	const float dx_ij[3] = {(pi->x[0] - pj->x[0])*cosmo->a,
				(pi->x[1] - pj->x[1])*cosmo->a,
				(pi->x[2] - pj->x[2])*cosmo->a};
	const float r2_ij = dx_ij[0] * dx_ij[0] + dx_ij[1] * dx_ij[1] + dx_ij[2] * dx_ij[2];

	/* Compute the potentials */
	runner_iact_grav_pp_full(r2_ij, h_i*h_i, 1.f/h_i, 1.f/(h_i*h_i*h_i), pi->mass, &dummy, &pot_ij_i);
	runner_iact_grav_pp_full(r2_ij, h_j*h_j, 1.f/h_j, 1.f/(h_j*h_j*h_j), pj->mass, &dummy, &pot_ij_j);

	/* Add the contribution i-j to the gravitional energy */
	E_grav += 0.25*(pj->mass*pot_ij_i + pi->mass*pot_ij_j);
      }
    } /* End of second neighbour loop */
  } /* End of first neighbour loop */

  E_grav *= phys_const->const_newton_G;

  f = min(2.f*E_rot/fabs(E_grav), 1);

  return f;
}

/**
 * @brief Update the properties of a sink particles by swallowing
 * a gas particle in the regulated accretion scheme.
 *
 * @param sp The #sink to update.
 * @param p The #part that is swallowed.
 * @param xp The #xpart that is swallowed.
 * @param cosmo The current cosmological model.
 * @param delta_mass The mass removed from the #part during the regulated accretion.
 */
__attribute__((always_inline)) INLINE static void sink_swallow_part_regulated_accretion(
    struct sink* sp, struct part* p, struct xpart* xp,
    const struct cosmology* cosmo, const double delta_mass) {

  /* Notice that if a gas is entirely swallowed, we have delta_mass = gas_mass,
     since we have not set the gas_mass to 0. Also, sink_data.swallow_id is >= 0.
     However, if the gas is not entirely swallowed, delta_mass != gas_mass and
     sink_data.swallow_id < 0. Also, in this case, the orginal mass of the gas is
     gas_mass + delta_mass. */

  /* Get the current dynamical masses */
  const float gas_mass = hydro_get_mass(p);
  const float sink_mass = sp->mass; /* mass before update */

  /* store the mass of the sink part i */
  const float sink_mass_old = sp->mass;

  /* Increase the dynamical mass of the sink. */
  sp->mass += delta_mass;
  sp->gpart->mass += delta_mass;

  /* Physical velocity difference between the particles */
  const float dv[3] = {(sp->v[0] - p->v[0]) * cosmo->a_inv,
                       (sp->v[1] - p->v[1]) * cosmo->a_inv,
                       (sp->v[2] - p->v[2]) * cosmo->a_inv};

  /* Physical distance between the particles */
  const float dx[3] = {(sp->x[0] - p->x[0]) * cosmo->a,
                       (sp->x[1] - p->x[1]) * cosmo->a,
                       (sp->x[2] - p->x[2]) * cosmo->a};

  /* Collect the swallowed angular momentum */
  sp->swallowed_angular_momentum[0] +=
      delta_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
  sp->swallowed_angular_momentum[1] +=
      delta_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
  sp->swallowed_angular_momentum[2] +=
      delta_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

  /* Update the sink momentum */
  const float sink_mom[3] = {sink_mass * sp->v[0] + delta_mass * p->v[0],
                             sink_mass * sp->v[1] + delta_mass * p->v[1],
                             sink_mass * sp->v[2] + delta_mass * p->v[2]};

  sp->v[0] = sink_mom[0] / sp->mass;
  sp->v[1] = sink_mom[1] / sp->mass;
  sp->v[2] = sink_mom[2] / sp->mass;
  sp->gpart->v_full[0] = sp->v[0];
  sp->gpart->v_full[1] = sp->v[1];
  sp->gpart->v_full[2] = sp->v[2];

  /* If the sink swallowed entirely a gas particle */
  if (p->sink_data.swallow_id >= 0) {
    sp->number_of_gas_swallows++;
    sp->number_of_direct_gas_swallows++;
    /* Update the sink metal masses fraction */
    chemistry_add_part_to_sink(sp, p, sink_mass_old);
  } else {
    const float p_mass_orig = gas_mass + delta_mass ;
    /* Update the sink and also gas metal masses */
    chemistry_transfer_part_to_sink(sp, p, delta_mass, delta_mass / p_mass_orig);
  }

  /* Shall I update the number of neighbours ? This will be updated in the next
     build of the array... If i update it here, I should pay attention to the
     for loops with sink->N_neighbour and the neighbour_array size. */
  
#ifdef SWIFT_DEBUG_CHECKS
  const float dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  message("sink %lld swallow gas particle %lld. "
	  "(Mass = %e, "
	  "Delta_v = [%f, %f, %f] U_V, "
	  "Delta_x = [%f, %f, %f] U_L, "
	  "Delta_v_rad = %f)", sp->id, p->id, sp->mass, -dv[0], -dv[1], -dv[2], -dx[0], -dx[1], -dx[2],
      (dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]) / dr);
#endif
}

/**
 * @brief Compute an angular momentum feedback and apply it to the neighbour
 * and to the sink.
 *
 * This function is only used for the regulated accretion scheme.
 *
 * TODO: check momentum and angular momentum conservation with the
 *       feedback.
 *
 * @param e The #engine.
 * @param si The #sink.
 * @param dt_sink Timestep of the sink si.
 * @param t_disk Disk accretion timescale.
 */
INLINE static void sink_angular_momentum_feedback(struct engine* e, struct sink* restrict si, const double dt_sink, const double t_disk) {

  struct space *s = e->s;
  const struct cosmology* cosmo = e->cosmology;

  sink_neighbour_array* neighbour_array = si->neighbour_array;

  /* Compute angular momentum of si. Notice that a factor a and a^-1 cancel
     each other.  */
  double sink_angular_momentum[3] = {si->x[1] * si->v[2] - si->x[2] * si->v[1],
				     si->x[2] * si->v[0] - si->x[0] * si->v[2],
				     si->x[0] * si->v[1] - si->x[1] * si->v[0]};

  /* Sink angular momentum = "intrisic" angular momentum + swallowed angular
     momentum */
  sink_angular_momentum[0] += si->swallowed_angular_momentum[0];
  sink_angular_momentum[1] += si->swallowed_angular_momentum[1];
  sink_angular_momentum[2] += si->swallowed_angular_momentum[2];

  double sink_angular_momentum_norm = sqrt(sink_angular_momentum[0]*sink_angular_momentum[0] + sink_angular_momentum[1]*sink_angular_momentum[1] + sink_angular_momentum[2]*sink_angular_momentum[2]) ;

  /* Compute the feedback angular momentum norm */
  /* Notice that we use t_disk here, not t_acc */
  double delta_angular_momentum_norm_sink = sink_angular_momentum_norm * (1.0 - exp(- dt_sink/t_disk));

  /* Accumulator of the particle feedback velocities and and angular momentum */
  double delta_v_sink[3] = { 0.0 , 0.0, 0.0 } ;
  double delta_angular_momentum_sink[3] = { 0.0 , 0.0, 0.0 } ;

  /* Il faut itérer en deux fois parce que pour delta_v_j on a une somme sous
     le denominateur. Donc une première boulce calcule delta_v_j et une deuxième
     applique delta_v_j aux part et calcule la partie du sink. */

  double delta_v_j_denominator_vect[3] = {0.0, 0.0, 0.0};

  /* Neighbour loop ; compute the denominator of delta_v_j */
  for (size_t j = 0 ; j < neighbour_array->size ; ++j) {
    struct part *pj = neighbour_array->part_neighbours[j];

    /* Ignore inhibited particles (they have already been removed!) */
    if (part_is_inhibited(pj, e)) {
      continue;
    }

    /* If this part has already been entirely swallowed, skip it */
    if (pj->sink_data.swallow_id >= 0){
      continue;
    }

    /* Physical distance */
    const double dx[3] = {(pj->x[0] - si->x[0])*cosmo->a,
			  (pj->x[1] - si->x[1])*cosmo->a,
			  (pj->x[2] - si->x[2])*cosmo->a};

    const double dx_cross_L_sink[3] = {dx[1]*sink_angular_momentum[2] - dx[2]*sink_angular_momentum[1],
				       dx[2]*sink_angular_momentum[0] - dx[0]*sink_angular_momentum[2],
				       dx[0]*sink_angular_momentum[1] - dx[1]*sink_angular_momentum[0]} ;

    const double dx_cross_L_sink_cross_dx[3] = {dx_cross_L_sink[1]*dx[2] - dx_cross_L_sink[2]*dx[1],
						dx_cross_L_sink[2]*dx[0] - dx_cross_L_sink[0]*dx[2],
						dx_cross_L_sink[0]*dx[1] - dx_cross_L_sink[1]*dx[0]};

    delta_v_j_denominator_vect[0] += pj->mass * dx_cross_L_sink_cross_dx[0];
    delta_v_j_denominator_vect[1] += pj->mass * dx_cross_L_sink_cross_dx[1];
    delta_v_j_denominator_vect[2] += pj->mass * dx_cross_L_sink_cross_dx[2];

  } /* End of neighbour loop */

  const double delta_v_j_denominator = sqrt( delta_v_j_denominator_vect[0]* delta_v_j_denominator_vect[0] + delta_v_j_denominator_vect[1]*delta_v_j_denominator_vect[1] + delta_v_j_denominator_vect[2]*delta_v_j_denominator_vect[2] );

  /* Now, we can apply the velocity feedback kick to the neighbours */
  /* Neighbour loop  */
  for (size_t j = 0 ; j < neighbour_array->size ; ++j) {
    struct part *pj = neighbour_array->part_neighbours[j];

    /* Ignore inhibited particles (they have already been removed!) */
    if (part_is_inhibited(pj, e)) {
      continue;
    }

    /* If this part has already been entirely swallowed, skip it */
    if (pj->sink_data.swallow_id >= 0){
      continue;
    }

    /* Physical distance */
    const double dx[3] = {(pj->x[0] - si->x[0])*cosmo->a,
			  (pj->x[1] - si->x[1])*cosmo->a,
			  (pj->x[2] - si->x[2])*cosmo->a};

    const double L_sink_cross_dx[3] = {sink_angular_momentum[1]*dx[2] - sink_angular_momentum[2]*dx[1],
				       sink_angular_momentum[2]*dx[0] - sink_angular_momentum[0]*dx[2],
				       sink_angular_momentum[0]*dx[1] - sink_angular_momentum[1]*dx[0]} ;

    /* Finally compute the velocity feedback */
    const double delta_v_j[3] = {delta_angular_momentum_norm_sink*L_sink_cross_dx[0]/delta_v_j_denominator,
				 delta_angular_momentum_norm_sink*L_sink_cross_dx[1]/delta_v_j_denominator,
				 delta_angular_momentum_norm_sink*L_sink_cross_dx[2]/delta_v_j_denominator};


    /* Lock space here */
    /* Update the part and its gpart */
    pj->v[0] += delta_v_j[0];
    pj->v[1] += delta_v_j[1];
    pj->v[2] += delta_v_j[2];
    pj->gpart->v_full[0] += delta_v_j[0];
    pj->gpart->v_full[1] += delta_v_j[1];
    pj->gpart->v_full[2] += delta_v_j[2];

    /* Update the sink kick */
    delta_v_sink[0] -= pj->mass*delta_v_j[0];
    delta_v_sink[1] -= pj->mass*delta_v_j[1];
    delta_v_sink[2] -= pj->mass*delta_v_j[2];

    /* I'm not sure what to do with this. I'm supposed to apply it to the
       sink angular momentum, but if the velocity changed, the angular momentum
       will change accordingly, since its not an attribute. Nope, because the
       formula takes the relative distances between the sink and the particle
       j. The sink momentum is computed only with its postion, so I need to
       remove this delta_L to the swallowed angular momentum. */
    delta_angular_momentum_sink[0] -= pj->mass*(dx[1]*delta_v_j[2] - dx[2]*delta_v_j[1]);
    delta_angular_momentum_sink[1] -= pj->mass*(dx[0]*delta_v_j[2] - dx[2]*delta_v_j[0]);
    delta_angular_momentum_sink[2] -= pj->mass*(dx[0]*delta_v_j[1] - dx[1]*delta_v_j[0]);

  } /* End of neighbour loop */

    /* Divide the sink velocity kick by the sink mass */
  delta_v_sink[0] /= si->mass;
  delta_v_sink[1] /= si->mass;
  delta_v_sink[2] /= si->mass;

  /* Apply the sink velocity kick to the sink. Do not forget the gpart! */
  si->v[0] += delta_v_sink[0];
  si->v[1] += delta_v_sink[1];
  si->v[2] += delta_v_sink[2];
  si->gpart->v_full[0] += delta_v_sink[0];
  si->gpart->v_full[1] += delta_v_sink[1];
  si->gpart->v_full[2] += delta_v_sink[2];

  /* Apply the sink angular momentum kick to the sink. */
  si->swallowed_angular_momentum[0] += delta_angular_momentum_sink[0];
  si->swallowed_angular_momentum[1] += delta_angular_momentum_sink[1];
  si->swallowed_angular_momentum[2] += delta_angular_momentum_sink[2];
}

#endif /* SWIFT_GEAR_SINK_H */
