/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Local includes */
#include "chemistry.h"
#include "cooling.h"
#include "minmax.h"
#include "random.h"
#include "sink_part.h"
#include "sink_properties.h"
#include "feedback.h"


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


  float random_number= random_unit_interval_part_ID_and_loop_idx(sink->id, rloop, e->ti_current,random_number_sink_formation);
        
  struct feedback_props*  feedback_props = e->feedback_props;
  
  /* Pick the correct table. (if only one table, threshold is < 0) */
  
  const float metal =
      chemistry_get_sink_total_iron_mass_fraction_for_feedback(sink);
  const float threshold = feedback_props->metallicity_max_first_stars;

  
  struct stellar_model* model;
  double minimal_discrete_mass;
  
  if(metal >= threshold) {
    model = &feedback_props->stellar_model;
    minimal_discrete_mass = sink_props->minimal_discrete_mass_first_stars;
  } else {
    model = &feedback_props->stellar_model_first_stars;
    minimal_discrete_mass = sink_props->minimal_discrete_mass;
  }  
  
  const struct initial_mass_function *imf=&model->imf;
  
  
  if (random_number < imf->sink_Pc) {
    // we are dealing with the continous part of the IMF
    sink->target_mass = imf->sink_stellar_particle_mass;
    sink->target_type = 0;
  } else {
    // we are dealing with the discrete part of the IMF  
    random_number = random_unit_interval_part_ID_and_loop_idx(sink->id, rloop+1, e->ti_current,random_number_sink_formation);
    double m = random_sample_power_law(minimal_discrete_mass, imf->mass_max, imf->exp[imf->n_parts-1], random_number);
    sink->target_mass = m;
    sink->target_type = 1;
  }
  
}        




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
    struct sink* sp, const struct sink_props* sink_props) {

  sp->r_cut = sink_props->cut_off_radius;
  sp->time_bin = 0;

  sp->number_of_gas_swallows = 0;
  sp->number_of_direct_gas_swallows = 0;
  sp->number_of_sink_swallows = 0;
  sp->number_of_direct_sink_swallows = 0;
  sp->swallowed_angular_momentum[0] = 0.f;
  sp->swallowed_angular_momentum[1] = 0.f;
  sp->swallowed_angular_momentum[2] = 0.f;

  sink_mark_sink_as_not_swallowed(&sp->merger_data);
    
}

/**
 * @brief Initialisation of particle data before the hydro density loop.
 * Note: during initalisation (space_init)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sink_init_part(
    struct part* restrict p) {

  struct sink_part_data* cpd = &p->sink_data;

  cpd->can_form_sink = 1;
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

  const float temperature_max = sink_props->maximal_temperature;
  const float temperature = cooling_get_temperature(phys_const, hydro_props, us,
                                                    cosmo, cooling, p, xp);

  const float density_threashold = sink_props->density_threashold;
  const float density = hydro_get_physical_density(p, cosmo);

  if (density > density_threashold && temperature < temperature_max) {
    //message("forming a sink particle ! %lld", p->id);
    return 1;
  }

  return 0;
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
  sink->r_cut = e->sink_properties->cut_off_radius;

  /* Flag it as not swallowed */
  sink_mark_sink_as_not_swallowed(&sink->merger_data);
  
  /* Additional initialisation */
  
  /* setup the target mass for sink star formation */
  sink_update_target_mass(sink, sink_props, e, 0);
  
  
  /* Copy the chemistry properties */
  chemistry_copy_sink_properties(p, xp, sink);  
  
  
  /* test the mass distribution */
  //for (int i=0;i<1000000;i++)
  //  {
  //    sink_update_target_mass(sink, sink_props, e, i);
  //    message("%g %d",sink->target_mass,sink->target_type);
  //  }
  //exit(-1);  
  
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

  message("sink %lld swallow gas particle %lld", sp->id,p->id);

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

  /*
  const float dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  message(
      "sink %lld swallowing gas particle %lld "
      "(Delta_v = [%f, %f, %f] U_V, "
      "Delta_x = [%f, %f, %f] U_L, "
      "Delta_v_rad = %f)",
      sp->id, p->id, -dv[0], -dv[1], -dv[2], -dx[0], -dx[1], -dx[2],
      (dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]) / dr);
  */

  /* Update the sink metal masses fraction */
  chemistry_add_part_to_sink(sp, p, msp_old);

  /* This sink swallowed a gas particle */
  sp->number_of_gas_swallows++;
  sp->number_of_direct_gas_swallows++;
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


  message("sink %lld swallow sink particle %lld", spi->id,spj->id);


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
  spi->number_of_sink_swallows+=spj->number_of_sink_swallows;
  spi->number_of_gas_swallows+=spj->number_of_gas_swallows;

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

  if(sink->mass > sink->target_mass*phys_const->const_solar_mass) 
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
INLINE static void sink_star_formation_separate_particles(const struct engine* e,
							  struct sink* si,
							  struct spart* sp) {
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
  sp->mass =  sink->target_mass*phys_const->const_solar_mass;   
  
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
  //sp->tracers_data = sink->tracers_data;

  /* Move over the splitting data */
  //sp->split_data = sink->split_data;

  ///* Store the birth density in the star particle */
  //sp->sf_data.birth_density = hydro_get_physical_density(p, cosmo);

  ///* Store the birth temperature*/
  //sp->sf_data.birth_temperature = cooling_get_temperature(
  //    phys_const, hydro_props, us, cosmo, cooling, p, xp);


  /* Copy the chemistry properties */
  chemistry_copy_sink_properties_to_star(sink, sp);
  

  /* Copy the progenitor id */
  sp->sf_data.progenitor_id = sink->id;




}

#endif /* SWIFT_GEAR_SINK_H */
