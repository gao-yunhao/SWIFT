/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_ENTROPY_FLOOR_EAGLE_H
#define SWIFT_ENTROPY_FLOOR_EAGLE_H

#include "adiabatic_index.h"
#include "cosmology.h"
#include "hydro_properties.h"
#include "parser.h"
#include "units.h"

/**
 * @file src/entropy_floor/EAGLE/entropy_floor.h
 * @brief Entropy floor used in the EAGLE model
 */

/**
 * @brief Properties of the entropy floor in the EAGLE model.
 */
struct entropy_floor_properties {

  /*! Density threshold for the Jeans floor in Hydrogen atoms per cubic cm */
  float Jeans_density_threshold_H_p_cm3;

  /*! Density threshold for the Jeans floor in internal units */
  float Jeans_density_threshold;

  /*! Inverse of the density threshold for the Jeans floor in internal units */
  float Jeans_density_threshold_inv;

  /*! Over-density threshold for the Jeans floor */
  float Jeans_over_density_threshold;

  /*! Slope of the Jeans floor power-law */
  float Jeans_gamma_effective;

  /*! Temperature of the Jeans floor at the density threshold in Kelvin */
  float Jeans_temperature_norm_K;

  /*! Temperature of the Jeans floor at the density thresh. in internal units */
  float Jeans_temperature_norm;

  /*! Pressure of the Jeans floor at the density thresh. in internal units */
  float Jeans_pressure_norm;

  /*! Density threshold for the Cool floor in Hydrogen atoms per cubic cm */
  float Cool_density_threshold_H_p_cm3;

  /*! Density threshold for the Cool floor in internal units */
  float Cool_density_threshold;

  /*! Inverse of the density threshold for the Cool floor in internal units */
  float Cool_density_threshold_inv;

  /*! Over-density threshold for the Cool floor */
  float Cool_over_density_threshold;

  /*! Slope of the Cool floor power-law */
  float Cool_gamma_effective;

  /*! Temperature of the Cool floor at the density threshold in Kelvin */
  float Cool_temperature_norm_K;

  /*! Temperature of the Cool floor at the density thresh. in internal units */
  float Cool_temperature_norm;

  /*! Pressure of the Cool floor at the density thresh. in internal units */
  float Cool_pressure_norm;
};

/**
 * @brief Compute the entropy floor of a given #part.
 *
 * Note that the particle is not updated!!
 *
 * @param p The #part.
 * @param cosmo The cosmological model.
 * @param props The properties of the entropy floor.
 */
float entropy_floor(const struct part *p, const struct cosmology *cosmo,
                    const struct entropy_floor_properties *props) {

  /* Physical density in internal units */
  const float rho = hydro_get_physical_density(p, cosmo);

  /* Critical density at this redshift.
   * Recall that this is 0 in a non-cosmological run */
  const float rho_crit = cosmo->critical_density;

  float pressure = 0.;

  /* Are we in the regime of the Jeans equation of state? */
  if ((rho >= rho_crit * props->Jeans_over_density_threshold) &&
      (rho >= props->Jeans_density_threshold)) {

    const float pressure_Jeans = props->Jeans_pressure_norm *
                                 powf(rho * props->Jeans_density_threshold_inv,
                                      props->Jeans_gamma_effective);

    pressure = max(pressure, pressure_Jeans);
  }

  /* Are we in the regime of the Cool equation of state? */
  if ((rho >= rho_crit * props->Cool_over_density_threshold) &&
      (rho >= props->Cool_density_threshold)) {

    const float pressure_Cool = props->Cool_pressure_norm *
                                powf(rho * props->Cool_density_threshold_inv,
                                     props->Cool_gamma_effective);

    pressure = max(pressure, pressure_Cool);
  }

  /* Convert to an entropy.
   * (Recall that the entropy is the same in co-moving and phycial frames) */
  return gas_entropy_from_pressure(rho, pressure);
}

/**
 * @brief Initialise the entropy floor by reading the parameters and converting
 * to internal units.
 *
 * @param params The YAML parameter file.
 * @param us The system of units used internally.
 * @param phys_cont The physical constants.
 * @param props The entropy floor properties to fill.
 */
void entropy_floor_init(struct swift_params *params,
                        const struct unit_system *us,
                        const struct phys_const *phys_const,
                        const struct hydro_props *hydro_props,
                        struct entropy_floor_properties *props) {

  /* Read the parameters in the units they are set */
  props->Jeans_density_threshold_H_p_cm3 = parser_get_param_float(
      params, "EAGLEEntropyFloor:Jeans_density_threshold_H_p_cm3");
  props->Jeans_over_density_threshold = parser_get_param_float(
      params, "EAGLEEntropyFloor:Jeans_over_density_threshold");
  props->Jeans_temperature_norm_K = parser_get_param_float(
      params, "EAGLEEntropyFloor:Jeans_temperature_norm_K");
  props->Jeans_gamma_effective =
      parser_get_param_float(params, "EAGLEEntropyFloor:Jeans_gamma_effective");

  props->Cool_density_threshold_H_p_cm3 = parser_get_param_float(
      params, "EAGLEEntropyFloor:Cool_density_threshold_H_p_cm3");
  props->Cool_over_density_threshold = parser_get_param_float(
      params, "EAGLEEntropyFloor:Cool_over_density_threshold");
  props->Cool_temperature_norm_K = parser_get_param_float(
      params, "EAGLEEntropyFloor:Cool_temperature_norm_K");
  props->Cool_gamma_effective =
      parser_get_param_float(params, "EAGLEEntropyFloor:Cool_gamma_effective");

  /* Initial Hydrogen abundance (mass fraction) */
  const float X_H = hydro_props->hydrogen_mass_fraction;

  /* Now convert to internal units */
  props->Jeans_temperature_norm =
      props->Jeans_temperature_norm_K /
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  props->Jeans_density_threshold =
      props->Jeans_density_threshold /
      units_cgs_conversion_factor(us, UNIT_CONV_VOLUME) *
      phys_const->const_proton_mass / X_H;

  props->Cool_temperature_norm =
      props->Cool_temperature_norm_K /
      units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  props->Cool_density_threshold =
      props->Cool_density_threshold /
      units_cgs_conversion_factor(us, UNIT_CONV_VOLUME) *
      phys_const->const_proton_mass / X_H;

  /* We assume neutral gas */
  const float mean_molecular_weight = hydro_props->mu_neutral;

  /* Get the common terms */
  props->Jeans_density_threshold_inv = 1.f / props->Jeans_density_threshold;
  props->Cool_density_threshold_inv = 1.f / props->Cool_density_threshold;

  /* P_norm = (k_B * T) / (m_p * mu) * rho_threshold */
  props->Jeans_pressure_norm =
      (phys_const->const_boltzmann_k * props->Jeans_temperature_norm) /
      (phys_const->const_proton_mass * mean_molecular_weight) *
      props->Jeans_density_threshold;

  props->Cool_pressure_norm =
      (phys_const->const_boltzmann_k * props->Cool_temperature_norm) /
      (phys_const->const_proton_mass * mean_molecular_weight) *
      props->Cool_density_threshold;
}

/**
 * @brief Write an entropy floor struct to the given FILE as a stream of bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
void entropy_floor_struct_dump(const struct entropy_floor_properties *props,
                               FILE *stream) {

  restart_write_blocks((void *)props, sizeof(struct entropy_floor_properties),
                       1, stream, "entropy floor", "entropy floor properties");
}

/**
 * @brief Restore a entropy floor struct from the given FILE as a stream of
 * bytes.
 *
 * @param props the struct
 * @param stream the file stream
 */
void entropy_floor_struct_restore(struct entropy_floor_properties *props,
                                  FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct entropy_floor_properties), 1,
                      stream, NULL, "entropy floor properties");
}

#endif /* SWIFT_ENTROPY_FLOOR_EAGLE_H */
