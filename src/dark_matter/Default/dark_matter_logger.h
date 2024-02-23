/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Camila Correa (camila.correa@cea.fr)
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
 *******************************************************************************/
#ifndef SWIFT_DEFAULT_DARK_MATTER_LOGGER_H
#define SWIFT_DEFAULT_DARK_MATTER_LOGGER_H

/* Local includes */
#include "cell.h"
#include "dark_matter_logger_struct.h"
#include "units.h"


/**
 * @brief log a SIDM event
 *
 * @param dmi the dmpart of the SIDM event
 * @param cosmo the cosmology struct
 * @param num_events number of events per time step
 */
INLINE static void dark_matter_log_num_events(struct sidm_history *sidm_history, const int num_events) {}


/**
 * @brief log a SIDM event
 *
 * @param dmi the dmpart of the SIDM event
 * @param cosmo the cosmology struct
 * @param num_events number of events per time step
 */
INLINE static void dark_matter_log_total_kinetic_energy(
           struct sidm_history *sidm_history,
           const double energy_before, const double energy_after) {}

/**
 * @brief add a star formation history struct to an other star formation history
 * struct
 *
 * @param sf_add the star formation struct which we want to add to the star
 * formation history
 * @param sf_update the star formation structure which we want to update
 */
INLINE static void dark_matter_logger_add(
             struct sidm_history *sh_update,
             const struct sidm_history *sh_add) {}

/**
 * @brief add a star formation history struct to the engine star formation
 * history accumulator struct
 *
 * @param sf_add the star formation accumulator struct which we want to add to
 * the star formation history
 * @param sf_update the star formation structure which we want to update
 */
INLINE static void dark_matter_logger_add_to_accumulator(
        struct sidm_history_accumulator *sh_update,
        const struct sidm_history *sh_add, long long n_parts_active) {}


/**
 * @brief Initialize the SIDM history structure in the #engine
 *
 * @param sh The pointer to the sidm history structure
 */
INLINE static void dark_matter_logger_accumulator_init(
       struct sidm_history_accumulator *sh) {}

/**
 * @brief Initialize the SIDM history structure in the #engine
 *
 * @param sfh The pointer to the SIDM history structure
 */
INLINE static void dark_matter_logger_init(struct sidm_history *sh) {}

/**
 * @brief Write the final status to a file
 *
 * @param fp The file to write to.
 * @param time the simulation time (time since Big Bang) in internal units.
 * @param a the scale factor.
 * @param z the redshift.
 * @param sh the #sidm_history struct.
 * @param step The time-step of the simulation.
 */
INLINE static void dark_matter_write_to_log_file(
           FILE *fp, const double time, const double a, const double z,
           struct sidm_history_accumulator sh, const int step) {}

/**
 * @brief Initialize the SIDM logger file
 *
 * @param fp the file pointer
 * @param us The current internal system of units.
 * @param phys_const Physical constants in internal units
 */
INLINE static void dark_matter_logger_init_log_file(
               FILE *fp, const struct unit_system *restrict us,
               const struct phys_const *phys_const) {}

#ifdef WITH_MPI
/**
 * @brief Do the MPI communication for the SIDM events logger
 *
 * @param e the engine we are running
 */
INLINE static void dark_matter_logger_MPI_Reduce(const struct engine *e) {}
#endif


#endif
