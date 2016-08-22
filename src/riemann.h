/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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
#ifndef SWIFT_RIEMANN_H
#define SWIFT_RIEMANN_H

/* gives us const_hydro_gamma and tells us which floating point type to use */
#include "const.h"
#include "error.h"
#include "float.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#if defined(RIEMANN_SOLVER_EXACT)

#define RIEMANN_SOLVER_IMPLEMENTATION "Exact Riemann solver (Toro 2009)"
#include "riemann/riemann_exact.h"

#elif defined(RIEMANN_SOLVER_TRRS)

#define RIEMANN_SOLVER_IMPLEMENTATION \
  "Two Rarefaction Riemann Solver (Toro 2009)"
#include "riemann/riemann_trrs.h"

#elif defined(RIEMANN_SOLVER_HLLC)

#define RIEMANN_SOLVER_IMPLEMENTATION \
  "Harten-Lax-van Leer-Contact Riemann solver (Toro 2009)"
#include "riemann/riemann_hllc.h"

#else

#error "Error: no Riemann solver selected!"

#endif

#endif /* SWIFT_RIEMANN_H */
