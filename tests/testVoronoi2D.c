/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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

#include "hydro/Shadowswift/voronoi2d_algorithm.h"

/* Number of cells used to test the 2D interaction algorithm */
#define TESTVORONOI2D_NUMCELL 100

/* declare the global variables needed to store the simulation box size */
VORONOI_DECLARE_GLOBAL_VARIABLES()

int main() {

  /* initialize simulation box */
  float anchor[3] = {-0.5f, -0.5f, -0.5f};
  float side[3] = {2.0f, 2.0f, 2.0f};
  voronoi_set_box(anchor, side);

  /* test initialization and finalization algorithms */
  {
    struct voronoi_cell cell;
    double x[3] = {0.5, 0.5, 0.5};

    voronoi_cell_init(&cell, x);

    float maxradius = voronoi_cell_finalize(&cell);

    assert(maxradius == 2.0f * sqrtf(2.0f));

    assert(cell.volume == 4.0f);

    assert(cell.centroid[0] == 0.5f);
    assert(cell.centroid[1] == 0.5f);
  }

  /* test interaction algorithm */
  {
    /* create 100 cells with random positions in [0,1]x[0,1] */
    struct voronoi_cell cells[TESTVORONOI2D_NUMCELL];
    double x[2];
    float dx[2];
    int i, j;
    float Atot;
    struct voronoi_cell *cell_i, *cell_j;

    for (i = 0; i < TESTVORONOI2D_NUMCELL; ++i) {
      x[0] = ((double)rand()) / ((double)RAND_MAX);
      x[1] = ((double)rand()) / ((double)RAND_MAX);
      voronoi_cell_init(&cells[i], x);
#ifdef VORONOI_VERBOSE
      message("cell[%i]: %g %g", i, x[0], x[1]);
#endif
    }

    /* interact the cells (with periodic boundaries) */
    for (i = 0; i < TESTVORONOI2D_NUMCELL; ++i) {
      cell_i = &cells[i];
      for (j = 0; j < TESTVORONOI2D_NUMCELL; ++j) {
        if (i != j) {
          cell_j = &cells[j];
          dx[0] = cell_i->x[0] - cell_j->x[0];
          dx[1] = cell_i->x[1] - cell_j->x[1];
          /* periodic boundaries */
          if (dx[0] >= 0.5f) {
            dx[0] -= 1.0f;
          }
          if (dx[0] < -0.5f) {
            dx[0] += 1.0f;
          }
          if (dx[1] >= 0.5f) {
            dx[1] -= 1.0f;
          }
          if (dx[1] < -0.5f) {
            dx[1] += 1.0f;
          }
#ifdef VORONOI_VERBOSE
          message("Cell %i before:", i);
          voronoi_print_cell(&cells[i]);
          message("Interacting cell %i with cell %i (%g %g, %g %g", i, j,
                  cells[i].x[0], cells[i].x[1], cells[j].x[0], cells[j].x[1]);
#endif
          voronoi_cell_interact(cell_i, dx, j);
        }
      }
    }

    /* print the cells to the stdout */
    for (i = 0; i < TESTVORONOI2D_NUMCELL; ++i) {
      voronoi_print_cell(&cells[i]);
      voronoi_cell_finalize(&cells[i]);
      Atot += cells[i].volume;
    }

    assert(fabs(Atot - 1.0f) < 1.e-6);
  }

  return 0;
}
