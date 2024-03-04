/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2024 Will J. Roper (w.roper@sussex.ac.uk).
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

/* Standard headers. */
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "cell.h"
#include "parser.h"
#include "space.h"
#include "zoom_region/zoom_init.h"

void make_mock_space(struct space *s) {

  /* Define the members we need for the test. */
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;
  s->nr_gparts = 18;

  /* Allocate memory for the gparts. */
  struct gpart *gparts =
      (struct gpart *)malloc(s->nr_gparts * sizeof(struct gpart));
  if (gparts == NULL) {
    error("Failed to allocate memory for gparts");
  }
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  /* We need the engine to be NULL for the logic. */
  s->e = NULL;

  /* Define the corners of the region */
  double cube_corners[8][3] = {
      {590, 590, 590}, {590, 515, 590}, {515, 590, 590}, {515, 515, 590},
      {590, 590, 515}, {590, 515, 515}, {515, 590, 515}, {515, 515, 515}};

  /* Loop over the gparts and set up baxckground and zoom particles. */
  for (size_t i = 0; i < s->nr_gparts; i++) {
    gparts[i].mass = 1.0;

    /* Handle background and zoom region particles differently. */
    if (i < 10) {
      /* Set background particles to be evenly spaced. */
      gparts[i].x[0] = s->dim[0] / s->nr_gparts * i;
      gparts[i].x[1] = s->dim[1] / s->nr_gparts * i;
      gparts[i].x[2] = s->dim[2] / s->nr_gparts * i;
      gparts[i].type = swift_type_dark_matter_background;

    } else {
      /* Set zoom region particles to be at the corners of the region. */
      gparts[i].x[0] = cube_corners[i - 10][0];
      gparts[i].x[1] = cube_corners[i - 10][1];
      gparts[i].x[2] = cube_corners[i - 10][2];
      gparts[i].type = swift_type_dark_matter;
    }

    s->gparts = gparts;
  }
}

void make_mock_cells(struct space *s) {
  /* Get the zoom properties */
  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Allocate cells. */
  s->cells_top = (struct cell *)malloc(s->nr_cells * sizeof(struct cell));
  if (s->cells_top == NULL) {
    error("Failed to allocate memory for cells");
  }
  bzero(s->cells_top, s->nr_cells * sizeof(struct cell));

  /* Get some zoom region properties */
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const int buffer_offset = zoom_props->buffer_cell_offset;
  const double zoom_bounds[3] = {zoom_props->region_lower_bounds[0],
                                 zoom_props->region_lower_bounds[1],
                                 zoom_props->region_lower_bounds[2]};
  const double buffer_bounds[3] = {zoom_props->buffer_lower_bounds[0],
                                   zoom_props->buffer_lower_bounds[1],
                                   zoom_props->buffer_lower_bounds[2]};

  /* Loop over zoom cells and set locations and initial values */
  for (int i = 0; i < zoom_props->cdim[0]; i++) {
    for (int j = 0; j < zoom_props->cdim[1]; j++) {
      for (int k = 0; k < zoom_props->cdim[2]; k++) {
        const size_t cid = cell_getid(zoom_props->cdim, i, j, k);
        struct cell *c = &s->cells_top[cid];
        c->loc[0] = (i * zoom_props->width[0]) + zoom_bounds[0];
        c->loc[1] = (j * zoom_props->width[1]) + zoom_bounds[1];
        c->loc[2] = (k * zoom_props->width[2]) + zoom_bounds[2];
        c->width[0] = zoom_props->width[0];
        c->width[1] = zoom_props->width[1];
        c->width[2] = zoom_props->width[2];

        c->type = cell_type_zoom;
        c->subtype = cell_subtype_regular;
      }
    }
  }

  /* Loop over natural cells and set locations and initial values */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        const size_t cid = cell_getid_offset(s->cdim, bkg_cell_offset, i, j, k);
        struct cell *c = &s->cells_top[cid];
        c->loc[0] = i * s->width[0];
        c->loc[1] = j * s->width[1];
        c->loc[2] = k * s->width[2];
        c->width[0] = s->width[0];
        c->width[1] = s->width[1];
        c->width[2] = s->width[2];
        c->type = cell_type_bkg;
        c->subtype = cell_subtype_regular;
      }
    }
  }

  /* Loop over buffer cells and set locations and initial values */
  for (int i = 0; i < zoom_props->buffer_cdim[0]; i++) {
    for (int j = 0; j < zoom_props->buffer_cdim[1]; j++) {
      for (int k = 0; k < zoom_props->buffer_cdim[2]; k++) {
        const size_t cid =
            cell_getid_offset(zoom_props->buffer_cdim, buffer_offset, i, j, k);
        struct cell *c = &s->cells_top[cid];
        c->loc[0] = (i * zoom_props->buffer_width[0]) + buffer_bounds[0];
        c->loc[1] = (j * zoom_props->buffer_width[1]) + buffer_bounds[1];
        c->loc[2] = (k * zoom_props->buffer_width[2]) + buffer_bounds[2];
        c->width[0] = zoom_props->buffer_width[0];
        c->width[1] = zoom_props->buffer_width[1];
        c->width[2] = zoom_props->buffer_width[2];
        c->type = cell_type_buffer;
        c->subtype = cell_subtype_regular;
      }
    }
  }

  /* Label void cells. */
  for (int cid = zoom_props->buffer_cell_offset;
       cid < zoom_props->buffer_cell_offset + zoom_props->nr_buffer_cells;
       cid++) {

    /* Get the cell */
    struct cell *c = &s->cells_top[cid];

    /* Get the middle of the cell. */
    double mid[3] = {c->loc[0] + 0.5 * c->width[0],
                     c->loc[1] + 0.5 * c->width[1],
                     c->loc[2] + 0.5 * c->width[2]};

    /* Label this cell if it contains the zoom region. */
    if ((mid[0] > s->zoom_props->region_lower_bounds[0]) &&
        (mid[0] < s->zoom_props->region_upper_bounds[0]) &&
        (mid[1] > s->zoom_props->region_lower_bounds[1]) &&
        (mid[1] < s->zoom_props->region_upper_bounds[1]) &&
        (mid[2] > s->zoom_props->region_lower_bounds[2]) &&
        (mid[2] < s->zoom_props->region_upper_bounds[2])) {
      c->subtype = cell_type_void;
    }
  }

  /* Label the empty cells. */
  for (int cid = zoom_props->bkg_cell_offset;
       cid < zoom_props->bkg_cell_offset + zoom_props->nr_bkg_cells; cid++) {

    /* Get this cell. */
    struct cell *c = &s->cells_top[cid];

    /* Get the middle of the cell. */
    double mid[3] = {c->loc[0] + 0.5 * c->width[0],
                     c->loc[1] + 0.5 * c->width[1],
                     c->loc[2] + 0.5 * c->width[2]};

    /* Assign the cell type. */
    if ((mid[0] > s->zoom_props->buffer_lower_bounds[0]) &&
        (mid[0] < s->zoom_props->buffer_upper_bounds[0]) &&
        (mid[1] > s->zoom_props->buffer_lower_bounds[1]) &&
        (mid[1] < s->zoom_props->buffer_upper_bounds[1]) &&
        (mid[2] > s->zoom_props->buffer_lower_bounds[2]) &&
        (mid[2] < s->zoom_props->buffer_upper_bounds[2])) {
      c->subtype = cell_type_empty;
    }
  }
}

int main(int argc, char *argv[]) {

  const char *input_file = argv[1];

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Create a space structure. */
  struct space *s = malloc(sizeof(struct space));
  if (s == NULL) {
    error("Failed to allocate memory for space");
  }
  bzero(s, sizeof(struct space));
  make_mock_space(s);

  /* Flag that we are running a zoom. */
  s->with_zoom_region = 1;

  /* Run the zoom_init function. */
  zoom_region_init(&param_file, s, 0);

  /* Make the cells. */
  make_mock_cells(s);

  /* Get a zoom, buffer and background cell using cell_getid functions. */
  const size_t zoom_cell_id = cell_getid_from_pos(s, 500.0, 500.0, 500.0);
  const size_t buffer_cell_id = cell_getid_from_pos(s, 450.0, 550., 450.0);
  const size_t bkg_cell_id = cell_getid_from_pos(s, 10.0, 990.0, 800.0);

  /* Test that the write type of cell was returned. */
  assert(s->cells_top[zoom_cell_id].type == cell_type_zoom);
  assert(s->cells_top[buffer_cell_id].type == cell_type_buffer);
  assert(s->cells_top[bkg_cell_id].type == cell_type_bkg);

  free(s->cells_top);
  free(s->zoom_props);
  free(s);

  return 0;
}
