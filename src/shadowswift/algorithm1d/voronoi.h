//
// Created by yuyttenh on 9/12/22.
//

#ifndef SWIFTSIM_VORONOI_H
#define SWIFTSIM_VORONOI_H

#include "../queues.h"
#include "./delaunay.h"

/**
 * @brief Voronoi interface.
 *
 * An interface is a connection between two neighbouring Voronoi cells. It is
 * completely defined by the indices of the generators that generate the two
 * neighbouring cells, a surface area and a midpoint position.
 */
struct voronoi_pair {
  /*! idx of the particle on the right of this pair in its respective swift
   * cell. Since the left particle is always local this is also the index of the
   * corresponding cell in this voronoi tesselation. */
  int left_idx;

  /*! idx of the particle on the right of this pair in its respective swift cell
   * if that cell is the same as the cell holding this Voronoi tesselation (i.e.
   * the particle is local) or in the super cell of its respective swift cell if
   * that swift cell is foreign. For local particles, this is also the index of
   * the corresponding cell in this voronoi tesselation. */
  int right_idx;

  /*! Real sid of this pair (boundary faces are stored under sid 27) */
  int sid;

  /*! Surface area of the interface. */
  double surface_area;

  /*! Midpoint of the interface. */
  double midpoint[3];
};

/**
 * @brief Voronoi grid.
 *
 * The grid stores a copy of the coordinates of the grid generators, the
 * coordinates of the grid vertices and the edge connections that make up the
 * grid. For every generator, it stores the number of vertices for the cell
 * generated by it, and the offset of the cell edges in the edge array.
 */
struct voronoi {

  /*! @brief Voronoi cell pairs. We store these per (SWIFT) cell, i.e. pairs[0]
   *  contains all pairs that are completely contained within this cell, while
   *  pairs[1] corresponds to pairs crossing the boundary between this cell and
   *  the cell with coordinates that are lower in all coordinate directions (the
   *  cell to the left, front, bottom, sid=0), and so on. SID 27 is reserved
   *  for boundary condition particles (e.g. reflective boundary conditions). */
  struct voronoi_pair *pairs[28];

  /*! @brief Current number of pairs per cell index. */
  int pair_count[28];

  /*! @brief Allocated number of pairs per cell index. */
  int pair_size[28];

  /*! @brief cell pair connection. Queue of 2-tuples containing the index of
   * the pair and the sid of the pair */
  struct int2_lifo_queue cell_pair_connections;
};

static inline int voronoi_add_pair(struct voronoi *v, const struct delaunay *d,
                                   int gen_idx_in_d, int ngb_gen_idx_in_d,
                                   struct part *parts, double ax);

/**
 * @brief Initialise the Voronoi grid based on the given Delaunay tessellation.
 *
 * This function allocates the memory for the Voronoi grid arrays and creates
 * the grid in linear time by
 *  1. Computing the grid vertices as the midpoints of Delaunay lines.
 *  2. Looping over all vertices and for each vertex creating the two faces for
 *     that vertex.
 *
 * During the second step, the geometrical properties (cell centroid, volume
 * and face midpoint) are computed as well.
 *
 * @param v Voronoi grid.
 * @param d Delaunay tessellation (read-only).
 * @param parts Array of parts of the SWIFT cell for which we are building the
 * voronoi tesselation.
 * @param part_is_active Flags indicating whether the particle is active.
 * @param count The number of hydro particles in #parts
 */
static inline void voronoi_build(struct voronoi *v, struct delaunay *d,
                                 struct part *parts) {
  voronoi_assert(d->vertex_end > 0);

  /* loop over the lines in the Delaunay tessellation and compute their
     midpoints. These happen to be the vertices of the Voronoi grid (because
     they are the points of equal distance to 2 generators) */
  double *vertices = (double *)malloc((d->line_index - 2) * sizeof(double));
  for (int i = 0; i < d->line_index - 2; ++i) {
    struct line *l = &d->lines[i + 2];
    int v0 = l->vertices[0];
    int v1 = l->vertices[1];

    /* if the line is not linked to a non-ghost, non-dummy vertex,
     * corresponding to an active particle, it is not a grid vertex and we can
     * skip it. */
    if ((v0 >= d->vertex_end || v0 < d->vertex_start) &&
        (v1 >= d->vertex_end || v1 < d->vertex_start))
      continue;

    if (v0 < d->vertex_start) {
      /* This could mean that a neighbouring cell of this grids cell is empty!
       * Or that we did not add all the necessary ghost vertices to the delaunay
       * tesselation. */
      error(
          "Vertex is part of line with Dummy vertex! This could mean that "
          "one of the neighbouring cells is empty.");
    }
    double v0x = d->vertices[v0];

    if (v1 < d->vertex_start) {
      error(
          "Vertex is part of line with Dummy vertex! This could mean that "
          "one of the neighbouring cells is empty.");
    }
    double v1x = d->vertices[v1];

    /* Compute midpoint */
    vertices[i] = 0.5 * (v0x + v1x);
  }

  /* Now loop over all the active, local generators and compute their faces */
  for (int gen_idx_in_d = d->vertex_start; gen_idx_in_d < d->vertex_end; ++gen_idx_in_d) {

    /* Get the corresponding particle idx */
    int p_idx = d->vertex_part_idx[gen_idx_in_d];

    struct part *p = &parts[p_idx];

    /* Get both lines and neighbouring generators connected to this generator */
    int l_right = d->vertex_line[gen_idx_in_d];
    int n_right_in_d = d->lines[l_right].vertices[1];
    int l_left = d->lines[l_right].neighbours[0];
    int n_left_in_d = d->lines[l_left].vertices[0];

    /* Compute volume, centroid and faces */
    double vor_vert0 = vertices[l_left - 2];
    double vor_vert1 = vertices[l_right - 2];
    double generator = p->x[0];

    p->geometry.pair_connections_offset = v->cell_pair_connections.index;
    voronoi_add_pair(v, d, gen_idx_in_d, n_left_in_d, parts,
                     vor_vert0);
    voronoi_add_pair(v, d, gen_idx_in_d, n_right_in_d, parts,
                     vor_vert1);
    p->geometry.volume = (float)(vor_vert1 - vor_vert0);
    p->geometry.centroid[0] = (float)(0.5 * (vor_vert0 + vor_vert1) - generator);
    p->geometry.centroid[1] = 0.f;
    p->geometry.centroid[2] = 0.f;
    p->geometry.nface = 2;
    p->geometry.min_face_dist =
        (float)fmin(generator - vor_vert0, vor_vert1 - generator);
  }

  /* Be clean */
  free(vertices);
}

inline static void voronoi_reset(struct voronoi *restrict v,
                                 int number_of_cells, double dmin) {
  /* reset indices for the voronoi pairs (faces). */
  for (int i = 0; i < 28; ++i) {
    v->pair_count[i] = 0;
  }

  /* Reset the cell_pair connections */
  int2_lifo_queue_reset(&v->cell_pair_connections);
}

inline static struct voronoi *voronoi_malloc(int number_of_cells, double dmin) {

  struct voronoi *v = malloc(sizeof(struct voronoi));

  /* Allocate memory for the voronoi pairs (faces). */
  for (int i = 0; i < 28; ++i) {
    v->pairs[i] = (struct voronoi_pair *)swift_malloc(
        "voronoi", 10 * sizeof(struct voronoi_pair));
    v->pair_count[i] = 0;
    v->pair_size[i] = 10;
  }

  /* Allocate memory for the cell_pair connections */
  int2_lifo_queue_init(&v->cell_pair_connections, 2 * number_of_cells);

  return v;
}

/**
 * @brief Free up all memory used by the Voronoi grid.
 *
 * @param v Voronoi grid.
 */
static inline void voronoi_destroy(struct voronoi *restrict v) {
  if (v->cell_pair_connections.values != NULL) {
    int2_lifo_queue_destroy(&v->cell_pair_connections);
  }

  for (int i = 0; i < 28; i++) {
    if (v->pairs[i] != NULL) {
      swift_free("voronoi", v->pairs[i]);
      v->pairs[i] = NULL;
      v->pair_count[i] = 0;
      v->pair_size[i] = 0;
    }
  }

  free(v);
}

/**
 * @brief Add a two particle pair to the grid.
 *
 * This function also adds the correct tuple to the cell_pair_connections queue.
 *
 * The grid connectivity is stored per cell sid: sid=0 corresponds to particle
 * pairs encountered during a self task (both particles are within the local
 * cell), while sid=1-26 correspond to particle interactions for which the right
 * neighbour is part of one of the 26 neighbouring cells.
 *
 * For each pair, we compute and store all the quantities required to compute
 * fluxes between the Voronoi cells: the surface area and midpoint of the
 * interface.
 *
 * @param v Voronoi grid.
 * @param d Delaunay tesselation, dual of this Voronoi grid.
 * @param gen_idx_in_d The index in the delaunay tesselation of the left vertex
 * of the new pair.
 * @param ngb_gen_idx_in_d The index in the delaunay tesselation of the right
 * vertex of the new pair.
 * @param part_is_active Array of flags indicating whether the corresponding
 * hydro parts are active (we only construct voronoi cells for active
 * particles).
 * @param ax vertex of the interface.
 * @return 1 if the face was valid (non degenerate surface area), else 0.
 */
static inline int voronoi_add_pair(struct voronoi *v, const struct delaunay *d,
                                   int gen_idx_in_d, int ngb_gen_idx_in_d,
                                   struct part *parts, double ax) {
  int sid;
  int left_part_idx = d->vertex_part_idx[gen_idx_in_d];
  int right_part_idx = d->vertex_part_idx[ngb_gen_idx_in_d];

  /* Pair between local active particles? */
  if (ngb_gen_idx_in_d < d->vertex_end) {
    /* Already processed this pair? */
    if (ngb_gen_idx_in_d < gen_idx_in_d) {
      /* Find the existing pair and add it to the cell_pair_connections. */
      struct part *ngb = &parts[right_part_idx];
      for (int i = 0; i < ngb->geometry.nface; i++) {
        int2 connection =
            v->cell_pair_connections
                .values[ngb->geometry.pair_connections_offset + i];
        if (v->pairs[connection._1][connection._0].right_idx == left_part_idx) {
          int2_lifo_queue_push(&v->cell_pair_connections, connection);
          return 1;
        }
      }
      /* If no pair is found, the face must have been degenerate, nothing left
       * to do. */
      return 0;
    }
    sid = 13;
  } else {
    sid = d->ghost_cell_sids[ngb_gen_idx_in_d - d->vertex_end];
  }

  /* Boundary particle? */
  int actual_sid = sid;
  if (sid & 1 << 5) {
    actual_sid &= ~(1 << 5);
    /* We store all boundary faces under fictive sid 27 */
    sid = 27;
  }

  /* Need to reallocate? */
  if (v->pair_count[sid] == v->pair_size[sid]) {
    v->pair_size[sid] <<= 1;
    v->pairs[sid] = (struct voronoi_pair *)swift_realloc(
        "voronoi", v->pairs[sid],
        v->pair_size[sid] * sizeof(struct voronoi_pair));
  }

  struct voronoi_pair *this_pair = &v->pairs[sid][v->pair_count[sid]];
  this_pair->surface_area = 1.;
  this_pair->midpoint[0] = ax;
  this_pair->midpoint[1] = 0.;
  this_pair->midpoint[2] = 0.;
  this_pair->left_idx = left_part_idx;
  this_pair->right_idx = right_part_idx;
  this_pair->sid = actual_sid;

  /* Add cell_pair_connection */
  int2 connection = {._0 = v->pair_count[sid], ._1 = sid};
  int2_lifo_queue_push(&v->cell_pair_connections, connection);
  ++v->pair_count[sid];
  return 1;
}

/**
 * @brief Write the Voronoi grid information to the given file.
 *
 * The output depends on the configuration. The maximal output contains 3
 * different types of output lines:
 *  - "G\tgx\tgy: x and y position of a single grid generator (optional).
 *  - "C\tcx\tcy\tV\tnface": centroid position, volume and (optionally) number
 *    of faces for a single Voronoi cell.
 *  - "F\tax\tay\tbx\tby\tleft\tngb\tright\tA\tmx\tmy": edge positions
 *    (optional), left and right generator index (and ngb cell index), surface
 *    area and midpoint position for a single two-pair interface.
 *
 * @param v Voronoi grid.
 * @param file File to write to.
 */
static inline void voronoi_write_grid(const struct voronoi *restrict v,
                                      const struct part *parts, int count,
                                      FILE *file, size_t *offset) {}

#endif  // SWIFTSIM_VORONOI_H
