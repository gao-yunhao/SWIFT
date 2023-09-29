#include "bvh.h"

#include <stdlib.h>
#include <string.h>

/**
 * @file Non static function implementations from bvh.h.
 **/

/**
 * @brief Finds a particle from this bvh that contains a given candidate
 * position in its search radius.
 *
 * Particles that are further away than a safety radius `r` from the candidate
 * are discarded and not considered a hit.
 *
 * @param bvh The #flat_bvh to search
 * @param node_id The node to start searching from (initially the root).
 * @param parts The #part stored in this bvh.
 * @param x, y, z The candidate position
 * @param r2 The square of the safety radius.
 * @returns The index of the "hit", i.e.: the particle that contains the
 * candidate position in it's search radius and is itself contained in the
 * safety radius of the candidate position. If no hit is found, -1 is returned.
 */
int flat_bvh_hit_rec(const struct flat_bvh *bvh, int node_id,
                     struct part *parts, double x, double y, double z,
                     double r2) {
  struct flat_bvh_node *node = &bvh->nodes[node_id];

  /* Anything to do here? */
  if (!bbox_contains(&node->bbox, x, y, z)) return -1;

  /* Are we in a leaf? */
  if (node->is_leaf) {

    /* Check individual particles */
    for (int i = 0; i < node->data[BVH_DATA_SIZE]; i++) {
      int pid = node->data[i];
      const struct part *p = &parts[pid];
      double hit_r2 = p->geometry.search_radius * p->geometry.search_radius;
      double dx[3] = {p->x[0] - x, p->x[1] - y, p->x[2] - z};
      double dx2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (dx2 <= r2 && dx2 < hit_r2) return pid;
    }

    /* no hit */
    return -1;
  }

  /* Else, recurse */
  /* Check if hit with central part of this node */
  int pid = node->data[BVH_DATA_SIZE];
  const struct part *p = &parts[pid];
  double hit_r2 = p->geometry.search_radius * p->geometry.search_radius;
  double dx[3] = {p->x[0] - x, p->x[1] - y, p->x[2] - z};
  double dx2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  if (dx2 <= r2 && dx2 < hit_r2) return pid;

  /* Else check subtrees */
  int hit = flat_bvh_hit_rec(bvh, node->children.left, parts, x, y, z, r2);
  if (hit == -1) {
    hit = flat_bvh_hit_rec(bvh, node->children.right, parts, x, y, z, r2);
  }
  return hit;
}

/**
 * @brief Recursively construct a BVH from hydro particles.
 *
 * @param bvh The BVH under construction
 * @param node_id The index of the root of the current subtree under
 * construction.
 * @param parts Array of particles used for construction.
 * @param pid The indices of particles to be added to the current subtree.
 * @param count The length of the `pid` array.
 **/
void flat_bvh_populate_rec(struct flat_bvh *bvh, int node_id,
                           const struct part *parts, int *pid, int count) {

  struct flat_bvh_node *node = &bvh->nodes[node_id];

  /* Construct leaf node? */
  if (count <= BVH_DATA_SIZE) {
    /* Copy data */
    node->data[BVH_DATA_SIZE] = count;
    for (int i = 0; i < count; i++) node->data[i] = pid[i];

    /* Set remaining fields for leaf */
    node->bbox = bbox_from_parts(parts, pid, count);
    node->is_leaf = 1;
    return;
  }

  /* Else, build bvh recursively */

  /* Determine split_direction (direction with the largest spread) */
  double anchor[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
  double opposite[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  for (int i = 0; i < count; i++) {
    const struct part *p = &parts[pid[i]];
    for (int k = 0; k < 3; k++) {
      anchor[k] = fmin(anchor[k], p->x[k]);
      opposite[k] = fmax(opposite[k], p->x[k]);
    }
  }
  double width[3] = {opposite[0] - anchor[0], opposite[1] - anchor[1],
                     opposite[2] - anchor[2]};
  enum direction split_direction;
  if (width[0] >= width[1] && width[0] >= width[2]) {
    split_direction = X_axis;
  } else if (width[1] >= width[0] && width[1] >= width[2]) {
    split_direction = Y_axis;
  } else {
    split_direction = Z_axis;
  }

  /* Now flip particles around, until they are divided into  two groups along
   * the split direction (those below and above the midpoint) */
  double midpoint = 0.5 * (anchor[split_direction] + opposite[split_direction]);
  int head = 0, tail = count - 1;
  while (head < tail) {
    /* Find the first part that belongs in the second group */
    while (parts[pid[head]].x[split_direction] <= midpoint && head < tail) {
      head++;
    }
    /* Find the last part that belongs in the first group */
    while (parts[pid[tail]].x[split_direction] > midpoint && head < tail) {
      tail--;
    }
    if (head < tail) {
      /* Swap the elements at head and tail */
      int tmp = pid[head];
      pid[head] = pid[tail];
      pid[tail] = tmp;
      head++;
      tail--;
    }
  }

  /* Finally find central particle of this node. Take it from the largest
   * subtree */
  int pid_central;
  if (head > (count - head)) {
    /* Search for the right-most particle in the left subtree */
    int idx_max = 0;
    float v_max = parts[pid[0]].x[split_direction];
    for (int i = 1; i < head; i++) {
      float v = parts[pid[i]].x[split_direction];
      if (v > v_max) {
        idx_max = i;
        v_max = v;
      }
    }
    /* Save it and flip the end of this subtrees array in its place */
    pid_central = pid[idx_max];
    pid[idx_max] = pid[head - 1];
    pid[head - 1] = pid_central;
  } else {
    /* Search for the left-most particle in the right subtree */
    int idx_min = head;
    float v_min = parts[pid[head]].x[split_direction];
    for (int i = head + 1; i < count; i++) {
      float v = parts[pid[i]].x[split_direction];
      if (v < v_min) {
        idx_min = i;
        v_min = v;
      }
    }
    /* Save it and flip the end of this subtrees array in its place */
    pid_central = pid[idx_min];
    pid[idx_min] = pid[head];
    pid[head] = pid_central;
    /* Also increase head to indicate that the start of the right subtrees array
     * has shifted */
    head++;
  }

  if (head >= count) error("Failed to partition elements in bvh construction!");

  /* Populate the left and right subtrees */
  int left_id = flat_bvh_new_node(bvh);
  flat_bvh_populate_rec(bvh, left_id, parts, pid, head - 1);
  int right_id = flat_bvh_new_node(bvh);
  flat_bvh_populate_rec(bvh, right_id, parts, &pid[head], count - head);

  /* Nodes might have been moved */
  node = &bvh->nodes[node_id];
  /* Update fields of this node */
  node->children.left = left_id;
  node->children.right = right_id;
  node->data[BVH_DATA_SIZE] = pid_central;
  node->is_leaf = 0;
  bbox_wrap(&bvh->nodes[left_id].bbox, &bvh->nodes[right_id].bbox, &node->bbox);
#ifdef SWIFT_DEBUG_CHECKS
  if (node->is_leaf) {
    assert(node->data[BVH_DATA_SIZE] <= BVH_DATA_SIZE);
  } else {
    assert(node->children.left < bvh->count && 0 <= node->children.left);
    assert(node->children.right < bvh->count && 0 <= node->children.right);
  }
#endif
}

/**
 * @brief Recursively construct a BVH from hydro particles.
 *
 * This method splits nodes by their midpoint, which produces a less balanced
 * BVH, but is faster to compute.
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The current BVH (subtree) under construction
 * @param parts Array of particles used for construction.
 * @param pid The indices of particles to be added to the current BVH.
 * @param count The length of the `pid` array.
 **/
void bvh_populate_rec_midpoint(struct BVH *bvh, const struct part *parts,
                               int *pid, int count) {
  if (count <= 3) {
    /* Set unused fields of this bvh */
    bvh->left = NULL;
    bvh->right = NULL;

    /* Set used fields for leaf */
    size_t size = count * sizeof(int);
    bvh->data = malloc(size);
    memcpy(bvh->data, pid, size);
    bvh->count = count;
    bvh->bbox = bbox_from_parts(parts, pid, count);
    return;
  }

  /* Determine split_direction (direction with the largest spread) */
  double anchor[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
  double opposite[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
  for (int i = 0; i < count; i++) {
    const struct part *p = &parts[pid[i]];
    for (int k = 0; k < 3; k++) {
      anchor[k] = fmin(anchor[k], p->x[k]);
      opposite[k] = fmax(opposite[k], p->x[k]);
    }
  }
  double width[3] = {opposite[0] - anchor[0], opposite[1] - anchor[1],
                     opposite[2] - anchor[2]};
  enum direction split_direction;
  if (width[0] >= width[1] && width[0] >= width[2]) {
    split_direction = X_axis;
  } else if (width[1] >= width[0] && width[1] >= width[2]) {
    split_direction = Y_axis;
  } else {
    split_direction = Z_axis;
  }

  /* Now flip particles around, until they are divided into  two groups along
   * the split direction (those below and above the midpoint) */
  double midpoint = 0.5 * (anchor[split_direction] + opposite[split_direction]);
  int head = 0, tail = count - 1;
  while (head < tail) {
    /* Find the first part that belongs in the second group */
    while (parts[pid[head]].x[split_direction] <= midpoint && head < tail) {
      head++;
    }
    /* Find the last part that belongs in the first group */
    while (parts[pid[tail]].x[split_direction] > midpoint && head < tail) {
      tail--;
    }
    if (head < tail) {
      /* Swap the elements at head and tail */
      int tmp = pid[head];
      pid[head] = pid[tail];
      pid[tail] = tmp;
      head++;
      tail--;
    }
  }

  if (head >= count) error("Failed to partition elements in bvh construction!");

  /* Populate the left and right subtrees */
  bvh->left = malloc(sizeof(struct BVH));
  bvh_populate_rec_midpoint(bvh->left, parts, pid, head);
  bvh->right = malloc(sizeof(struct BVH));
  bvh_populate_rec_midpoint(bvh->right, parts, &pid[head], count - head);

  /* Set the other fields of this bvh */
  bbox_wrap(&bvh->left->bbox, &bvh->right->bbox, &bvh->bbox);
  bvh->data = NULL;
  bvh->count = 0;
}

/**
 * @brief Recursively construct a BVH from hydro particles.
 *
 * This method splits nodes by their median position along the direction of the
 * maximal width. This produces a well balanced BVH, but is slower to compute.
 *
 * This method uses the less efficient BVH struct.
 *
 * @param bvh The current BVH (subtree) under construction
 * @param parts Array of particles used for construction.
 * @param coords Array of 3 arrays containing the x, y and z coordinates of the
 * particles to be added respectively.
 * @param pid The indices of particles to be added to the current BVH.
 * @param count The length of the `pid` array.
 **/
void bvh_populate_rec(struct BVH *bvh, const struct part *parts,
                      double **coords, int *restrict pid, int count) {
  if (count <= 4) {
    /* Set unused fields of this bvh */
    bvh->left = NULL;
    bvh->right = NULL;

    /* Set used fields for leaf */
    size_t size = count * sizeof(int);
    bvh->data = malloc(size);
    memcpy(bvh->data, pid, size);
    bvh->count = count;
    bvh->bbox = bbox_from_parts(parts, pid, count);
    return;
  }

  /* Determine split_direction (direction with the largest spread) */
  double min_x = DBL_MAX;
  double max_x = -DBL_MAX;
  double min_y = DBL_MAX;
  double max_y = -DBL_MAX;
  double min_z = DBL_MAX;
  double max_z = -DBL_MAX;

  for (int i = 0; i < count; i++) {
    const struct part *p = &parts[pid[i]];
    min_x = min(min_x, p->x[0]);
    max_x = max(max_x, p->x[0]);
    min_y = min(min_y, p->x[1]);
    max_y = max(max_y, p->x[1]);
    min_z = min(min_z, p->x[2]);
    max_z = max(max_z, p->x[2]);
  }

  enum direction split_direction;
  if (max_x - min_x >= max_y - min_y && max_x - min_x >= max_z - min_z) {
    split_direction = X_axis;
  } else if (max_y - min_y >= max_x - min_x && max_y - min_y >= max_z - min_z) {
    split_direction = Y_axis;
  } else {
    split_direction = Z_axis;
  }

  /* Sort particles along splitting direction and apply median splitting */
  qsort_r(pid, count, sizeof(int), &cmp, coords[split_direction]);
  int median_idx = count / 2 + 1;

  /* Populate the left and right subtree of this bvh */
  bvh->left = malloc(sizeof(struct BVH));
  bvh_populate_rec(bvh->left, parts, coords, pid, median_idx);
  bvh->right = malloc(sizeof(struct BVH));
  bvh_populate_rec(bvh->right, parts, coords, &pid[median_idx],
                   count - median_idx);

  /* Set the other fields of this bvh */
  bbox_wrap(&bvh->left->bbox, &bvh->right->bbox, &bvh->bbox);
  bvh->data = NULL;
  bvh->count = 0;
}
