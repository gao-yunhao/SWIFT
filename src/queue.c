/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "queue.h"

/* Local headers. */
#include "atomic.h"
#include "error.h"
#include "memswap.h"

/**
 * @brief Push the task at the given index up the heap until it is either at the
 * top or smaller than its parent.
 *
 * @param q The task #queue.
 * @param ind The index of the task to be sifted-down in the queue.
 *
 * @return The new index of the entry.
 */
int queue_bubble_up(struct queue *q, int ind) {
  /* Set some pointers we will use often. */
  struct queue_entry *entries = q->entries;
  const float w = entries[ind].weight;

  /* While we are not yet at the top of the heap... */
  while (ind > 0) {
    /* Check if the parent is larger and bail if not.. */
    const int parent = (ind - 1) / 2;
    if (w < entries[parent].weight) break;

    /* Parent is not larger, so swap. */
    memswap(&entries[ind], &entries[parent], sizeof(struct queue_entry));
    ind = parent;
  }

  return ind;
} //  对权重进行堆排序？返回值是插入堆结束后元素的位置；

/**
 * @brief Push the task at the given index down the heap until both its children
 * have a smaller weight.
 *
 * @param q The task #queue.
 * @param ind The index of the task to be sifted-down in the queue.
 *
 * @return The new index of the entry.
 */
int queue_sift_down(struct queue *q, int ind) {
  /* Set some pointers we will use often. */
  struct queue_entry *entries = q->entries;
  const int qcount = q->count;
  const float w = entries[ind].weight;

  /* While we still have at least one child... */
  while (1) {
    /* Check if we still have children. */
    int child = 2 * ind + 1;
    if (child >= qcount) break;

    /* Which of both children is the largest? */
    if (child + 1 < qcount && entries[child + 1].weight > entries[child].weight)
      child += 1;

    /* Do we want to swap with the largest child? */
    if (entries[child].weight > w) {
      memswap(&entries[ind], &entries[child], sizeof(struct queue_entry));
      ind = child;
    } else
      break;
  }

  return ind;
} //  上面那个函数的逆过程

/**
 * @brief Enqueue all tasks in the incoming DEQ.
 *
 * @param q The #queue, assumed to be locked.
 */
void queue_get_incoming(struct queue *q) {

  struct queue_entry *entries = q->entries;

  /* Loop over the incoming DEQ. */
  while (1) {

    /* Is there a next element? */
    const int ind = q->first_incoming % queue_incoming_size;
    if (q->tid_incoming[ind] < 0) break;

    /* Get the next offset off the DEQ. */
    const int offset = atomic_swap(&q->tid_incoming[ind], -1); //  将第二个参数赋值给第一个参数，并输出第一个参数原来的值；
    atomic_inc(&q->first_incoming); //  将参数加一（原子操作）

    /* Does the queue need to be grown? */
    if (q->count == q->size) {
      struct queue_entry *temp;
      q->size *= queue_sizegrow;
      if ((temp = (struct queue_entry *)malloc(sizeof(struct queue_entry) *
                                               q->size)) == NULL)
        error("Failed to allocate new indices.");
      memcpy(temp, entries, sizeof(struct queue_entry) * q->count);
      free(entries);
      q->entries = entries = temp;
    }

    /* Drop the task at the end of the queue. */
    entries[q->count].tid = offset;
    entries[q->count].weight = q->tasks[offset].weight;
    q->count += 1;
    atomic_dec(&q->count_incoming); //  将参数减一（原子操作）

    /* Re-heap by bubbling up the new (last) element. */
    queue_bubble_up(q, q->count - 1);

#ifdef SWIFT_DEBUG_CHECK
    /* Check the queue's consistency. */
    for (int k = 1; k < q->count; k++)
      if (entries[(k - 1) / 2].weight < entries[k].weight)
        error("Queue heap is disordered.");
#endif
  }
} //  将incoming的任务id和权重全部入队（实际上是放入堆中）

/**
 * @brief Insert a used tasks into the given queue.
 *
 * @param q The #queue.
 * @param t The #task.
 */
void queue_insert(struct queue *q, struct task *t) {
  /* Get an index in the DEQ. */
  const int ind = atomic_inc(&q->last_incoming) % queue_incoming_size; //  如果incoming队列没有装满，就直接从last_incoming的下一个位置开始插入；

  /* Spin until the new offset can be stored. */
  while (atomic_cas(&q->tid_incoming[ind], -1, t - q->tasks) != -1) { //  如果前两个参数相等，就把第三个参数的值赋给第一个参数；始终返回第一个参数的旧值；
                                                                      //  t - q->tasks：指针相减，返回这两个指针之间元素个数；
    /* Try to get the queue lock, non-blocking, ensures that at
       least somebody is working on this queue. */
    if (lock_trylock(&q->lock) == 0) { //  等待其他线程处理完他们（在这个队列上）的任务再执行下面的代码，等待过程会继续运行这个循环；

      /* Clean up the incoming DEQ. */
      queue_get_incoming(q);           //  等待结束，清空incoming队列，清空后从此前的last_incoming的下一个位置开始插入；

      /* Release the queue lock. */
      if (lock_unlock(&q->lock) != 0) {//  我猜是释放锁，让其他（可能存在的）被锁住的线程开始运行；
        error("Unlocking the qlock failed.\n");
      }
    }
  }

  /* Increase the incoming count. */
  atomic_inc(&q->count_incoming);
} //  将一个任务插入到incoming的缓冲区中；

/**
 * @brief Initialize the given queue.
 *
 * @param q The #queue.
 * @param tasks List of tasks to which the queue indices refer to.
 */
void queue_init(struct queue *q, struct task *tasks) {

  /* Allocate the task list if needed. */
  q->size = queue_sizeinit;
  if ((q->entries = (struct queue_entry *)malloc(sizeof(struct queue_entry) *
                                                 q->size)) == NULL)
    error("Failed to allocate queue entries.");

  /* Set the tasks pointer. */
  q->tasks = tasks;

  /* Init counters. */
  q->count = 0;

  /* Init the queue lock. */
  if (lock_init(&q->lock) != 0) error("Failed to init queue lock.");

  /* Init the incoming DEQ. */
  if ((q->tid_incoming = (int *)malloc(sizeof(int) * queue_incoming_size)) ==
      NULL)
    error("Failed to allocate queue incoming buffer.");
  for (int k = 0; k < queue_incoming_size; k++) {
    q->tid_incoming[k] = -1;
  }
  q->first_incoming = 0;
  q->last_incoming = 0;
  q->count_incoming = 0;
} //  使用这个函数时的语句：for (int k = 0; k < nr_queues; k++) queue_init(&s->queues[k], NULL); 也就是说初始化时存储的task是空指针；

/**
 * @brief Get a task free of dependencies and conflicts.
 *
 * @param q The task #queue.
 * @param prev The previous #task extracted from this #queue.
 * @param blocking Block until access to the queue is granted.
 */
struct task *queue_gettask(struct queue *q, const struct task *prev, //  prev参数似乎没有用到？
                           int blocking) {

  swift_lock_type *qlock = &q->lock;
  struct task *res = NULL;

  /* Grab the task lock. */
  if (blocking) {
    if (lock_lock(qlock) != 0) error("Locking the qlock failed.\n");
  } else {
    if (lock_trylock(qlock) != 0) return NULL; //  如果当前队列被其他线程占用，那么直接查找任务失败？
  }

  /* Fill any tasks from the incoming DEQ. */
  queue_get_incoming(q);

  /* If there are no tasks, leave immediately. */
  if (q->count == 0) {
    lock_unlock_blind(qlock);
    return NULL;
  }

  /* Set some pointers we will use often. */
  struct queue_entry *entries = q->entries;
  struct task *qtasks = q->tasks;
  const int old_qcount = q->count;

  /* Loop over the queue entries. */
  int ind;
  for (ind = 0; ind < old_qcount; ind++) {

    /* Try to lock the next task. */
    if (task_lock(&qtasks[entries[ind].tid])) break; //  判断当前cell是否具备执行各种task的条件；

    /* Should we de-prioritize this task? */

    // MATTHIEU: We now have more than 64 tasks so the bit-wise
    // operation here is problematic.
    // However, the mask was such that this condition is always true anyway.
    if (1 /* (1ULL << qtasks[entries[ind].tid].type) & */
        /* queue_lock_fail_reweight_mask */) {

      /* Scale the task's weight. */
      entries[ind].weight *= queue_lock_fail_reweight_factor;

      /* Send it down the binary heap. */
      if (queue_sift_down(q, ind) != ind) ind -= 1;
    }
  } //  我怀疑这个循环有可能会变成死循环，比如初始的时候task的权重，1，0.8，此后交替×0.5并且ind-=1，那么如果260行一直没有break，那么ind就会一直等于1，这个循环就会一直执行；
  /* Did we get a task? */ //  我猜设计代码的时候可能认为当sift down几次之后队列元素的位置就不会发生变化了，这样ind就不会减一，循环就可以正常运行；
  if (ind < old_qcount) {

    /* Another one bites the dust. */
    const int qcount = q->count -= 1;

    /* Get a pointer on the task that we want to return. */
    res = &qtasks[entries[ind].tid];

    /* Swap this task with the last task and re-heap. */
    if (ind < qcount) { //  将选出来的task出队；
      entries[ind] = entries[qcount];
      ind = queue_bubble_up(q, ind);
      ind = queue_sift_down(q, ind);
    }

  } else
    res = NULL;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the queue's consistency. */
  for (int k = 1; k < q->count; k++)
    if (entries[(k - 1) / 2].weight < entries[k].weight)
      error("Queue heap is disordered.");
#endif

  /* Release the task lock. */
  if (lock_unlock(qlock) != 0) error("Unlocking the qlock failed.\n");

  /* Take the money and run. */
  return res;
} //  从这个队列中按优先级取出一个任务；

void queue_clean(struct queue *q) {

  free(q->entries);
  free(q->tid_incoming);
}

/**
 * @brief Dump a formatted list of tasks in the queue to the given file stream.
 *
 * @param nodeID the node id of this rank.
 * @param index a number for this queue, added to the output.
 * @param file the FILE stream, should opened for write.
 * @param q The task #queue.
 */
void queue_dump(int nodeID, int index, FILE *file, struct queue *q) {

  swift_lock_type *qlock = &q->lock;

  /* Grab the queue lock. */
  if (lock_lock(qlock) != 0) error("Locking the qlock failed.\n");

  /* Fill any tasks from the incoming DEQ. */
  queue_get_incoming(q);

  /* Loop over the queue entries. */
  for (int k = 0; k < q->count; k++) {
    struct task *t = &q->tasks[q->entries[k].tid];

    fprintf(file, "%d %d %d %s %s %.2f\n", nodeID, index, k,
            taskID_names[t->type], subtaskID_names[t->subtype], t->weight);
  }

  /* Release the task lock. */
  if (lock_unlock(qlock) != 0) error("Unlocking the qlock failed.\n");
} //  输出任务名称到某个文件中？
