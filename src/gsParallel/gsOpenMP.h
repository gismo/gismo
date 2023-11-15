/** @file gsOpenMP.h

    @brief OpenMP stub routines to be used when omp.h is not available

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#pragma once

#ifdef _OPENMP

#if _OPENMP >= 202111
#define GISMO_HAS_OPENMP_52 1
#else
#define GISMO_HAS_OPENMP_52 0
#endif

#if _OPENMP >= 202011
#define GISMO_HAS_OPENMP_51 1
#else
#define GISMO_HAS_OPENMP_51 0
#endif

#if _OPENMP >= 201811
#define GISMO_HAS_OPENMP_50 1
#else
#define GISMO_HAS_OPENMP_50 0
#endif

#if _OPENMP >= 201511
#define GISMO_HAS_OPENMP_45 1
#else
#define GISMO_HAS_OPENMP_45 0
#endif

#if _OPENMP >= 201307
#define GISMO_HAS_OPENMP_40 1
#else
#define GISMO_HAS_OPENMP_40 0
#endif

#if _OPENMP >= 201107
#define GISMO_HAS_OPENMP_31 1
#else
#define GISMO_HAS_OPENMP_31 0
#endif

#if _OPENMP >= 200805
#define GISMO_HAS_OPENMP_30 1
#else
#define GISMO_HAS_OPENMP_30 0
#endif

#if _OPENMP >= 200505
#define GISMO_HAS_OPENMP_25 1
#else
#define GISMO_HAS_OPENMP_25 0
#endif

#include <omp.h>

#else

#define GISMO_HAS_OPENMP_52 0
#define GISMO_HAS_OPENMP_51 0
#define GISMO_HAS_OPENMP_50 0
#define GISMO_HAS_OPENMP_45 0
#define GISMO_HAS_OPENMP_40 0
#define GISMO_HAS_OPENMP_31 0
#define GISMO_HAS_OPENMP_30 0
#define GISMO_HAS_OPENMP_25 0

#include <gsCore/gsForwardDeclarations.h>

void GISMO_EXPORT omp_set_num_threads(int num_threads);

int GISMO_EXPORT omp_get_num_threads(void);

int GISMO_EXPORT omp_get_max_threads(void);

int GISMO_EXPORT omp_get_thread_num(void);

int GISMO_EXPORT omp_get_num_procs(void);

int GISMO_EXPORT omp_in_parallel(void);

void GISMO_EXPORT omp_set_dynamic(int dynamic_threads);

int GISMO_EXPORT omp_get_dynamic(void);

int GISMO_EXPORT omp_get_cancellation(void);

void GISMO_EXPORT omp_set_nested(int nested);

int GISMO_EXPORT omp_get_nested(void);

typedef enum omp_sched_t {
    omp_sched_static  = 1,
    omp_sched_dynamic = 2,
    omp_sched_guided  = 3,
    omp_sched_auto    = 4,
    omp_sched_monotonic = 0x80000000
} omp_sched_t;

void GISMO_EXPORT omp_set_schedule(omp_sched_t kind, int chunk_size);

void GISMO_EXPORT omp_get_schedule(omp_sched_t *kind, int *chunk_size);

int GISMO_EXPORT omp_get_thread_limit(void);

void GISMO_EXPORT omp_set_max_active_levels(int max_active_levels);

int GISMO_EXPORT omp_get_max_active_levels(void);

int GISMO_EXPORT omp_get_level(void);

int GISMO_EXPORT omp_get_ancestor_thread_num(int level);

int GISMO_EXPORT omp_get_team_size(int level);

int GISMO_EXPORT omp_get_active_level(void);

int GISMO_EXPORT omp_in_final(void);

typedef enum omp_proc_bind_t {
    omp_proc_bind_false = 0,
    omp_proc_bind_true = 1,
    omp_proc_bind_master = 2,
    omp_proc_bind_close = 3,
    omp_proc_bind_spread = 4
} omp_proc_bind_t;

omp_proc_bind_t omp_get_proc_bind(void);

int GISMO_EXPORT omp_get_num_places(void);

int GISMO_EXPORT omp_get_place_num_procs(int place_num);

void GISMO_EXPORT omp_get_place_proc_ids(int place_num, int *ids);

int GISMO_EXPORT omp_get_place_num(void);

int GISMO_EXPORT omp_get_partition_num_places(void);

void GISMO_EXPORT omp_get_partition_place_nums(int *place_nums);

void GISMO_EXPORT omp_set_default_device(int device_num);

int GISMO_EXPORT omp_get_default_device(void);

int GISMO_EXPORT omp_get_num_devices(void);

int GISMO_EXPORT omp_get_num_teams(void);

int GISMO_EXPORT omp_get_team_num(void);

int GISMO_EXPORT omp_is_initial_device(void);

int GISMO_EXPORT omp_get_initial_device(void);

int GISMO_EXPORT omp_get_max_task_priority(void);

typedef struct omp_lock_t {
    int lock;
} omp_lock_t;

enum { OMP_UNLOCKED = -1, OMP_INIT, OMP_LOCKED };

void GISMO_EXPORT omp_init_lock(omp_lock_t *arg);

typedef enum omp_sync_hint_t {
    omp_sync_hint_none           = 0,
    omp_lock_hint_none           = omp_sync_hint_none,
    omp_sync_hint_uncontended    = 1,
    omp_lock_hint_uncontended    = omp_sync_hint_uncontended,
    omp_sync_hint_contended      = (1<<1),
    omp_lock_hint_contended      = omp_sync_hint_contended,
    omp_sync_hint_nonspeculative = (1<<2),
    omp_lock_hint_nonspeculative = omp_sync_hint_nonspeculative,
    omp_sync_hint_speculative    = (1<<3),
    omp_lock_hint_speculative    = omp_sync_hint_speculative,
    kmp_lock_hint_hle            = (1<<16),
    kmp_lock_hint_rtm            = (1<<17),
    kmp_lock_hint_adaptive       = (1<<18)
} omp_sync_hint_t;

typedef omp_sync_hint_t omp_lock_hint_t;

void GISMO_EXPORT omp_init_lock_with_hint(omp_lock_t *arg, omp_lock_hint_t hint);

void GISMO_EXPORT omp_destroy_lock(omp_lock_t *arg);

void GISMO_EXPORT omp_set_lock(omp_lock_t *arg);

void GISMO_EXPORT omp_unset_lock(omp_lock_t *arg);

int GISMO_EXPORT omp_test_lock(omp_lock_t *arg);

typedef struct omp_nest_lock_t {
    int owner;
    int count;
} omp_nest_lock_t;

enum { OMP_NOOWNER = -1, OMP_MASTER = 0 };

void GISMO_EXPORT omp_init_nest_lock(omp_nest_lock_t *arg);

void GISMO_EXPORT omp_init_nest_lock_with_hint(omp_nest_lock_t *arg,
                                               omp_lock_hint_t hint);

void GISMO_EXPORT omp_destroy_nest_lock(omp_nest_lock_t *arg);

void GISMO_EXPORT omp_set_nest_lock(omp_nest_lock_t *arg);

void GISMO_EXPORT omp_unset_nest_lock(omp_nest_lock_t *arg);

int GISMO_EXPORT omp_test_nest_lock(omp_nest_lock_t *arg);

double omp_get_wtime(void);

double omp_get_wtick(void);

void * omp_target_alloc(size_t size, int device_num);

void GISMO_EXPORT omp_target_free(void *device_ptr, int device_num);

int GISMO_EXPORT omp_target_is_present(void *ptr, int device_num);

int GISMO_EXPORT omp_target_memcpy(void *dst, void *src, size_t length,
                                   size_t dst_offset, size_t src_offset,
                                   int dst_device, int src_device);

int GISMO_EXPORT omp_target_memcpy_rect(void *dst, void *src,
                                        size_t element_size,
                                        int num_dims,
                                        const size_t *volume,
                                        const size_t *dst_offsets,
                                        const size_t *src_offsets,
                                        const size_t *dst_dimensions,
                                        const size_t *src_dimensions,
                                        int dst_device_num, int src_device_num);

int GISMO_EXPORT omp_target_associate_ptr(void *host_ptr, void *device_ptr,
                                          size_t size, size_t device_offset,
                                          int device_num);

int GISMO_EXPORT omp_target_disassociate_ptr(void *ptr, int device_num);
#endif // _OPENMP
