/** @file gsOpenMP.cpp

    @brief Implementation of OpenMP stub routines to be used when libomp is not available

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#if !defined(_OPENMP)

#include <gsCore/gsExport.h>
#include <gsParallel/gsOpenMP.h>

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void omp_set_num_threads(int num_threads)
{}

int  omp_get_num_threads(void)
{
    return 1;
}

int  omp_get_max_threads(void)
{
    return 1;
}

int  omp_get_thread_num(void)
{
    return 0;
}

int  omp_get_num_procs(void)
{
    return 1;
}

int  omp_in_parallel(void)
{
    return 0;
}

void omp_set_dynamic(int dynamic_threads)
{}

int  omp_get_dynamic(void)
{
    return 0;
}

int  omp_get_cancellation(void)
{
    return 0;
}

void omp_set_nested(int nested)
{}

int  omp_get_nested(void)
{
    return 0;
}

void omp_set_schedule(omp_sched_t kind, int chunk_size)
{}

void omp_get_schedule(omp_sched_t *kind, int *chunk_size)
{
    *kind = omp_sched_static;
    *chunk_size = 0;
}

int omp_get_thread_limit(void)
{
    return 1;
}

void omp_set_max_active_levels(int max_active_levels)
{}

int omp_get_max_active_levels(void)
{
    return 0;
}

int omp_get_level(void)
{
    return 0;
}

int omp_get_ancestor_thread_num(int level)
{
    return level == 0 ? 0 : -1;
}

int omp_get_team_size(int level)
{
    return level == 0 ? 1 : -1;
}

int omp_get_active_level(void)
{
    return 0;
}

int omp_in_final(void)
{
    return 1;
}

omp_proc_bind_t omp_get_proc_bind(void)
{
    return omp_proc_bind_false;
}

int omp_get_num_places(void)
{
    return 0;
}

int omp_get_place_num_procs(int place_num)
{
    return 0;
}

void  omp_get_place_proc_ids(int place_num, int *ids)
{}

int omp_get_place_num(void)
{
    return -1;
}

int omp_get_partition_num_places(void)
{
    return 0;
}

void omp_get_partition_place_nums(int *place_nums)
{}

void omp_set_default_device(int device_num)
{}

int omp_get_default_device(void)
{
    return 0;
}

int omp_get_num_devices(void)
{
    return 0;
}

int omp_get_num_teams(void)
{
    return 1;
}

int omp_get_team_num(void)
{
    return 0;
}

int omp_is_initial_device(void)
{
    return 1;
}

int omp_get_initial_device(void)
{
    return -10;
}

int omp_get_max_task_priority(void)
{
    return 0;
}

void omp_init_lock(omp_lock_t *arg)
{
    arg->lock = OMP_UNLOCKED;
}

void omp_init_lock_with_hint(omp_lock_t *arg, omp_lock_hint_t hint)
{
    omp_init_lock(arg);
}

void omp_destroy_lock(omp_lock_t *arg)
{
    arg->lock = OMP_INIT;
}

void omp_set_lock(omp_lock_t *arg)
{
    if (arg->lock == OMP_UNLOCKED)
    {
        arg->lock = OMP_LOCKED;
    }
    else if (arg->lock == OMP_LOCKED)
    {
        fprintf(stderr, "error: deadlock in using lock variable\n");
        exit(1);
    }
    else
    {
        exit(1);
    }
}

void omp_unset_lock(omp_lock_t *arg)
{
    if (arg->lock == OMP_LOCKED)
    {
        arg->lock = OMP_UNLOCKED;
    }
    else if (arg->lock == OMP_UNLOCKED)
    {
        fprintf(stderr, "error: lock not set\n");
        exit(1);
    }
    else
    {
        fprintf(stderr, "error: lock not initialized\n");
        exit(1);
    }
}

int omp_test_lock(omp_lock_t *arg)
{
    if (arg->lock == OMP_UNLOCKED)
    {
        arg->lock = OMP_LOCKED;
        return 1;
    }
    else if (arg->lock == OMP_LOCKED)
    {
        return 0;
    }
    else {
        fprintf(stderr, "error: lock not initialized\n");
        exit(1);
    }
}

void omp_init_nest_lock(omp_nest_lock_t *arg)
{
    arg->owner = OMP_NOOWNER;
    arg->count = 0;
}

void omp_init_nest_lock_with_hint(omp_nest_lock_t *arg,
                                  omp_lock_hint_t hint)
{
    omp_init_nest_lock(arg);
}

void omp_destroy_nest_lock(omp_nest_lock_t *arg)
{
    arg->owner = OMP_NOOWNER;
    arg->count = OMP_UNLOCKED;
}

void omp_set_nest_lock(omp_nest_lock_t *arg)
{
    if (arg->owner == OMP_MASTER && arg->count >= 1)
    {
        arg->count++;
    }
    else if (arg->owner == OMP_NOOWNER && arg->count == 0)
    {
        arg->owner = OMP_MASTER;
        arg->count = 1;
    }
    else
    {
        fprintf(stderr, "error: lock corrupted or not initialized\n");
        exit(1);
    }
}

void omp_unset_nest_lock(omp_nest_lock_t *arg)
{
    if (arg->owner == OMP_MASTER && arg->count >= 1)
    {
        arg->count--;
        if (arg->count == 0)
        {
            arg->owner = OMP_NOOWNER;
        }
    }
    else if (arg->owner == OMP_NOOWNER && arg->count == 0)
    {
        fprintf(stderr, "error: lock not set\n");
        exit(1);
    }
    else
    {
        fprintf(stderr, "error: lock corrupted or not initialized\n");
        exit(1);
    }
}

int omp_test_nest_lock(omp_nest_lock_t *arg)
{
    omp_set_nest_lock(arg);
    return arg->count;
}

double omp_get_wtime(void)
{
    /* This function does not provide a working
     * wallclock timer. Replace it with a version
     * customized for the target machine.
     */
    return 0.0;
}

double omp_get_wtick(void)
{
    /* This function does not provide a working
     * clock tick function. Replace it with
     * a version customized for the target machine.
     */
    return 365. * 86400.;
}

void * omp_target_alloc(size_t size, int device_num)
{
    if (device_num != -10)
        return NULL;
    return malloc(size);
}

void omp_target_free(void *device_ptr, int device_num)
{
    free(device_ptr);
}

int omp_target_is_present(void *ptr, int device_num)
{
    return 1;
}

int omp_target_memcpy(void *dst, void *src, size_t length,
                      size_t dst_offset, size_t src_offset,
                      int dst_device, int src_device)
{
    // only the default device is valid in a stub
    if (dst_device != -10 || src_device != -10
        || ! dst || ! src )
        return EINVAL;
    memcpy((char *)dst + dst_offset,
           (char *)src + src_offset,
           length);
    return 0;
}

int omp_target_memcpy_rect(void *dst, void *src,
                           size_t element_size,
                           int num_dims,
                           const size_t *volume,
                           const size_t *dst_offsets,
                           const size_t *src_offsets,
                           const size_t *dst_dimensions,
                           const size_t *src_dimensions,
                           int dst_device_num, int src_device_num)
{
    int ret=0;
    // Both null, return number of dimensions supported,
    // this stub supports an arbitrary number
    if (dst == NULL && src == NULL) return INT_MAX;
  
    if (!volume || !dst_offsets || !src_offsets
        || !dst_dimensions || !src_dimensions
        || num_dims < 1 ) {
        ret = EINVAL;
        goto done;
    }
    if (num_dims == 1) {
        ret = omp_target_memcpy(dst, src,
                                element_size * volume[0],
                                dst_offsets[0] * element_size,
                                src_offsets[0] * element_size,
                                dst_device_num, src_device_num);
        if(ret) goto done;
    } else {
        size_t dst_slice_size = element_size;
        size_t src_slice_size = element_size;
        for (int i=1; i < num_dims; i++) {
            dst_slice_size *= dst_dimensions[i];
            src_slice_size *= src_dimensions[i];
        }
        size_t dst_off = dst_offsets[0] * dst_slice_size;
        size_t src_off = src_offsets[0] * src_slice_size;
        for (size_t i=0; i < volume[0]; i++) {
            ret = omp_target_memcpy_rect(
                (char *)dst + dst_off + dst_slice_size*i,
                (char *)src + src_off + src_slice_size*i,
                element_size,
                num_dims - 1,
                volume + 1,
                dst_offsets + 1,
                src_offsets + 1,
                dst_dimensions + 1,
                src_dimensions + 1,
                dst_device_num,
                src_device_num);
            if (ret) goto done;
        }
    }
done:
    return ret;
}

int omp_target_associate_ptr(void *host_ptr, void *device_ptr,
                             size_t size, size_t device_offset,
                             int device_num)
{
    // No association is possible because all host pointers
    // are considered present
    return EINVAL;
}

int omp_target_disassociate_ptr(void *ptr, int device_num)
{
    return EINVAL;
}
#endif // !defined(_OPENMP)
