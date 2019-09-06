#ifndef LEVELS_CXX_HPP
#define LEVELS_CXX_HPP

#include <stddef.h>
#include "flcl-cxx.hpp"

typedef struct _levels_t {
    int32_t n;
    int32_t ntop;
    int32_t chunk_count;
    int32_t *ltop;
    int32_t *chunk_start;
    int32_t *chunk_end;
    int32_t *chunk_ids;
} levels_t;

typedef struct _c_levels_t {
    int32_t n;
    flcl_ndarray_t nd_ltop;
    flcl_ndarray_t nd_chunk_start;
    flcl_ndarray_t nd_chunk_end;
    flcl_ndarray_t nd_chunk_ids;
} c_levels_t;

#endif // LEVELS_CXX_HPP