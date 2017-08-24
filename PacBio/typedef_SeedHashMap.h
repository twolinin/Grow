#ifndef TYPEDEF_SEEDHASHMAP_H
#define TYPEDEF_SEEDHASHMAP_H

#include "tbb/concurrent_hash_map.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "PacBioScaffoldProcess.h"


typedef tbb::concurrent_hash_map<int64_t, SeedSequenceInfo > SeedHashMap;