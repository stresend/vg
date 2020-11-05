/**
 * \file vectorized_banded_aligner.hpp
 *
 * Defines a vectorized banded global alignment algorithm
 *
 */
#ifndef VG_VECTORIZED_BANDED_ALIGNER_HPP_INCLUDED
#define VG_VECTORIZED_BANDED_ALIGNER_HPP_INCLUDED

#include<limits>
#include<cstdint>
#include<immintrin.h>
#include<algorithm>

namespace vg {

using namespace std;

struct SWVector {
    __m128i match;
    __m128i insert_row;
    __m128i insert_col;
};

struct BandedVectorMatrix {
    
    // returns whether a coordinate in the full matrix is inside
    // the band
    bool is_inside_band(int i, int j);
    
    // returns the index of the vector that contains a coordinate
    // from the full matrix
    int vector_index(int i, int j);
    
    // returns the index within the vector conaining these
    // coordinates that corresponds to this exact coordinate
    int index_within_vector(int i, int j);
    
    // returns coordinates in the full matrix that correspond
    // to an entry in a given vector
    pair<int, int> coordinate(int vec_idx, int idx_within_vec);
    
    // returns the index of the vector that corresponds to the
    // entries directly above the indicated vector in the full matrix
    int vector_above(int vec_idx);
    
    // returns the index of the vector that corresponds to the
    // entries directly below the indicated vector in the full matrix
    int vector_below(int vec_idx);
    
    // returns the index of the vector that corresponds to the
    // entries up and to the left of the indicated vector in the
    // full matrix
    int vector_diagonal(int vec_idx);
   
    // returns number of vectors in band
    int get_band_size();
    
    SWVector* vectors;
    int64_t first_diag;
    int64_t num_diags;
    int64_t num_cols;
};

void init_BandedVectorMatrix(BandedVectorMatrix& matrix, int64_t first_diag,
                             int64_t num_diags, int64_t num_cols);

}

#endif
