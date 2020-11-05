/**
 * \file vectorized_banded_aligner.cpp
 *
 * Implements a the vectorized banded global aligner
 *
 */

#define init_BandedVectorMatrix_debug

#include "vectorized_banded_aligner.hpp"
#include<iostream>
using namespace std;

//namespace vg {
    /* Note on vector indexing:
     *    0 index reserved for a "buffer" vector, calls to vector_above with vec_idx < num_cols returns this index
     *    1->v indices are vectors that are filled during dynamic programming
     *    v+1 index is another "buffer" vector, call to vect_below on vec_idx > v-num_cols returns this index
     */ 
    bool BandedVectorMatrix::is_inside_band(int i, int j){
        return i >= first_diag+j && i < first_diag+j+num_diags;
    }

    int BandedVectorMatrix::vector_index(int i, int j){
        return index_within_vector(i-j-first_diag, j);
    }

    int BandedVectorMatrix::index_within_vector(int i, int j){
        return ((i+1)/8)*(num_cols+1)+j%num_cols+1;
    }

    pair<int, int> BandedVectorMatrix::coordinate(int vec_idx, int idx_within_vec){
        int matrix_j = vec_idx%num_cols;
        int matrix_i =  idx_within_vec +first_diag - 1 + matrix_j + 8*(vec_idx/num_cols);
        return pair<int, int>(matrix_i, matrix_j);
    }


    int BandedVectorMatrix::vector_above(int vec_idx){
        return max((int64_t)0, vec_idx - num_cols);
    }

    int BandedVectorMatrix::vector_below(int vec_idx){
        return min(num_cols+(num_cols+1)*((num_diags-1)/8)+1, vec_idx+num_cols+((vec_idx)/num_cols));
    }

    int BandedVectorMatrix::vector_diagonal(int vec_idx){
        return vec_idx -1;
    }

    int BandedVectorMatrix::get_band_size(){
        return num_cols+(num_cols+1)*((num_diags-1)/8)+1;
    }
    /* Notes on vector initialization:
     *    Here's an example of a banded graph turned vectors(let 1's be the buffer area, 0's be the dp region)
     * 0001111111        111111111111
     * 0000111111        111000000001
     * 0000011111        110000000001
     * 1000001111        100000000001
     * 1100000111        100000000001
     * 1110000011   ->   100000000001
     * 1111000001        100000000001
     * 1111100000        
     * 1111110000        
     * 1111111000
     * 
     * Notice: this implemetation adds a single row of buffer to the top of the vectorized band.
     * also notice that an extra buffer column is inserted before the first column 
     * As well as another buffer column at the end of the matrix, this column will act as the vectors below
     */
    void init_BandedVectorMatrix(BandedVectorMatrix& matrix, int64_t first_diag,
                                 int64_t num_diags, int64_t num_cols){
        matrix.num_cols = num_cols;
        matrix.first_diag = first_diag;
        matrix.num_diags = num_diags;

        void* banded_ad;
        if(posix_memalign(&banded_ad, 16, sizeof(SWVector)*matrix.get_band_size())){
            perror("posix_memalign fail for BandedVectorMatrix");
        }
        matrix.vectors = (SWVector*) banded_ad;

        short init_mask[8];
        _mm_store_si128((__m128i*)init_mask, _mm_setzero_si128());
        for(int idx = 0; idx < -first_diag +1 && idx < 8; idx++){
            init_mask[idx] = 1;
        }
# ifdef init_BandedVectorMatrix_debug
        cout << "init_mask before iterations: ";
        for(int i = 0; i < 8; i++){
            cout << init_mask[i] << " ";
        }
        cout << endl;
#endif

        // index zero is all buffer
        matrix.vectors[0].match = _mm_set1_epi16(1);
        matrix.vectors[0].insert_row = _mm_set1_epi16(1);
        matrix.vectors[0].insert_col = _mm_set1_epi16(1);

        // initialize top left corner of rectangularized vector matrix
        for(int idx =1; idx < -first_diag + 1; idx++){
            matrix.vectors[idx].match = _mm_load_si128((__m128i*)init_mask);
            matrix.vectors[idx].insert_row = _mm_load_si128((__m128i*)init_mask);
            matrix.vectors[idx].insert_col = _mm_load_si128((__m128i*)init_mask);
            init_mask[-first_diag-idx+1] = 0;
#ifdef init_BandedVectorMatrix_debug
            cout << "init_mask during iteration " << idx << ": ";
            for(int i = 0; i < 8; i++){
                cout << init_mask[i] << " ";
            }
            cout << endl;
#endif
        }
        // initialize first row of vectors to init_mask
        for(int idx = -first_diag +1; idx < matrix.num_cols +1; idx++){
            matrix.vectors[idx].match = _mm_load_si128((__m128i*)init_mask);
            matrix.vectors[idx].insert_row = _mm_load_si128((__m128i*)init_mask);
            matrix.vectors[idx].insert_col = _mm_load_si128((__m128i*)init_mask);
        )
        // initialize all vectors as zeroes
        for(int idx = -first_diag; idx < matrix.get_band_size(); idx++){
            matrix.vectors[idx].match = _mm_setzero_si128();
            matrix.vectors[idx].insert_row = _mm_setzero_si128();
            matrix.vectors[idx].insert_col = _mm_setzero_si128();
        }
        // initialize all side vectors as buffer
        for(int idx = matrix.num_cols+2; idx < matrix.get_band_size(); idx += num_cols+1){
            matrix.vectors[idx].match = _mm_set1_epi16(1);
            matrix.vectors[idx].insert_row = _mm_set1_epi16(1);
            matrix.vectors[idx].insert_col = _mm_set1_epi16(1);
        }
    }
//}
