/**
 * \file vectorized_banded_aligner.cpp
 *
 * Implements the vectorized banded global aligner
 *
 */

#define init_BandedVectorMatrix_debug

#include "vectorized_banded_aligner.hpp"
#include<iostream>
#include<limits>

using namespace std;

namespace vg {


const size_t BandedVectorHeap::init_block_size = (1 << 20); // 1 MB

BandedVectorAligner* BandedVectorAligner::init(Alignment& alignment, HandleGraph& g, int64_t band_padding, 
                                               bool permissive_banding, bool adjust_for_base_quality) {

    BandedVectorHeap* heap = BandedVectorHeap::init();
    BandedVectorAligner* aligner = (BandedVectorAligner*) heap->irreversible_alloc(sizeof(BandedVectorAligner));
    aligner->heap = heap;
    aligner->graph = g;
    aligner->alignment = alignment;

    //convert vector made from topological sort into array
    vector<handle_t> topo_order = algorithms::lazier_topological_order(&g);
    aligner->number_of_nodes = topo_order.size();
    aligner->topological_order = (handle_t*) heap->alloc(sizeof(handle_t)*aligner->number_of_nodes);
    for(int i = 0; i < aligner->number_of_nodes; i++){
	    aligner->topological_order[i] = topo_order[i];
    }
    return aligner;
}

void BandedVectorAligner::destroy() {
    heap->destroy();
}

BandedVectorHeap* BandedVectorHeap::init() {
    static_assert(init_block_size > sizeof(BandedVectorHeap), "heap must fit within the initial block");
    BandedVectorMemoryBlock* init_block = BandedVectorMemoryBlock::init(init_block_size);
    BandedVectorHeap* heap = (BandedVectorHeap*) init_block->irreversible_alloc(sizeof(BandedVectorHeap));
    heap->head = heap->curr = init_block;
    return heap;
}

void BandedVectorHeap::reset() {
    BandedVectorMemoryBlock* prev = nullptr;
    BandedVectorMemoryBlock* resetter = head;
    while (prev != curr) {
        resetter->reset();
        prev = resetter;
        resetter = resetter->next_block();
    }
    curr = head;
}

void* BandedVectorHeap::alloc(size_t _size) {
    curr = curr->find_space(_size);
    return curr->alloc(_size);
}

void* BandedVectorHeap::irreversible_alloc(size_t _size) {
    curr = curr->find_space(_size);
    return curr->irreversible_alloc(_size);
}


void BandedVectorHeap::destroy() {

    BandedVectorMemoryBlock* destroyer = head;
    while (destroyer) {
        BandedVectorMemoryBlock* next = destroyer->next_block();
        destroyer->destroy();
        destroyer = next;
    }
}

BandedVectorHeap::BandedVectorMemoryBlock* BandedVectorHeap::BandedVectorMemoryBlock::init(size_t _size) {

    // add room for the block itself
    size_t block_size = round_to_align(sizeof(BandedVectorMemoryBlock));
    _size = round_to_align(_size) + block_size;

    void* mem;
    int return_code = posix_memalign(&mem, sizeof(__m128i), _size);
    if (return_code) {
        cerr << "error: unable to posix_memalign " << _size << " bytes with alignment " << sizeof(__m128i) << endl;
        exit(1);
    }

    BandedVectorMemoryBlock* block = (BandedVectorMemoryBlock*) mem;

    // the pointer to be free'd
    block->block_begin = (uint8_t*) mem;
    // where the allocations can occur
    block->bottom = block->curr = block->block_begin + block_size;
    block->top = block->block_begin + _size;

    return block;
}


// from https://stackoverflow.com/questions/3407012/c-rounding-up-to-the-nearest-multiple-of-a-number
size_t BandedVectorHeap::BandedVectorMemoryBlock::round_to_align(size_t _size) {
    return (_size + sizeof(__m128i) - 1) & -sizeof(__m128i);
}

size_t BandedVectorHeap::BandedVectorMemoryBlock::remaining() const {
    return top - curr;
}

size_t BandedVectorHeap::BandedVectorMemoryBlock::size() const {
    return top - block_begin;
}

void* BandedVectorHeap::BandedVectorMemoryBlock::alloc(size_t _size) {
    _size = round_to_align(_size);
    void* ptr = curr;
    curr += _size;
    return ptr;
}

void* BandedVectorHeap::BandedVectorMemoryBlock::irreversible_alloc(size_t _size) {
    void* ptr = alloc(_size);
    bottom = curr;
    return ptr;
}

void BandedVectorHeap::BandedVectorMemoryBlock::reset() {
    curr = bottom;
}

BandedVectorHeap::BandedVectorMemoryBlock* BandedVectorHeap::BandedVectorMemoryBlock::find_space(size_t _size) {
    BandedVectorMemoryBlock* prev = nullptr;
    BandedVectorMemoryBlock* block = this;
    while (block && block->remaining() < _size) {
        prev = block;
        block = block->next_block();
    }
    if (!block) {
        block = init(max(_size, 2 * size()));
        prev->next = block;
    }
    return block;
}

BandedVectorHeap::BandedVectorMemoryBlock* BandedVectorHeap::BandedVectorMemoryBlock::next_block() {
    return next;
}

void BandedVectorHeap::BandedVectorMemoryBlock::destroy() {
    free(block_begin);
}


/* Note on vector indexing:
 *    0 index reserved for a "buffer" vector, calls to vector_above with vec_idx < num_cols returns this index
 *    1->v indices are vectors that are filled during dynamic programming
 *    v+1 index is another "buffer" vector, call to vect_below on vec_idx > v-num_cols returns this index
 */
bool BandedVectorMatrix::is_inside_band(int i, int j){
    return i >= first_diag+j && i < first_diag+j+num_diags;
}

int BandedVectorMatrix::vector_index(int i, int j){
    return ((i-j-first_diag+1)/8)*(num_cols+1)+j%num_cols+1;
}

int BandedVectorMatrix::index_within_vector(int i, int j){
    return (i-first_diag-j+1)%8;
}

pair<int, int> BandedVectorMatrix::coordinate(int vec_idx, int idx_within_vec){
    int matrix_j = vec_idx%num_cols;
    int matrix_i =  idx_within_vec+first_diag-1+matrix_j+8*(vec_idx/num_cols);
    return pair<int, int>(matrix_i, matrix_j);
}


int BandedVectorMatrix::vector_above(int vec_idx){
    return max<int64_t>(0, vec_idx - num_cols);
}

int BandedVectorMatrix::vector_below(int vec_idx){
    return min<int64_t>(get_band_size(), vec_idx+num_cols+1)%get_band_size();
}

int BandedVectorMatrix::vector_diagonal(int vec_idx){
    return vec_idx-1;
}

int BandedVectorMatrix::get_band_size(){
    return num_cols+(num_cols+1)*((num_diags-1)/8)+1;
}

//dynamic programming section
void BandedVectorMatrix::fill_matrix(HandleGraph& graph, int8_t* score_mat, int8_t* nt_table, 
		                     int8_t gap_open, int8_t gap_extend, bool qual_adjusted){
    //initialize gap_open and gap_extend vectors -- would it be best to keep in BandedVectorMatrix or just 
    __m128i gap_extend_vector = _mm_set1_epi16(gap_extend);
    __m128i gap_open_vector = _mm_set1_epi16(gap_open);


    //initialize extra gap extension vectors -- thanks dozeu
    __m128i gap_extend1 = _mm_load_si128(&gap_extend_vector);//I can see this being useful for readability, may not keep
    __m128i gap_extend2 = _mm_add_epi16(gap_extend1, gap_extend1);
    __m128i gap_extend4 = _mm_slli_epi16(gap_extend1, 2);
    __m128i gap_extend8 = _mm_slli_epi16(gap_extend1, 3);

    //outer loop goes through each column, inner loop goes through each vector in each column
    for(int j = 1; j < num_cols+1; j++){
        for(int idx = j; idx<get_band_size(); idx+=num_cols+1){
            //determine score vector -- need to implement
            __m128i score_vector;

	    //max of insert_row
	    __m128i left_insert_row = _mm_alignr_epi8(
			    vectors[vector_diagonal(idx)].insert_row, 
			    vectors[vector_below(vector_diagonal(idx))].insert_row, 
			    2);
	    __m128i left_match      = _mm_alignr_epi8(
			    vectors[vector_diagonal(idx)].match,
			    vectors[vector_below(vector_diagonal(idx))].match,
			    2);
	    vectors[idx].insert_row = _mm_max_epi16(_mm_sub_epi16(left_insert_row, gap_extend_vector),
			                            _mm_sub_epi16(left_match, gap_open_vector));
	    //max of insert_col -- thanks dozeu
	    
	    //max of match -- has to be done last
	    __m128i diagonal_match = vectors[vector_diagonal(idx)].match;
	    vectors[idx].match = _mm_max_epi16(
			    vectors[idx].insert_row, 
			    _mm_max_epi16(
				    vectors[idx].insert_col, 
                                    _mm_adds_epi16(score_vector, diagonal_match)));
	}
    }

}

void BandedVectorMatrix::print_full_matrix(){
    //TODO: print full matrix
}

void BandedVectorMatrix::print_rectangularized_band(){
    //TODO: add ability to select which matrix to print
    const string& read = alignment.sequence();
    cerr << "Print rectangualrized band of match"<<endl;
    for(int j = 0; j < num_cols; j++){
        for(int i = 0; i < num_cols - first_diag - num_diags; i++){
             alignas(16) short temp[8]; 
	     _mm_store_si128((__m128i*)temp, vectors[vector_index(i, j)].match);
	     cerr << temp[index_within_vector(i,j)];
	}
	cerr << endl;
    }
}

/* Notes on vector initialization:
 *    Here's an example of a banded graph turned vectors(let 1's be the buffer area, 0's be the dp region)
 * 0001111111        11111111111
 * 0000111111        11100000000
 * 0000011111        11000000000
 * 1000001111        10000000000
 * 1100000111        10000000000
 * 1110000011   ->   10000000000
 * 1111000001        10000000000
 * 1111100000
 * 1111110000
 * 1111111000
 *
 * Notice: this implemetation adds a single row of buffer to the top of the vectorized band.
 * also notice that an extra buffer column is inserted before the first column
 * 
 */
void init_BandedVectorMatrix(BandedVectorMatrix& matrix, BandedVectorHeap* heap, 
		             Alignment& alignment, handle_t node, int8_t* score_mat, 
			     int64_t first_diag, int64_t num_diags, int64_t num_cols){

    matrix.num_cols = num_cols;
    matrix.first_diag = first_diag;
    matrix.num_diags = num_diags;
    matrix.alignment = alignment;
    matrix.node = node;

    //initialize copy of score_mat that is stored in each matrix, find max and determine buffer value
    int8_t max_num = 0;
    matrix.score_mat = (int8_t*) heap->alloc(sizeof(int8_t) * 32);//I'm pretty sure score_mat is 32 ints, this may change
    for(int idx = 0; idx < 32; idx++){
	max_num = max<int8_t>(max_num, score_mat[idx]);
        matrix.score_mat[idx] = score_mat[idx];
    }
    short matrix_buffer = numeric_limits<int8_t>::max() - max_num;
    
    //alloc vectors
    matrix.vectors = (SWVector*) heap->alloc(sizeof(SWVector)*matrix.get_band_size());
    //initialize init_arr to represent corner of matrix
    alignas(16) short init_arr[8];
    _mm_store_si128((__m128i*)init_arr, _mm_setzero_si128());
    for(int idx = 0; idx < -first_diag +1 && idx < 8; idx++){
        init_arr[idx] = matrix_buffer;
    }

# ifdef init_BandedVectorMatrix_debug
    cout << "band_size test: " << matrix.get_band_size() << endl;
    cout << "init_arr before iterations: ";
    for(int i = 0; i < 8; i++){
        cout << init_arr[i] << " ";
    }
    cout << endl;
#endif

    // initialize top left corner of rectangularized vector matrix
    for(int idx = 1; idx<-first_diag+1; idx++){
        matrix.vectors[idx].match = _mm_load_si128((__m128i*)init_arr);
        matrix.vectors[idx].insert_row = _mm_load_si128((__m128i*)init_arr);
        matrix.vectors[idx].insert_col = _mm_load_si128((__m128i*)init_arr);
        init_arr[-first_diag-idx+1] = 0;
#ifdef init_BandedVectorMatrix_debug
        cout << "init_arr during iteration " << idx << ": ";
        for(int i = 0; i < 8; i++){
            cout << init_arr[i] << " ";
        }
        cout << endl;
#endif
    }

#ifdef init_BandedVectorMatrix_debug
    cout << "Init_arr initialization indices: ";
#endif
    // initialize first row of vectors to init_mask
    for(int idx = -first_diag+1; idx < matrix.num_cols+1 && idx < matrix.get_band_size(); idx++){
        matrix.vectors[idx].match = _mm_load_si128((__m128i*)init_arr);
        matrix.vectors[idx].insert_row = _mm_load_si128((__m128i*)init_arr);
        matrix.vectors[idx].insert_col = _mm_load_si128((__m128i*)init_arr);
#ifdef init_BandedVectorMatrix_debug
        cout << idx << ", ";
#endif
    }

#ifdef init_BandedVectorMatrix_debug
    cout << endl << "Zero_vector indices: ";
#endif
    // initialize all vectors as zeroes
    for(int idx = matrix.num_cols+1; idx < matrix.get_band_size(); idx++){
        matrix.vectors[idx].match = _mm_setzero_si128();
        matrix.vectors[idx].insert_row = _mm_setzero_si128();
        matrix.vectors[idx].insert_col = _mm_setzero_si128();
#ifdef init_BandedVectorMatrix_debug
        cout << idx << ", ";
#endif
    }
#ifdef init_BandedVectorMatrix_debug
    cout << endl << "buffer_vector indices: ";
#endif
    // initialize all side vectors as buffer
    for(int idx = 0; idx < matrix.get_band_size(); idx += num_cols+1){
        matrix.vectors[idx].match = _mm_set1_epi16(matrix_buffer);
        matrix.vectors[idx].insert_row = _mm_set1_epi16(matrix_buffer);
        matrix.vectors[idx].insert_col = _mm_set1_epi16(matrix_buffer);
#ifdef init_BandedVectorMatrix_debug
        cout << idx << ", ";
#endif
    }
#ifdef init_BandedVectorMatrix_debug
    cout << endl << "End debug init_BandedVectorMatrix"<< endl;
#endif
}

}
