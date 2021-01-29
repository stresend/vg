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


BandedVectorAligner* BandedVectorAligner::init(int8_t* score_mat, int8_t* nt_table, int16_t gap_open, int16_t gap_extend, bool adjust_for_base_quality) {
    //create heap, save in aligner
    BandedVectorHeap* heap = BandedVectorHeap::init();
    BandedVectorAligner* aligner = (BandedVectorAligner*) heap->irreversible_alloc(sizeof(BandedVectorAligner));
    
    aligner->heap = heap;
    aligner->score_mat = (int8_t*) heap->irreversible_alloc(sizeof(int8_t)*16);//16 elements in score_mat
    for(int idx = 0; idx < 16; idx++){
        aligner->score_mat[idx] = score_mat[idx];
    }
    aligner->nt_table = (int8_t*) heap->irreversible_alloc(sizeof(int8_t)*256);//256 elements in nt_table
    for(int idx = 0; idx < 256; idx++){
        aligner->nt_table[idx] = nt_table[idx];
    }

    aligner->adjust_for_base_quality =  adjust_for_base_quality;
    return aligner;
}

void BandedVectorAligner::new_instance(AlignerInstance* instance, HandleGraph* graph, Alignment* alignment, 
		                vector<handle_t>& topological_order, int16_t gap_open, int16_t gap_extend){

    // reset heap before creating new instance
    heap->reset();

    // save instance on heap and then save
    instance = (AlignerInstance*) heap->alloc(sizeof(AlignerInstance));
    instance->graph = graph;
    instance->alignment = alignment;
    instance->gap_open = gap_open;
    instance->gap_extend = gap_extend;

    for(int i = 0; i < topological_order.size(); i++){
        instance->node_id_to_idx[graph->get_id(topological_order[i])] = i;
    }
}

void BandedVectorAligner::align_instance(){
    for(int i = 0; i < current_instance->num_nodes; i++){
        current_instance->matrices[i].fill_matrix(heap, score_mat, nt_table, adjust_for_base_quality);
    }
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
    int matrix_i = idx_within_vec+first_diag-1+matrix_j+8*(vec_idx/num_cols);
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
void BandedVectorMatrix::fill_matrix(BandedVectorHeap* heap, int8_t* score_mat, int8_t* nt_table, bool qual_adjusted){
    
    int query_size = num_diags + 1 + (16*num_diags)-(num_diags+1)%16;// replace with round up function for readability
    int8_t* query = (int8_t*) heap->alloc(sizeof(int8_t)*query_size);
    query_forward(query, score_mat, nt_table, node_seq[0], 0);
    update_first_column(query);
    //outer loop goes through each column, inner loop goes through each vector in each column
    for(int j = 2; j < num_cols+1; j++){
        //determine score vector
        query_forward(query, score_mat, nt_table, node_seq[j-1], j-1);    
        int query_idx = 0;
        for(int idx = j, band_size = get_band_size(), step = num_cols+1; idx<band_size; idx+=step){
            __m128i left_insert_col = _mm_alignr_epi8(
			        vectors[vector_diagonal(idx)].insert_col, 
			        vectors[vector_below(vector_diagonal(idx))].insert_col, 
			        14);
	        __m128i left_match      = vectors[vector_diagonal(idx)].match;
			query_idx = update_vector(left_insert_col, left_match, query, query_idx, idx);
        }
    }

}

//fills out query. Needs to be vectorized!!
void BandedVectorMatrix::query_forward(int8_t* query, int8_t* score_mat, int8_t* nt_table, char node_char, int col_num){
    int query_size = num_diags + 1 + (16*num_diags)-(num_diags+1)%16;     
    const string& read = alignment->sequence();
    
    //placing read into paddded_read
    //TODO: this step needs to be vectorized
    
    //first step adding whatever padding is needed at the start
    int query_idx = 0;
    for(int temp = first_diag + col_num; temp < 0; temp++){
       query[query_idx++] = 0x00;
    }
    //next step, start inserting relevant read information
    for(int temp = max<int>(first_diag + col_num, 0), read_idx = first_diag+col_num + temp; read_idx < read.size() && temp < query_size; temp++){
        query[query_idx++] = score_mat[5 * nt_table[node_char] + nt_table[read[read_idx++]]];
    }

    //final step, fill up tail with buffer
    while(query_idx < query_size){
        query[query_idx++] = 0x00;
    }   
}

//dynamic programming -- thanks dozeu
int BandedVectorMatrix::update_vector(__m128i& left_insert_col, __m128i& left_match, int8_t* query, int query_idx, int idx){

    //TODO: All of this should probably be stored somewhere, I just gotta figure out where to put it
    __m128i gap_open_vector = _mm_set1_epi16(gap_open);
    __m128i min_vector = _mm_set1_epi16(numeric_limits<int8_t>::min());
    //initialize extra gap extension vectors -- thanks dozeu
    __m128i gap_extend1 = _mm_set1_epi16(gap_extend);//I can see this being useful for readability, may not keep
    __m128i gap_extend2 = _mm_add_epi16(gap_extend1, gap_extend1);
    __m128i gap_extend4 = _mm_slli_epi16(gap_extend1, 2);

	__m128i score_vector = _mm_cvtepi8_epi16(_mm_loadu_si128((__m128i*) &query[8*query_idx++]));
	__m128i max_insert_col = _mm_max_epi16(
		   _mm_subs_epi16(
			   left_insert_col, 
			   gap_extend1),
									
		   _mm_subs_epi16(
			   left_match, 
			   gap_open_vector
			   )
		   );

	__m128i max_match = _mm_max_epi16(
			max_insert_col, 
			_mm_adds_epi16(
				score_vector, 
				_mm_alignr_epi8(
					vectors[idx].match, 
					vectors[vector_above(idx)].match, 
					14)
				)
			);
	__m128i max_insert_row = _mm_subs_epi16(
			_mm_max_epi16(
				_mm_subs_epi16(
					_mm_alignr_epi8(
						max_match,
						vectors[vector_above(idx)].match,
						14), 
					gap_open_vector), 
				_mm_alignr_epi8(
					min_vector, 
					vectors[idx].insert_col, 
					14)),
			gap_extend1);
	max_insert_row = _mm_max_epi16(
			max_insert_row, 
			_mm_subs_epi16(
				_mm_alignr_epi8(
					max_insert_row, 
					min_vector, 
					14), 
				gap_extend1
				)
			);
	max_insert_row = _mm_max_epi16(
			max_insert_row, 
			_mm_subs_epi16(
				_mm_alignr_epi8(
					max_insert_row, 
					min_vector, 
					12), 
				gap_extend2)
			);
	max_insert_row = _mm_max_epi16(
			max_insert_row,
			_mm_subs_epi16(
				_mm_alignr_epi8(
					max_insert_row, 
					min_vector,
					8), 
				gap_extend4)
			);
	
	max_match = _mm_max_epi16(max_match, max_insert_row);

	vectors[idx].match = max_match;
	vectors[idx].insert_row = max_insert_row;
	vectors[idx].insert_col = max_insert_col;

    return query_idx;
}

void BandedVectorMatrix::update_first_column(int8_t* query){
    if(is_source){
        query_forward(query, score_mat, nt_table, node_seq[0], 1);
        for(int i = 0; i < number_of_seeds; i++){
            int query_idx = 0;
            for(int idx = 1, band_size = get_band_size(), step = num_cols+1; idx<band_size; idx+=step){
                __m128i left_insert_col = _mm_alignr_epi8(
                        seeds[i]->get_vector_insert_col(seeds[i],idx - 1), 
			            seeds[i]->get_vector_insert_col(seeds[i],idx + 7),// 8 rows below previous vectors
			            14);
	            __m128i left_match      = _mm_alignr_epi8(
			            seeds[i]->get_vector_match(seeds[i], idx-1),
			            seeds[i]->get_vector_match(seeds[i], idx+7),
			            14);
		    	query_idx = update_vector(left_insert_col, left_match, query, query_idx, idx);
            }
        }   
    }else{        
        query_forward(query, score_mat, nt_table, node_seq[0], 1);	    
        int query_idx = 0;
	    
        for(int idx = 1, band_size = get_band_size(), step = num_cols+1; idx<band_size; idx+=step){
            __m128i left_insert_col = _mm_alignr_epi8(
			        vectors[vector_diagonal(idx)].insert_col, 
			        vectors[vector_below(vector_diagonal(idx))].insert_col, 
			        14);
	        __m128i left_match      = vectors[vector_diagonal(idx)].match;
			query_idx = update_vector(left_insert_col, left_match, query, query_idx, idx);
        }
    }        
}

// this function and get_vector_insert_col are very similar, may consider combining to reduce duplicate code
__m128i BandedVectorMatrix::get_vector_match(BandedVectorMatrix* seed, int y){
    alignas(16) short vector_arr[8];
    int vec_idx = 0;
    while(!seed->is_inside_band(seed->num_cols - 1, y) && vec_idx < 8){
        //I would normally have limit<short>::min()-highest_score as buffer, without recalculating highest_score, this is fastest
        vector_arr[vec_idx++] = numeric_limits<short>::min()+numeric_limits<int8_t>::max();
        y++;
    }
    
    while(seed->is_inside_band(seed->num_cols-1, y)&& vec_idx < 8){
        alignas(16) short temp_arr[8];
        _mm_store_si128((__m128i*)temp_arr, seed->vectors[seed->vector_index(y, seed->num_cols-1)].match);
        vector_arr[vec_idx++] = temp_arr[seed->index_within_vector(seed->num_cols-1, y++)];
    }
    while(!seed->is_inside_band(seed->num_cols - 1, y)&& vec_idx < 8){
        //I would normally have limit<short>::min()-highest_score as buffer, without recalculating highest_score, this is fastest
        vector_arr[vec_idx++] = numeric_limits<short>::min()+numeric_limits<int8_t>::max();
        y++;
    }

    return _mm_load_si128((__m128i*)vector_arr);
}

__m128i BandedVectorMatrix::get_vector_insert_col(BandedVectorMatrix* seed, int y){
    alignas(16) short vector_arr[8];
    int vec_idx = 0;
    while(!seed->is_inside_band(seed->num_cols - 1, y) && vec_idx < 8){
        //I would normally have limit<short>::min()-highest_score as buffer, without recalculating highest_score, this is fastest
        vector_arr[vec_idx++] = numeric_limits<short>::min()+numeric_limits<int8_t>::max();
        y++;
    }
    
    while(seed->is_inside_band(seed->num_cols-1, y)&& vec_idx < 8){
        alignas(16) short temp_arr[8];
        _mm_store_si128((__m128i*)temp_arr, seed->vectors[seed->vector_index(y,seed->num_cols-1)].insert_col);
        vector_arr[vec_idx++] = temp_arr[seed->index_within_vector(seed->num_cols-1, y++)];
    }
    while(!seed->is_inside_band(seed->num_cols - 1, y)&& vec_idx < 8){
        //I would normally have limit<short>::min()-highest_score as buffer, without recalculating highest_score, this is fastest
        vector_arr[vec_idx++] = numeric_limits<short>::min()+numeric_limits<int8_t>::max();
        y++;
    }

    return _mm_load_si128((__m128i*)vector_arr);
}

void BandedVectorMatrix::print_full_matrix(){
    //TODO: print full matrix
}

void BandedVectorMatrix::print_rectangularized_band(){
    //TODO: add ability to select which matrix to print
    const string& read = alignment->sequence();
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
void init_BandedVectorMatrix(BandedVectorMatrix& matrix, BandedVectorHeap* heap, AlignerInstance* instance,  
        int8_t* score_mat, int8_t* nt_table, handle_t node, int64_t first_diag, int64_t num_diags,
        int64_t num_cols){

    matrix.num_cols = num_cols;
    matrix.first_diag = first_diag;
    matrix.num_diags = num_diags;
    matrix.alignment = instance->alignment;
    matrix.node = node;
    matrix.node_seq = instance->graph->get_sequence(node);
    matrix.score_mat = score_mat;
    matrix.nt_table = nt_table;
    matrix.gap_extend = instance->gap_extend;
    matrix.gap_open = instance->gap_open;

    short matrix_buffer = numeric_limits<short>::min() + numeric_limits<int8_t>::max();
    
    //alloc vectors
    matrix.vectors = (SWVector*) heap->alloc(sizeof(SWVector)*matrix.get_band_size());
    //initialize init_arr to represent corner of matrix
    alignas(16) short init_arr[8];
    _mm_store_si128((__m128i*)init_arr, _mm_setzero_si128());
    for(int idx = 0; idx < -first_diag +1 && idx < 8; idx++){
        init_arr[idx] = matrix_buffer;
    }


    //fill out seeds, I run follow_edges twice to count and then alloc
    int num_seeds = 0;
    instance->graph->follow_edges(node, true, [&](const handle_t& prev) {
                num_seeds++;
            });
    matrix.seeds = (BandedVectorMatrix**)heap->alloc(sizeof(BandedVectorMatrix*)*num_seeds);
    matrix.number_of_seeds = num_seeds;
    num_seeds = 0;
    instance->graph->follow_edges(node, true, [&](const handle_t& prev) {
                matrix.seeds[num_seeds] = &(instance->matrices[instance->node_id_to_idx[instance->graph->get_id(prev)]]);
                num_seeds++;
            });
     
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
    cout << endl << "End debug init_BandedVectorMatrix"<< endl;
#endif
       
}

}
