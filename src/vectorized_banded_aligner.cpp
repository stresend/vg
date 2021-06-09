/**
 * \file vectorized_banded_aligner.cpp
 *
 * Implements the vectorized banded global aligner
 *
 */

//#define init_BandedVectorMatrix_debug
//#define BandedVectorAligner_debug

#include "vectorized_banded_aligner.hpp"
#include<iostream>
#include<limits>

using namespace std;

namespace vg {


const size_t BandedVectorHeap::init_block_size = (1 << 20); // 1 MB


BandedVectorAligner* BandedVectorAligner::init(const int8_t* score_mat, int8_t* nt_table, int16_t gap_open, int16_t gap_extend, bool adjust_for_base_quality) {
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
    aligner->gap_extend = gap_extend;
    aligner->gap_open = gap_open;
    return aligner;
}

void BandedVectorAligner::new_instance(const HandleGraph& graph, Alignment& alignment, int32_t band_padding){

    // reset heap before creating new instance
    heap->reset();

    // save instance on heap and then save
    current_instance = (AlignerInstance*) heap->alloc(sizeof(AlignerInstance));
	current_instance->topological_order = handlealgs::lazier_topological_order(&graph);
	current_instance->source_nodes     = handlealgs::head_nodes(&graph);
	current_instance->sink_nodes         = handlealgs::tail_nodes(&graph);
    current_instance->matrices = (BandedVectorMatrix*) heap->alloc(current_instance->topological_order.size() * sizeof(BandedVectorMatrix));
    current_instance->graph = &graph;
    current_instance->alignment = &alignment;

    for(int i = 0; i < current_instance->topological_order.size(); i++){
        current_instance->handle_to_idx.insert(pair<handle_t, int>(current_instance->topological_order[i], i));
    }

    //initialize matrices
    //first, find banded path
    vector<bool> node_masked;
    vector<pair<int, int>> band_ends;
    find_banded_paths(current_instance, node_masked, band_ends);
    vector<int64_t> seq_lens_out;
	shortest_seq_paths(seq_lens_out);
    for(int i = 0; i < current_instance->topological_order.size(); i++){
        int num_diags = 1 + band_ends[i].second - band_ends[i].first;
        int num_cols = graph.get_length(current_instance->topological_order[i]);
        init_BandedVectorMatrix(current_instance->matrices[i], heap, current_instance, current_instance->topological_order[i], band_ends[i].first, num_diags, num_cols, seq_lens_out[i]);
    }
}

void BandedVectorAligner::align_instance(){
    for(int i = 0, num_nodes = current_instance->topological_order.size(); 
        i < num_nodes; i++){
        current_instance->matrices[i].fill_matrix(current_instance, heap, score_mat, nt_table, gap_open, gap_extend, adjust_for_base_quality);
    }
	traceback();
}

// borrowed from BandedGlobalAligner
void BandedVectorAligner::path_lengths_to_sinks(AlignerInstance* instance, vector<int>& shortest_path_to_sink,
                                                         vector<int>& longest_path_to_sink) {
#ifdef BandedVectorAligner_debug
    cerr << "[BandedVectorAligner_debug::path_lengths_to_sinks]: finding longest and shortest paths to sink node" << endl;
#endif
    
    // find the longest path from the right side of each matrix to the end of the graph
    // set initial values
    longest_path_to_sink.resize(instance->topological_order.size(), 0);
    shortest_path_to_sink.resize(instance->topological_order.size(), numeric_limits<int>::max());
    
    // set base case (longest path already set to 0)
    for (const handle_t& handle : instance->sink_nodes) {
        shortest_path_to_sink[instance->handle_to_idx.at(handle)] = 0;
    }
    
    // iterate in reverse order
    for (int64_t i = instance->topological_order.size() - 1; i >= 0; i--) {
        
        int64_t node_seq_len = instance->graph->get_length(instance->topological_order[i]);
        // compute longest path through this node to right side of incoming matrices
        int64_t longest_path_length = longest_path_to_sink[i] + node_seq_len;
        int64_t shortest_path_length = shortest_path_to_sink[i] + node_seq_len;
        
        instance->graph->follow_edges(instance->topological_order[i], true, [&](const handle_t& prev) {
            
            int64_t prev_idx = instance->handle_to_idx.at(prev);
            
            if (longest_path_to_sink[prev_idx] < longest_path_length) {
                
#ifdef BandedVectorAligner_debug
                cerr << "[BandedVectorAligner_debug::path_lengths_to_sinks]: path through " << instance->graph.get_id(prev) << " of length " << longest_path_length << " to node at index " << prev_idx << " is longer than current longest path " << longest_path_to_sink[prev_idx] << ", updating it now" << endl;
#endif
                
                longest_path_to_sink[prev_idx] = longest_path_length;
            }
            if (shortest_path_to_sink[prev_idx] > shortest_path_length) {
                
#ifdef BandedVectorAligner_debug
                cerr << "[BandedVectorAligner_debug::path_lengths_to_sinks]: path through " << instance->graph.get_id(prev) << " of length " << shortest_path_length << " to node at index " << prev_idx << " is shorter than current shortest path " << shortest_path_to_sink[prev_idx] << ", updating it now" << endl;
#endif
                
                shortest_path_to_sink[prev_idx] = shortest_path_length;
            }
        });
    }
}



void BandedVectorAligner::traceback() {
    
    int64_t read_length = current_instance->alignment->sequence().length();
    int32_t empty_score = read_length > 0 ? -gap_open - (read_length - 1) * gap_extend : 0;
    
        
    Alignment* alignment = current_instance->alignment;

            
    // what node does the alignment start at
    int64_t node_id;
    matrix_t mat = Match;
	int64_t i = read_length, j;
	int16_t max_score = numeric_limits<int16_t>::min();
	
	for(auto sink : current_instance->sink_nodes){
		int16_t matrix_max = current_instance->matrices[current_instance->handle_to_idx.at(sink)].get_score_at(mat, current_instance->graph->get_length(sink), read_length);
		if(max_score < matrix_max){
			max_score = matrix_max;
			node_id = current_instance->graph->get_id(sink);
			j = current_instance->graph->get_length(sink);
		}
		
	}
    // we only start in a gap if the sequence if there is nothing to align
    bool in_lead_gap = alignment->sequence().empty();
            
            
            
    // do traceback
   VBBuilder builder(*alignment);
            
    while (node_id != 0) {
        int64_t node_idx = current_instance->handle_to_idx[current_instance->graph->get_handle(node_id)];
        // trace through the matrix
        current_instance->matrices[node_idx].traceback(*current_instance->graph, builder, i, j, mat, in_lead_gap, score_mat,
                                                    nt_table, gap_open, gap_extend);
        // trace over edges
        current_instance->matrices[node_idx].traceback_over_edge(*current_instance->graph,   *current_instance->alignment, builder,  i, j, mat, in_lead_gap,
                                                               node_id, score_mat, nt_table, gap_open, gap_extend);
    }
            
    // construct the alignment path
    builder.finalize_alignment();
            
    // add score to alignment
    alignment->set_score(max_score);

}
    




// some code borrowed from BandedGlobalAligner
// fills vectors with whether nodes are masked by the band width, and the band ends of each node
void BandedVectorAligner::find_banded_paths(AlignerInstance* instance,
                                                     vector<bool>& node_masked,
                                                     vector<pair<int, int>>& band_ends) {
    
    // keeps track of which nodes cannot reach the bottom corner within the band
    node_masked.resize(instance->topological_order.size(), false);
    
    // the bottom and top indices of the band in the rightmost column of each node's matrix
    band_ends.resize(instance->topological_order.size(), make_pair(numeric_limits<int>::max(),
                                                         numeric_limits<int>::min()));
    
    // find the longest and shortest path from each node to any sink
    vector<int> shortest_path_to_sink;
    vector<int> longest_path_to_sink;
    path_lengths_to_sinks(instance, shortest_path_to_sink, longest_path_to_sink);
    

	// initialize with wide enough bands that every source can hit every connected sink
	for (const handle_t& init_node : instance->source_nodes) {
		int64_t init_node_idx = instance->handle_to_idx.at(init_node);
		int64_t init_node_seq_len = instance->graph->get_length(init_node);
		band_ends[init_node_idx].first = min<int>(-instance->band_padding,
											 instance->alignment->sequence().size()
											 - init_node_seq_len + longest_path_to_sink[init_node_idx] - instance->band_padding);
		band_ends[init_node_idx].second = max<int>(instance->band_padding,
											  instance->alignment->sequence().size()
											  - init_node_seq_len + shortest_path_to_sink[init_node_idx] + instance->band_padding);
#ifdef BandedVectorAligner_debug
		cerr << "[BandedVectorAligner::find_banded_paths]: initializing band path end at node " << instance->graph.get_id(init_node) << " at index " << init_node_idx << " to top " << band_ends[init_node_idx].first << ", and bottom " << band_ends[init_node_idx].second << " from shortest and longest paths of length " << shortest_path_to_sink[init_node_idx] << " and " << longest_path_to_sink[init_node_idx] << " compared to read length " << instance->alignment->sequence().size() << " with padding " << instance->band_padding << endl;
#endif
	}
        
    
    // iterate through the rest of the nodes in topological order
    for (int64_t i = 0; i < instance->topological_order.size(); i++) {
        const handle_t& node = instance->topological_order[i];
        int64_t node_seq_len = instance->graph->get_length(node);
        
        int64_t extended_band_top = band_ends[i].first + node_seq_len;
        int64_t extended_band_bottom = band_ends[i].second + node_seq_len;
        
#ifdef BandedVectorAligner_debug
        cerr << "[BandedVectorAligner::find_banded_paths]: following edges out of node " << instance->graph.get_id(node) << " at index " << i << " with sequence " << instance->graph.get_sequence(node) << ", band of " << band_ends[i].first << ", " << band_ends[i].second << " extending to " << extended_band_top << ", " << extended_band_bottom << endl;
#endif
        // can alignments from this node reach the bottom right corner within the band?
        if (extended_band_top + shortest_path_to_sink[i] > instance->alignment->sequence().size()
            || extended_band_bottom + longest_path_to_sink[i] < instance->alignment->sequence().size()) {
            
            node_masked[i] = true;
            
#ifdef BandedVectorAligner_debug
            cerr << "[BandedVectorAligner::find_banded_paths]: cannot complete alignment to read of length " << alignment.sequence().size() << " along shortest path " << shortest_path_to_sink[i] << " or longest path " << longest_path_to_sink[i] << ", which reach range " << extended_band_top + shortest_path_to_sink[i] << ", " << extended_band_bottom + longest_path_to_sink[i] << endl;
#endif
            continue;
        }
        
       instance->graph->follow_edges(node, false, [&](const handle_t& next) {
            
            int64_t node_out_idx = instance->handle_to_idx.at(next);
            
#ifdef BandedVectorAligner_debug
            cerr << "[BandedVectorAligner::find_banded_paths]: extending band to node at index " << node_out_idx << endl;
#endif
            
            if (extended_band_top < band_ends[node_out_idx].first) {
#ifdef BandedVectorAligner_debug
                cerr << "[BandedVectorAligner::find_banded_paths]: updating top band limit from  " << band_ends[node_out_idx].first << " to " << extended_band_top << endl;
#endif
                band_ends[node_out_idx].first = extended_band_top;
            }
            
            if (extended_band_bottom > band_ends[node_out_idx].second) {
#ifdef BandedVectorAligner_debug
                cerr << "[BandedVectorAligner::find_banded_paths]: updating bottom band limit from  " << band_ends[node_out_idx].second << " to " << extended_band_bottom << endl;
#endif
                band_ends[node_out_idx].second = extended_band_bottom;
            }
        });
    }
}

// returns the shortest sequence from any source node to each node

void BandedVectorAligner::shortest_seq_paths(vector<int64_t>& seq_lens_out) {
    
    // initialize vector with min identity to store sequence lengths
    seq_lens_out.resize(current_instance->topological_order.size(), numeric_limits<int64_t>::max());
    
    // base cases
    for (const handle_t& handle : current_instance->source_nodes) {
        seq_lens_out[current_instance->handle_to_idx[handle]] = 0;
    }
    
    // dynamic programming to calculate sequence lengths for rest of nodes
    for (size_t i = 0; i < current_instance->topological_order.size(); i++) {
        int64_t seq_len = current_instance->graph->get_length(current_instance->topological_order[i]) + seq_lens_out[i];

        current_instance->graph->follow_edges(current_instance->topological_order[i], false, [&](const handle_t& handle) {
            int64_t target_idx = current_instance->handle_to_idx.at(handle);
            // find the shortest sequence that can reach the top left corner of the matrix
            if (seq_len < seq_lens_out[target_idx]) {
                seq_lens_out[target_idx] = seq_len;
            }
        });
    }
}

void BandedVectorAligner::destroy() {
    heap->destroy();
}

VBBuilder::VBBuilder(Alignment& alignment) :
                                                   alignment(alignment),
                                                   matrix_state(Match),
                                                   matching(false),
                                                   current_node_id(0),
                                                   edit_length(0),
                                                   edit_read_end_idx(0){

}

void VBBuilder::update_state(const HandleGraph& graph, matrix_t matrix,
                                                           const handle_t& node, int64_t read_idx, int64_t node_idx,
                                                           bool empty_node_seq) {
#ifdef BandedVectorAligner_traceback_debug
    cerr << "[VBBuilder::update_state] beginning " << (empty_node_seq ? "" : "non-") << "empty state update for read index " << read_idx << ", node seq index " << node_idx << endl;
#endif
    if (graph.get_id(node) != current_node_id) {
#ifdef BandedVectorAligner_traceback_debug
        cerr << "[VBBuilder::update_state] at new node " << graph.get_id(node) << " previously " << current_node_id << endl;
#endif
        // conclude current mapping and proceed to next node
        finish_current_node();
        current_node_id = graph.get_id(node);
        current_node_sequence = graph.get_sequence(node);
        matrix_state = matrix;
        if (matrix_state == Match) {
            matching = (alignment.sequence()[read_idx] == current_node_sequence[node_idx]);
        }
        edit_length = !empty_node_seq;
        edit_read_end_idx = read_idx;
    }
    else if (matrix != matrix_state) {
#ifdef BandedVectorAligner_traceback_debug
        cerr << "[VBBuilder::update_state] transitioning into a new matrix" << endl;
#endif
        // transitioned into another matrix, finish current edit and begin new one
        if (!empty_node_seq) {
            finish_current_edit();
        }
        matrix_state = matrix;
        if (matrix_state == Match) {
            matching = (alignment.sequence()[read_idx] == current_node_sequence[node_idx]);
        }
        edit_length = 1;
        edit_read_end_idx = read_idx;
    }
    else if (matrix == Match &&
             (alignment.sequence()[read_idx] == current_node_sequence[node_idx]) != matching) {
#ifdef BandedVectorAligner_traceback_debug
        cerr << "[VBBuilder::update_state] switching between match and mismatch" << endl;
#endif
        // switch from match to mismatch state or vice versa
        finish_current_edit();
        matching = !matching;
        edit_length = 1;
        edit_read_end_idx = read_idx;
    }
    else {
        // same edit, extend length
        edit_length++;
    }
    
#ifdef BandedVectorAligner_traceback_debug
    cerr << "[VBBuilder::update_state] finished updating state, matrix is " << (matrix_state == Match ? "match" : (matrix_state == InsertRow ? "insert row" : "insert column" )) << ", is matching? " << (matching ? "yes" : "no") << ", edit length " << edit_length << ", edit end index (on read) " << edit_read_end_idx << ", current node " << current_node_id << endl;
#endif
}

void VBBuilder::finish_current_edit() {
    
#ifdef BandedVectorAligner_traceback_debug
    cerr << "[VBBuilder::finish_current_edit] finishing edit" << endl;
#endif
    
    mapping_edits.emplace_front();
    
    switch (matrix_state) {
        case Match:
            mapping_edits.front().set_from_length(edit_length);
            mapping_edits.front().set_to_length(edit_length);
            
            if (!matching) {
                mapping_edits.front().set_sequence(alignment.sequence().substr(edit_read_end_idx - edit_length + 1,
                                                                               edit_length));
            }
            
            break;
            
        case InsertRow:
            mapping_edits.front().set_from_length(0);
            mapping_edits.front().set_to_length(edit_length);
            mapping_edits.front().set_sequence(alignment.sequence().substr(edit_read_end_idx - edit_length + 1,
                                                                           edit_length));
            break;
            
        case InsertCol:
            mapping_edits.front().set_from_length(edit_length);
            mapping_edits.front().set_to_length(0);
            break;
            
        default:
            cerr << "error:[BandedGlobalAligner] unrecognized matrix type" << endl;
            assert(0);
            break;
    }
#ifdef BandedVectorAligner_traceback_debug
    cerr << "[VBBuilder::finish_current_edit] edit: " << pb2json(mapping_edits.front()) << endl;
#endif
}

void VBBuilder::finish_current_node() {

    // sentinel for first iteration
    if (current_node_id == 0) {
#ifdef BandedVectorAligner_traceback_debug
        cerr << "[VBBuilder::finish_current_node] at beginning of traceback, not creating a mapping" << endl;
#endif
        return;
    }
    
    finish_current_edit();
    
#ifdef BandedVectorAligner_traceback_debug
    cerr << "[VBBuilder::finish_current_node] finishing mapping for node " << current_node_id << endl;
#endif
    
    node_mappings.emplace_front();
    for (Edit edit : mapping_edits) {
        *(node_mappings.front().add_edit()) = edit;
    }
    
#ifdef BandedVectorAligner_traceback_debug
    cerr << "[VBBuilder::finish_current_node] mapping: " << pb2json(node_mappings.front()) << endl;
#endif
    
    mapping_edits.clear();
    
    (*(node_mappings.front().mutable_position())).set_node_id(current_node_id);
    // note: global alignment always starts at beginning of node, default offset 0 is correct
}

void VBBuilder::finalize_alignment() {
    
    finish_current_node();
    
#ifdef debug_banded_aligner_traceback
    cerr << "[VBBuilder::finalize_alignment] finalizing alignment" << endl;
#endif
    
    alignment.clear_path();
    Path* path = alignment.mutable_path();
    
    int32_t mapping_rank = 1;
    for (Mapping& mapping : node_mappings) {
        mapping.set_rank(mapping_rank);
        *(path->add_mapping()) = mapping;
        mapping_rank++;
    }
    
    node_mappings.clear();
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BABuilder::finalize_alignment] alignment: " << pb2json(alignment) << endl;
#endif
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
void BandedVectorMatrix::fill_matrix(AlignerInstance* instance, BandedVectorHeap* heap, int8_t* score_mat, int8_t* nt_table, int16_t gap_open, int16_t gap_extend, bool qual_adjusted){
    
    int query_size = num_diags + 1 + (16*num_diags)-(num_diags+1)%16;// replace with round up function for readability
    int8_t* query = (int8_t*) heap->alloc(sizeof(int8_t)*query_size);
    query_forward(query, score_mat, nt_table, node_seq[0], 0);
    update_first_column(score_mat, nt_table, gap_open, gap_extend, query);
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
			query_idx = update_vector(left_insert_col, left_match, gap_open, gap_extend, query, query_idx, idx);
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
int BandedVectorMatrix::update_vector(__m128i& left_insert_col, __m128i& left_match, int16_t gap_open, int16_t gap_extend, int8_t* query, int query_idx, int idx){

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

void BandedVectorMatrix::update_first_column(int8_t* score_mat, int8_t* nt_table, int16_t gap_open, int16_t gap_extend, int8_t* query){
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
		    	query_idx = update_vector(left_insert_col, left_match, gap_open, gap_extend, query, query_idx, idx);
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
			query_idx = update_vector(left_insert_col, left_match, gap_open, gap_extend, query, query_idx, idx);
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
int16_t BandedVectorMatrix::get_score_at(matrix_t mat, int i, int j){
	alignas(16) int16_t temp_vector[8];
	switch(mat){
		case Match: {
			 _mm_store_si128((__m128i*)temp_vector, vectors[vector_index(i, j)].match);

			break;
		}
		case InsertRow: {
			 _mm_store_si128((__m128i*)temp_vector, vectors[vector_index(i, j)].insert_row);
			break;
		}
		case InsertCol:{
			 _mm_store_si128((__m128i*)temp_vector, vectors[vector_index(i, j)].insert_col);
			break;
		}
	}
	return temp_vector[index_within_vector(i, j)];
}

void BandedVectorMatrix::traceback(const HandleGraph& graph, VBBuilder& builder,
                                                       int64_t& i, int64_t& j, matrix_t& mat, bool& in_lead_gap,
                                                       const int8_t* score_mat, const int8_t* nt_table,
                                                       const int8_t gap_open,  const int8_t gap_extend) {
    
#ifdef debug_banded_aligner_traceback
    cerr << "[BandedVectorMatrix::traceback] starting traceback back through node " << graph.get_id(node) << " from rectangular coordinates (" << i << ", " << j << "), currently " << (in_lead_gap ? "" : "not ") << "in a lead gap" << endl;
#endif
    
   const string& read = alignment->sequence();
    
   string node_seq = graph.get_sequence(node);
   int64_t ncols = node_seq.size();
   int64_t node_id = graph.get_id(node);
    
    int16_t min_inf = numeric_limits<int16_t>::min();
    int16_t curr_score;
    int16_t source_score;
    int16_t score_diff;
	
    // do node traceback unless we are in the lead gap implied at the edge of the DP matrix or we
    // are already at a node boundary trying to get across
    while ((j > 0 || mat == InsertRow) && !in_lead_gap) {
        
#ifdef debug_banded_aligner_traceback
     cerr << "[BAMatrix::traceback] traceback coordinates (" << i << ", " << j << "), current matrix is " << (mat == Match ? "match" : (mat == InsertCol ? "insert column" : "insert row")) << endl;
#endif
        
     // add traceback step to alignment
    builder.update_state(graph, mat, node, i + first_diag + j, j);
        
            
#ifdef debug_banded_aligner_traceback
         cerr << "[BandedVectorMatrix::traceback] taking inside matrix deflection to " << (mat == Match ? "match" : (mat == InsertCol ? "insert column" : "insert row")) << endl;
#endif
            
        continue;
     }
        
    // find optimal traceback
    bool found_trace = false;
    switch (mat) {
        case Match:
        {
            if (i + j == -first_diag) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedVectorMatrix::traceback] next cell is outside matrix, opening implied lead gap" << endl;
#endif
                // at top of matrix, move into implied lead gap along top edge
                mat = InsertCol;
                j--;
                in_lead_gap = true;
                // no where else to go, so break out of switch statement without checking alts
                break;
            }
			
            curr_score =  get_score_at(Match, i, j);

                
            int16_t match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + first_diag + j]]];
                
#ifdef debug_banded_aligner_traceback
            cerr << "[BandedVectorMatrix::traceback] transitioning from match, current score " << (int) curr_score << " match/mismatch score " << (int) match_score << " from node char " << j << " (" << node_seq[j] << ") and read char " << i + top_diag + j << " (" << read[i + top_diag + j] << ")" << endl;
#endif
                
            source_score = get_score_at(Match, i, j-1);
            score_diff = curr_score - (source_score + match_score);
            if (score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedVectorMatrix::traceback] found next cell in match matrix with score " << source_score  << endl;
#endif
                mat = Match;
                found_trace = true;
            }
                
            source_score = get_score_at(InsertRow, i, j-1);
            if (source_score > min_inf) {
                score_diff = curr_score - (source_score + match_score);
                if (!found_trace && score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BandedVectorMatrix::traceback] found next cell in insert row matrix with score " << (int) source_score << endl;
#endif
                    mat = InsertRow;
                    found_trace = true;
                }
            }
                
            source_score = get_score_at(InsertCol, i, j-1);
            if (source_score > min_inf) {
                score_diff = curr_score - (source_score + match_score);
                if (!found_trace && score_diff == 0) {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BandedVectorMatrix::traceback] found next cell in insert column matrix with score " << (int) source_score<< endl;
#endif
                    mat = InsertCol;
                    found_trace = true;
                }

            }
                
            if (!found_trace) {
                cerr << "error:[BandedVectorAligner] traceback stuck in match matrix interior" << endl;
                assert(0);
            }
                
            j--;
                
            break;
        }
                
        case InsertRow:
        {
            if (i == 0) {
                // along top of band
                cerr << "error:[BandedVectorAligner] traceback attempted to leave band from top" << endl;
                assert(0);
            }
                
            if (i + j == -first_diag) {
                // at top of matrix, move into implied lead gap along top edge
                in_lead_gap = true;
                // no where else to go, so break out of switch statement without checking alts
                i--;
                break;
            }
            curr_score = get_score_at(InsertRow, i, j);
             
            source_score = get_score_at(Match, i-1, j);
            score_diff = curr_score - (source_score - gap_open);
            if (score_diff == 0) {
                mat = Match;
                found_trace = true;
            }

            source_score = get_score_at(InsertRow, i-1, j);
            if (source_score > min_inf) {
                score_diff = curr_score - (source_score - gap_extend);
                if (!found_trace && score_diff == 0) {
                    mat = InsertRow;
                    found_trace = true;
                }
            }

            source_score = get_score_at(InsertCol, i-1, j);
            if (source_score > min_inf) {
                score_diff = curr_score - (source_score - gap_open);
                if (!found_trace && score_diff == 0) {
                    mat = InsertCol;
                    found_trace = true;
                }
            }
               
            if (!found_trace) {
                cerr << "error:[BandedVectorAligner] traceback stuck in insert row matrix interior" << endl;
                assert(0);
            }
               
            i--;
                
            break;
        }
            
        case InsertCol:
        {
            if (i == get_band_size() - 1) {
                // along bottom of band
                cerr << "error:[BandedVectorAligner] traceback attempted to leave band from bottom" << endl;
                assert(0);
            }

            curr_score =  get_score_at(InsertCol, i, j);

            source_score = get_score_at(Match, i+1, j-1);
            score_diff = curr_score - (source_score - gap_open);
            if (score_diff == 0) {
                mat = Match;
                found_trace = true;
            }
               
            source_score = get_score_at(InsertRow, i+1, j-1);
            if (source_score > min_inf) {
                score_diff = curr_score - (source_score - gap_open);
                if (!found_trace && score_diff == 0) {
                    mat = InsertRow;
                    found_trace = true;
                }
            }
               
            source_score = get_score_at(InsertCol, i+1, j-1);
            if (source_score > min_inf) {
                score_diff = curr_score - (source_score - gap_extend);
                if (!found_trace && score_diff == 0) {
                    mat = InsertCol;
                    found_trace = true;
                }
            }
                
            if (!found_trace) {
                cerr << "error:[BandedVectorAligner] traceback stuck in insert column matrix interior" << endl;
                assert(0);
            }
                
            i++;
            j--;
              
            break;
        }
      
        default:
        {
            cerr << "error:[BandedVectorAligner] unrecognized matrix type given to traceback" << endl;
            assert(0);
            break;
        }
    }

    
    if (in_lead_gap) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BandedVectorMatrix::traceback] running through node sequence in a lead gap" << endl;
#endif
        // add lead column gaps until reaching edge of node
        mat = InsertCol;
        while (j > 0) {
            builder.update_state(graph, mat, node, -1, j);
            j--;
            i++;
        }
    }
}


void BandedVectorMatrix::traceback_over_edge(const HandleGraph& graph, Alignment& alignment,
                                                                 VBBuilder& builder,
                                                                 int64_t& i, int64_t& j, matrix_t& mat,
                                                                 bool& in_lead_gap, int64_t& node_id,
                                                                 const int8_t* score_mat, const int8_t* nt_table,
                                                                 const int8_t gap_open, const int8_t gap_extend) {
    
    // begin POA across the boundary
    
    bool treat_as_source = false;
    
    BandedVectorMatrix* traceback_seed = nullptr;
    int64_t traceback_seed_row = std::numeric_limits<int64_t>::min();
    int64_t traceback_seed_col = std::numeric_limits<int64_t>::min();
    matrix_t traceback_mat = Match;
    
    const string& read = alignment.sequence();
    
    int64_t idx, next_idx;
    int16_t score_diff;
    int16_t curr_score;
    int16_t source_score;
    
    int64_t curr_diag = first_diag + i;
    string node_seq = graph.get_sequence(node);
    int64_t ncols = graph.get_length(node);
    int16_t min_inf = numeric_limits<int16_t>::min();
    
    bool found_trace = false;
    // if we traverse through nodes with no sequence, we need to keep track of which ones
    vector<BandedVectorMatrix*> empty_intermediate_nodes;
    vector<BandedVectorMatrix*> empty_source_path;
    
    // a queue of seeds and their empty predecessors
    vector<pair<BandedVectorMatrix*, vector<BandedVectorMatrix*>>> seed_queue;
    for (int idx = 0; idx < number_of_seeds; idx++) {
        seed_queue.emplace_back(seeds[i], vector<BandedVectorMatrix*>());
    }
    
    if (in_lead_gap) {
#ifdef debug_banded_aligner_traceback
        cerr << "[BandedVectorMatrix::traceback_over_edge] at boundary, following seed backward from a lead gap" << endl;
#endif
        // we are in the implied lead gap along the top of the matrix
        // take the shortest path back to the origin of the global alignment
        
        // add final read deletion of node
        builder.update_state(graph, mat, node, -1, j);
        
        while (!seed_queue.empty()) {
            
            auto seed_record = seed_queue.back();
            BandedVectorMatrix* seed = seed_record.first;
            seed_queue.pop_back();
            
            // is seed masked?
            if (seed == nullptr) {
                continue;
            }
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BandedVectorMatrix::traceback_over_edge] checking seed node " << graph.get_id(seed->node) << endl;
#endif
            
            // if this node is empty, add its predecessors to the queue
            if (graph.get_length(seed->node) == 0) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedVectorMatrix::traceback_over_edge] seed node " << graph.get_id(seed->node) << " is empty, checking predecessors" << endl;
#endif
                // record that this seed comes before its predecessors in the traceback
                seed_record.second.push_back(seed);
                if (seed->number_of_seeds==0) {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BandedVectorMatrix::traceback_over_edge] empty seed node " << graph.get_id(seed->node) << " is a source" << endl;
#endif
                    treat_as_source = true;
                    empty_source_path = seed_record.second;
                }
                
                for (int idx = 0; idx < seed->number_of_seeds; idx++){
                    seed_queue.push_back(make_pair(seed->seeds[i], seed_record.second));
                }
                continue;
            }
            
            score_diff = gap_extend * (seed->cumulative_seq_len + graph.get_length(seed->node) - cumulative_seq_len);
            if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedVectorMatrix::traceback_over_edge] found a lead gap traceback to node " << graph.get_id(seed->node) << endl;
#endif
                traceback_seed = seed;
                empty_intermediate_nodes = seed_record.second;
                found_trace = true;
            }
        }
        
        if (traceback_seed) {
            // where in the matrix is this?
            int64_t seed_ncols = graph.get_length(traceback_seed->node);
            int64_t seed_extended_top_diag = traceback_seed->first_diag + seed_ncols;
            traceback_seed_row = first_diag - seed_extended_top_diag + i + 1;
            traceback_seed_col = seed_ncols - 1;
        }
    }
    else {
        
        builder.update_state(graph, mat, node, curr_diag, 0);
        
        int16_t match_score;
        switch (mat) {
            case Match:
            {
                curr_score = get_score_at(Match, i*ncols, 0);
                match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + first_diag]]];

                break;
            }
                
            case InsertCol:
            {
                curr_score = get_score_at(InsertCol, i*ncols, 0);
                break;
            }
                
            case InsertRow:
            {
                cerr << "error:[BandedVectorAligner] traceback attempted to open row gap over node boundary" << endl;
                assert(0);
                break;
            }
                
            default:
            {
                cerr << "error:[BandedVectorAligner] unrecognized matrix type given to traceback" << endl;
                assert(0);
                break;
            }
        }
        
#ifdef debug_banded_aligner_traceback
        cerr << "[BandedVectorMatrix::traceback_over_edge] at boundary, in node " << graph.get_id(node) << " following seed backward from " << (mat == Match ? "match" : "insert column") << " matrix with score " << (int) curr_score << endl;
#endif
        
        // matches stay on same diagonal, column insertions move over one diagonal
        // note that the indexing is relative to THIS matrix, not the seed
        
        // check traceback goes to each seed matrix
        while (!seed_queue.empty()) {
            auto seed_record = seed_queue.back();
            BandedVectorMatrix* seed = seed_record.first;
            seed_queue.pop_back();
            
            // is the matrix masked?
            if (seed == nullptr) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedVectorMatrix::traceback_over_edge] seed is masked" << endl;
#endif
                continue;
            }
            
            if (graph.get_length(seed->node) == 0) {
                
                for (int idx = 0; idx < seed->number_of_seeds; idx++) {
                    seed_queue.push_back(make_pair(seed->seeds[i], seed_record.second));
                    seed_queue.back().second.push_back(seed);
                }
                
                if (seed->number_of_seeds == 0) {
                    treat_as_source = true;
                    // keep track of the path through empty nodes to a source
                    empty_source_path = seed_record.second;
                    empty_source_path.push_back(seed);
                }
                continue;
            }
            
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BandedVectorMatrix::traceback_over_edge] checking seed node " << graph.get_id(seed->node) << endl;
#endif
            
            int64_t seed_node_id = graph.get_id(seed->node);
            int64_t seed_ncols = graph.get_length(seed->node);
            
            // the diagonals in the current matrix that this seed extends to
            int64_t seed_extended_top_diag = seed->first_diag + seed_ncols;
            int64_t seed_extended_bottom_diag = seed->first_diag +num_diags+ seed_ncols;
            
            // does the traceback diagonal extend backward to this matrix?
            if (curr_diag > seed_extended_bottom_diag - (mat == InsertCol) // col inserts hit 1 less of band
                || curr_diag < seed_extended_top_diag) {
#ifdef debug_banded_aligner_traceback
                cerr << "[BandedVectorMatrix::traceback_over_edge] seed extended diags are top: " << seed_extended_top_diag << ", bottom: " << seed_extended_bottom_diag << " and curr mat is " << (mat == InsertCol ? "" : "not") << " insert column, so we cannot extend from this seed to the current diag " << curr_diag << endl;
#endif
                continue;
            }
            
            int64_t seed_col = seed_ncols - 1;
            int64_t seed_row = -(seed_extended_top_diag - first_diag) + i + (mat == InsertCol);
            next_idx = seed_row * seed_ncols + seed_col;
            
#ifdef debug_banded_aligner_traceback
            cerr << "[BandedVectorMatrix::traceback_over_edge] checking seed rectangular coordinates (" << seed_row << ", " << seed_col << "), with indices calculated from current diagonal " << curr_diag << " (top diag " << top_diag << " + offset " << i << "), seed top diagonal " << seed->top_diag << ", seed seq length " << seed_ncols << " with insert column offset " << (mat == InsertCol) << endl;
#endif
            
            switch (mat) {
                case Match:
                {
#ifdef debug_banded_aligner_traceback
                    cerr << "[BandedVectorMatrix::traceback_over_edge] poa backwards from match, seed extended top diag " << seed_extended_top_diag << endl;
#endif
                    // does match lead into a lead row gap?
                    if (seed->first_diag + seed_row + seed_col == -1) {
#ifdef debug_banded_aligner_traceback
                        cerr << "[BandedVectorMatrix::traceback_over_edge] traceback points to a lead column gap of length " << seed->cumulative_seq_len + seed_ncols << " with score " << (int) -gap_open - (seed->cumulative_seq_len + seed_ncols - 1) * gap_extend << " extending to score here of " << (int) curr_score << " with match score " << (int) match_score << endl;
#endif
                        // score of implied column gap
                        source_score = -gap_open - (seed->cumulative_seq_len + seed_ncols - 1) * gap_extend;
                        score_diff = curr_score - (source_score + match_score);
                        if (score_diff == 0 && !found_trace) {
                            traceback_mat = InsertCol;
                            traceback_seed = seed;
                            traceback_seed_row = seed_row;
                            traceback_seed_col = seed_col;
                            in_lead_gap = true;
                            found_trace = true;
                            empty_intermediate_nodes = seed_record.second;
#ifdef debug_banded_aligner_traceback
                            cerr << "[BandedVectorMatrix::traceback_over_edge] hit found in lead gap with score " << -gap_open - (seed->cumulative_seq_len + seed_ncols - 1) * gap_extend << endl;
#endif
                        }
                        
                        // don't check any of the matrices because they will have garbage values in this position or seg fault
                        break;
                    }
                    
                    source_score = seed->get_score_at(Match, seed_row*seed_ncols, seed_col);
                    // don't need to check edge condition because match does not have min inf
                    score_diff = curr_score - (source_score + match_score);
                    if (score_diff == 0 && !found_trace) {
                        traceback_mat = Match;
                        traceback_seed = seed;
                        traceback_seed_row = seed_row;
                        traceback_seed_col = seed_col;
                        found_trace = true;
                        empty_intermediate_nodes = seed_record.second;

                    }
                    
                    source_score = seed->get_score_at(InsertCol, seed_row*seed_ncols, seed_col);
                    // check edge condition
                    if (source_score > min_inf) {
                        score_diff = curr_score - (source_score + match_score);
                        if (score_diff == 0 && !found_trace) {

                            traceback_mat = InsertCol;
                            traceback_seed = seed;
                            traceback_seed_row = seed_row;
                            traceback_seed_col = seed_col;
                            found_trace = true;
                            empty_intermediate_nodes = seed_record.second;
                        }
                    }
                    
                    source_score = seed->get_score_at(InsertRow, seed_row*seed_ncols, seed_col);
                    // check edge condition
                    if (source_score > min_inf) {
                        score_diff = curr_score - (source_score + match_score);
                        if (score_diff == 0 && !found_trace) {

                            traceback_mat = InsertRow;
                            traceback_seed = seed;
                            traceback_seed_row = seed_row;
                            traceback_seed_col = seed_col;
                            found_trace = true;
                            empty_intermediate_nodes = seed_record.second;
                        }
                    }
                    
                    break;
                }
                    
                case InsertCol:
                {
                    source_score = seed->get_score_at(Match, seed_row*seed_ncols, seed_col);
                    // don't need to check edge condition because match does not have min inf
                    score_diff = curr_score - (source_score - gap_open);
                    if (score_diff == 0 && !found_trace) {

                        traceback_mat = Match;
                        traceback_seed = seed;
                        traceback_seed_row = seed_row;
                        traceback_seed_col = seed_col;
                        found_trace = true;
                        empty_intermediate_nodes = seed_record.second;
                    }
                    
                    source_score = seed->get_score_at(InsertCol, seed_row*seed_ncols, seed_col);
                    // check edge condition
                    if (source_score > min_inf) {
                        score_diff = curr_score - (source_score - gap_extend);
                        if (score_diff == 0 && !found_trace) {
                            traceback_mat = InsertCol;
                            traceback_seed = seed;
                            traceback_seed_row = seed_row;
                            traceback_seed_col = seed_col;
                            found_trace = true;
                            empty_intermediate_nodes = seed_record.second;
                        }
                    }
                    
                    source_score = seed->get_score_at(InsertRow, seed_row*seed_ncols, seed_col);
                    // check edge condition
                    if (source_score > min_inf) {
                        score_diff = curr_score - (source_score - gap_open);
                        if (score_diff == 0 && !found_trace) {
                            traceback_mat = InsertRow;
                            traceback_seed = seed;
                            traceback_seed_row = seed_row;
                            traceback_seed_col = seed_col;
                            found_trace = true;
                            empty_intermediate_nodes = seed_record.second;
                        }

                    }
                    
                    break;
                }
                    
                case InsertRow:
                {
                    cerr << "error:[BandedVectorAligner] illegal matrix type for moving across node boundary" << endl;
                    assert(0);
                    break;
                }
                    
                default:
                {
                    cerr << "error:[BandedVectorAligner] unrecognized matrix type given to traceback" << endl;
                    assert(0);
                    break;
                }
            }
        }
    }
    
    bool found_source_trace = false;
    
    if (treat_as_source) {
        // this is a source node, or it is connected to one by a zero-length path
        
#ifdef debug_banded_aligner_traceback
        cerr << "[BandedVectorMatrix::traceback_over_edge] at beginning of first node in alignment" << endl;
#endif
        if (in_lead_gap && !found_trace) {
            // this will always be the shortest gap
            found_source_trace = true;
            j--;
            i++;
        }
        else {
            switch (mat) {
                case Match:
                {
                    int16_t match_score;

                    match_score = score_mat[5 * nt_table[node_seq[j]] + nt_table[read[i + first_diag + j]]];
                    
                    source_score = curr_diag > 0 ? -gap_open - (curr_diag - 1) * gap_extend : 0;
                    score_diff = curr_score - (source_score + match_score);
                    if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                        cerr << "[BandedVectorMatrix::traceback_over_edge] alignment starts with match, adding read char " << first_diag + i << ": " << read[top_diag + i] << endl;
#endif
                        mat = InsertRow;
                        found_source_trace = true;
                    }
                    j--;
                    break;
                }
                    
                case InsertCol:
                {
                    source_score = -gap_open - (i + first_diag) * gap_extend;
                    score_diff = curr_score - (source_score - gap_open);
                    if (score_diff == 0 && !found_trace) {
#ifdef debug_banded_aligner_traceback
                        cerr << "[BandedVectorMatrix::traceback_over_edge] alignment starts with column gap" << endl;
#endif
                        mat = InsertRow;
                        found_source_trace = true;
                    }
                    j--;
                    i++;
                    break;
                }
                    
                default:
                {
                    cerr << "error:[BandedVectorAligner] invalid matrix type for final traceback column" << endl;
                    assert(0);
                    break;
                }
            }
        }
    }
    
    
    if (found_source_trace) {
        // if we traversed any empty nodes before finding the traceback, add empty updates for them
        for (BandedVectorMatrix* seed_path_node : empty_source_path) {
            builder.update_state(graph, InsertCol, seed_path_node->node, i + first_diag, 0, true);
        }
        
        // add any lead row gaps necessary
        while (first_diag + i > 0) {
#ifdef debug_banded_aligner_traceback
            cerr << "[BandedVectorMatrix::traceback_over_edge] initial row gaps are present, adding read char " << top_diag + i - 1 << ": " << read[top_diag + i - 1] << endl;
#endif
            i--;
            const handle_t& end_node = empty_source_path.empty() ? node : empty_source_path.back()->node;
            builder.update_state(graph, InsertRow, end_node, i + first_diag, -1, graph.get_length(end_node) == 0);
        }
        
        // set the node ID to 0 to indicate completion
        node_id = 0;
        return;
    }
    else {
#ifdef debug_banded_aligner_traceback
        cerr << "[BAMatrix::traceback_over_edge] traversed " << empty_intermediate_nodes.size() << " empty nodes before finding trace" << endl;
#endif
        // if we traversed any empty nodes before finding the traceback, add empty updates for them
        for (BandedVectorMatrix* intermediate_node : empty_intermediate_nodes) {
            builder.update_state(graph, mat, intermediate_node->node, i + first_diag, 0, true);
        }
        
    }
    
    if (!found_trace) {
        cerr << "error:[BandedVectorAligner] traceback stuck at node boundary" << endl;
        assert(0);
    }
    
    // set the traceback values for the next matrix
    i = traceback_seed_row;
    j = traceback_seed_col;
    mat = traceback_mat;
    node_id = graph.get_id(traceback_seed->node);
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
        handle_t node, int64_t first_diag, int64_t num_diags,
        int64_t num_cols, int64_t cu_seq_len){

    matrix.num_cols = num_cols;
    matrix.first_diag = first_diag;
    matrix.num_diags = num_diags;
    matrix.node = node;
    matrix.node_seq = instance->graph->get_sequence(node);
	matrix.cumulative_seq_len = cu_seq_len;

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
                matrix.seeds[num_seeds] = &(instance->matrices[instance->handle_to_idx.at(prev)]);
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
