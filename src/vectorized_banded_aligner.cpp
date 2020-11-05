/**
 * \file vectorized_banded_aligner.cpp
 *
 * Implements a the vectorized banded global aligner
 *
 */

#include "vectorized_banded_aligner.hpp"

namespace vg {

const size_t BandedVectorHeap::init_block_size = (1 << 20); // 1 MB

BandedVectorAligner* BandedVectorAligner::init() {
    
    BandedVectorHeap* heap = BandedVectorHeap::init();
    BandedVectorAligner* aligner = (BandedVectorAligner*) heap->irreversible_alloc(sizeof(BandedVectorAligner));
    aligner->heap = heap;
    return aligner;
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



}
