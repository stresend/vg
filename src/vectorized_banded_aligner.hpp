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
#include<iostream>

namespace vg {

using namespace std;

// forward declaration
struct BandedVectorHeap;

/*
 * Core aligner object for vectorized banded global alignment
 */
struct BandedVectorAligner {
    
    // constructor
    static BandedVectorAligner* init();
    // destructor
    void destroy();
    
    // TODO: design interface
    
private:
    
    // TODO: design members
    
    BandedVectorHeap* heap;
    
    
public:
    // delete all constructors, assignents operators, and the destructor
    BandedVectorAligner() = delete;
    ~BandedVectorAligner() = delete;
    BandedVectorAligner(const BandedVectorAligner& other) = delete;
    BandedVectorAligner(BandedVectorAligner&& other) = delete;
    BandedVectorAligner& operator=(BandedVectorAligner&& other) = delete;
    BandedVectorAligner& operator=(const BandedVectorAligner& other) = delete;
};

/*
 * A heap that allocates 16-byte aligned memory pointers and allows
 * collective freeing of many pointers simultaneously
 */
struct BandedVectorHeap {
    
public:
    
    // constructor
    static BandedVectorHeap* init();
    
    // invalidate all existing pointers from this heap, and free up
    // the memory to be reallocated
    void reset();
    
    // allocate an aligned block of the size indicated.
    void* alloc(size_t size);
    
    // allocate an aligned block of the size indicated, which will not
    // be free'd when calling reset(). may make any allocation since the
    // last time reset() was called permanant as well.
    void* irreversible_alloc(size_t size);
    
    // destructor
    void destroy();
    
private:
    
    static const size_t init_block_size;
    
    struct BandedVectorMemoryBlock {
    public:
        static BandedVectorMemoryBlock* init(size_t size);
        
        size_t remaining() const;
        size_t size() const;
        void* alloc(size_t size);
        void* irreversible_alloc(size_t size);
        void reset();
        BandedVectorMemoryBlock* next_block();
        BandedVectorMemoryBlock* find_space(size_t size);
        
        void destroy();
        
    private:
        static size_t round_to_align(size_t size);
        
        uint8_t* block_begin;
        uint8_t* bottom;
        uint8_t* top;
        uint8_t* curr;
        BandedVectorMemoryBlock* next;
    };
    
    BandedVectorMemoryBlock* head;
    BandedVectorMemoryBlock* curr;
    
public:
    // delete all constructors, assignents operators, and the destructor
    BandedVectorHeap() = delete;
    ~BandedVectorHeap() = delete;
    BandedVectorHeap(const BandedVectorHeap& other) = delete;
    BandedVectorHeap(BandedVectorHeap&& other) = delete;
    BandedVectorHeap& operator=(BandedVectorHeap&& other) = delete;
    BandedVectorHeap& operator=(const BandedVectorHeap& other) = delete;
};



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
