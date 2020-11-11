/// \file vectorized_banded_aligner
/// 
/// Unit tests for vectorized banded aligner
///

#include <iostream>
#include "../vectorized_banded_aligner.cpp"
#include "catch.hpp"

namespace vg{
namespace unittest{

TEST_CASE("vectorized matrix of size 0", "[vector-banded]"){
    BandedVectorHeap* heap=BandedVectorHeap::init();
    BandedVectorMatrix M_test;
    init_BandedVectorMatrix(heap,M_test, 0,0,0);
    alignas(16) short test_arr[8];
    _mm_store_si128((__m128i*)test_arr, M_test.vectors[0].match);

    REQUIRE(M_test.get_band_size() == 1);
    REQUIRE(M_test.vector_index(0, 0) == 1);
    REQUIRE(M_test.index_within_vector(0,0) ==1);
    REQUIRE(test_arr[0]==1);
    REQUIRE(test_arr[1]==1);
    REQUIRE(test_arr[2]==1);
    REQUIRE(test_arr[3]==1);
    REQUIRE(test_arr[4]==1);
    REQUIRE(test_arr[5]==1);
    REQUIRE(test_arr[6]==1);
    REQUIRE(test_arr[7]==1);
    heap->destroy();
}

TEST_CASE("vectorized matrix of size 6x1(num_diags x num_cols)", "[vector-banded]"){
    BandedVectorHeap* heap = BandedVectorHeap::init();
    BandedVectorMatrix test;
    init_BandedVectorMatrix(heap, test, 0, 6, 1);
    alignas(16) short test_arr[8];
    _mm_store_si128((__m128i*)test_arr, test.vectors[1].match);

    REQUIRE(test.get_band_size() == 2);
    REQUIRE(test.vector_index(5,0)==1);
    REQUIRE(test.index_within_vector(0,0)==1);
    REQUIRE(test_arr[0]==1);
    REQUIRE(test_arr[1]==0);
    REQUIRE(test_arr[2]==0);
    REQUIRE(test_arr[3]==0);
    REQUIRE(test_arr[4]==0);
    REQUIRE(test_arr[5]==0);
    REQUIRE(test_arr[6]==0);
    REQUIRE(test_arr[7]==0);
    heap->destroy();
}

TEST_CASE("vectorized matrix of size 6x10(num_diags x num_cols)", "[vector-banded]"){
    BandedVectorHeap* heap = BandedVectorHeap::init();
    BandedVectorMatrix test;
    init_BandedVectorMatrix(heap, test, 0, 6, 10);
    alignas(16) short test_arr[8];
    _mm_store_si128((__m128i*)test_arr, test.vectors[1].match);

    REQUIRE(test.get_band_size() == 11);
    REQUIRE(test.vector_index(8,6)==7);
    REQUIRE(test.index_within_vector(8,6)==3);
    REQUIRE(test_arr[0]==1);
    REQUIRE(test_arr[1]==0);
    REQUIRE(test_arr[2]==0);
    REQUIRE(test_arr[3]==0);
    REQUIRE(test_arr[4]==0);
    REQUIRE(test_arr[5]==0);
    REQUIRE(test_arr[6]==0);
    REQUIRE(test_arr[7]==0);
    heap->destroy();

}

TEST_CASE("vectorized matrix of size 40x20(num_diags x num_cols)","[vector-banded]"){
    BandedVectorHeap* heap = BandedVectorHeap::init();
    BandedVectorMatrix test;
    init_BandedVectorMatrix(heap, test, 0, 40, 20);
    alignas(16) short test_arr[8];
    _mm_store_si128((__m128i*)test_arr, test.vectors[21].match);

    REQUIRE(test.get_band_size()==105);
    REQUIRE(test.vector_index(10,10)==11);
    REQUIRE(test.index_within_vector(10,10)==1);

    //vector at index 21 should be all buffer
    REQUIRE(test_arr[0]==1);
    REQUIRE(test_arr[7]==1);

    //vector at index 22 should be all zeroes
    _mm_store_si128((__m128i*)test_arr, test.vectors[22].match);
    REQUIRE(test_arr[0]==0);
    REQUIRE(test_arr[7]==0);

    heap->destroy();
}

TEST_CASE("vectorized matrix with non-zero first_diag", "[vector-banded]"){
    BandedVectorHeap* heap = BandedVectorHeap::init();
    BandedVectorMatrix test;
    init_BandedVectorMatrix(heap, test, -4, 20, 20);
    alignas(16) short test_arr[8];
    _mm_store_si128((__m128i*)test_arr, test.vectors[1].match);

    REQUIRE(test.vector_index(4,0)==22);
    REQUIRE(test.index_within_vector(2,0)==7);

    //vector at index 1 should have taper
    REQUIRE(test_arr[0]==1);
    REQUIRE(test_arr[1]==1);
    REQUIRE(test_arr[2]==1);
    REQUIRE(test_arr[3]==1);
    REQUIRE(test_arr[4]==1);
    REQUIRE(test_arr[5]==0);
    REQUIRE(test_arr[6]==0);
    REQUIRE(test_arr[7]==0);

    _mm_store_si128((__m128i*)test_arr, test.vectors[2].match);
    //vector at index 2 should taper
    REQUIRE(test_arr[0]==1);
    REQUIRE(test_arr[1]==1);
    REQUIRE(test_arr[2]==1);
    REQUIRE(test_arr[3]==1);
    REQUIRE(test_arr[4]==0);
    REQUIRE(test_arr[5]==0);
    REQUIRE(test_arr[6]==0);
    REQUIRE(test_arr[7]==0);

    heap->destroy();
}

}
}
