/// \file minimizer_mapper.cpp
///  
/// unit tests for the minimizer mapper

#include <iostream>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "../minimizer_mapper.hpp"
#include "../build_index.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

// We define a child class to expose protected stuff for testing
class TestMinimizerMapper : public MinimizerMapper {
public:
    TestMinimizerMapper(
        gbwtgraph::GBWTGraph gbwt_graph,
        gbwtgraph::DefaultMinimizerIndex minimizer_index,
        MinimumDistanceIndex distance_index,
        PathPositionHandleGraph* handle_graph) 
        : MinimizerMapper(gbwt_graph, minimizer_index, distance_index, handle_graph){};
    using MinimizerMapper::MinimizerMapper;
    using MinimizerMapper::score_extension_group;
    using MinimizerMapper::Minimizer;
    using MinimizerMapper::fragment_length_distr;
};

TEST_CASE("MinimizerMapper::score_extension_group works", "[giraffe][mapping]") {

    // Define an Alignment to score extensions against.
    Alignment aln;
    aln.set_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    
    // Define a bunch of fake GaplessExtensions to chain up and score
    vector<GaplessExtension> to_score;
    
    SECTION("Score of no gapless extensions is 0") {
        REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 0);
    }
    
    SECTION("One-base extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 1;
        to_score.back().read_interval.second = 2;
        to_score.back().score = 1;
        
        SECTION("Score of a 1-base gapless extension is 1") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 1);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 2;
        to_score.back().read_interval.second = 3;
        to_score.back().score = 1;
        
        SECTION("Score of two adjacent 1-base gapless extensions is 2") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 2);
        }
        
    }
    
    SECTION("Longer extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 1;
        to_score.back().read_interval.second = 10;
        to_score.back().score = 9;
        
        SECTION("Score of one 9-base extension is 9") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 9);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 10;
        to_score.back().read_interval.second = 19;
        to_score.back().score = 9;
        
        SECTION("Score of two 9-base extensions abutting is 18") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 18);
        }
        
        to_score.back().read_interval.first++;
        to_score.back().read_interval.second++;
        
        SECTION("Score of two 9-base extensions separated by a 1-base gap is 12") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 12);
        }
        
        to_score.back().read_interval.first++;
        to_score.back().read_interval.second++;
        
        SECTION("Score of two 9-base extensions separated by a 2-base gap is 11") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 11);
        }
        
        to_score.back().read_interval.first -= 3;
        to_score.back().read_interval.second -= 3;
        
        SECTION("Score of two 9-base extensions with a 1-base ovealap is 12") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 11);
        }
    }
    
    SECTION("Many possibly overlapping extensions work") {
    
        to_score.emplace_back();
        to_score.back().read_interval.first = 0;
        to_score.back().read_interval.second = 1;
        to_score.back().score = 1;
    
        for (size_t i = 0; i < 35; i++) {
        
            to_score.emplace_back();
            to_score.back().read_interval.first = i + 1;
            to_score.back().read_interval.second = i + 1 + 30;
            to_score.back().score = 30;
        
        }
    
        
        
        SECTION("Score of one 1-base extensions and 2 30-base extensions is 61") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == 61);
        }
        
        to_score.emplace_back();
        to_score.back().read_interval.first = 28;
        to_score.back().read_interval.second = 28 + 45;
        to_score.back().score = 45;
        
        // Sort by read interval as is required
        std::sort(to_score.begin(), to_score.end(), [](const GaplessExtension& a, const GaplessExtension& b) {
            return (a.read_interval.first < b.read_interval.first) ||
                (a.read_interval.first == b.read_interval.first && a.read_interval.second < b.read_interval.second);
        });
        
        SECTION("Score of one 1-base extension, a 30-base extension, a backtrack of 4, and a 45-base extension is correct") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == (31 + 45 - 6 - 3));
        }
        
        
        for (size_t i = 3; i < 29; i++) {
             to_score.emplace_back();
            to_score.back().read_interval.first = i;
            to_score.back().read_interval.second = i + 15;
            to_score.back().score = 15;
        }
        
        // Sort by read interval as is required
        std::sort(to_score.begin(), to_score.end(), [](const GaplessExtension& a, const GaplessExtension& b) {
            return (a.read_interval.first < b.read_interval.first) ||
                (a.read_interval.first == b.read_interval.first && a.read_interval.second < b.read_interval.second);
        });
        
        SECTION("Score of one 1-base extension, a 30-base extension, a backtrack of 4, and a 45-base extension is not distracted") {
            REQUIRE(TestMinimizerMapper::score_extension_group(aln, to_score, 6, 1) == (31 + 45 - 6 - 3));
        }
    
    }
}

TEST_CASE("Fragment length distribution gets reasonable value", "[giraffe][mapping]") {

        vector<int64_t> distances { 27, 69, 76, 88, 107, 114, 119, 121, 124, 124, 125, 125, 125, 125, 126, 126, 126, 127, 127, 127, 127, 127, 127, 
            134, 134, 134, 134, 135, 135, 136, 136, 136, 136, 136, 136, 136, 136, 137, 137, 137, 138, 139, 139, 139, 140, 140, 140, 140, 140, 140, 
            144, 144, 144, 145, 145, 145, 145, 145, 145, 145, 145, 145, 146, 146, 146, 146, 146, 147, 147, 147, 147, 147, 148, 148, 148, 148, 148, 
            154, 154, 154, 154, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 156, 156, 156, 156, 156, 156, 156, 156, 
            160, 160, 160, 160, 160, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 162, 162, 162, 162, 163, 163, 163, 
            167, 167, 167, 167, 167, 167, 167, 168, 168, 168, 168, 169, 169, 169, 169, 169, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 171, 
            174, 174, 174, 174, 174, 175, 175, 175, 175, 175, 175, 176, 176, 176, 176, 176, 176, 176, 177, 178, 178, 178, 178, 178, 178, 178, 178,
            182, 182, 182, 182, 182, 182, 182, 182, 183, 183, 183, 183, 183, 183, 183, 183, 183, 184, 184, 184, 184, 184, 185, 185, 185, 185, 185,
            189, 189, 189, 189, 189, 189, 189, 189, 189, 189, 189, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 190, 191, 191, 191, 
            195, 195, 195, 195, 196, 196, 196, 196, 196, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 198, 198, 198, 198, 198, 198, 198, 
            202, 202, 202, 202, 202, 202, 202, 202, 202, 203, 203, 203, 203, 203, 203, 203, 203, 204, 204, 204, 204, 204, 204, 204, 204, 204, 
            207, 207, 208, 208, 208, 208, 208, 208, 208, 208, 208, 209, 209, 209, 209, 209, 209, 210, 210, 210, 210, 210, 210, 210, 211, 211, 
            214, 214, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 215, 216, 216, 216, 216, 216, 216, 216, 216, 216, 217, 217, 217, 217, 
            220, 220, 220, 220, 220, 220, 220, 220, 220, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 221, 222, 222, 222, 222, 222, 
            226, 226, 226, 226, 226, 226, 226, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 227, 228, 228, 228, 229, 229, 229, 229, 
            233, 233, 233, 233, 233, 233, 233, 233, 234, 234, 234, 234, 234, 234, 234, 234, 234, 234, 234, 235, 235, 235, 235, 235, 236, 236, 
            238, 238, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 239, 240, 240, 240, 240, 240, 240, 240, 241, 241, 
            244, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 245, 246, 246, 246, 246, 246, 246, 246, 246, 246, 246, 246, 247, 247, 247, 
            249, 250, 250, 250, 250, 250, 250, 250, 250, 251, 251, 251, 251, 251, 251, 251, 252, 252, 252, 252, 252, 252, 252, 252, 252, 252, 
            255, 255, 255, 255, 255, 255, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 256, 257, 257, 257, 257, 257, 257, 257, 257, 
            261, 261, 261, 261, 261, 262, 262, 262, 262, 262, 262, 262, 262, 262, 262, 262, 263, 263, 263, 263, 264, 264, 264, 264, 264, 264, 
            268, 268, 268, 268, 268, 269, 269, 269, 269, 269, 269, 269, 269, 269, 269, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 270, 
            274, 275, 275, 275, 276, 276, 276, 276, 276, 276, 276, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 278, 278, 278, 278, 
            281, 281, 281, 281, 281, 281, 281, 282, 282, 282, 282, 282, 282, 282, 282, 283, 283, 283, 283, 283, 283, 283, 283, 284, 284, 284, 
            288, 288, 288, 288, 288, 289, 289, 289, 289, 289, 289, 289, 290, 290, 290, 291, 291, 291, 291, 291, 291, 291, 291, 291, 292, 292, 
            295, 296, 296, 296, 296, 296, 296, 296, 296, 297, 297, 297, 297, 297, 297, 297, 298, 298, 298, 298, 298, 298, 298, 298, 298, 299,
            302, 302, 302, 302, 302, 303, 303, 303, 303, 303, 303, 303, 303, 304, 304, 304, 304, 304, 304, 304, 305, 305, 305, 305, 306, 306,
            312, 312, 312, 312, 312, 313, 313, 313, 313, 313, 314, 314, 314, 314, 314, 315, 315, 315, 315, 315, 315, 316, 316, 316, 316, 316,
            320, 320, 321, 321, 321, 321, 321, 321, 321, 322, 322, 322, 322, 322, 322, 323, 323, 323, 323, 323, 324, 324, 324, 324, 325, 325,
            332, 333, 333, 333, 333, 333, 333, 333, 333, 333, 334, 334, 334, 334, 335, 335, 335, 335, 335, 336, 336, 336, 336, 336, 337, 337,
            342, 343, 344, 344, 344, 344, 344, 344, 345, 345, 346, 346, 346, 346, 346, 347, 347, 347, 347, 348, 348, 348, 348, 348, 348, 348,
            356, 356, 357, 357, 358, 358, 359, 359, 359, 359, 361, 361, 361, 361, 362, 362, 362, 363, 363, 363, 363, 364, 364, 364, 364, 365,
            374, 374, 374, 375, 376, 376, 376, 377, 377, 377, 378, 378, 379, 379, 379, 379, 380, 380, 380, 380, 381, 381, 382, 382, 383, 384,
            397, 398, 398, 398, 399, 399, 399, 399, 399, 399, 400, 400, 402, 402, 402, 403, 403, 403, 403, 404, 404, 404, 404, 406, 407, 407,
            463, 465, 466, 466, 467, 470, 471, 473, 474, 479, 479, 480, 481, 482, 483, 485, 491, 493, 493, 495, 496, 497, 498, 501, 507, 511,
            512, 512, 515, 519, 521, 521, 521, 523, 524, 527, 531, 535, 537, 541, 542, 550, 556, 557, 557, 559, 561, 569, 572, 573, 575, 579,
            580, 585, 589, 622, 627, 647, 667, 715, 757, 780, 1302, 4224, 10520, 16912, 17605, 17773, 18141, 18234, 19908, 22806, 26071, 
            33167, 39460, 39642, 55666, 59773, 63297, 74729, 82293, 84261, 103051, 125968, 126638, 133620, 134000, 156120, 156834, 158566, 159945,
            163617, 168131, 170576, 186151, 187063, 196981, 199264, 205006, 211618, 214498, 232698, 241394, 242280, 253799, 254397, 257196, 261598,
            316979, 329457, 654834,
            332, 333, 333, 333, 333, 333, 333, 333, 333, 333, 334, 334, 334, 334, 335, 335, 335, 335, 335, 336, 336, 336, 336, 336, 337, 337,
            342, 343, 344, 344, 344, 344, 344, 344, 345, 345, 346, 346, 346, 346, 346, 347, 347, 347, 347, 348, 348, 348, 348, 348, 348, 348,
            356, 356, 357, 357, 358, 358, 359, 359, 359, 359, 361, 361, 361, 361, 362, 362, 362, 363, 363, 363, 363, 364, 364, 364, 364, 365,
            274, 275, 275, 275, 276, 276, 276, 276, 276, 276, 276, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 277, 278, 278, 278, 278 };


        //Make an empty minimizer mapper just to get the fragment length distr
        gbwtgraph::GBWTGraph gbwt_graph;
        gbwtgraph::DefaultMinimizerIndex minimizer_index;
        MinimumDistanceIndex distance_index;
        PathPositionHandleGraph* handle_graph;
        TestMinimizerMapper test_mapper (gbwt_graph, minimizer_index, distance_index, handle_graph);
        for (int64_t dist : distances) {
            if (dist <= test_mapper.max_fragment_length) {
                test_mapper.fragment_length_distr.register_fragment_length(dist);
            }
        }

        SECTION("Make sure we have a reasonable distribution") {
            //Distribution should ignore outliers
            REQUIRE(test_mapper.fragment_length_distr.std_dev() <= 400);
        }
}

class TestableMinimizerMapper : public MinimizerMapper {
public:
    using MinimizerMapper::Minimizer;
    using MinimizerMapper::faster_cap;
};

TEST_CASE("Mapping quality cap cannot be confused by excessive Gs", "[giraffe][mapping]") {
    string sequence;
    string quality;
    for (size_t i = 0; i < 150; i++) {
        sequence.push_back('G');
        quality.push_back((char)0x1E);
    }
    
    // Cover the read in 25bp cores with 10bp flanks on each side
    int core_width = 25;
    int flank_width = 10;
    vector<TestableMinimizerMapper::Minimizer> minimizers;
    // They are all going to be explored
    vector<size_t> minimizers_explored;
    
    string min_seq;
    for (int i = 0; i < core_width; i++) {
        min_seq.push_back('G');
    }
    auto encoded = gbwtgraph::DefaultMinimizerIndex::key_type::encode(min_seq);
    
    for (int core_start = 0; core_start + core_width < sequence.size(); core_start++) {
        minimizers_explored.push_back(minimizers.size());
        minimizers.emplace_back();
        TestableMinimizerMapper::Minimizer& m = minimizers.back();
        
        if (core_start <= flank_width) {
            // Partial left flank
            m.agglomeration_start = 0;
            m.agglomeration_length = core_width + flank_width + core_start;
        } else if (sequence.size() - core_start - core_width <= flank_width) {
            // Partial right flank
            m.agglomeration_start = core_start - flank_width;
            m.agglomeration_length = sequence.size() - m.agglomeration_start - 1;
        } else {
            // Full flanks
            m.agglomeration_start = core_start - flank_width;
            m.agglomeration_length = core_width + flank_width * 2;
        }
        
        // We need to set the key and its hash
        m.value.key = encoded;
        m.value.hash = m.value.key.hash();
        m.value.offset = core_start;
        m.value.is_reverse = false;
        
        m.hits = 229;
        // We knowe the occurrences won't be used.
        m.occs = nullptr;
        m.length = core_width;
        m.candidates_per_window = flank_width + 1;
        m.score = 1;
    }
    
    
    // Compute the MAPQ cap
    double cap = TestableMinimizerMapper::faster_cap(minimizers, minimizers_explored, sequence, quality);
    
    // The MAPQ cap should not be infinite.
    REQUIRE(!isinf(cap));
}





}

}

