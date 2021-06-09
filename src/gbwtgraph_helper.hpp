#ifndef VG_GBWTGRAPH_HELPER_HPP_INCLUDED
#define VG_GBWTGRAPH_HELPER_HPP_INCLUDED

/** \file 
 * Utility classes and functions for working with GBWTGraph.
 */

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

namespace vg {

//------------------------------------------------------------------------------

/*
    These are the proper ways of saving and loading GBWTGraph structures.
    Loading with `vg::io::VPKG::load_one` is also supported.
*/

/// Load GBWTGraph from the file.
/// NOTE: Call `graph.set_gbwt()` afterwards with the appropriate GBWT index.
void load_gbwtgraph(gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Load GBZ from the file.
void load_gbz(gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress = false);

/// Load GBZ from separate GBWT / GBWTGraph files.
void load_gbz(gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress = false);

/// Load GBWT and GBWTGraph from the GBZ file.
void load_gbz(gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Load a minimizer index from the file.
void load_minimizer(gbwtgraph::DefaultMinimizerIndex& index, const std::string& filename, bool show_progress = false);

/// Save GBWTGraph to the file.
void save_gbwtgraph(const gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Save GBZ to the file.
void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& filename, bool show_progress = false);

/// Save GBWT and GBWTGraph to the GBZ file.
void save_gbz(const gbwt::GBWT& index, gbwtgraph::GBWTGraph& graph, const std::string& filename, bool show_progress = false);

/// Save GBZ to separate GBWT / GBWTGraph files.
void save_gbz(const gbwtgraph::GBZ& gbz, const std::string& gbwt_name, const std::string& graph_name, bool show_progress = false);

/// Save a minimizer index to the file.
void save_minimizer(const gbwtgraph::DefaultMinimizerIndex& index, const std::string& filename, bool show_progress = false);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GBWTGRAPH_HELPER_HPP_INCLUDED
