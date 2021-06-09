#ifndef VG_GBWT_HELPER_HPP_INCLUDED
#define VG_GBWT_HELPER_HPP_INCLUDED

/** \file 
 * Utility classes and functions for working with GBWT.
 */

#include <vector>

#include "position.hpp"

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/fast_locate.h>
#include <handlegraph/mutable_path_handle_graph.hpp>

namespace vg {

std::vector<std::string> parseGenotypes(const std::string& vcf_line, size_t num_samples);

//------------------------------------------------------------------------------

/// Convert gbwt::node_type to handle_t.
inline handle_t gbwt_to_handle(const HandleGraph& graph, gbwt::node_type node) {
    return graph.get_handle(gbwt::Node::id(node), gbwt::Node::is_reverse(node));
}

/// Convert gbwt::node_type and an offset as size_t to pos_t.
inline pos_t gbwt_to_pos(gbwt::node_type node, size_t offset) {
    return make_pos_t(gbwt::Node::id(node), gbwt::Node::is_reverse(node), offset);
}

/// Convert handle_t to gbwt::node_type.
inline gbwt::node_type handle_to_gbwt(const HandleGraph& graph, handle_t handle) {
    return gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle));
}

/// Extract gbwt::node_type from pos_t.
inline gbwt::node_type pos_to_gbwt(pos_t pos) {
    return gbwt::Node::encode(id(pos), is_rev(pos));
}

/// Convert Mapping to gbwt::node_type.
inline gbwt::node_type mapping_to_gbwt(const Mapping& mapping) {
    return gbwt::Node::encode(mapping.position().node_id(), mapping.position().is_reverse());
}

/// Convert Path to a GBWT path.
inline gbwt::vector_type path_to_gbwt(const Path& path) {
    gbwt::vector_type result(path.mapping_size());
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = mapping_to_gbwt(path.mapping(i));
    }
    return result;
}

/// Extract a path as a GBWT path.
gbwt::vector_type extract_as_gbwt_path(const PathHandleGraph& graph, const std::string& path_name);

// Find all predecessor nodes of the path, ignoring self-loops.
gbwt::vector_type path_predecessors(const PathHandleGraph& graph, const std::string& path_name);

//------------------------------------------------------------------------------

// GBWT construction helpers.

/// Determine the node width in bits for the GBWT nodes based on the given graph.
gbwt::size_type gbwt_node_width(const HandleGraph& graph);

/// Finish GBWT construction and optionally print the metadata.
void finish_gbwt_constuction(gbwt::GBWTBuilder& builder,
    const std::vector<std::string>& sample_names,
    const std::vector<std::string>& contig_names,
    size_t haplotype_count, bool print_metadata,
    const std::string& header = "GBWT");

//------------------------------------------------------------------------------

/*
    These are the proper ways of saving and loading GBWT structures.
    Loading with `vg::io::VPKG::load_one` is also supported.
*/

/// Load a compressed GBWT from the file.
void load_gbwt(gbwt::GBWT& index, const std::string& filename, bool show_progress = false);

/// Load a dynamic GBWT from the file.
void load_gbwt(gbwt::DynamicGBWT& index, const std::string& filename, bool show_progress = false);

/// Load an r-index from the file.
void load_r_index(gbwt::FastLocate& index, const std::string& filename, bool show_progress = false);

/// Save a compressed GBWT to the file.
void save_gbwt(const gbwt::GBWT& index, const std::string& filename, bool show_progress = false);

/// Save a dynamic GBWT to the file.
void save_gbwt(const gbwt::DynamicGBWT& index, const std::string& filename, bool show_progress = false);

/// Save an r-index to the file.
void save_r_index(const gbwt::FastLocate& index, const std::string& filename, bool show_progress = false);

//------------------------------------------------------------------------------

/**
 * Helper class that stores either a GBWT or a DynamicGBWT and loads them from a file
 * or converts between them when necessary.
 */
struct GBWTHandler {
    enum index_type { index_none, index_compressed, index_dynamic };

    /// Compressed GBWT.
    gbwt::GBWT compressed;

    /// Dynamic GBWT.
    gbwt::DynamicGBWT dynamic;

    /// Which index is in use.
    index_type in_use = index_none;

    /// The in-memory indexes are backed by this file.
    std::string filename;

    /// Print progress information to stderr when loading/converting indexes.
    bool show_progress = false;

    /// Switch to a compressed GBWT, converting it from the dynamic GBWT or reading it
    /// from a file if necessary.
    void use_compressed();

    /// Switch to a dynamic GBWT, converting it from the compressed GBWT or reading it
    /// from a file if necessary.
    void use_dynamic();

    /// Start using this compressed GBWT. Clears the index used as the argument.
    /// The GBWT is not backed by a file.
    void use(gbwt::GBWT& new_index);

    /// Start using this dynamic GBWT. Clears the index used as the argument.
    /// The GBWT is not backed by a file.
    void use(gbwt::DynamicGBWT& new_index);

    /// The GBWT is no longer backed by a file.
    void unbacked();

    /// Serialize the in-memory index to this file and start using it as the backing file.
    void serialize(const std::string& new_filename);

    /// Clear the in-memory index.
    void clear();
};

//------------------------------------------------------------------------------

/// Return the list of thread ids / gbwt path ids for the given sample.
std::vector<gbwt::size_type> threads_for_sample(const gbwt::GBWT& gbwt_index, const std::string& sample_name);

/// Return the list of thread ids / gbwt path ids for the given contig.
std::vector<gbwt::size_type> threads_for_contig(const gbwt::GBWT& gbwt_index, const std::string& contig_name);

/// Insert a GBWT thread into the graph and return its name. Returns an empty string on failure.
/// If a path name is specified and not empty, that name will be used for the inserted path.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
std::string insert_gbwt_path(MutablePathHandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id, std::string path_name = "");

/// Extract a GBWT thread as a path in the given graph.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
Path extract_gbwt_path(const HandleGraph& graph, const gbwt::GBWT& gbwt_index, gbwt::size_type id);

/// Get a string representation of a thread name stored in GBWT metadata.
/// When short_name is true, return a short name made of just the sample and contig
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
std::string thread_name(const gbwt::GBWT& gbwt_index, gbwt::size_type id, bool short_name = false);

/// Get a sample name of a thread stored in GBWT metadata.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
std::string thread_sample(const gbwt::GBWT& gbwt_index, gbwt::size_type id);

/// Get phase of a thread stored in GBWT metadata.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
int thread_phase(const gbwt::GBWT& gbwt_index, gbwt::size_type id);

/// Get count of a thread stored in GBWT metadata.
/// NOTE: id is a gbwt path id, not a gbwt sequence id.
gbwt::size_type thread_count(const gbwt::GBWT& gbwt_index, gbwt::size_type id);

//------------------------------------------------------------------------------

/// Transform the paths into a GBWT index. Primarily for testing.
gbwt::GBWT get_gbwt(const std::vector<gbwt::vector_type>& paths);

//------------------------------------------------------------------------------

/// Load a translation file (created with vg gbwt --translation) and return a mapping
/// original segment ids to a list of chopped node ids
unordered_map<nid_t, vector<nid_t>> load_translation_map(ifstream& input_stream);

/// Load a translation file (created with vg gbwt --translation) and return a backwards mapping
/// of chopped node to original segment position (id,offset pair)
/// NOTE: hopefully this is just a short-term hack, and we get a general interface baked into
//        the handlegraphs themselves
unordered_map<nid_t, pair<nid_t, size_t>> load_translation_back_map(HandleGraph& graph, ifstream& input_stream);

//------------------------------------------------------------------------------

} // namespace vg

#endif // VG_GBWT_HELPER_HPP_INCLUDED
