/** \file validate_main.cpp
 *
 * Defines the "vg validate" subcommand
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../alignment.hpp"
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_validate(char** argv) {
    cerr << "usage: " << argv[0] << " validate [options] [graph]" << endl
         << "Validate the graph." << endl
         << endl
         << "options:" << endl
         << "    default: check all aspects of the graph, if options are specified do only those" << endl
         << "    -o, --orphans   verify that all nodes have edges" << endl
         << "    -a, --gam FILE  verify that edits in the alignment fit on nodes in the graph" << endl;
}

int main_validate(int argc, char** argv) {

    if (argc <= 2) {
        help_validate(argv);
        return 1;
    }

    bool check_orphans = false;
    string gam_path;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"orphans", no_argument, 0, 'o'},
            {"gam", required_argument, 0, 'a'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hoa:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'o':
                check_orphans = true;
                break;

            case 'a':
                gam_path = optarg;
                break;

            case 'h':
            case '?':
                help_validate(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    // load the graph
    unique_ptr<PathHandleGraph> graph;
    string graph_filename = get_input_file_name(optind, argc, argv);
    graph = vg::io::VPKG::load_one<PathHandleGraph>(graph_filename);

    // validate the alignment if given
    bool valid_aln = true;
    if (!gam_path.empty()) {
        get_input_file(gam_path, [&](istream& in) {
                vg::io::for_each<Alignment>(in, [&](Alignment& aln) {
                        if (!alignment_is_valid(aln, graph.get())) {
                            valid_aln = false;
                        }
                    });
            });
    }

    // VG's a little less structured, so try its own logic
    bool valid_graph = true;
    VG* vg_graph = dynamic_cast<VG*>(graph.get());
    if (vg_graph != nullptr) {
        if (!vg_graph->is_valid(true, true, check_orphans, true)) {
            valid_graph = false;
        }
    }

    if (valid_graph) {
        // I don't think this is possible with any libbdsg implementations, but check edges just in case
        graph->for_each_edge([&](const edge_t& edge) {
                if (!graph->has_node(graph->get_id(edge.first))) {
                    cerr << "graph invalid: source node missing for edge " 
                         << graph->get_id(edge.first) << ":" << graph->get_is_reverse(edge.first) << " -> "
                         << graph->get_id(edge.second) << ":" << graph->get_is_reverse(edge.second) << endl;
                valid_graph = false;
                }
                if (!graph->has_node(graph->get_id(edge.second))) {
                    cerr << "graph invalid: sink node missing for edge " 
                         << graph->get_id(edge.first) << ":" << graph->get_is_reverse(edge.first) << " -> "
                         << graph->get_id(edge.second) << ":" << graph->get_is_reverse(edge.second) << endl;
                    valid_graph = false;
                }
            });

        graph->for_each_path_handle([&](path_handle_t path_handle) {
                size_t i = 0;
                handle_t prev;
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        handle_t handle = graph->get_handle_of_step(step_handle);
                        if (i > 0) {
                            if (!graph->has_edge(prev, handle)) {
                                cerr << "graph invalid: missing edge between " << (i-1) << "th step ("
                                     << graph->get_id(prev) << ":" << graph->get_is_reverse(prev) << ") and "
                                     << i << "th step (" << graph->get_id(handle) << ":" << graph->get_is_reverse(handle)
                                     << ") of path " << graph->get_path_name(path_handle) << endl;
                                valid_graph = false;
                            }
                            if (!graph->has_edge(graph->flip(handle), graph->flip(prev))) {
                                cerr << "graph invalid: missing edge between " << (i) << "th step ("
                                     << graph->get_id(handle) << ":" << !graph->get_is_reverse(handle) << ") and "
                                     << (i-1) << "th step (" << graph->get_id(prev) << ":" << graph->get_is_reverse(prev)
                                     << ") of path " << graph->get_path_name(path_handle) << endl;
                                valid_graph = false;
                            }
                        }
                        ++i;
                        prev = handle;
                    });
            });

        if (check_orphans) {
            graph->for_each_handle([&](handle_t handle) {
                    if (graph->get_degree(handle, true) + graph->get_degree(handle, false) == 0) {
                        cerr << "graph invalid: orphan node found: " << graph->get_id(handle) << endl;
                        valid_graph = false;
                    }
                });
            
        }
    }

    if (!gam_path.empty()) {
        cerr << "alignment: " << (valid_aln ? "valid" : "invalid") << endl;
    }
    cerr << "graph: " << (valid_graph ? "valid" : "invalid") << endl;

    return valid_aln && valid_graph ? 0 : 1;
}

// Register subcommand
static Subcommand vg_validate("validate", "validate the semantics of a graph or gam", DEVELOPMENT, main_validate);

