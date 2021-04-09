/** \file gbwt_main.cpp
 *
 * Defines the "vg gbwt" subcommand for building, merging, and manipulating GBWT indexes
 * and GBWTGraphs.
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <set>

#include "subcommand.hpp"
#include "../gbwt_helper.hpp"
#include "../haplotype_indexer.hpp"
#include "../path.hpp"
#include "../region.hpp"

#include <vg/io/vpkg.hpp>

#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gfa.h>
#include <gbwtgraph/path_cover.h>

using namespace vg;

struct GBWTConfig {
    enum build_mode { build_none, build_vcf, build_gfa, build_paths, build_alignments };
    enum merge_mode { merge_none, merge_insert, merge_fast, merge_parallel };
    enum path_cover_mode { path_cover_none, path_cover_augment, path_cover_local, path_cover_greedy };

    // Requirements and modes.
    bool produces_one_gbwt = false; // Steps 1-4 eventually produce one input GBWT regardless of the number of input args.
    bool produces_sequence_source = false; // Produces a `SequenceSource` that can be used in place of a graph in step 5.
    build_mode build = build_none;
    merge_mode merge = merge_none;
    path_cover_mode path_cover = path_cover_none;
    bool metadata_mode = false, thread_mode = false;

    // Input GBWT construction.
    HaplotypeIndexer haplotype_indexer;
    bool gam_format = false, inputs_as_jobs = false, parse_only = false;
    size_t build_jobs = default_build_jobs();

    // GFA parsing.
    gbwtgraph::GFAParsingParameters gfa_parameters;

    // Parallel merging.
    gbwt::MergeParameters merge_parameters;

    // Other parameters and flags.
    bool show_progress = false;
    bool count_threads = false;
    bool metadata = false, contigs = false, haplotypes = false, samples = false, list_names = false, thread_names = false;
    size_t num_paths = default_num_paths(), context_length = default_context_length();
    bool num_paths_set = false;
    size_t search_threads = omp_get_max_threads();

    // Input data.
    std::vector<std::string> input_filenames;
    std::string gbwt_name; // Single input GBWT to load.
    std::string graph_name;

    // File/sample names.
    std::string gbwt_output; // Output GBWT.
    std::string thread_output; // Threads in SDSL format.
    std::string graph_output; // Output GBWTGraph.
    std::string segment_translation; // Segment to node translation output.
    std::string r_index_name; // Output r-index.
    std::set<std::string> to_remove; // Sample names to remove.

    GBWTConfig() {
        this->merge_parameters.setMergeJobs(default_merge_jobs());
    }

    static size_t default_build_jobs() {
        return std::max(static_cast<size_t>(1), static_cast<size_t>(omp_get_max_threads() / 2));
    }

    static constexpr size_t default_num_paths() {
        return gbwtgraph::PATH_COVER_DEFAULT_N;
    }

    static constexpr size_t default_num_paths_local() {
        return gbwtgraph::LOCAL_HAPLOTYPES_DEFAULT_N;
    }

    static constexpr size_t default_context_length() {
        return gbwtgraph::PATH_COVER_DEFAULT_K;
    }

    static size_t default_merge_jobs() {
        return std::min(static_cast<size_t>(gbwt::MergeParameters::MERGE_JOBS), std::max(static_cast<size_t>(1), static_cast<size_t>(omp_get_max_threads() / 2)));
    }
};

struct GraphHandler {
    enum graph_type { graph_none, graph_path, graph_source };

    std::unique_ptr<PathHandleGraph> path_graph = nullptr;
    std::unique_ptr<gbwtgraph::SequenceSource> sequence_source = nullptr;
    graph_type in_use = graph_none;

    // Load the `PathHandleGraph` specified in the config and release other graphs.
    // No effect if the handler already contains a `PathHandleGraph`.
    void get_graph(const GBWTConfig& config);

    // Take the ownership of the provided `SequenceSource` and store it in the handler.
    // Releases other graphs.
    void use(std::unique_ptr<gbwtgraph::SequenceSource>& source);

    void clear();

    // If the handler contains a `SequenceSource`, serialize it according to the config.
    void serialize_segment_translation(const GBWTConfig& config) const;
};

//----------------------------------------------------------------------------

GBWTConfig parse_gbwt_config(int argc, char** argv);
void validate_gbwt_config(GBWTConfig& config);

void step_1_build_gbwts(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config);
void step_2_merge_gbwts(GBWTHandler& gbwts, GBWTConfig& config);
void step_3_remove_samples(GBWTHandler& gbwts, GBWTConfig& config);
void step_4_path_cover(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config);
void step_5_gbwtgraph(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config);
void step_6_r_index(GBWTHandler& gbwts, GBWTConfig& config);
void step_7_metadata(GBWTHandler& gbwts, GBWTConfig& config);
void step_8_threads(GBWTHandler& gbwts, GBWTConfig& config);

void report_time_memory(const std::string& what, double start_time, const GBWTConfig& config);
void print_metadata(std::ostream& out, const GBWTHandler& gbwts);

//----------------------------------------------------------------------------

int main_gbwt(int argc, char** argv) {
    GBWTConfig config = parse_gbwt_config(argc, argv);
    validate_gbwt_config(config);

    // Let GBWT operate silently and use the same temporary directory as vg.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::TempFile::setDirectory(temp_file::get_dir());

    // This is the data we are using.
    GBWTHandler gbwts;
    gbwts.filename = config.gbwt_name;
    gbwts.show_progress = config.show_progress;
    GraphHandler graphs;

    // Input GBWT construction.
    if (config.build != GBWTConfig::build_none) {
        step_1_build_gbwts(gbwts, graphs, config);
    }

    // Merge multiple input GBWTs.
    if (config.merge != GBWTConfig::merge_none) {
        step_2_merge_gbwts(gbwts, config);
    }

    // Remove samples from the GBWT.
    if (!config.to_remove.empty()) {
        step_3_remove_samples(gbwts, config);
    }

    // Path cover construction.
    if (config.path_cover != GBWTConfig::path_cover_none) {
        step_4_path_cover(gbwts, graphs, config);
    }

    // Now we can serialize the GBWT.
    if (!config.gbwt_output.empty()) {
        double start = gbwt::readTimer();
        gbwts.serialize(config.gbwt_output);
        graphs.serialize_segment_translation(config);
        report_time_memory("GBWT serialized", start, config);
    }

    // GBWTGraph construction.
    if (!config.graph_output.empty()) {
        step_5_gbwtgraph(gbwts, graphs, config);
    }

    // We no longer need the graph.
    graphs.clear();

    // R-index construction.
    if (!config.r_index_name.empty()) {
        step_6_r_index(gbwts, config);
    }

    // Metadata options.
    if (config.metadata_mode) {
        step_7_metadata(gbwts, config);
    }

    // Thread options.
    if (config.thread_mode) {
        step_8_threads(gbwts, config);
    }

    return 0;
}

//----------------------------------------------------------------------------

void help_gbwt(char** argv) {
    std::cerr << "usage: " << argv[0] << " gbwt [options] [args]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Manipulate GBWTs. Input GBWTs are loaded from input args or built in earlier steps." << std::endl;
    std::cerr << std::endl;
    std::cerr << "General options:" << std::endl;
    std::cerr << "    -x, --xg-name FILE      read the graph from FILE" << std::endl;
    std::cerr << "    -o, --output FILE       write output GBWT to FILE" << std::endl;
    std::cerr << "    -d, --temp-dir DIR      use directory DIR for temporary files" << std::endl;
    std::cerr << "    -p, --progress          show progress and statistics" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GBWT construction parameters (for steps 1 and 4):" << std::endl;
    std::cerr << "        --buffer-size N     GBWT construction buffer size in millions of nodes (default " << (gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION) << ")" << std::endl;
    std::cerr << "        --id-interval N     store path ids at one out of N positions (default " << gbwt::DynamicGBWT::SAMPLE_INTERVAL << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Search parameters (for -b and -r):" << std::endl;
    std::cerr << "        --num-threads N     use N parallel search threads (default " << omp_get_max_threads() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 1: GBWT construction (requires -o and one of { -v, -G, -E, A }):" << std::endl;
    std::cerr << "    -v, --vcf-input         index the haplotypes in the VCF files specified in input args in parallel" << std::endl;
    std::cerr << "                            (inputs must be over different contigs; requires -x, implies -f)" << std::endl;
    std::cerr << "        --preset X          use preset X (available: 1000gp)" << std::endl;
    std::cerr << "        --num-jobs N        use at most N parallel build jobs (default " << GBWTConfig::default_build_jobs() << ")" << std::endl;
    std::cerr << "        --inputs-as-jobs    create one build job for each input instead of using first-fit heuristic" << std::endl;
    std::cerr << "        --parse-only        store the VCF parses without building GBWTs" << std::endl;
    std::cerr << "                            (use -o for the file name prefix; skips subsequent steps)" << std::endl;
    std::cerr << "        --ignore-missing    do not warn when variants are missing from the graph" << std::endl;
    std::cerr << "        --actual-phasing    do not interpret unphased homozygous genotypes as phased" << std::endl;
    std::cerr << "        --force-phasing     replace unphased genotypes with randomly phased ones" << std::endl;
    std::cerr << "        --discard-overlaps  skip overlapping alternate alleles if the overlap cannot be resolved" << std::endl;
    std::cerr << "                            instead of creating a phase break" << std::endl;
    std::cerr << "        --batch-size N      index the haplotypes in batches of N samples (default 200)" << std::endl; // FIXME source for the default
    std::cerr << "        --sample-range X-Y  index samples X to Y (inclusive, 0-based)" << std::endl;
    std::cerr << "        --rename V=P        VCF contig V matches path P in the graph (may repeat)" << std::endl;
    std::cerr << "        --vcf-variants      variants in the graph use VCF contig names instead of path names" << std::endl;
    std::cerr << "        --vcf-region C:X-Y  restrict VCF contig C to coordinates X to Y (inclusive, 1-based; may repeat)" << std::endl;
    std::cerr << "        --exclude-sample X  do not index the sample with name X (faster than -R; may repeat)" << std::endl;
    std::cerr << "    -G, --gfa-input         index the walks or paths in the GFA file (one input arg)" << std::endl;
    std::cerr << "        --max-node N        chop long segments into nodes of at most N bp (default " << gbwtgraph::MAX_NODE_LENGTH << ", use 0 to disable)" << std::endl;
    std::cerr << "        --path-regex X      parse metadata from path names using regex X (default " << gbwtgraph::GFAParsingParameters::DEFAULT_REGEX << ")" << std::endl;
    std::cerr << "        --path-fields X     map regex submatches to these fields (default " << gbwtgraph::GFAParsingParameters::DEFAULT_FIELDS << ")" << std::endl;
    std::cerr << "        --translation FILE  write the segment to node translation table to FILE" << std::endl;
    std::cerr << "    -E, --index-paths       index the embedded non-alt paths in the graph (requires -x, no input args)" << std::endl;
    std::cerr << "        --paths-as-samples  each path becomes a sample instead of a contig in the metadata" << std::endl;
    std::cerr << "    -A, --alignment-input   index the alignments in the GAF files specified in input args (requires -x)" << std::endl;
    std::cerr << "        --gam-format        the input files are in GAM format instead of GAF format" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 2: Merge multiple input GBWTs (requires -o):" << std::endl;
    std::cerr << "    -m, --merge             use the insertion algorithm" << std::endl;
    std::cerr << "    -f, --fast              fast merging algorithm (node ids must not overlap)" << std::endl;
    std::cerr << "    -b, --parallel          use the parallel algorithm" << std::endl;
    std::cerr << "        --chunk-size N      search in chunks of N sequences (default " << gbwt::MergeParameters::CHUNK_SIZE << ")" << std::endl;
    std::cerr << "        --pos-buffer N      use N MiB position buffers for each search thread (default " << gbwt::MergeParameters::POS_BUFFER_SIZE << ")" << std::endl;
    std::cerr << "        --thread-buffer N   use N MiB thread buffers for each search thread (default " << gbwt::MergeParameters::THREAD_BUFFER_SIZE << ")" << std::endl;
    std::cerr << "        --merge-buffers N   merge 2^N thread buffers into one file per merge job (default " << gbwt::MergeParameters::MERGE_BUFFERS << ")" << std::endl;
    std::cerr << "        --merge-jobs N      run N parallel merge jobs (default " << GBWTConfig::default_merge_jobs() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 3: Remove samples (requires -o and one input GBWT):" << std::endl;
    std::cerr << "    -R, --remove-sample X   remove the sample with name X from the index (may repeat)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 4: Path cover GBWT construction (requires -o, -x, and one of { -a, -l, -P }):" << std::endl;
    std::cerr << "    -a, --augment-gbwt      add a path cover of missing components (one input GBWT)" << std::endl;
    std::cerr << "    -l, --local-haplotypes  sample local haplotypes (one input GBWT)" << std::endl;
    std::cerr << "    -P, --path-cover        build a greedy path cover (no input GBWTs)" << std::endl;
    std::cerr << "    -n, --num-paths N       find N paths per component (default " << GBWTConfig::default_num_paths_local() << " for -l, " << GBWTConfig::default_num_paths() << " otherwise)" << std::endl;
    std::cerr << "    -k, --context-length N  use N-node contexts (default " << GBWTConfig::default_context_length() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 5: GBWTGraph construction (requires one of { -x, -G } and one input GBWT):" << std::endl;
    std::cerr << "    -g, --graph-name FILE   build GBWTGraph and store it in FILE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 6: R-index construction (one input GBWT):" << std::endl;
    std::cerr << "    -r, --r-index FILE      build an r-index and store it in FILE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 7: Metadata (one input GBWT):" << std::endl;
    std::cerr << "    -M, --metadata          print basic metadata" << std::endl;
    std::cerr << "    -C, --contigs           print the number of contigs" << std::endl;
    std::cerr << "    -H, --haplotypes        print the number of haplotypes" << std::endl;
    std::cerr << "    -S, --samples           print the number of samples" << std::endl;
    std::cerr << "    -L, --list-names        list contig/sample names (use with -C or -S)" << std::endl;
    std::cerr << "    -T, --thread-names      list thread names" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 8: Threads (one input GBWT):" << std::endl;
    std::cerr << "    -c, --count-threads     print the number of threads" << std::endl;
    std::cerr << "    -e, --extract FILE      extract threads in SDSL format to FILE" << std::endl;
    std::cerr << std::endl;
}

//----------------------------------------------------------------------------

void use_preset(std::string preset_name, GBWTConfig& config) {
    for (char& c : preset_name) {
        c = std::tolower(c);
    }
    if (preset_name == "1000gp") {
        config.haplotype_indexer.gbwt_buffer_size = 200;
        config.haplotype_indexer.samples_in_batch = 100;
        config.haplotype_indexer.force_phasing = true;
        config.haplotype_indexer.discard_overlaps = true;
    } else {
        std::cerr << "error: [vg gbwt] invalid preset: " << preset_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

GBWTConfig parse_gbwt_config(int argc, char** argv) {
    if (argc == 2) {
        help_gbwt(argv);
        std::exit(EXIT_FAILURE);
    }

    // Long options with no corresponding short options.
    constexpr int OPT_BUFFER_SIZE = 1000;
    constexpr int OPT_ID_INTERVAL = 1001;
    constexpr int OPT_NUM_THREADS = 1002;
    constexpr int OPT_PRESET = 1100;
    constexpr int OPT_NUM_JOBS = 1101;
    constexpr int OPT_INPUTS_AS_JOBS = 1102;
    constexpr int OPT_PARSE_ONLY = 1103;
    constexpr int OPT_IGNORE_MISSING = 1104;
    constexpr int OPT_ACTUAL_PHASING = 1105;
    constexpr int OPT_FORCE_PHASING = 1106;
    constexpr int OPT_DISCARD_OVERLAPS = 1107;
    constexpr int OPT_BATCH_SIZE = 1108;
    constexpr int OPT_SAMPLE_RANGE = 1109;
    constexpr int OPT_RENAME = 1110;
    constexpr int OPT_VCF_VARIANTS = 1111;
    constexpr int OPT_VCF_REGION = 1112;
    constexpr int OPT_EXCLUDE_SAMPLE = 1113;
    constexpr int OPT_MAX_NODE = 1114;
    constexpr int OPT_PATH_REGEX = 1115;
    constexpr int OPT_PATH_FIELDS = 1116;
    constexpr int OPT_TRANSLATION = 1117;
    constexpr int OPT_PATHS_AS_SAMPLES = 1118;
    constexpr int OPT_GAM_FORMAT = 1119;
    constexpr int OPT_CHUNK_SIZE = 1200;
    constexpr int OPT_POS_BUFFER = 1201;
    constexpr int OPT_THREAD_BUFFER = 1202;
    constexpr int OPT_MERGE_BUFFERS = 1203;
    constexpr int OPT_MERGE_JOBS = 1204;

    static struct option long_options[] =
    {
        // General
        { "xg-name", required_argument, 0, 'x' },
        { "output", required_argument, 0, 'o' },
        { "temp-dir", required_argument, 0, 'd' },
        { "progress",  no_argument, 0, 'p' },

        // GBWT construction parameters
        { "buffer-size", required_argument, 0, OPT_BUFFER_SIZE },
        { "id-interval", required_argument, 0, OPT_ID_INTERVAL },

        // Search parameters
        { "num-threads", required_argument, 0, OPT_NUM_THREADS },

        // Input GBWT construction: VCF
        { "vcf-input", no_argument, 0, 'v' },
        { "preset", required_argument, 0, OPT_PRESET },
        { "num-jobs", required_argument, 0, OPT_NUM_JOBS },
        { "inputs-as-jobs", no_argument, 0, OPT_INPUTS_AS_JOBS },
        { "parse-only", no_argument, 0, OPT_PARSE_ONLY },
        { "ignore-missing", no_argument, 0, OPT_IGNORE_MISSING },
        { "actual-phasing", no_argument, 0, OPT_ACTUAL_PHASING },
        { "force-phasing", no_argument, 0, OPT_FORCE_PHASING },
        { "discard-overlaps", no_argument, 0, OPT_DISCARD_OVERLAPS },
        { "batch-size", required_argument, 0, OPT_BATCH_SIZE },
        { "sample-range", required_argument, 0, OPT_SAMPLE_RANGE },
        { "rename", required_argument, 0, OPT_RENAME },
        { "vcf-variants", no_argument, 0, OPT_VCF_VARIANTS },
        { "vcf-region", required_argument, 0, OPT_VCF_REGION },
        { "exclude-sample", required_argument, 0, OPT_EXCLUDE_SAMPLE },

        // Input GBWT construction: GFA
        { "gfa-input", no_argument, 0, 'G' },
        { "max-node", required_argument, 0, OPT_MAX_NODE },
        { "path-regex", required_argument, 0, OPT_PATH_REGEX },
        { "path-fields", required_argument, 0, OPT_PATH_FIELDS },
        { "translation", required_argument, 0, OPT_TRANSLATION },

        // Input GBWT construction: paths
        { "index-paths", no_argument, 0, 'E' },
        { "paths-as-samples", no_argument, 0, OPT_PATHS_AS_SAMPLES },

        // Input GBWT construction: GAF/GAM
        { "alignment-input", no_argument, 0, 'A' },
        { "gam-format", no_argument, 0, OPT_GAM_FORMAT },

        // Merging
        { "merge", no_argument, 0, 'm' },
        { "fast", no_argument, 0, 'f' },
        { "parallel", no_argument, 0, 'b' },
        { "chunk-size", required_argument, 0, OPT_CHUNK_SIZE },
        { "pos-buffer", required_argument, 0, OPT_POS_BUFFER },
        { "thread-buffer", required_argument, 0, OPT_THREAD_BUFFER },
        { "merge-buffers", required_argument, 0, OPT_MERGE_BUFFERS },
        { "merge-jobs", required_argument, 0, OPT_MERGE_JOBS },

        // Remove sample
        { "remove-sample", required_argument, 0, 'R' },

        // Path cover
        { "augment-gbwt", no_argument, 0, 'a' },
        { "local-haplotypes", no_argument, 0, 'l' },
        { "path-cover", no_argument, 0, 'P' },
        { "num-paths", required_argument, 0, 'n' },
        { "context-length", required_argument, 0, 'k' },

        // GBWTGraph
        { "graph-name", required_argument, 0, 'g' },

        // R-index
        { "r-index", required_argument, 0, 'r' },

        // Metadata
        { "metadata", no_argument, 0, 'M' },
        { "contigs", no_argument, 0, 'C' },
        { "haplotypes", no_argument, 0, 'H' },
        { "samples", no_argument, 0, 'S' },
        { "list-names", no_argument, 0, 'L' },
        { "thread-names", no_argument, 0, 'T' },

        // Threads
        { "count-threads", no_argument, 0, 'c' },
        { "extract", required_argument, 0, 'e' },

        { "help", no_argument, 0, 'h' },
        { 0, 0, 0, 0 }
    };

    int c;
    optind = 2; // force optind past command positional argument
    GBWTConfig config;
    while (true) {
        int option_index = 0;
        c = getopt_long(argc, argv, "x:o:d:pvGEAmfbR:alPn:k:g:r:MCHSLTce:h?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'x':
            config.graph_name = optarg;
            break;
        case 'o':
            config.gbwt_output = optarg;
            break;
        case 'd':
            temp_file::set_dir(optarg);
            break;
        case 'p':
            config.show_progress = true;
            break;

        // GBWT construction parameters
        case OPT_BUFFER_SIZE:
            config.haplotype_indexer.gbwt_buffer_size = std::max(parse<size_t>(optarg), 1ul);
            config.gfa_parameters.automatic_batch_size = false; // User-defined buffer size overrides heuristics.
            break;
        case OPT_ID_INTERVAL:
            config.haplotype_indexer.id_interval = parse<size_t>(optarg);
            break;

        // Search parameters
        case OPT_NUM_THREADS:
            config.search_threads = std::max(parse<size_t>(optarg), 1ul);
            break;

        // Input GBWT construction: VCF
        case 'v':
            assert(config.build == GBWTConfig::build_none);
            config.build = GBWTConfig::build_vcf;
            config.produces_one_gbwt = true;
            break;
        case OPT_PRESET:
            use_preset(optarg, config);
            break;
        case OPT_NUM_JOBS:
            config.build_jobs = parse<size_t>(optarg);
            break;
        case OPT_INPUTS_AS_JOBS:
            config.inputs_as_jobs = true;
            break;
        case OPT_PARSE_ONLY:
            config.parse_only = true;
            break;
        case OPT_IGNORE_MISSING:
            config.haplotype_indexer.warn_on_missing_variants = false;
            break;
        case OPT_ACTUAL_PHASING:
            config.haplotype_indexer.phase_homozygous = false;
            break;
        case OPT_FORCE_PHASING:
            config.haplotype_indexer.force_phasing = true;
            break;
        case OPT_DISCARD_OVERLAPS:
            config.haplotype_indexer.discard_overlaps = true;
            break;
        case OPT_BATCH_SIZE:
            config.haplotype_indexer.samples_in_batch = std::max(parse<size_t>(optarg), 1ul);
            break;
        case OPT_SAMPLE_RANGE:
            {
                // Parse first-last
                string range(optarg);
                size_t found = range.find("-");
                if(found == std::string::npos || found == 0 || found + 1 == range.size()) {
                    cerr << "error: [vg gbwt] could not parse range " << range << endl;
                    std::exit(EXIT_FAILURE);
                }
                config.haplotype_indexer.sample_range.first = parse<size_t>(range.substr(0, found));
                config.haplotype_indexer.sample_range.second = parse<size_t>(range.substr(found + 1)) + 1;
            }
            break;
        case OPT_RENAME:
            {
                // Parse old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error: [vg gbwt] could not parse rename " << key_value << endl;
                    std::exit(EXIT_FAILURE);
                }
                // Parse out the two parts
                string vcf_contig = key_value.substr(0, found);
                string graph_contig = key_value.substr(found + 1);
                // Add the name mapping
                config.haplotype_indexer.path_to_vcf[graph_contig] = vcf_contig;
            }
            break;
        case OPT_VCF_VARIANTS:
            config.haplotype_indexer.rename_variants = false;
            break;
        case OPT_VCF_REGION:
            {
                // Parse contig:first-last
                std::string region(optarg);
                Region parsed;
                parse_region(region, parsed);
                if (parsed.start <= 0 || parsed.end <= 0) {
                    // We need both range bounds, and we can't accept 0 since input is 1-based.
                    cerr << "error: [vg gbwt] could not parse 1-based region " << region << endl;
                }
                // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                config.haplotype_indexer.regions[parsed.seq] = std::make_pair((size_t) (parsed.start - 1), (size_t) parsed.end);
            }
            break;
        case OPT_EXCLUDE_SAMPLE:
            config.haplotype_indexer.excluded_samples.insert(optarg);
            break;

        // Input GBWT construction: GFA
        case 'G':
            assert(config.build == GBWTConfig::build_none);
            config.build = GBWTConfig::build_gfa;
            config.produces_one_gbwt = true;
            config.produces_sequence_source = true;
            break;
        case OPT_MAX_NODE:
            config.gfa_parameters.max_node_length = parse<size_t>(optarg);
            break;
        case OPT_PATH_REGEX:
            config.gfa_parameters.path_name_regex = optarg;
            break;
        case OPT_PATH_FIELDS:
            config.gfa_parameters.path_name_fields = optarg;
            break;
        case OPT_TRANSLATION:
            config.segment_translation = optarg;
            break;

        // Input GBWT construction: Paths
        case 'E':
            assert(config.build == GBWTConfig::build_none);
            config.build = GBWTConfig::build_paths;
            config.produces_one_gbwt = true;
            break;
        case OPT_PATHS_AS_SAMPLES:
            config.haplotype_indexer.paths_as_samples = true;
            break;

        // Input GBWT construction: GAF/GAM
        case 'A':
            assert(config.build == GBWTConfig::build_none);
            config.build = GBWTConfig::build_alignments;
            config.produces_one_gbwt = true;
            break;
        case OPT_GAM_FORMAT:
            config.gam_format = true;
            break;

        // Merging
        case 'm':
            config.merge = GBWTConfig::merge_insert;
            config.produces_one_gbwt = true;
            break;
        case 'f':
            config.merge = GBWTConfig::merge_fast;
            config.produces_one_gbwt = true;
            break;
        case 'b':
            config.merge = GBWTConfig::merge_parallel;
            config.produces_one_gbwt = true;
            break;
        case OPT_CHUNK_SIZE:
            config.merge_parameters.setChunkSize(parse<size_t>(optarg));
            break;
        case OPT_POS_BUFFER:
            config.merge_parameters.setPosBufferSize(parse<size_t>(optarg));
            break;
        case OPT_THREAD_BUFFER:
            config.merge_parameters.setThreadBufferSize(parse<size_t>(optarg));
            break;
        case OPT_MERGE_BUFFERS:
            config.merge_parameters.setMergeBuffers(parse<size_t>(optarg));
            break;
        case OPT_MERGE_JOBS:
            config.merge_parameters.setMergeJobs(parse<size_t>(optarg));
            break;

        // Remove sample
        case 'R':
            config.to_remove.insert(optarg);
            break;

        // Path cover
        case 'a':
            assert(config.path_cover == GBWTConfig::path_cover_none);
            config.path_cover = GBWTConfig::path_cover_augment;
            break;
        case 'l':
            assert(config.path_cover == GBWTConfig::path_cover_none);
            config.path_cover = GBWTConfig::path_cover_local;
            if (!config.num_paths_set) {
                config.num_paths = GBWTConfig::default_num_paths_local();
            }
            break;
        case 'P':
            assert(config.path_cover == GBWTConfig::path_cover_none);
            config.path_cover = GBWTConfig::path_cover_greedy;
            config.produces_one_gbwt = true;
            break;
        case 'n':
            config.num_paths = parse<size_t>(optarg);
            config.num_paths_set = true;
            break;
        case 'k':
            config.context_length = parse<size_t>(optarg);
            break;

        // GBWTGraph
        case 'g':
            config.graph_output = optarg;
            break;

        // Build r-index
        case 'r':
            config.r_index_name = optarg;
            break;

        // Metadata
        case 'M':
            config.metadata = true;
            config.metadata_mode = true;
            break;
        case 'C':
            config.contigs = true;
            config.metadata_mode = true;
            break;
        case 'H':
            config.haplotypes = true;
            config.metadata_mode = true;
            break;
        case 'S':
            config.samples = true;
            config.metadata_mode = true;
            break;
        case 'L':
            config.list_names = true;
            config.metadata_mode = true;
            break;
        case 'T':
            config.thread_names = true;
            config.metadata_mode = true;
            break;

        // Threads
        case 'c':
            config.count_threads = true;
            config.thread_mode = true;
            break;
        case 'e':
            config.thread_output = optarg;
            config.thread_mode = true;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_gbwt(argv);
            std::exit(EXIT_FAILURE);
            break;

        default:
            std::exit(EXIT_FAILURE);
        }
    }

    // The remaining args are input args.
    for (int arg = optind; arg < argc; arg++) {
        config.input_filenames.push_back(argv[arg]);
    }
    config.gbwt_name = (config.input_filenames.empty() ? "" : config.input_filenames.front());

    // Copy information from primary fields to redundant fields.
    config.haplotype_indexer.show_progress = config.show_progress;
    config.gfa_parameters.show_progress = config.show_progress;
    config.gfa_parameters.batch_size = config.haplotype_indexer.gbwt_buffer_size * gbwt::MILLION;
    config.gfa_parameters.sample_interval = config.haplotype_indexer.id_interval;

    return config;
}

//----------------------------------------------------------------------------

void validate_gbwt_config(GBWTConfig& config) {
    if (config.build != GBWTConfig::build_none) {
        if (config.gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt] GBWT construction requires output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (config.build == GBWTConfig::build_vcf) {
            if (config.graph_name.empty() || config.input_filenames.empty()) {
                std::cerr << "error: [vg gbwt] GBWT construction from VCF files requires a graph and input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (config.parse_only) {
                config.haplotype_indexer.batch_file_prefix = config.gbwt_output;
            }
        } else if (config.build == GBWTConfig::build_gfa) {
            if (!config.graph_name.empty()) {
                std::cerr << "error: [vg gbwt] GBWT construction from GFA does not use a graph" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (config.input_filenames.size() != 1) {
                std::cerr << "error: [vg gbwt] GBWT construction from GFA requires one input arg" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        } else if (config.build == GBWTConfig::build_alignments) {
            if (config.graph_name.empty() || config.input_filenames.empty()) {
                std::cerr << "error: [vg gbwt] GBWT construction from alignments requires a graph and input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        } else if (config.build == GBWTConfig::build_paths) {
            if (config.graph_name.empty()) {
                std::cerr << "error: [vg gbwt] GBWT construction from embedded paths requires a graph" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (!config.input_filenames.empty()) {
                std::cerr << "error: [vg gbwt] GBWT construction from embedded paths does not use input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }

    if (config.merge != GBWTConfig::merge_none) {
        if (config.input_filenames.size() < 2 || config.gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt] merging requires multiple input GBWTs and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    if (!config.to_remove.empty()) {
        if (!(config.input_filenames.size() == 1 || config.merge != GBWTConfig::merge_none) || config.gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt] removing a sample requires one input GBWT and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    if (config.path_cover != GBWTConfig::path_cover_none) {
        if (config.gbwt_output.empty() || config.graph_name.empty()) {
            std::cerr << "error: [vg gbwt] path cover options require output GBWT and a graph" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (config.path_cover == GBWTConfig::path_cover_greedy && !config.input_filenames.empty()) {
            std::cerr << "error: [vg gbwt] greedy path cover does not use input GBWTs" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if ((config.path_cover == GBWTConfig::path_cover_local || config.path_cover == GBWTConfig::path_cover_augment) && !(config.input_filenames.size() == 1 || config.merge != GBWTConfig::merge_none)) {
            std::cerr << "error: [vg gbwt] path cover options -a and -l require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (config.num_paths == 0) {
            std::cerr << "error: [vg gbwt] number of paths must be non-zero for path cover" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (config.context_length < gbwtgraph::PATH_COVER_MIN_K) {
            std::cerr << "error: [vg gbwt] context length must be at least " << gbwtgraph::PATH_COVER_MIN_K << " for path cover" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    if (!config.graph_output.empty()) {
        if ((config.graph_name.empty() && !config.produces_sequence_source) || !(config.input_filenames.size() == 1 || config.produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt] GBWTGraph construction requires a graph and one input GBWT or GFA input" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    if (!config.r_index_name.empty()) {
        if (!(config.input_filenames.size() == 1 || config.produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt] r-index construction requires one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    if (config.metadata_mode) {
        if (!(config.input_filenames.size() == 1 || config.produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt] metadata operations require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    if (config.thread_mode) {
        if (!(config.input_filenames.size() == 1 || config.produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt] thread operations require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
}

//----------------------------------------------------------------------------

struct job_type {
    std::string filename;
    std::vector<path_handle_t> paths;
    size_t size;

    // Large jobs first.
    bool operator<(const job_type& another) const {
        return (this->size > another.size);
    }

    typedef std::pair<path_handle_t, size_t> path_type;

    void insert(path_type path) {
        this->paths.push_back(path.first);
        this->size += path.second;
    }
};

std::vector<job_type> determine_jobs(std::unique_ptr<PathHandleGraph>& graph, const GBWTConfig& config) {

    std::vector<job_type> result;

    // Determine the non-alt paths.
    std::vector<job_type::path_type> paths;
    size_t max_length = 0;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        if (!Paths::is_alt(graph->get_path_name(path_handle))) {
            paths.emplace_back(path_handle, graph->get_step_count(path_handle));
            max_length = std::max(max_length, paths.back().second);
        }
    });
    if (paths.empty()) {
        return result;
    }

    struct vcf_paths {
        size_t file;
        std::vector<job_type::path_type> paths;

        // In descending order by length.
        void sort_paths() {
            std::sort(this->paths.begin(), this->paths.end(), [](job_type::path_type a, job_type::path_type b) -> bool {
                return (a.second > b.second);
            });
        }

        void insert(job_type::path_type path) {
            this->paths.push_back(path);
        }

        bool empty() const {
            return this->paths.empty();
        }
    };

    // Initialize the files.
    std::vector<vcf_paths> paths_by_file;
    for (size_t i = 0; i < config.input_filenames.size(); i++) {
        paths_by_file.push_back({ i, {} });
    }

    // Determine which VCF file contains each path.
    std::vector<size_t> path_found_in(paths.size(), config.input_filenames.size());
    for (size_t i = 0; i < config.input_filenames.size(); i++) {
        std::string filename = config.input_filenames[i];
        vcflib::VariantCallFile variant_file;
        variant_file.parseSamples = false;
        variant_file.open(filename);
        if (!variant_file.is_open()) {
            std::cerr << "error: [vg gbwt] could not open VCF file " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }
        for (size_t j = 0; j < paths.size(); j++) {
            std::string contig_name = graph->get_path_name(paths[j].first);
            if (config.haplotype_indexer.path_to_vcf.find(contig_name) != config.haplotype_indexer.path_to_vcf.end()) {
                contig_name = config.haplotype_indexer.path_to_vcf.at(contig_name);
            }
            variant_file.setRegion(contig_name);
            vcflib::Variant var(variant_file);
            if (!(variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == contig_name)) {
                continue;
            }
            if (path_found_in[j] < config.input_filenames.size()) {
                std::cerr << "error: [vg gbwt] contig " << contig_name << " found in files " << config.input_filenames[path_found_in[j]] << " and " << filename << std::endl;
                std::exit(EXIT_FAILURE);
            }
            paths_by_file[i].insert(paths[j]);
            path_found_in[j] = i;
        }
    }

    // Special case: Each input file is a single job.
    if (config.inputs_as_jobs) {
        for (const vcf_paths& curr : paths_by_file) {
            if (curr.empty()) {
                continue;
            }
            job_type job({ config.input_filenames[curr.file], {}, 0 });
            for (auto path : curr.paths) {
                job.insert(path);
            }
            result.push_back(job);
        }
        return result;
    }

    // First-fit heuristic: Create jobs of size at most max_length from each file.
    for (size_t i = 0; i < paths_by_file.size(); i++) {
        paths_by_file[i].sort_paths();
    }
    for (const vcf_paths& curr : paths_by_file) {
        if (curr.empty()) {
            continue;
        }
        std::vector<job_type> jobs;
        for (auto path : curr.paths) {
            bool inserted = false;
            for (size_t i = 0; i < jobs.size(); i++) {
                if (jobs[i].size + path.second <= max_length) {
                    jobs[i].insert(path);
                    inserted = true;
                    break;
                }
            }
            if (!inserted) {
                jobs.push_back({ config.input_filenames[curr.file], {}, 0 });
                jobs.back().insert(path);
            }
        }
        result.insert(result.end(), jobs.begin(), jobs.end());
    }

    // Sort the jobs in descending order by size.
    std::sort(result.begin(), result.end());
    return result;
}

void use_or_save(std::unique_ptr<gbwt::DynamicGBWT>& index, GBWTHandler& gbwts, std::vector<std::string>& filenames, size_t i, bool show_progress) {
    if (filenames.size() == 1) {
        gbwts.use(*index);
    } else {
        std::string temp = temp_file::create("gbwt-" + std::to_string(i) + "-");
        if (show_progress) {
            #pragma omp critical
            {
                std::cerr << "Job " << i << ": Saving the GBWT to " << temp << std::endl;
            }
        }
        vg::io::VPKG::save(*index, temp);
        filenames[i] = temp;
    }
}

void step_1_build_gbwts(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        std::cerr << "Building input GBWTs" << std::endl;
    }
    gbwts.unbacked(); // We will build a new GBWT.
    if (config.build != GBWTConfig::build_gfa) {
        graphs.get_graph(config);
    }

    if (config.build == GBWTConfig::build_vcf) {
        if (config.show_progress) {
            std::cerr << "Input type: VCF" << std::endl;
        }
        omp_set_num_threads(config.build_jobs);
        // Process each VCF contig corresponding to a non-alt path.
        std::vector<job_type> jobs = determine_jobs(graphs.path_graph, config);
        if (jobs.size() > 1 && config.merge == GBWTConfig::merge_none) {
            config.merge = GBWTConfig::merge_fast;
        }
        std::vector<std::vector<std::string>> vcf_parses(jobs.size());
        if (config.show_progress) {
            std::cerr << "Parsing " << jobs.size() << " VCF files using up to " << config.build_jobs << " parallel jobs" << std::endl;
        }
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t i = 0; i < jobs.size(); i++) {
            std::string job_name = "Job " + std::to_string(i);
            if (config.show_progress) {
                #pragma omp critical
                {
                    std::cerr << job_name << ": File " << jobs[i].filename << ", paths {";
                    for (path_handle_t handle : jobs[i].paths) {
                        std::cerr << " " << graphs.path_graph->get_path_name(handle);
                    }
                    std::cerr << " }" << std::endl;
                }
            }
            vcf_parses[i] = config.haplotype_indexer.parse_vcf(jobs[i].filename, *(graphs.path_graph), jobs[i].paths, job_name);
        }
        graphs.clear(); // Delete the graph to save memory.
        if (!config.parse_only) {
            std::vector<std::string> gbwt_files(vcf_parses.size(), "");
            if (config.show_progress) {
                std::cerr << "Building " << vcf_parses.size() << " GBWTs using up to " << config.build_jobs << " parallel jobs" << std::endl;
            }
            #pragma omp parallel for schedule(dynamic, 1)
            for (size_t i = 0; i < vcf_parses.size(); i++) {
                std::string job_name = "Job " + std::to_string(i);
                std::unique_ptr<gbwt::DynamicGBWT> parsed = config.haplotype_indexer.build_gbwt(vcf_parses[i], job_name);
                use_or_save(parsed, gbwts, gbwt_files, i, config.show_progress);
            }
            if (vcf_parses.size() > 1) {
                config.input_filenames = gbwt_files; // Use the temporary GBWTs as inputs.
            }
        }
    } else if (config.build == GBWTConfig::build_gfa) {
        if(config.show_progress) {
            std::cerr << "Input type: GFA" << std::endl;
        }
        auto result = gbwtgraph::gfa_to_gbwt(config.gbwt_name, config.gfa_parameters);
        if (result.first.get() == nullptr || result.second.get() == nullptr) {
            std::cerr << "error: [vg gbwt] GBWT construction from GFA failed" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        gbwts.use(*(result.first));
        graphs.use(result.second);
    } else if (config.build == GBWTConfig::build_paths) {
        if(config.show_progress) {
            std::cerr << "Input type: embedded paths" << std::endl;
        }
        std::unique_ptr<gbwt::DynamicGBWT> temp = config.haplotype_indexer.build_gbwt(*(graphs.path_graph));
        gbwts.use(*temp);
    } else if (config.build == GBWTConfig::build_alignments) {
        if (config.show_progress) {
            std::cerr << "Input type: " << (config.gam_format ? "GAM" : "GAF") << std::endl;
        }
        std::unique_ptr<gbwt::DynamicGBWT> temp = config.haplotype_indexer.build_gbwt(*(graphs.path_graph), config.input_filenames, (config.gam_format ? "GAM" : "GAF"));
        gbwts.use(*temp);
    }

    report_time_memory("GBWTs built", start, config);
    if (config.parse_only) {
        std::exit(EXIT_SUCCESS); // VCF parsing does not produce GBWTs to continue with.
    }
}
//----------------------------------------------------------------------------

void step_2_merge_gbwts(GBWTHandler& gbwts, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        std::string algo_name;
        if (config.merge == GBWTConfig::merge_fast) {
            algo_name = "fast";
        } else if (config.merge == GBWTConfig::merge_insert) {
            algo_name = "insertion";
        } else if (config.merge == GBWTConfig::merge_parallel) {
            algo_name = "parallel";
        }
        std::cerr << "Merging " << config.input_filenames.size() << " input GBWTs (" << algo_name << " algorithm)" << std::endl;
    }

    if (config.merge == GBWTConfig::merge_fast) {
        std::vector<gbwt::GBWT> indexes(config.input_filenames.size());
        for (size_t i = 0; i < config.input_filenames.size(); i++) {
            load_gbwt(config.input_filenames[i], indexes[i], config.show_progress);
        }
        if (config.show_progress) {
            std::cerr << "Merging the GBWTs" << std::endl;
        }
        gbwt::GBWT merged(indexes);
        gbwts.use(merged);
    } else if (config.merge == GBWTConfig::merge_insert) {
        gbwts.use_dynamic();
        for (size_t i = 1; i < config.input_filenames.size(); i++) {
            gbwt::GBWT next;
            load_gbwt(config.input_filenames[i], next, config.show_progress);
            if (next.size() > 2 * gbwts.dynamic.size()) {
                std::cerr << "warning: [vg gbwt] merging " << config.input_filenames[i] << " into a substantially smaller index" << std::endl;
                std::cerr << "warning: [vg gbwt] merging would be faster in another order" << std::endl;
            }
            if (config.show_progress) {
                std::cerr << "Inserting " << next.sequences() << " sequences of total length " << next.size() << std::endl;
            }
            gbwts.dynamic.merge(next);
        }
    } else if (config.merge == GBWTConfig::merge_parallel) {
        gbwts.use_dynamic();
        omp_set_num_threads(config.search_threads);
        for (size_t i = 1; i < config.input_filenames.size(); i++) {
            gbwt::DynamicGBWT next;
            load_gbwt(config.input_filenames[i], next, config.show_progress);
            if (next.size() > 2 * gbwts.dynamic.size()) {
                std::cerr << "warning: [vg gbwt] merging " << config.input_filenames[i] << " into a substantially smaller index" << std::endl;
                std::cerr << "warning: [vg gbwt] merging would be faster in another order" << std::endl;
            }
            if (config.show_progress) {
                std::cerr << "Inserting " << next.sequences() << " sequences of total length " << next.size() << std::endl;
            }
            gbwts.dynamic.merge(next, config.merge_parameters);
        }
    }
    gbwts.unbacked(); // We modified the GBWT.

    if (config.show_progress) {
        print_metadata(std::cerr, gbwts);
        report_time_memory("GBWTs merged", start, config);
    }
}

//----------------------------------------------------------------------------

void step_3_remove_samples(GBWTHandler& gbwts, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        std::cerr << "Removing " << config.to_remove.size() << " sample(s) from the index" << std::endl;
    }

    gbwts.use_dynamic();
    if (!(gbwts.dynamic.hasMetadata() && gbwts.dynamic.metadata.hasPathNames() && gbwts.dynamic.metadata.hasSampleNames())) {
        std::cerr << "error: [vg gbwt] the index does not contain metadata with thread and sample names" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::set<gbwt::size_type> sample_ids;
    for (const std::string& sample_name : config.to_remove) {
        gbwt::size_type sample_id = gbwts.dynamic.metadata.sample(sample_name);
        if (sample_id >= gbwts.dynamic.metadata.samples()) {
            std::cerr << "warning: [vg gbwt] the index does not contain sample " << sample_name << std::endl;
        } else {
            sample_ids.insert(sample_id);
        }
    }
    std::vector<gbwt::size_type> path_ids;
    for (gbwt::size_type sample_id : sample_ids) {
        std::vector<gbwt::size_type> current_paths = gbwts.dynamic.metadata.removeSample(sample_id);
        path_ids.insert(path_ids.end(), current_paths.begin(), current_paths.end());
    }
    if (path_ids.empty()) {
        std::cerr << "warning: [vg gbwt] no threads associated with the removed samples" << std::endl;
    } else {
        if (config.show_progress) {
            std::cerr << "Removing " << path_ids.size() << " threads" << std::endl;
        }
        size_t foo = gbwts.dynamic.remove(path_ids);
    }
    gbwts.unbacked(); // We modified the GBWT.

    report_time_memory("Samples removed", start, config);
}

//----------------------------------------------------------------------------

void step_4_path_cover(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        std::cerr << "Finding a " << config.num_paths << "-path cover with context length " << config.context_length << std::endl;
    }

    graphs.get_graph(config);
    if (config.path_cover == GBWTConfig::path_cover_greedy) {
        if (config.show_progress) {
            std::cerr << "Algorithm: greedy" << std::endl;
        }
        gbwt::GBWT cover = gbwtgraph::path_cover_gbwt(*(graphs.path_graph), config.num_paths, config.context_length, config.haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, config.haplotype_indexer.id_interval, config.show_progress);
        gbwts.use(cover);
    } else if (config.path_cover == GBWTConfig::path_cover_augment) {
        if (config.show_progress) {
            std::cerr << "Algorithm: augment" << std::endl;
        }
        gbwts.use_dynamic();
        gbwtgraph::augment_gbwt(*(graphs.path_graph), gbwts.dynamic, config.num_paths, config.context_length, config.haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, config.haplotype_indexer.id_interval, config.show_progress);
    } else {
        if (config.show_progress) {
            std::cerr << "Algorithm: local haplotypes" << std::endl;
        }
        gbwts.use_compressed();
        gbwt::GBWT cover = gbwtgraph::local_haplotypes(*(graphs.path_graph), gbwts.compressed, config.num_paths, config.context_length, config.haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, config.haplotype_indexer.id_interval, config.show_progress);
        gbwts.use(cover);
    }
    gbwts.unbacked(); // We modified the GBWT.

    report_time_memory("Path cover built", start, config);
}

//----------------------------------------------------------------------------

void step_5_gbwtgraph(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        std::cerr << "Building GBWTGraph" << std::endl;
    }

    gbwts.use_compressed();
    if (config.show_progress) {
        std::cerr << "Starting the construction" << std::endl;
    }
    gbwtgraph::GBWTGraph graph;
    if (graphs.in_use == GraphHandler::graph_source) {
        graph = gbwtgraph::GBWTGraph(gbwts.compressed, *(graphs.sequence_source));
    } else {
        graphs.get_graph(config);
        graph = gbwtgraph::GBWTGraph(gbwts.compressed, *(graphs.path_graph));
    }
    if (config.show_progress) {
        std::cerr << "Serializing GBWTGraph to " << config.graph_output << std::endl;
    }
    vg::io::VPKG::save(graph, config.graph_output);

    report_time_memory("GBWTGraph built", start, config);
}

//----------------------------------------------------------------------------

void step_6_r_index(GBWTHandler& gbwts, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        std::cerr << "Building r-index" << std::endl;
    }

    omp_set_num_threads(config.search_threads);
    gbwts.use_compressed();
    if (config.show_progress) {
        std::cerr << "Starting the construction" << std::endl;
    }
    gbwt::FastLocate r_index(gbwts.compressed);
    if (config.show_progress) {
        std::cerr << "Serializing the r-index to " << config.r_index_name << std::endl;
    }
    vg::io::VPKG::save(r_index, config.r_index_name);

    report_time_memory("R-index built", start, config);
}

//----------------------------------------------------------------------------

void step_7_metadata(GBWTHandler& gbwts, GBWTConfig& config) {
    gbwts.use_compressed();
    if (!gbwts.compressed.hasMetadata()) {
        std::cerr << "error: [vg gbwt] the GBWT does not contain metadata" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (config.metadata) {
        print_metadata(std::cout, gbwts);
    }

    if (config.contigs) {
        if (config.list_names) {
            if (gbwts.compressed.metadata.hasContigNames()) {
                for (size_t i = 0; i < gbwts.compressed.metadata.contigs(); i++) {
                    std::cout << gbwts.compressed.metadata.contig(i) << std::endl;
                }
            } else {
                std::cerr << "error: [vg gbwt] the metadata does not contain contig names" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        } else {
            std::cout << gbwts.compressed.metadata.contigs() << std::endl;
        }
    }

    if (config.haplotypes) {
        std::cout << gbwts.compressed.metadata.haplotypes() << std::endl;
    }

    if (config.samples) {
        if (config.list_names) {
            if (gbwts.compressed.metadata.hasSampleNames()) {
                for (size_t i = 0; i < gbwts.compressed.metadata.samples(); i++) {
                    std::cout << gbwts.compressed.metadata.sample(i) << std::endl;
                }
            } else {
                std::cerr << "error: [vg gbwt] the metadata does not contain sample names" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        } else {
            std::cout << gbwts.compressed.metadata.samples() << std::endl;
        }
    }

    if (config.thread_names) {
        if (gbwts.compressed.metadata.hasPathNames()) {
            for (size_t i = 0; i < gbwts.compressed.metadata.paths(); i++) {
                std::cout << thread_name(gbwts.compressed, i) << std::endl;
            }
        } else {
            std::cerr << "error: [vg gbwt] the metadata does not contain thread names" << std::endl;
        }
    }
}

//----------------------------------------------------------------------------

void step_8_threads(GBWTHandler& gbwts, GBWTConfig& config) {
    // Extract threads in SDSL format.
    if (!config.thread_output.empty()) {
        double start = gbwt::readTimer();
        if (config.show_progress) {
            std::cerr << "Extracting threads to " << config.thread_output << std::endl;
        }
        gbwts.use_compressed();
        if (config.show_progress) {
            std::cerr << "Starting the extraction" << std::endl;
        }
        gbwt::size_type node_width = gbwt::bit_length(gbwts.compressed.sigma() - 1);
        gbwt::text_buffer_type out(config.thread_output, std::ios::out, gbwt::MEGABYTE, node_width);
        for (gbwt::size_type id = 0; id < gbwts.compressed.sequences(); id += 2) { // Ignore reverse complements.
            gbwt::vector_type sequence = gbwts.compressed.extract(id);
            for (auto node : sequence) {
                out.push_back(node);
            }
            out.push_back(gbwt::ENDMARKER);
        }
        out.close();
        report_time_memory("Threads extracted", start, config);
    }

    // There are two sequences for each thread.
    if (config.count_threads) {
        gbwts.use_compressed();
        std::cout << (gbwts.compressed.sequences() / 2) << std::endl;
    }
}

//----------------------------------------------------------------------------

void GraphHandler::get_graph(const GBWTConfig& config) {
    if (this->in_use == graph_path) {
        return;
    } else {
        if (config.show_progress) {
            std::cerr << "Loading input graph from " << config.graph_name << std::endl;
        }
        this->sequence_source.reset();
        this->path_graph = vg::io::VPKG::load_one<PathHandleGraph>(config.graph_name);
        if (this->path_graph == nullptr) {
            std::cerr << "error: [vg gbwt] could not load graph " << config.graph_name << std::endl;
            std::exit(EXIT_FAILURE);
        }
        this->in_use = graph_path;
    }
}

void GraphHandler::use(std::unique_ptr<gbwtgraph::SequenceSource>& source) {
    this->path_graph.reset();
    this->sequence_source = std::move(source);
    this->in_use = graph_source;
}

void GraphHandler::clear() {
    this->path_graph.reset();
    this->sequence_source.reset();
    this->in_use = graph_none;
}

void GraphHandler::serialize_segment_translation(const GBWTConfig& config) const {
    if (this->in_use != graph_source || config.segment_translation.empty()) {
        return;
    }
    if (config.show_progress) {
        std::cerr << "Serializing segment to node translation to " << config.segment_translation << std::endl;
    }

    std::ofstream out(config.segment_translation, std::ios_base::binary);
    if (this->sequence_source->uses_translation()) {
        auto& translation = this->sequence_source->segment_translation;
        for (auto iter = translation.begin(); iter != translation.end(); ++iter) {
            out << "T\t" << iter->first << "\t" << iter->second.first;
            for(nid_t i = iter->second.first + 1; i < iter->second.second; i++) {
            out << "," << i;
            }
            out << "\n";
        }
    }
    out.close();
}

//----------------------------------------------------------------------------

void report_time_memory(const std::string& what, double start_time, const GBWTConfig& config) {
    if (config.show_progress) {
        double seconds = gbwt::readTimer() - start_time;
        std::cerr << what << " in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
        std::cerr << std::endl;
    }
}

void print_metadata(std::ostream& out, const GBWTHandler& gbwts) {
    if (gbwts.in_use == GBWTHandler::index_compressed) {
        gbwt::operator<<(out, gbwts.compressed.metadata) << std::endl;
    } else if (gbwts.in_use == GBWTHandler::index_dynamic) {
        gbwt::operator<<(out, gbwts.dynamic.metadata) << std::endl;
    }
}

//----------------------------------------------------------------------------

// Register subcommand
static vg::subcommand::Subcommand vg_gbwt("gbwt", "build and manipulate GBWTs", vg::subcommand::TOOLKIT, main_gbwt);
