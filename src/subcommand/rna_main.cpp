/** \file rna_main.cpp
 *
 * Defines the "vg rna" subcommand.
 */

#include <unistd.h>
#include <getopt.h>
#include <chrono>

#include "subcommand.hpp"

#include "../transcriptome.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include "../gbwt_helper.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;


void help_rna(char** argv) {
    cerr << "\nusage: " << argv[0] << " rna [options] graph.[vg|pg|hg|og] > splice_graph.[vg|pg|hg|og]" << endl

         << "\nGeneral options:" << endl

         << "    -t, --threads INT          number of compute threads to use [1]" << endl
         << "    -p, --progress             show progress" << endl
         << "    -h, --help                 print help message" << endl

         << "\nInput options:" << endl

         << "    -n, --transcripts FILE     transcript file(s) in gtf/gff format; may repeat" << endl
         << "    -m, --introns FILE         intron file(s) in bed format; may repeat " << endl
         << "    -y, --feature-type NAME    parse only this feature type in the gtf/gff (parse all if empty) [exon]" << endl
         << "    -s, --transcript-tag NAME  use this attribute tag in the gtf/gff file(s) as id [transcript_id]" << endl
         << "    -l, --haplotypes FILE      project transcripts onto haplotypes in GBWT index file" << endl

         << "\nConstruction options:" << endl

         << "    -e, --use-all-paths        project transcripts onto all embedded haplotype paths" << endl
         << "    -c, --do-not-collapse      do not collapse identical transcripts across haplotypes" << endl
         << "    -d, --remove-non-gene      remove intergenic and intronic regions (deletes reference paths)" << endl
         << "    -o, --do-not-sort          do not topological sort and compact splice graph (default if -l is used)" << endl
         << "    -r, --add-ref-paths        add reference transcripts as embedded paths in the splice graph" << endl
         << "    -a, --add-non-ref-paths    add non-reference transcripts as embedded paths in the splice graph" << endl

         << "\nOutput options:" << endl

         << "    -u, --out-ref-paths        output reference transcripts in GBWT, fasta and info" << endl
         << "    -b, --write-gbwt FILE      write transcripts as threads to GBWT index file" << endl
         << "    -g, --gbwt-bidirectional   add transcripts as bidirectional threads to GBWT index" << endl
         << "    -w, --gbwt-id-interval     transcript id interval in GBWT index [32]" << endl
         << "    -f, --write-fasta FILE     write transcripts as sequences to fasta file" << endl
         << "    -i, --write-info FILE      write transcript origin info to tsv file" << endl

         << endl;
}

int32_t main_rna(int32_t argc, char** argv) {

    if (argc == 2) {
        help_rna(argv);
        return 1;
    }
    
    vector<string> transcript_filenames;
    vector<string> intron_filenames;
    string feature_type = "exon";
    string transcript_tag = "transcript_id";
    string haplotypes_filename;
    bool use_all_paths = false;
    bool collapse_transcript_paths = true;
    bool remove_non_transcribed = false;
    bool sort_collapse_graph = true;
    bool add_reference_transcript_paths = false;
    bool add_non_reference_transcript_paths = false;
    bool output_reference_transcript_paths = false;
    string gbwt_out_filename = "";
    bool gbwt_add_bidirectional = false;
    int32_t gbwt_id_interval = 32;
    string fasta_out_filename = "";
    string info_out_filename = "";
    int32_t num_threads = 1;
    bool show_progress = false;

    int32_t c;
    optind = 2;

    while (true) {
        static struct option long_options[] =
            {
                {"transcripts",  no_argument, 0, 'n'},
                {"introns",  no_argument, 0, 'm'},
                {"feature-type",  no_argument, 0, 'y'},
                {"transcript-tag",  no_argument, 0, 's'},
                {"haplotypes",  no_argument, 0, 'l'},
                {"use-all-paths",  no_argument, 0, 'e'},
                {"do-not-collapse",  no_argument, 0, 'c'},
                {"remove-non-gene",  no_argument, 0, 'd'},
                {"do-not-sort",  no_argument, 0, 'o'},
                {"add-ref-paths",  no_argument, 0, 'r'},
                {"add-non-ref-paths",  no_argument, 0, 'a'},
                {"out-ref-paths",  no_argument, 0, 'u'},           
                {"write-gbwt",  no_argument, 0, 'b'},
                {"gbwt-bidirectional",  no_argument, 0, 'g'},
                {"gbwt-id-interval",  no_argument, 0, 'w'},
                {"write-fasta",  no_argument, 0, 'f'},
                {"write-info",  no_argument, 0, 'i'},
                {"threads",  no_argument, 0, 't'},
                {"progress",  no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int32_t option_index = 0;
        c = getopt_long(argc, argv, "n:m:y:s:l:ercdoraub:gw:f:i:t:ph?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'n':
            transcript_filenames.push_back(optarg);
            break;

        case 'm':
            intron_filenames.push_back(optarg);
            break;

        case 'y':
            feature_type = optarg;
            break;

        case 's':
            transcript_tag = optarg;
            break;

        case 'l':
            haplotypes_filename = optarg;
            break;

        case 'e':
            use_all_paths = true;
            break;

        case 'c':
            collapse_transcript_paths = false;
            break;

        case 'd':
            remove_non_transcribed = true;
            break;

        case 'o':
            sort_collapse_graph = false;
            break;

        case 'r':
            add_reference_transcript_paths = true;
            break;

        case 'a':
            add_non_reference_transcript_paths = true;
            break;

        case 'u':
            output_reference_transcript_paths = true;
            break;

        case 'b':
            gbwt_out_filename = optarg;
            break;

        case 'g':
            gbwt_add_bidirectional = true;
            break;

        case 'w':
            gbwt_id_interval = stoi(optarg);
            break;

        case 'f':
            fasta_out_filename = optarg;
            break;

        case 'i':
            info_out_filename = optarg;
            break;

        case 't':
            num_threads = stoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;

        case 'h':
        case '?':
            help_rna(argv);
            exit(1);
            break;

        default:
            abort();
        }
    }

    if (argc < optind + 1) {
        help_rna(argv);
        return 1;
    }

    if (transcript_filenames.empty() && intron_filenames.empty()) {

        cerr << "[vg rna] ERROR: No transcripts or introns were given. Use --transcripts FILE and/or --introns FILE." << endl;
        return 1;       
    }

    if (remove_non_transcribed && !add_reference_transcript_paths && !add_non_reference_transcript_paths) {

        cerr << "[vg rna] WARNING: Reference paths are deleted when removing intergenic and intronic regions. Consider adding transcripts as embedded paths using --add-ref-paths and/or --add-non-ref-paths." << endl;
    }

    double time_parsing_start = gcsa::readTimer();
    if (show_progress) { cerr << "[vg rna] Parsing graph file ..." << endl; }

    unique_ptr<MutablePathDeletableHandleGraph> splice_graph(nullptr);

    // Load variation graph.
    string splice_graph_filename = get_input_file_name(optind, argc, argv);
    splice_graph = move(vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(splice_graph_filename));

    if (splice_graph == nullptr) {
        cerr << "[transcriptome] ERROR: Could not load graph." << endl;
        exit(1);
    }

    // Construct transcriptome and parse graph.
    Transcriptome transcriptome(move(splice_graph));
    assert(splice_graph == nullptr);

    unique_ptr<gbwt::GBWT> haplotype_index;

    if (!haplotypes_filename.empty()) {

        // Load haplotype GBWT index.
        if (show_progress) { cerr << "[vg rna] Parsing haplotype GBWT index file ..." << endl; }
        haplotype_index = vg::io::VPKG::load_one<gbwt::GBWT>(haplotypes_filename);
        assert(haplotype_index->bidirectional());

        sort_collapse_graph = false;

    } else {

        // Construct empty GBWT index if no is given. 
        haplotype_index = unique_ptr<gbwt::GBWT>(new gbwt::GBWT());
    }

    if (show_progress) { cerr << "[vg rna] Graph " << ((!haplotype_index->empty()) ? "and GBWT index " : "") << "parsed in " << gcsa::readTimer() - time_parsing_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };

    transcriptome.num_threads = num_threads;
    transcriptome.feature_type = feature_type;
    transcriptome.transcript_tag = transcript_tag;
    transcriptome.use_all_paths = use_all_paths;
    transcriptome.use_reference_paths = (add_reference_transcript_paths || output_reference_transcript_paths);
    transcriptome.collapse_transcript_paths = collapse_transcript_paths;


    double time_splice_start = gcsa::readTimer();
    if (show_progress) { cerr << "[vg rna] Adding novel exon boundaries and splice-junctions to graph ..." << endl; }

    int32_t num_introns_added = 0;

    // Add introns as novel splice-junctions to graph.
    for (auto & filename: intron_filenames) {

        get_input_file(filename, [&](istream& transcript_stream) {
            num_introns_added += transcriptome.add_intron_splice_junctions(transcript_stream, haplotype_index);
        });
    }

    int32_t num_transcripts_added = 0;

    // Add transcripts as novel exon boundaries and splice-junctions to graph.
    for (auto & filename: transcript_filenames) {

        get_input_file(filename, [&](istream& transcript_stream) {
            num_transcripts_added += transcriptome.add_transcript_splice_junctions(transcript_stream, haplotype_index);
        });
    }

    if (show_progress) { cerr << "[vg rna] " << num_introns_added << " introns and " << num_transcripts_added <<  " transcripts parsed, and graph augmented " << ((transcriptome.splice_graph_node_updated()) ? "" : "(no novel exon boundaries) ") << "in " << gcsa::readTimer() - time_splice_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };


    if (sort_collapse_graph) {

        double time_sort_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Topological sorting and compacting splice graph ..." << endl; }
        
        transcriptome.compact_ordered();
        
        if (show_progress) { cerr << "[vg rna] Splice graph sorted and compacted in " << gcsa::readTimer() - time_sort_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    if (!haplotype_index->empty() || transcriptome.use_all_paths || transcriptome.use_reference_paths) {

        double time_project_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Projecting haplotype-specfic transcripts ..." << endl; }

        int32_t num_transcripts_projected = 0;

        // Add transcripts to transcriptome by projecting them onto embedded paths 
        // in a graph and/or haplotypes in a GBWT index.
        for (auto & filename: transcript_filenames) {

            get_input_file(filename, [&](istream& transcript_stream) {
                num_transcripts_projected += transcriptome.add_transcripts(transcript_stream, *haplotype_index);
            });
        }

        if (show_progress) { cerr << "[vg rna] " << num_transcripts_projected <<  " haplotype-specfic transcripts projected " << "in " << gcsa::readTimer() - time_project_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }

    if (remove_non_transcribed) {

        double time_remove_start = gcsa::readTimer();
        if (show_progress) { cerr << "[vg rna] Removing non-transcribed regions ..." << endl; }

        transcriptome.remove_non_transcribed(!(add_reference_transcript_paths || add_non_reference_transcript_paths));

        if (show_progress) { cerr << "[vg rna] Regions removed in " << gcsa::readTimer() - time_remove_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    if (add_reference_transcript_paths || add_non_reference_transcript_paths) {

        double time_add_start = gcsa::readTimer();

        if (add_reference_transcript_paths && add_non_reference_transcript_paths) {

            if (show_progress) { cerr << "[vg rna] Adding transcripts as embedded paths in the splice graph ..." << endl; }

        } else {

            if (show_progress) { cerr << "[vg rna] Adding " << ((add_reference_transcript_paths) ? "reference" : "non-reference") << " transcripts as embedded paths in the splice graph ..." << endl; }
        }

        auto num_embedded_paths = transcriptome.embed_transcript_paths(add_reference_transcript_paths, add_non_reference_transcript_paths);

        if (show_progress) { cerr << "[vg rna] " << num_embedded_paths << " paths added in " << gcsa::readTimer() - time_add_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };
    }


    double time_writing_start = gcsa::readTimer();
    int32_t num_written_transcripts = 0;

    // Construct and write GBWT index of transcript paths in transcriptome.
    if (!gbwt_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing transcripts as " << ((gbwt_add_bidirectional) ? "bidirectional " : "") << "threads to GBWT index file ..." << endl; }

        // Silence GBWT index construction. 
        gbwt::Verbosity::set(gbwt::Verbosity::SILENT); 
        gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(transcriptome.splice_graph().max_node_id(), true)), gbwt::DynamicGBWT::INSERT_BATCH_SIZE, gbwt_id_interval);

        auto num_added_threads = transcriptome.add_transcripts_to_gbwt(&gbwt_builder, output_reference_transcript_paths, gbwt_add_bidirectional);

        // Finish contruction and recode index.
        gbwt_builder.finish();

        vg::io::VPKG::save(gbwt_builder.index, gbwt_out_filename);

        assert(num_written_transcripts == num_added_threads || num_written_transcripts == 0);
        num_written_transcripts = num_added_threads;
    }

    // Write transcript path sequences in transcriptome to fasta file.
    if (!fasta_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing transcripts as sequences to fasta file ..." << endl; }

        ofstream fasta_ostream;
        fasta_ostream.open(fasta_out_filename);
        auto num_written_sequences = transcriptome.write_sequences(&fasta_ostream, output_reference_transcript_paths);
        fasta_ostream.close();

        assert(num_written_transcripts == num_written_sequences || num_written_transcripts == 0);
        num_written_transcripts = num_written_sequences;
    }    

    // Write origin info on transcripts in transcriptome to tsv file.
    if (!info_out_filename.empty()) {

        if (show_progress) { cerr << "[vg rna] Writing origin info on transcripts to tsv file ..." << endl; }

        ofstream info_ostream;
        info_ostream.open(info_out_filename);
        auto num_written_info = transcriptome.write_info(&info_ostream, *haplotype_index, output_reference_transcript_paths);
        info_ostream.close();

        assert(num_written_transcripts == num_written_info || num_written_transcripts == 0);
        num_written_transcripts = num_written_info;
    }    

    if (show_progress) { cerr << "[vg rna] Writing splice graph to stdout ..." << endl; }

    // Write splice graph to stdout 
    transcriptome.write_splice_graph(&cout);

    if (show_progress) { cerr << "[vg rna] Splice graph " << ((!gbwt_out_filename.empty() || !fasta_out_filename.empty() || !info_out_filename.empty()) ? "and " + to_string(num_written_transcripts) + " transcripts " : "") << "written in " << gcsa::readTimer() - time_writing_start << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl; };

    return 0;
}

// Register subcommand
static Subcommand vg_rna("rna", "construct spliced variation graphs and transcript paths", PIPELINE, 3, main_rna);

