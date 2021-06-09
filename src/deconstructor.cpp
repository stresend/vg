#include "deconstructor.hpp"
#include "traversal_finder.hpp"
#include <gbwtgraph/gbwtgraph.h>

//#define debug

using namespace std;


namespace vg {
Deconstructor::Deconstructor() : VCFOutputCaller("") {
}
Deconstructor::~Deconstructor(){
    for (auto& c : gbwt_pos_caches) {
        delete c;
    }
}

/**
 * Takes in a vector of snarltraversals
 * returns their sequences as a vector<string>
 * returns a boolean hasRef
 * if a reference path is present, hasRef is set to true and the first
 * string in the vector is the reference allele
 * otherwise, hasRef is set to false and all strings are alt alleles.
 */
vector<int> Deconstructor::get_alleles(vcflib::Variant& v, const vector<SnarlTraversal>& travs, int ref_path_idx,
                                       char prev_char, bool use_start) {

    assert(ref_path_idx >=0 && ref_path_idx < travs.size());

    // map strings to allele numbers
    // (we are using the traversal finder in such a way that duplicate alleles can get returned
    // in order to be able to preserve the path names)
    map<string, int> allele_idx;
    size_t cur_alt = 1;

    // go from traversals number (offset in travs) to allele number
    vector<int> trav_to_allele(travs.size());

    // compute the allele as a string
    auto trav_to_string = [&](const SnarlTraversal& trav) {
        string allele;
        // we skip the snarl endpoints
        for (int j = 1; j < trav.visit_size() - 1; ++j) {
            const string& node_sequence = graph->get_sequence(graph->get_handle(trav.visit(j).node_id()));
            allele += trav.visit(j).backward() ? reverse_complement(node_sequence) : node_sequence;
        }
        return toUppercase(allele);
    };
       

    // set the reference allele
    string ref_allele = trav_to_string(travs.at(ref_path_idx));
    allele_idx[ref_allele] = 0;
    trav_to_allele[ref_path_idx] = 0;
    bool substitution = true;
        
    // set the other alleles (they can end up as 0 alleles too if their strings match the reference)
    for (int i = 0; i < travs.size(); ++i) {
        if (i != ref_path_idx) {
            string allele = trav_to_string(travs[i]);
            auto ai_it = allele_idx.find(allele);
            if (ai_it == allele_idx.end()) {
                // make a new allele for this string
                allele_idx[allele] = cur_alt;
                trav_to_allele.at(i) = cur_alt;
                ++cur_alt;
                substitution = substitution && allele.size() == ref_allele.size();
            } else {
                // allele string has been seen, map this traversal to it
                trav_to_allele.at(i) = ai_it->second;
            }
        }
    }

    // fill in the variant
    v.alleles.resize(allele_idx.size());
    assert(allele_idx.size() > 0);
    v.alt.resize(allele_idx.size() - 1);

    for (auto ai_pair : allele_idx) {
        string allele_string = ai_pair.first;
        if (!use_start) {
            reverse_complement_in_place(allele_string);
        }
        if (!substitution) {
            allele_string = string(1, prev_char) + allele_string;
        }
        v.alleles[ai_pair.second] = allele_string;
        if (ai_pair.second > 0) {
            v.alt[ai_pair.second - 1] = allele_string;
        } else {
            v.ref = allele_string;
        }
    }

    // shift our variant back if it's an indel
    if (!substitution) {
        assert(v.position >= 2);
        --v.position;
    }

    v.updateAlleleIndexes();

    return trav_to_allele;
}

void Deconstructor::get_genotypes(vcflib::Variant& v, const vector<string>& names,
                                  const vector<int>& trav_to_allele,
                                  const vector<gbwt::size_type>& trav_thread_ids) {
    assert(names.size() == trav_to_allele.size());
    // set up our variant fields
    v.format.push_back("GT");
    if (path_to_sample && path_restricted) {
        v.format.push_back("PI");
    }

    // get a list of traversals for every vcf sample
    // (this will be 1:1 unless we're using the path_to_sample name map)
    map<string, vector<int> > sample_to_traversals;
    // phasing information from the gbwt where applicable
    vector<int> gbwt_phases(trav_to_allele.size(), -1);
    for (int i = 0; i < names.size(); ++i) {
        string sample_name;
        if (trav_thread_ids[i] != numeric_limits<gbwt::size_type>::max()) {
            sample_name = thread_sample(gbwt_trav_finder->get_gbwt(), gbwt::Path::id(trav_thread_ids[i]));
            gbwt_phases[i] = thread_phase(gbwt_trav_finder->get_gbwt(), gbwt::Path::id(trav_thread_ids[i]));
        }
        else if (path_to_sample && path_to_sample->count(names[i])) {
            sample_name = path_to_sample->find(names[i])->second;
        } else {
            sample_name = names[i];
        }
        if (sample_names.count(sample_name)) {
            sample_to_traversals[sample_name].push_back(i);
        }
    }

    // write out the genotype for each sample
    // if we're mapping a vg path name to its prefix for the sample name, we stick some information about the full
    // path name in the PI part of format
    set<string> conflicts;
    for (const auto& sample_name : sample_names) {
        if (sample_to_traversals.count(sample_name)) {
            const vector<int>& travs = sample_to_traversals[sample_name];
            assert(!travs.empty());
            vector<int> chosen_travs;
            bool conflict;
            std::tie(chosen_travs, conflict) = choose_traversals(sample_name, travs, trav_to_allele, names, gbwt_phases);
            if (conflict) {
                conflicts.insert(sample_name);
            }            
            string genotype;
            for (int i = 0; i < chosen_travs.size(); ++i) {
                if (i > 0) {
                    genotype += gbwt_trav_finder.get() ? "|" : "/";
                }
                genotype += chosen_travs[i] != -1 ? std::to_string(trav_to_allele[chosen_travs[i]]) : ".";
            }
            v.samples[sample_name]["GT"] = {genotype};
            if (path_to_sample) {
                for (auto trav : travs) {
                    v.samples[sample_name]["PI"].push_back(names[trav] + "=" + std::to_string(trav_to_allele[trav]));
                }
            }
        } else {
            string blank_gt = ".";
            if (gbwt_sample_to_phase_range.count(sample_name)) {
                auto& phase_range = gbwt_sample_to_phase_range[sample_name];
                for (int phase = phase_range.first + 1; phase <= phase_range.second; ++phase) {
                    blank_gt += "|.";
                }
            }
            v.samples[sample_name]["GT"] = {blank_gt};
            if (path_to_sample && path_restricted) {
                v.samples[sample_name]["PI"] = {blank_gt};
            }
        }
    }
    for (auto& conflict_sample : conflicts) {
        v.info["CONFLICT"].push_back(conflict_sample);
    }
}

pair<vector<int>, bool> Deconstructor::choose_traversals(const string& sample_name,
                                                         const vector<int>& travs, const vector<int>& trav_to_allele,
                                                         const vector<string>& trav_to_name,
                                                         const vector<int>& gbwt_phases) {

    assert(trav_to_name.size() == trav_to_allele.size());
    assert(gbwt_phases.size() == trav_to_name.size());    
    assert(!travs.empty());
    // count the number of times each allele comes up in a traversal
    vector<int> allele_frequencies(*max_element(trav_to_allele.begin(), trav_to_allele.end()) + 1, 0);
    for (auto trav : travs) {
        // we always want to choose alt over ref when possible in sorting logic below, so
        // cap ref frequency at 1
        int allele = trav_to_allele.at(trav);
        if (allele > 0 || allele_frequencies[allele] == 0) {
            ++allele_frequencies[allele];
        }        
    }
    // sort on frquency
    function<bool(int, int)> comp = [&] (int trav1, int trav2) {
        if (allele_frequencies[trav_to_allele[trav1]] < allele_frequencies[trav_to_allele[trav2]]) {
            return true;
        } else if (allele_frequencies[trav_to_allele[trav1]] == allele_frequencies[trav_to_allele[trav2]]) {
            // prefer non-ref when possible
            if (trav_to_allele[trav1] < trav_to_allele[trav2]) {
                return true;
            } else if (trav_to_allele[trav1] == trav_to_allele[trav2]) {
                return trav_to_name[trav1] < trav_to_name[trav2];
            }
        } 
        return false;
    };
    vector<int> sorted_travs = travs;
    std::sort(sorted_travs.begin(), sorted_travs.end(), comp);
    // find the <ploidy> most frequent traversals
    vector<int> most_frequent_travs;

    // try to pull out unique phases if available
    bool has_phasing = gbwt_sample_to_phase_range.count(sample_name) &&
        std::any_of(gbwt_phases.begin(), gbwt_phases.end(), [](int i) { return i >= 0; });
    bool phasing_conflict = false;
    int sample_ploidy = ploidy;
    int min_phase = 1;
    int max_phase = ploidy;
    if (has_phasing) {
        // override ploidy with information about all phases found in input
        std::tie(min_phase, max_phase) = gbwt_sample_to_phase_range.at(sample_name);
        // shift left by 1 unless min phase is 0
        sample_ploidy = min_phase == 0 ? max_phase + 1 : max_phase;
        assert(sample_ploidy > 0);
        
        set<int> used_phases;
        for (int i = sorted_travs.size() - 1; i >= 0 && most_frequent_travs.size() < sample_ploidy; --i) {
            int phase = gbwt_phases.at(sorted_travs[i]);
            if (!used_phases.count(phase)) {
                most_frequent_travs.push_back(sorted_travs[i]);
                used_phases.insert(phase);
            } else {
                phasing_conflict = true;
            }
        }
    } else {
        for (int i = sorted_travs.size() - 1; i >= 0 && most_frequent_travs.size() < sample_ploidy; --i) {
            most_frequent_travs.push_back(sorted_travs[i]);
        }
    }

    // sort by phase
    if (has_phasing) {
        std::sort(most_frequent_travs.begin(), most_frequent_travs.end(),
                  [&](int t1, int t2) {return gbwt_phases.at(t1) < gbwt_phases.at(t2);});
        if (max_phase > 0) {
            // pad out by phase
            assert(gbwt_phases.at(most_frequent_travs.back()) <= max_phase);
            assert(max_phase < 1000);
            // we normally expect to have phases 1,2,3, ...
            // in this case, we shift them all back, otherwise leave 0-based
            int offset = min_phase != 0 ? -1 : 0;
            vector<int> padded_travs(max_phase + 1 + offset, -1);
            for (auto ft : most_frequent_travs) {
                int phase = gbwt_phases.at(ft) + offset;
                padded_travs.at(phase) = ft;
            }
            swap(padded_travs, most_frequent_travs);
        }
    }
    // check if there's a conflict
    size_t zero_count = std::count(allele_frequencies.begin(), allele_frequencies.end(), 0);
    bool conflict = phasing_conflict || allele_frequencies.size() - zero_count > sample_ploidy;
    return make_pair(most_frequent_travs, conflict);
}
    
bool Deconstructor::deconstruct_site(const Snarl* snarl) {

    auto contents = snarl_manager->shallow_contents(snarl, *graph, false);
    if (contents.first.empty()) {
        // Nothing but the boundary nodes in this snarl
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping empty site " << pb2json(*snarl) << endl;
#endif
        return false;
    }
#ifdef debug
#pragma omp crtiical (cerr)
    cerr << "Computing traversals of site " << pb2json(*snarl) << endl;
#endif

    // find every traversal that runs through a path in the graph
    pair<vector<SnarlTraversal>, vector<pair<step_handle_t, step_handle_t> > > path_travs;
    path_travs = path_trav_finder->find_path_traversals(*snarl);
    vector<string> path_trav_names;

    // pick out the traversal corresponding to a reference path, breaking ties consistently
    string ref_trav_name;
    for (int i = 0; i < path_travs.first.size(); ++i) {
        string path_trav_name = graph->get_path_name(graph->get_path_handle_of_step(path_travs.second[i].first));
#ifdef debug
#pragma omp critical (cerr)
        {
            cerr << "Traversal " << i << ": name=" << path_trav_name << ", size=" << path_travs.first[i].visit_size()
                 << ", start=" << graph->get_position_of_step(path_travs.second[i].first)
                 << ", end=" << graph->get_position_of_step(path_travs.second[i].second) << endl
                 << " trav=" << pb2json(path_travs.first[i]) << endl;
        }
#endif
        tuple<bool, string, size_t, size_t> subpath_parse = Paths::parse_subpath_name(path_trav_name);
        if (get<0>(subpath_parse)) {
            path_trav_name = get<1>(subpath_parse);
        }
        if (ref_paths.count(path_trav_name) &&
            (ref_trav_name.empty() || path_trav_name < ref_trav_name)) {
            ref_trav_name = path_trav_name;
        }
        path_trav_names.push_back(path_trav_name);
    }
    
    // remember all the reference traversals (there can be more than one only in the case of a
    // cycle in the reference path
    vector<int> ref_travs;
    // hacky subpath support -- gets added to variant on output
    vector<int64_t> ref_offsets;
    if (!ref_trav_name.empty()) {
        for (int i = 0; i < path_travs.first.size(); ++i) {
            string path_trav_name = graph->get_path_name(graph->get_path_handle_of_step(path_travs.second[i].first));
            tuple<bool, string, size_t, size_t> subpath_parse = Paths::parse_subpath_name(path_trav_name);
            int64_t sub_offset = 0;
            if (get<0>(subpath_parse)) {
                path_trav_name = get<1>(subpath_parse);
                sub_offset = (int64_t)get<2>(subpath_parse);
            }
            if (path_trav_name == ref_trav_name) {
                ref_travs.push_back(i);
                ref_offsets.push_back(sub_offset);
            }
        }
    }

    // add in the gbwt traversals if we can
    size_t first_gbwt_trav_idx = path_travs.first.size();
    vector<gbwt::size_type> trav_thread_ids(first_gbwt_trav_idx, numeric_limits<gbwt::size_type>::max());
    vector<int64_t> gbwt_trav_offsets;
    if (gbwt_trav_finder.get() != nullptr) {
        pair<vector<SnarlTraversal>, vector<gbwt::size_type>> thread_travs = gbwt_trav_finder->find_path_traversals(*snarl);
        for (int i = 0; i < thread_travs.first.size(); ++i) {
            string gbwt_sample_name = thread_sample(gbwt_trav_finder->get_gbwt(), gbwt::Path::id(thread_travs.second[i]));
            // we count on convention of reference as embedded path above, so ignore it here
            // todo: would be nice to be more flexible...
            if (gbwt_sample_name != gbwtgraph::REFERENCE_PATH_SAMPLE_NAME) {
                string name = thread_name(gbwt_trav_finder->get_gbwt(), gbwt::Path::id(thread_travs.second[i]), true);
                
                path_trav_names.push_back(name);
                path_travs.first.push_back(thread_travs.first[i]);
                // dummy handles so we can use the same code as the named path traversals above
                path_travs.second.push_back(make_pair(step_handle_t(), step_handle_t()));
                // but we keep the thread id for later
                trav_thread_ids.push_back(thread_travs.second[i]);
                // keep the offset (which is stored in the contig field)
                gbwt_trav_offsets.push_back(thread_count(gbwt_trav_finder->get_gbwt(), gbwt::Path::id(thread_travs.second[i])));
            }
        }
    }

    // if there's no reference traversal, go fishing in the gbwt
    if (ref_travs.empty() && gbwt_trav_finder.get()) {
        int gbwt_ref_trav = -1;
        int64_t gbwt_ref_offset = 0;
        for (int i = first_gbwt_trav_idx; i < path_travs.first.size(); ++i) {
            if (ref_paths.count(path_trav_names[i]) &&
                (gbwt_ref_trav < 0 || path_trav_names[i] < path_trav_names[gbwt_ref_trav])) {
                gbwt_ref_trav = i;
                gbwt_ref_offset = gbwt_trav_offsets.at(i - first_gbwt_trav_idx);
            } 
        }
        if (gbwt_ref_trav >= 0) {
            ref_travs.push_back(gbwt_ref_trav);
            ref_offsets.push_back(gbwt_ref_offset);
            string& name = path_trav_names[gbwt_ref_trav];
            assert(name.compare(0, 8, "_thread_") != 0);
            ref_trav_name = name;
        }
    }
                                                          
    // there's no reference path through the snarl, so we can't make a variant
    // (todo: should we try to detect this before computing traversals?)
    if (ref_travs.empty()) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping site becuase no reference traversal was found " << pb2json(*snarl) << endl;
#endif
        return false;
    }

    if (!path_restricted && gbwt_trav_finder.get() == nullptr) {
        // add in the exhaustive traversals
        vector<SnarlTraversal> additional_travs;
                        
        // exhaustive traversal can't do all snarls
        if (snarl->type() != ULTRABUBBLE) {
            return false;
        }
        if (!check_max_nodes(snarl)) {
#pragma omp critical (cerr)
            cerr << "Warning: Skipping site because it is too complex for exhaustive traversal enumeration: " << pb2json(*snarl) << endl << "         Consider using -e to traverse embedded paths" << endl;
            return false;
        }
        additional_travs = explicit_exhaustive_traversals(snarl);
         
        // happens when there was a nested non-ultrabubble snarl
        if (additional_travs.empty()) {
            return false;
        }
        path_travs.first.insert(path_travs.first.end(), additional_travs.begin(), additional_travs.end());
        for (int i = 0; i < additional_travs.size(); ++i) {
            // dummy names so we can use the same code as the named path traversals above
            path_trav_names.push_back(" >>" + std::to_string(i));
            // dummy handles so we can use the same code as the named path traversals above
            path_travs.second.push_back(make_pair(step_handle_t(), step_handle_t()));
        }

    }
    
    // there's not alt path through the snarl, so we can't make an interesting variant
    if (path_travs.first.size() < 2) {
#ifdef debug
#pragma omp critical (cerr)
        cerr << "Skipping site because to alt traversal was found " << pb2json(*snarl) << endl;
#endif
        return false;
    }

    // we write a variant for every reference traversal
    for (size_t i = 0; i < ref_travs.size(); ++i) {
        auto& ref_trav_idx = ref_travs[i];
        auto& ref_trav_offset = ref_offsets[i];

        const SnarlTraversal& ref_trav = path_travs.first[ref_trav_idx];
        
        vcflib::Variant v;
        v.quality = 60;

        // write variant's sequenceName (VCF contig)
        v.sequenceName = ref_trav_name;

        // Map our snarl endpoints to oriented positions in the embedded path in the graph
        handle_t first_path_handle;
        size_t first_path_pos;
        bool use_start;
        if (ref_trav_idx < first_gbwt_trav_idx) {
            step_handle_t start_step = path_travs.second[ref_trav_idx].first;
            step_handle_t end_step = path_travs.second[ref_trav_idx].second;
            handle_t start_handle = graph->get_handle_of_step(start_step);
            handle_t end_handle = graph->get_handle_of_step(end_step);
            size_t start_pos = graph->get_position_of_step(start_step);
            size_t end_pos = graph->get_position_of_step(end_step);
            use_start = start_pos < end_pos;
            first_path_handle = use_start ? start_handle : end_handle;
            first_path_pos = use_start ? start_pos : end_pos;
        } else {
            std::tie(use_start, first_path_handle, first_path_pos) = get_gbwt_path_position(ref_trav, trav_thread_ids[ref_trav_idx]);
#ifdef debug
            cerr << "got " << use_start << " " << graph->get_id(first_path_handle) << ":" << graph->get_is_reverse(first_path_handle)
                 << " " << first_path_pos << " from gbwt for " << ref_trav_name << endl;
#endif
        }
        // Get the first visit of our snarl traversal
        const Visit& first_trav_visit = use_start ? ref_trav.visit(0) : ref_trav.visit(ref_trav.visit_size() - 1);

        char prev_char;
        if ((use_start && first_trav_visit.backward() == graph->get_is_reverse(first_path_handle)) ||
            (!use_start && first_trav_visit.backward() != graph->get_is_reverse(first_path_handle))) {
            // Our path and traversal have consistent orientation.  leave off the end of the start node going forward
            first_path_pos += graph->get_length(first_path_handle);
            prev_char = ::toupper(graph->get_sequence(first_path_handle)[graph->get_length(first_path_handle) - 1]);
        } else {
            // They are flipped: leave off the beginning of the start node going backward
            prev_char = reverse_complement(::toupper(graph->get_sequence(first_path_handle)[0]));
        }
        
        // shift from 0-based to 1-based for VCF
        first_path_pos += 1;

        v.position = first_path_pos + ref_trav_offset;

        v.id = snarl_name(snarl);
        
        // Convert the snarl traversals to strings and add them to the variant
        vector<int> trav_to_allele = get_alleles(v, path_travs.first, ref_trav_idx, prev_char, use_start);

        // Fill in the genotypes
        if (path_restricted || gbwt_trav_finder.get()) {
            get_genotypes(v, path_trav_names, trav_to_allele, trav_thread_ids);
        }

        // Fill in some snarl hierarchy information
        if (include_nested) {
            // would be nicer to do this constant time!
            size_t level = 0;
            for (const Snarl* cur = snarl; !snarl_manager->is_root(cur); cur = snarl_manager->parent_of(cur)) {
                ++level;
            }
            v.info["LV"].push_back(std::to_string(level));
            if (level > 0) {
                const Snarl* parent = snarl_manager->parent_of(snarl);
                string parent_id = snarl_name(parent);
                v.info["PS"].push_back(parent_id);
            } 
        }

        // we only bother printing out sites with at least 1 non-reference allele
        if (!std::all_of(trav_to_allele.begin(), trav_to_allele.end(), [](int i) { return i == 0; })) {
            if (path_restricted || gbwt_trav_finder.get()) {
                // run vcffixup to add some basic INFO like AC
                vcf_fixup(v);
            }
            add_variant(v);
        }
    }    
    return true;
}

/**
 * Convenience wrapper function for deconstruction of multiple paths.
 */
void Deconstructor::deconstruct(vector<string> ref_paths, const PathPositionHandleGraph* graph, SnarlManager* snarl_manager,
                                bool path_restricted_traversals, int ploidy, bool include_nested,
                                const unordered_map<string, string>* path_to_sample,
                                gbwt::GBWT* gbwt,
                                const unordered_map<nid_t, pair<nid_t, size_t>>* translation) {

    this->graph = graph;
    this->snarl_manager = snarl_manager;
    this->path_restricted = path_restricted_traversals;
    this->ploidy =ploidy;
    this->path_to_sample = path_to_sample;
    this->ref_paths = set<string>(ref_paths.begin(), ref_paths.end());
    this->include_nested = include_nested;
    this->translation = translation;
    assert(path_to_sample == nullptr || path_restricted || gbwt);
    if (gbwt) {
        this->gbwt_pos_caches.resize(get_thread_count(), nullptr);
        for (size_t i = 0; i < this->gbwt_pos_caches.size(); ++i) {
            if (this->gbwt_pos_caches[i] == nullptr) {
                this->gbwt_pos_caches[i] = new LRUCache<gbwt::size_type, shared_ptr<unordered_map<handle_t, size_t>>>(lru_size);
            }
        }
    }
    
    // Keep track of the non-reference paths in the graph.  They'll be our sample names
    sample_names.clear();
    graph->for_each_path_handle([&](const path_handle_t& path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (!this->ref_paths.count(path_name)) {
                // rely on the given map.  if a path isn't in it, it'll be ignored
                if (path_to_sample) {
                    if (path_to_sample->count(path_name)) {
                        sample_names.insert(path_to_sample->find(path_name)->second);
                    }
                    // if we have the map, we only consider paths there-in
                }
                else {
                    // no name mapping, just use every path as is
                    sample_names.insert(path_name);
                }
            }
        });
    if (gbwt) {
        // add in sample names from the gbwt
        for (size_t i = 0; i < gbwt->metadata.paths(); i++) {
            string sample_name = thread_sample(*gbwt, i);
            if (sample_name != gbwtgraph::REFERENCE_PATH_SAMPLE_NAME &&
                (path_to_sample == nullptr || path_to_sample->count(sample_name))) {
                sample_names.insert(thread_sample(*gbwt, i));
                int phase = thread_phase(*gbwt, i);
                if (!gbwt_sample_to_phase_range.count(sample_name)) {
                    gbwt_sample_to_phase_range[sample_name] = make_pair(phase, phase);
                } else {
                    pair<int, int>& phase_range = gbwt_sample_to_phase_range[sample_name];
                    phase_range.first = std::min(phase_range.first, phase);
                    phase_range.second = std::max(phase_range.second, phase);
                }
            }
        }
    }
    
    // print the VCF header
    stringstream stream;
    stream << "##fileformat=VCFv4.2" << endl;
    if (path_restricted || gbwt) {
        stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    }
    if (path_to_sample && path_restricted) {
        stream << "##FORMAT=<ID=PI,Number=.,Type=String,Description=\"Path information. Original vg path name for sample as well as its allele (can be many paths per sample)\">" << endl;
    }
    if (path_to_sample || gbwt) {
        stream << "##INFO=<ID=CONFLICT,Number=.,Type=String,Description=\"Sample names for which there are multiple paths in the graph with conflicting alleles";
        if (!gbwt) {
            stream << " (details in PI field)";
        }
        stream << "\">" << endl;
    }
    if (path_restricted || gbwt) {
        stream << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">" << endl;
        stream << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">" << endl;
        stream << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">" << endl;
        stream << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl;
    }
    if (include_nested) {
        stream << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">" << endl;
        stream << "##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">" << endl;
    }
    set<string> gbwt_ref_paths;
    for(auto& refpath : ref_paths) {
        if (graph->has_path(refpath)) {
            size_t path_len = 0;
            path_handle_t path_handle = graph->get_path_handle(refpath);
            for (handle_t handle : graph->scan_path(path_handle)) {
                path_len += graph->get_length(handle);
            }
            stream << "##contig=<ID=" << refpath << ",length=" << path_len << ">" << endl;
        } else {
            gbwt_ref_paths.insert(refpath);
        }       
    }
    if (!gbwt_ref_paths.empty()) {
        unordered_map<string, vector<gbwt::size_type>> gbwt_name_to_ids;
        for (size_t i = 0; i < gbwt->metadata.paths(); i++) {
            gbwt_name_to_ids[thread_name(*gbwt, i, true)].push_back(i);
        }
        for (const string& refpath : gbwt_ref_paths) {
            vector<gbwt::size_type>& thread_ids = gbwt_name_to_ids.at(refpath);
            size_t path_len = 0;
            for (gbwt::size_type thread_id : thread_ids) {
                size_t offset = thread_count(*gbwt, thread_id);
                size_t len = path_to_length(extract_gbwt_path(*graph, *gbwt, thread_id));
                path_len = std::max(path_len, offset + len);
            }
            stream << "##contig=<ID=" << refpath << ",length=" << path_len << ">" << endl;
        }
    }
    
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    if (path_restricted || gbwt) {
        for (auto& sample_name : sample_names) {
            stream << "\t" << sample_name;
        }
    }
    stream << endl;
    
    string hstr = stream.str();
    assert(output_vcf.openForOutput(hstr));
    cout << output_vcf.header << endl;

    // create the traversal finder
    map<string, const Alignment*> reads_by_name;
    path_trav_finder = unique_ptr<PathTraversalFinder>(new PathTraversalFinder(*graph,
                                                                               *snarl_manager));
    
    if (!path_restricted && !gbwt) {
        trav_finder = unique_ptr<TraversalFinder>(new ExhaustiveTraversalFinder(*graph,
                                                                                *snarl_manager,
                                                                                true));

    }
    
    if (gbwt != nullptr) {
        gbwt_trav_finder = unique_ptr<GBWTTraversalFinder>(new GBWTTraversalFinder(*graph, *gbwt));
    }
    
    // Do the top-level snarls in parallel
    snarl_manager->for_each_top_level_snarl_parallel([&](const Snarl* snarl) {
            vector<const Snarl*> todo(1, snarl);
            vector<const Snarl*> next;
            while (!todo.empty()) {
                for (auto next_snarl : todo) {
                    // if we can't make a variant from the snarl due to not finding
                    // paths through it, we try again on the children
                    // note: we may want to push the parallelism down a bit 
                    if (!deconstruct_site(next_snarl) || include_nested) {
                        const vector<const Snarl*>& children = snarl_manager->children_of(next_snarl);
                        next.insert(next.end(), children.begin(), children.end());
                    }
                }
                swap(todo, next);
                next.clear();
            }
        });
    
    // write variants in sorted order
    write_variants(cout);
}

bool Deconstructor::check_max_nodes(const Snarl* snarl)  {
    unordered_set<id_t> nodeset = snarl_manager->deep_contents(snarl, *graph, false).first;
    int node_count = 0;
    for (auto node_id : nodeset) {
        handle_t node = graph->get_handle(node_id);
        if (graph->get_degree(node, true) > 1 || graph->get_degree(node, false) > 1) {
            ++node_count;
            if (node_count > max_nodes_for_exhaustive) {
                return false;
            }
        }
    }
    return true;
};

vector<SnarlTraversal> Deconstructor::explicit_exhaustive_traversals(const Snarl* snarl){
    vector<SnarlTraversal> out_travs;
    bool ultra_all_the_way_down = true;
    function<void(const SnarlTraversal&, const Snarl&)> extend_trav =
        [&](const SnarlTraversal& trav, const Snarl& nested_snarl) {
        // exhaustive traversal finder is limited.  if we find something
        // that's not an ultrabubble, not much we can do
        if (nested_snarl.type() != ULTRABUBBLE) {
            ultra_all_the_way_down = false;
            return;
        }
        vector<SnarlTraversal> nested_travs = trav_finder->find_traversals(nested_snarl);
        for (auto& nested_trav : nested_travs) {
            SnarlTraversal extended_trav = trav;
            bool is_explicit = true;
            for (int i = 0; i < nested_trav.visit_size(); ++i) {
                if (nested_trav.visit(i).node_id() != 0) {
                    Visit* visit = extended_trav.add_visit();
                    *visit = nested_trav.visit(i);
                } else {
                    extend_trav(extended_trav, nested_trav.visit(i).snarl());
                    is_explicit = false;
                }
            }
            if (is_explicit) {
                out_travs.push_back(extended_trav);
            }
        }
    };
    SnarlTraversal trav;
    extend_trav(trav, *snarl);
    if (!ultra_all_the_way_down) {
        out_travs.clear();
    }        
    return out_travs;
}

string Deconstructor::snarl_name(const Snarl* snarl) {
    nid_t start_node = snarl->start().node_id();
    nid_t end_node = snarl->end().node_id();
    if (translation) {
        auto i = translation->find(start_node);
        if (i == translation->end()) {
            throw runtime_error("Error [vg deconstruct]: Unable to find node " + std::to_string(start_node) + " in translation file");
        }
        start_node = i->second.first;
        i = translation->find(end_node);
        if (i == translation->end()) {
            throw runtime_error("Error [vg deconstruct]: Unable to find node " + std::to_string(end_node) + " in translation file");
        }
        end_node = i->second.first;
    }
    return (snarl->start().backward() ? "<" : ">") + std::to_string(start_node) +
        (snarl->end().backward() ? "<" : ">") + std::to_string(end_node);
}

tuple<bool, handle_t, size_t> Deconstructor::get_gbwt_path_position(const SnarlTraversal& trav, const gbwt::size_type& thread) {

    // scan the whole thread in order to get the path positions -- there's no other way unless we build an index a priori
    const gbwt::GBWT& gbwt = gbwt_trav_finder->get_gbwt();
    bool thread_reversed = gbwt::Path::is_reverse(thread);
    gbwt::size_type sequence_id = gbwt::Path::encode(gbwt::Path::id(thread), false);
        ;
    // the handles we're looking out for when walking along the thread
    handle_t start_handle;
    handle_t end_handle;
    const Visit& v1 = trav.visit(0);
    const Visit& v2 = trav.visit(trav.visit_size() - 1);
    if (thread_reversed) {
        start_handle = graph->get_handle(v2.node_id(), !v2.backward());
        end_handle = graph->get_handle(v1.node_id(), !v1.backward());
    } else {
        start_handle = graph->get_handle(v1.node_id(), v1.backward());
        end_handle = graph->get_handle(v2.node_id(), v2.backward());
    }
    int64_t start_offset = -1;
    int64_t end_offset = -1;
    int64_t offset = 0;

    LRUCache<gbwt::size_type, shared_ptr<unordered_map<handle_t, size_t>>>* gbwt_pos_cache =
        gbwt_pos_caches[omp_get_thread_num()];
    
    pair<shared_ptr<unordered_map<handle_t, size_t>>, bool> cached = gbwt_pos_cache->retrieve(thread);
    if (cached.second) {
        start_offset = cached.first->at(start_handle);
    } else {
        shared_ptr<unordered_map<handle_t, size_t>> path_map = make_shared<unordered_map<handle_t, size_t>>();
        for (gbwt::edge_type pos = gbwt.start(sequence_id); pos.first != gbwt::ENDMARKER; pos = gbwt.LF(pos)) {
            handle_t handle = graph->get_handle(gbwt::Node::id(pos.first), gbwt::Node::is_reverse(pos.first));
            path_map->insert(make_pair(handle, offset));
            if (handle == start_handle) {
                start_offset = offset;
            }
            if (handle == end_handle && start_offset != -1) {
                end_offset = offset;
            }
            size_t len = graph->get_length(gbwt_to_handle(*graph, pos.first));
            offset += len;
        }
        assert(start_offset >= 0 && end_offset >= 0);
        gbwt_pos_cache->put(thread, path_map);
    }
  
    auto rval = make_tuple<bool, handle_t, size_t>((bool)!thread_reversed, (handle_t)start_handle, (size_t)start_offset);

    return rval;
}

}

