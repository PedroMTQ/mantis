import re

from unifunc.source import UniFunc

from mantis.cython_src.get_non_overlapping_hits import get_non_overlapping_hits


class Consensus():

    def __init__(self):
        self.unifunc = UniFunc()

    def get_ref_weight(self, ref):
        '''
        NCBI
        NOG
        Pfam-A
        kofam
        '''
        if re.search('NCBI[GT]', ref):
            corrected_ref = 'ncbi'
        elif re.search('NOG[GT]', ref):
            corrected_ref = 'nog'
        elif 'Pfam-A' in ref:
            corrected_ref = 'pfam'
        elif ref == 'kofam_merged':
            corrected_ref = 'kofam'
        else:
            corrected_ref = str(ref)

        if corrected_ref in self.mantis_ref_weights:
            return self.mantis_ref_weights[corrected_ref]
        else:
            return self.mantis_ref_weights['else']

    # this will add the descriptions from the respective GO IDs, but we only add them if we dont have too many GO IDs
    def add_from_go_obo(self, query_dict):
        for ref_file in query_dict:
            for ref_hit in query_dict[ref_file]:
                clean_query_identifiers = self.remove_unspecific_identifiers(
                    query_dict[ref_file][ref_hit]['identifiers'])
                temp_ref_hit_dict_identifiers = set(query_dict[ref_file][ref_hit]['identifiers'])
                temp_ref_hit_dict_description = set(query_dict[ref_file][ref_hit]['description'])
                for id_str in clean_query_identifiers:
                    if 'go:' in id_str[0:3]:
                        if id_str in self.go_dict:
                            extra_info = self.go_dict[id_str]
                            for go_id in extra_info['identifiers']:
                                temp_ref_hit_dict_identifiers.add(go_id)
                            for syn in extra_info['synonyms']:
                                temp_ref_hit_dict_description.add(syn)
                # we use the temp info to help matching annotations, but we dont want to add them to the final annotations
                query_dict[ref_file][ref_hit]['temp_identifiers'] = set(temp_ref_hit_dict_identifiers)
                for id_str in temp_ref_hit_dict_identifiers:
                    # some annotations have unspecific ec numbers, if so, we remove them from the similarity analysis
                    if 'enzyme_ec:' in id_str and '-' in id_str:
                        query_dict[ref_file][ref_hit]['temp_identifiers'].remove(id_str)
                query_dict[ref_file][ref_hit]['temp_description'] = temp_ref_hit_dict_description

    def read_interpreted_annotation(self, interpreted_annotation_tsv):
        with open(interpreted_annotation_tsv) as file:
            line = file.readline()
            line = file.readline()
            dict_annotations = {}
            while line:
                line = line.strip('\n')
                line = line.strip()
                line = line.split('\t')
                query, ref_file, ref_hit, ref_hit_accession, evalue, bitscore, direction, query_len, query_start, query_end, ref_start, ref_end, ref_len = line[0:13]
                annotations = line[14:]
                if query not in dict_annotations: dict_annotations[query] = {'query_len': int(query_len),
                                                                             'ref_files': {}}
                # since the same hmm file can have different hits for the same query
                if ref_file not in dict_annotations[query]['ref_files']:
                    dict_annotations[query]['ref_files'][ref_file] = {}
                dict_annotations[query]['ref_files'][ref_file][ref_hit] = {'identifiers': set(), 'description': set(),
                                                                           'evalue': float(evalue),
                                                                           'bitscore': float(bitscore),
                                                                           'direction': direction,
                                                                           'ref_hit_accession': ref_hit_accession,
                                                                           'query_start': int(query_start),
                                                                           'query_end': int(query_end),
                                                                           'ref_start': int(ref_start),
                                                                           'ref_end': int(ref_end),
                                                                           'ref_len': int(ref_len),
                                                                           }
                for annot in annotations:
                    if 'description:' in annot[0:len('description:')]:
                        annot_text = annot.split(':')[1:]
                        annot_text = ' '.join(annot_text)
                        dict_annotations[query]['ref_files'][ref_file][ref_hit]['description'].add(annot_text)
                    elif 'is_essential_gene:' in annot:
                        dict_annotations[query]['is_essential_gene'] = True
                    else:
                        dict_annotations[query]['ref_files'][ref_file][ref_hit]['identifiers'].add(annot)
                line = file.readline()
        for query in dict_annotations:
            self.add_from_go_obo(dict_annotations[query]['ref_files'])
        return dict_annotations

    def generate_gff_line_consensus(self, query,
                                    dict_annotations,
                                    is_essential,
                                    consensus_hits,
                                    total_hits,
                                    ref_files_consensus,
                                    ref_names_consensus):
        # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        # verified with http://genometools.org/cgi-bin/gff3validator.cgi
        for ref_file in dict_annotations:
            ref_annotations = dict_annotations[ref_file]
            for ref_id in ref_annotations:
                ref_id_annotations = dict_annotations[ref_file][ref_id]
                query_start = str(ref_id_annotations['query_start'])
                query_end = str(ref_id_annotations['query_end'])
                evalue = str(ref_id_annotations['evalue'])
                direction = ref_id_annotations['direction']
                ref_start = ref_id_annotations['ref_start']
                ref_end = ref_id_annotations['ref_end']
                ref_len = ref_id_annotations['ref_len']
                hit_accession = ref_id_annotations['ref_hit_accession']
                hit_identifiers = ref_id_annotations['identifiers']
                hit_description = ref_id_annotations['description']
                gff_line = '\t'.join([query, 'Mantis', 'CDS', query_start, query_end, evalue, direction, '0']) + '\t'
                if hit_accession != '-':
                    attributes = f'Name={ref_id};Target={ref_id} {ref_start} {ref_end};Alias={hit_accession}'
                else:
                    attributes = f'Name={ref_id};Target={ref_id} {ref_start} {ref_end}'
                notes = f'Note=ref_file:{ref_file},ref_len:{ref_len}'
                descriptions = []
                for d in hit_description:
                    descriptions.append(f'description:{d}')
                if descriptions:
                    notes += ',' + ','.join(descriptions)
                if is_essential:
                    notes += ',is_essential_gene:True'

                dbxref = []
                ontology_terms = []
                for i in hit_identifiers:
                    if not i.startswith('go:'):
                        dbxref.append(i)
                    else:
                        ontology_terms.append(i)
                all_annotations = None
                if dbxref and ontology_terms:
                    all_annotations = 'Dbxref=' + ','.join(dbxref) + ';' + 'Ontology_term=' + ','.join(ontology_terms)
                elif dbxref and not ontology_terms:
                    all_annotations = 'Dbxref=' + ','.join(dbxref)
                elif not dbxref and ontology_terms:
                    all_annotations = 'Ontology_term=' + ','.join(ontology_terms)
                if all_annotations:
                    gff_line += ';'.join([attributes, notes, all_annotations])
                else:
                    gff_line += ';'.join([attributes, notes])

                yield gff_line

    def yield_seq_regions(self, interpreted_annotation_tsv):
        already_added = set()
        with open(interpreted_annotation_tsv) as file:
            file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.strip()
                line = line.split('\t')
                query, ref_file, ref_hit, ref_hit_accession, evalue, bitscore, direction, query_len, query_start, query_end, ref_start, ref_end, ref_len = line[0:13]
                if query not in already_added:
                    already_added.add(query)
                    yield f'##sequence-region {query} 1 {query_len}\n'
                line = file.readline()

    def generate_consensus(self, interpreted_annotation_tsv, stdout_file_path=None):
        dict_annotations = self.read_interpreted_annotation(interpreted_annotation_tsv)
        self.query_match_hits = {}
        for query in dict_annotations:
            hits = []
            query_dict = dict_annotations[query]
            query_len = query_dict['query_len']
            ref_files = query_dict['ref_files']
            for ref_file in ref_files:
                for ref_hit in ref_files[ref_file]:
                    hits.append([ref_file, ref_hit, ref_files[ref_file][ref_hit]])
            if self.domain_algorithm == 'heuristic':
                best_hits = self.get_best_hits_approximation(list(hits), sorting_class='consensus',
                                                             sorting_type=self.sorting_type)
            elif self.domain_algorithm == 'bpo':
                best_hits = self.get_lowest_hit(list(hits), sorting_class='consensus', sorting_type=self.sorting_type)
            else:
                best_hits = self.get_best_cython_consensus(list(hits),
                                                           query_len,
                                                           query,
                                                           stdout_file_path=stdout_file_path)
            # expanding hits
            consensus_hits, total_hits, ref_files_consensus, ref_names_consensus = self.expand_best_combination(best_hits, dict_annotations[query]['ref_files'])
            if 'is_essential_gene' in dict_annotations[query]:
                is_essential = True
            else:
                is_essential = False
            consensus_line = self.generate_consensus_line(query,
                                                          dict_annotations[query]['ref_files'],
                                                          is_essential,
                                                          consensus_hits,
                                                          total_hits,
                                                          ref_files_consensus,
                                                          ref_names_consensus)
            gff_lines = []
            if self.output_gff:
                gff_lines = self.generate_gff_line_consensus(query,
                                                             dict_annotations[query]['ref_files'],
                                                             is_essential,
                                                             consensus_hits,
                                                             total_hits,
                                                             ref_files_consensus,
                                                             ref_names_consensus)
            yield consensus_line, gff_lines

    # same method as non_overlapping hits from Homology_processor , with some minor alterations to account for consensus
    # @timeit_function
    def get_best_cython_consensus(self, hits, query_len, query, stdout_file_path=None):
        try:
            best_hits = self.get_best_hits_Consensus(list(hits), query_len)
        except (TimeoutError, RecursionError):
            best_hits = self.get_best_hits_approximation(list(hits), sorting_class='consensus', sorting_type='evalue')
            stdout_file = open(stdout_file_path, 'a+')
            print(f'Query {query} was approximated during consensus generation', flush=True, file=stdout_file)
            stdout_file.close()
        return best_hits

    def query_hits_to_cython_Consensus(self, query_hits):
        conversion_dict = {}
        res = set()
        for hit_i in range(len(query_hits)):
            ref_file = query_hits[hit_i][0]
            ref_hit = query_hits[hit_i][1]
            hit_info = query_hits[hit_i][2]
            query_start, query_end = recalculate_coordinates(hit_info['query_start'],
                                                             hit_info['query_end'],
                                                             self.overlap_value)
            res.add(tuple([
                hit_i,
                query_start,
                query_end,
                ref_hit
            ]))
            conversion_dict[hit_i] = [ref_file, ref_hit, hit_info]
        return res, conversion_dict


    def get_min_max_alt_alg(self, query_hits):
        all_bitscore, all_evalue = [], []
        for hit in query_hits:
            ref_file, ref_hit, hit_info = hit
            if hit_info['evalue']:
                all_evalue.append(hit_info['evalue'])
            if hit_info['bitscore']:
                all_bitscore.append(hit_info['bitscore'])
        min_val_bitscore, max_val_bitscore = log10(min(all_bitscore)), log10(max(all_bitscore))
        if all_evalue:
            max_val_evalue, min_val_evalue = log10(min(all_evalue)), log10(max(all_evalue))
        else:
            min_val_evalue, max_val_evalue = 0, 0
        return min_val_evalue, max_val_evalue, min_val_bitscore, max_val_bitscore

    def add_scaled_values(self, query_hits):
        min_val_evalue, max_val_evalue, min_val_bitscore, max_val_bitscore = self.get_min_max_alt_alg(query_hits)
        for hit in query_hits:
            ref_file, ref_hit, hit_info = hit
            hit_weight = self.get_ref_weight(ref_file)
            hit_annotation_score = self.get_annotation_score(hit_info)
            if hit_info['evalue']:
                hit_info['scaled_evalue'] = min_max_scale(log10(hit_info['evalue']),
                                                          min_val_evalue,
                                                          max_val_evalue) + 0.01
            else:
                hit_info['scaled_evalue'] = 2
            hit_info['scaled_bitscore'] = min_max_scale(log10(hit_info['bitscore']),
                                                        min_val_bitscore,
                                                        max_val_bitscore) + 0.01

            hit_info['scaled_evalue'] *= (hit_annotation_score + hit_weight) / 2
            hit_info['scaled_bitscore'] *= (hit_annotation_score + hit_weight) / 2

    def get_annotation_score(self, hit_info):
        if hit_info['temp_identifiers'] and hit_info['temp_description']:
            annotation_score = 1
        elif hit_info['temp_identifiers']:
            annotation_score = 0.75
        elif hit_info['temp_description']:
            annotation_score = 0.5
        else:
            annotation_score = 0.25
        return annotation_score

    def get_intersecting_sources(self, combo, query_hits):
        res = 0
        for hit in combo:
            hit_source, hit_name, hit_info = hit
            for q_hit in query_hits:
                q_hit_source, q_hit_name, q_hit_info = q_hit
                if hit_source != q_hit_source:
                    if self.is_hit_match(hit_info, hit_name, hit_source, q_hit_info, q_hit_name, q_hit_source):
                        res += 1
        return (res + len(combo)) / len(query_hits)

    def get_combo_score_Consensus(self, combo, query_length, min_val, max_val, ref_sources):
        '''
        ref_sources= number of intersecting sources divided by number of hits from all sources
        metrics for consensus:
        intersecting sources = amount of different sources that intersect with the combo
        combo coverage  = percentage of the sequence covered by the combo
        combo evalue = summed scaled evalue of the combo. Is averaged
        combo weight = summed ref weight of the combo. Better sources should give better combinations of hits . Is averaged
        combo annotation score = score of the annotation (0.25 if just the ref name, 0.5 if just the description, 0.75 if just identifiers, 1 if identifiers and description
        '''
        average_value = 0
        average_hit_coverage = 0
        total_weight = 0
        annotation_score = 0
        hit_ranges = []

        for hit in combo:
            ref_file, ref_hit, hit_info = hit
            hit_start, hit_end, ref_start, ref_end = hit_info['query_start'], hit_info['query_end'], \
                                                     hit_info['ref_start'], hit_info['ref_end']
            ref_len = hit_info['ref_len']
            hit_evalue = hit_info['evalue']
            hit_bitscore = hit_info['bitscore']
            hit_coverage = (hit_end - hit_start + 1) / query_length
            # ref_coverage = (ref_end - ref_start + 1) / ref_len
            average_hit_coverage += hit_coverage
            # lower is best
            if self.sorting_type == 'evalue':
                if hit_evalue:
                    scaled_value = min_max_scale(log10(hit_evalue), min_val, max_val) + 0.01
                else:
                    scaled_value = 2
            # higher is best
            elif self.sorting_type == 'bitscore':
                if hit_bitscore:
                    scaled_value = min_max_scale(log10(hit_bitscore), min_val, max_val) + 0.01
                else:
                    scaled_value = 0
            average_value += scaled_value

            hit_ranges.append([hit_start, hit_end])
            total_weight += self.get_ref_weight(ref_file)
            annotation_score += self.get_annotation_score(hit_info)

        # sources score
        sources_score = (ref_sources + annotation_score / len(combo) + total_weight / len(combo)) / 3
        # sequence hit quality
        average_value = average_value / len(combo)
        # average coverage of all hits
        average_hit_coverage /= len(combo)
        # coverage of combination
        combination_coverage = get_combination_ranges(hit_ranges)
        combination_coverage /= query_length
        if self.best_combo_formula == 1:  # default
            combo_score = sources_score * average_value
        elif self.best_combo_formula == 2:
            combo_score = sources_score * average_value * combination_coverage * average_hit_coverage
        elif self.best_combo_formula == 3:
            combo_score = (sources_score * average_value * combination_coverage * average_hit_coverage) ** (1 / 4)
        elif self.best_combo_formula == 4:
            combo_score = sources_score * average_value * combination_coverage
        elif self.best_combo_formula == 5:
            combo_score = (sources_score * average_value * combination_coverage) ** (1 / 3)
        elif self.best_combo_formula == 6:
            combo_score = (sources_score * average_value + average_hit_coverage + combination_coverage) / 3
        elif self.best_combo_formula == 7:
            combo_score = (sources_score * 2 * average_value + average_hit_coverage + combination_coverage) / 4
        elif self.best_combo_formula == 8:
            combo_score = (sources_score * 2 * average_value + combination_coverage) / 3
        elif self.best_combo_formula == 9:
            combo_score = (sources_score * 3 * average_value + average_hit_coverage + combination_coverage) / 5
        elif self.best_combo_formula == 10:
            combo_score = (sources_score * 3 * average_value + combination_coverage) / 4
        elif self.best_combo_formula == 11:
            combo_score = average_value * sources_score * (combination_coverage) ** 2 * (average_hit_coverage) ** 2
        elif self.best_combo_formula == 12:
            combo_score = average_value * sources_score * (combination_coverage) ** 2

        return combo_score

    def get_best_hits_Consensus(self, query_hits, query_length):
        cython_hits, conversion_dict = self.query_hits_to_cython_Consensus(query_hits)
        cython_possible_combos = get_non_overlapping_hits(cython_hits, time_limit=self.time_limit)
        # sometimes this calculation is not feasible (when we have too many small hits and the query sequence is too long, the calculation of conbinations would take forever)
        if not cython_possible_combos: return None
        best_combo_score = 0
        best_combo = None
        min_val, max_val = self.get_min_max_dfs(cython_possible_combos, conversion_dict, sorting_class='consensus')
        for len_combo in cython_possible_combos:
            if len_combo:
                for cython_combo in cython_possible_combos[len_combo]:
                    combo = self.cython_to_query_hits(cython_combo, conversion_dict)
                    # we benefit the combo which intersect with other sources
                    ref_sources = self.get_intersecting_sources(combo, query_hits)
                    combo_score = self.get_combo_score_Consensus(combo, query_length, min_val, max_val, ref_sources)
                    if best_combo is None:
                        best_combo = combo
                        best_combo_score = combo_score
                    elif combo_score > best_combo_score:
                        best_combo_score = combo_score
                        best_combo = combo
        return best_combo

    # @timeit_function
    def expand_best_combination(self, best_hits, query_dict):
        hits_merged = set()
        ref_files_consensus = []
        ref_names_consensus = []
        for best_hit in best_hits:
            best_hit_file, best_hit_name, best_hit_info = best_hit
            hits_merged.add(best_hit_name)
            ref_files_consensus.append(best_hit_file)
            ref_names_consensus.append(best_hit_name)
            temp_best_hit_info = dict(best_hit_info)
            first_check = False
            # iterations might change already iterated over best_hit_info, this way we check if there is any info added, if so, we repeat the cycle one more otherwise we return
            while temp_best_hit_info != best_hit_info or not first_check:
                for ref_file in query_dict:
                    for hit_to_test in query_dict[ref_file]:
                        if self.is_hit_match(query_dict[ref_file][hit_to_test], hit_to_test, ref_file, best_hit_info,
                                             best_hit_name, best_hit_file):
                            self.add_to_hit(query_dict[ref_file][hit_to_test], best_hit_info)
                            hits_merged.add(hit_to_test)
                            ref_files_consensus.append(ref_file)
                            ref_names_consensus.append(hit_to_test)
                first_check = True
                temp_best_hit_info = dict(best_hit_info)
        # checking how many hits we managed to merge out of all the hits - coverage_consensus
        total_hits = 0
        for ref_file in query_dict:
            for i in query_dict[ref_file]:
                total_hits += 1
        return str(len(hits_merged)), str(total_hits), ref_files_consensus, ref_names_consensus

    def is_nlp_match(self, hit1_info_description, hit2_info_description):
        if self.no_unifunc:
            return False
        for hit1_d in hit1_info_description:
            for hit2_d in hit2_info_description:
                score = self.unifunc.get_similarity_score(hit1_d, hit2_d, only_return=True, verbose=False)
                if score > self.nlp_threshold:
                    return True
        return False

    # eggnog is not specific enough with the go identifiers, so if there is more than one, we dont use it for the consensus
    def remove_unspecific_identifiers(self, identifiers, specificity_threshold=10):
        res = set()
        dict_counter = {}
        for i in identifiers:
            link_type = i.split(':')[0]
            if link_type not in dict_counter: dict_counter[link_type] = 0
            dict_counter[link_type] += 1
        for i in identifiers:
            link_type = i.split(':')[0]
            if dict_counter[link_type] <= specificity_threshold:
                res.add(i)
        return res

    def is_overlap_coordinates(self, hit1, hit2):
        if hit1[1] - hit1[0] < hit2[1] - hit2[0]:
            smaller_hit = hit1
        else:
            smaller_hit = hit2
        min_len = smaller_hit[1] - smaller_hit[0] + 1
        hit1_set = set(range(hit1[0], hit1[1] + 1))
        hit2_set = set(range(hit2[0], hit2[1] + 1))
        intersection_hits = hit1_set.intersection(hit2_set)
        if len(intersection_hits) / min_len > self.minimum_consensus_overlap:
            return True
        return False

    # this makes mantis more efficient since it saves matches in memory, however it also makes it more ram intensive
    def is_hit_match(self, hit1_info, hit1_name, hit1_source, hit2_info, hit2_name, hit2_source):
        if hit1_source == hit2_source and hit1_name == hit2_name:
            return False
        if self.no_consensus_expansion:
            return False
        # cleaning up memory
        if float(psutil.virtual_memory().percent) > 95: self.query_match_hits = {}
        # to avoid adding duplicate entries
        sorted_hits = sorted([hit1_source, hit2_source])
        if hit1_source == sorted_hits[0]:
            first_hit_source = hit1_source
            first_hit_name = hit1_name
            second_hit_source = hit2_source
            second_hit_name = hit2_name
        else:
            first_hit_source = hit2_source
            first_hit_name = hit2_name
            second_hit_source = hit1_source
            second_hit_name = hit1_name
        dict_key = first_hit_source + '_' + first_hit_name + '__' + second_hit_source + '_' + second_hit_name
        sufficient_coord_overlap = self.is_overlap_coordinates([hit1_info['query_start'], hit1_info['query_end']],
                                                               [hit2_info['query_start'], hit2_info['query_end']])
        if dict_key in self.query_match_hits:
            if sufficient_coord_overlap:
                return self.query_match_hits[dict_key]
            else:
                return False
        # now to actually check
        hit1_info_description = hit1_info['temp_description']
        hit1_info_identifiers = hit1_info['temp_identifiers']
        hit2_info_description = hit2_info['temp_description']
        hit2_info_identifiers = hit2_info['temp_identifiers']
        temp_hit1_info_identifiers = self.remove_unspecific_identifiers(hit1_info_identifiers)
        temp_hit2_info_identifiers = self.remove_unspecific_identifiers(hit2_info_identifiers)

        if sufficient_coord_overlap:
            if temp_hit1_info_identifiers.intersection(temp_hit2_info_identifiers):
                self.query_match_hits[dict_key] = True
                return True
            elif self.is_nlp_match(hit1_info_description, hit2_info_description):
                self.query_match_hits[dict_key] = True
                return True
            else:
                self.query_match_hits[dict_key] = False
                return False
        else:
            return False

    def add_to_hit(self, hit_to_test_info, best_hit_info):
        best_hit_info['identifiers'].update(hit_to_test_info['identifiers'])
        best_hit_info['description'].update(hit_to_test_info['description'])

    # to remove bad descriptions
    def remove_trash_descriptions(self, all_descriptions):
        res = set()
        for d in all_descriptions:
            current_d = d.strip().lower()
            if d.strip().lower() not in [
                'enzyme',
                'domain',
                'protein',
                'unknown function',
                'domain of unknown function',
                'protein of unknown function',
                'uncharacterised protein family',
                'unknown protein family',
                'uncharacterized protein',
                'uncharacterised protein',
                'uncharacterised conserved protein',
                'uncharacterized conserved protein',
                'hypothetical protein',
            ]:
                if re.search(r'(protein|domain|domian|family|repeat|short repeats|region) (of|with) (unknown|unknwon) function(\s\(?[dp]uf\d{2,}\)?)?', current_d):
                    pass
                else:
                    res.add(d)
        return res

    # to remove redundant descriptions
    def remove_redundant_descriptions(self, all_descriptions):
        res = set()
        already_added = set()
        unspecific_tokens = ['enzyme', 'protein', 'domain']
        all_descriptions = sorted(all_descriptions)
        for d in all_descriptions:
            test = d.lower()
            for p in set(punctuation):
                test = test.replace(p, '')
            test = test.replace('\'', '')
            test = test.replace('\"', '')
            test = test.strip()
            test = test.split(' ')
            test = [i.strip() for i in test if i not in unspecific_tokens]
            test = ' '.join(test)

            if test not in already_added:
                res.add(d)
                already_added.add(test)
        res = sorted(res)
        return res

    def sort_ref_files_and_hits(self, ref_files_consensus, ref_names_consensus):
        res = {}
        for i in range(len(ref_files_consensus)):
            r_file = ref_files_consensus[i]
            if r_file not in res: res[r_file] = set()
            res[r_file].add(ref_names_consensus[i])
        ref_files = []
        ref_hits = []

        for r_file in sorted(res):
            for r_hit in sorted(res[r_file]):
                ref_files.append(r_file)
                ref_hits.append(r_hit)

        ref_files = ';'.join(ref_files)
        ref_hits = ';'.join(ref_hits)
        return ref_files, ref_hits

    def generate_consensus_line(self, query, query_dict, is_essential, consensus_hits, total_hits, ref_files_consensus,
                                ref_names_consensus):
        ref_files, ref_hits = self.sort_ref_files_and_hits(ref_files_consensus, ref_names_consensus)
        # consensus_coverage is a str with consensus sources/all sources, better as a str instead of float as its easier to understand
        row_start = [query, ref_files, ref_hits, consensus_hits, total_hits, '|']
        row_start = [str(i) for i in row_start]
        all_descriptions = set()
        all_identifiers = set()
        for ref_file in ref_files_consensus:
            for ref_hit_name in ref_names_consensus:
                if ref_hit_name in query_dict[ref_file]:
                    description = query_dict[ref_file][ref_hit_name]['description']
                    identifiers = query_dict[ref_file][ref_hit_name]['identifiers']
                    all_identifiers.update(identifiers)
                    all_descriptions.update(description)
        res = []
        sorted_identifiers = sorted(all_identifiers)
        for link in sorted_identifiers:
            if 'enzyme_ec' in link:
                if link not in res:
                    res.append(link)
        # so that description always comes in the end
        for link in sorted_identifiers:
            if 'kegg_map_lineage' not in link:
                if link not in res:
                    res.append(link)
        for link in sorted_identifiers:
            if 'kegg_map_lineage' in link:
                if link not in res:
                    res.append(link)
        if is_essential:
            res.append('is_essential_gene:True')
        clean_descriptions = self.remove_redundant_descriptions(all_descriptions)
        clean_descriptions = self.remove_trash_descriptions(clean_descriptions)
        for link in clean_descriptions:
            if 'description:' + link not in res:
                res.append('description:' + link)
        res = sorted(res)
        res = row_start + res
        return res

    def generate_consensus_output(self, interpreted_annotation_tsv, consensus_annotation_tsv, stdout_file_path=None):
        first_line = ['Query',
                      'Ref_Files',
                      'Ref_Hits',
                      'Consensus_hits',
                      'Total_hits',
                      '|',
                      'Links',
                      ]
        first_line = '\t'.join(first_line)
        output_file = open(consensus_annotation_tsv, 'w+')
        output_file.write(first_line + '\n')

        gff_file = None
        gff_file_path = consensus_annotation_tsv.replace('.tsv', '.gff')
        if self.output_gff:
            gff_file = open(gff_file_path, 'w+')
            gff_file.write('##gff-version 3' + '\n')

        consensus_annotation = self.generate_consensus(interpreted_annotation_tsv, stdout_file_path=stdout_file_path)
        if self.output_gff:
            seq_regions = self.yield_seq_regions(interpreted_annotation_tsv)
            for sr in seq_regions:
                gff_file.write(sr)

        for consensus_query, gff_lines in consensus_annotation:
            line = '\t'.join(consensus_query)
            output_file.write(line + '\n')
            if self.output_gff:
                for gff_l in gff_lines:
                    gff_file.write(gff_l + '\n')

        if gff_file:
            gff_file.close()
        output_file.close()


if __name__ == '__main__':
    m = Consensus()
