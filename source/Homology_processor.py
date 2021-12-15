try:
    from source.Assembler import *
except:
    from Assembler import *

try:
    from source.cython_src.get_non_overlapping_hits import get_non_overlapping_hits
except:
    try:
        from cython_src.get_non_overlapping_hits import get_non_overlapping_hits
    except:
        if not cython_compiled():
            compile_cython()
            try:
                from source.cython_src.get_non_overlapping_hits import get_non_overlapping_hits
            except:
                try:
                    from cython_src.get_non_overlapping_hits import get_non_overlapping_hits
                except:
                    kill_switch(CythonNotCompiled)
        # when cython's version is not compatible with the current python version, we need to recompile it
        else:
            compile_cython()
            try:
                from source.cython_src.get_non_overlapping_hits import get_non_overlapping_hits
            except:
                try:
                    from cython_src.get_non_overlapping_hits import get_non_overlapping_hits
                except:
                    kill_switch(CythonNotCompiled)

# This class will process the domtblout output and get the best hit for our queries, it is inherited by the MANTIS

class Homology_processor():

    def create_chunk_output_dirs(self, chunk_dir):
        to_create=['searchout']
        if self.keep_files:
            to_create.append('output_hmmer')
        for hmmer_dir in to_create:
            if not file_exists(chunk_dir + hmmer_dir):
                Path(chunk_dir + hmmer_dir).mkdir(parents=True, exist_ok=True)

    ######SPLITTING TARGET FASTA INTO CHUNKS######
    def chunk_dict_generator(self, protein_seqs, seq_chunks):
        res = {}
        for seq in seq_chunks:
            res[seq] = protein_seqs[seq]
        return res

    ######CHECKING HMMER JOBS######

    def get_taxon_chunks(self, taxon_id, domtblout_path, db):
        all_domtblout_with_chunks = os.listdir(domtblout_path)
        res = []
        for domtblout in all_domtblout_with_chunks:
            if db + str(taxon_id) in domtblout:
                res.append(domtblout)
        return res

    def taxon_annotation_finished(self, target_ref, output_folder, chunks_n):
        c = 0
        domtblout_folder = add_slash(add_slash(output_folder) + 'searchout')
        domtblout_files = os.listdir(domtblout_folder)
        for domtblout in domtblout_files:
            if target_ref in domtblout:
                if '_finished' in domtblout: c += 1
        if c == chunks_n:
            for domtblout in domtblout_files:
                if target_ref in domtblout:
                    if '_finished' in domtblout:
                        move_file(domtblout_folder + domtblout, domtblout_folder + domtblout.strip('_finished'))
            return True
        else:
            sleep(5)
            return False

    def save_temp_fasta_length(self, chunk_dir, hmm, count_seqs, db):
        with open(f'{chunk_dir}missing_annotation.{db}.length', 'a+') as file:
            file.write(hmm + '\t' + str(count_seqs) + '\n')

    ########To merge domtblout

    def group_output_chunks(self, searchout_path, all_searchout_with_chunks, chunk_suffix):
        # grouping the domtblout into the corresponding original 'hmm'
        res = {}
        for searchout in all_searchout_with_chunks:
            if 'chunk' in searchout:
                ref_name = searchout.split('_chunk_')[0]
                ref_searchout = ref_name + chunk_suffix
                if ref_searchout not in res: res[ref_searchout] = []
                res[ref_searchout].append(searchout_path + searchout)
        return res

    def merge_output_chunks(self, searchout_path, all_searchout_with_chunks, chunk_suffix, stdout_file=None):
        grouped_searchout = self.group_output_chunks(searchout_path, all_searchout_with_chunks, chunk_suffix)
        for searchout in grouped_searchout:
            concat_files(output_file=searchout_path + searchout, list_file_paths=grouped_searchout[searchout],
                         stdout_file=stdout_file)
            for to_delete in grouped_searchout[searchout]:
                os.remove(to_delete)



    def get_dmndout_line(self, line):
        res = line.strip('\n').split()
        if res: return res[0]

    def get_domtblout_line(self, line):
        if line[0] != '#':
            res = line.strip('\n').split()[0:22]
            if '---' not in line:
                if len(res) == 22:
                    return res[0]
            else:
                print('BAD LINE', line)


    def split_hits(self, searchout_path, domtblout, chunks_domtblout, current_chunk_dir, worker_count):
        # split the hits into different chunk allows for memory saving during hit processing
        split_searchout_path = searchout_path.replace('raw_searchout', 'searchout')
        list_of_files = [f'{split_searchout_path}{domtblout}_chunk_{i}' for i in range(worker_count)]
        Path(f'{current_chunk_dir}searchout').mkdir(parents=True, exist_ok=True)
        hit_to_file = {}
        file_yielder = yield_file(list_of_files)
        for chunk_searchout in chunks_domtblout:
            chunk_searchout_path = searchout_path + chunk_searchout
            if '.dmndout' in chunk_searchout: read_function=self.get_dmndout_line
            elif '.domtblout' in chunk_searchout: read_function=self.get_domtblout_line
            with open(chunk_searchout_path) as file:
                line = file.readline()
                while line:
                    query_name = read_function(line)
                    if query_name:
                        if query_name not in hit_to_file:
                            file_name = next(file_yielder)
                            hit_to_file[query_name] = file_name
                        out_file_path = hit_to_file[query_name]
                        with open(out_file_path, 'a+') as opened_file: opened_file.write(line)
                    line = file.readline()
            if not self.keep_files:
                os.remove(chunk_searchout_path)

    ########To calculate new evalue

    def get_temp_fasta_length(self, chunk_dir, domtblout, db):
        ref = domtblout.split('_')[0]
        with open(f'{chunk_dir}missing_annotation.{db}.length') as file:
            line = file.readline()
            while line:
                ref_key, ref_len = line.split('\t')
                if ref == ref_key: return int(ref_len)
                line = file.readline()

    def remove_temp_fasta_length(self, chunk_dir, db):
        if file_exists(f'{chunk_dir}missing_annotation.{db}.length'):
            os.remove(f'{chunk_dir}missing_annotation.{db}.length')


    def remove_annotated_queries(self, missing_queries, annotated_queries):
        for p in annotated_queries:
            if p in missing_queries:
                del missing_queries[p]




    def generate_temp_fasta(self, missing_queries, output_folder, db):
        temp_path = f'{output_folder}missing_annotations.{db}.tmp'
        if file_exists(temp_path):
            remove_temp_fasta(temp_path, db)
        with open(temp_path, 'w+') as file:
            for mq in missing_queries:
                chunks = [missing_queries[mq][x:x + 60] for x in range(0, len(missing_queries[mq]), 60)]
                chunk_str = '\n'.join(i for i in chunks)
                file.write('>' + mq + '\n' + chunk_str + '\n')
        return temp_path

    ########Processing HMMER hits

    def is_overlap(self, temp_queries, current_query):
        # the coordinates here already take into account the overlap value, so even if the y set is small or empty, it doesnt matter
        if not temp_queries or not current_query: return False
        y_start,y_end=recalculate_coordinates(current_query['hit_start'],
                                              current_query['hit_end'],
                                              self.overlap_value)
        y=set(range(y_start, y_end))

        for t in temp_queries:
            if t['hit_name'] == current_query['hit_name']:  return True
            x_start,x_end = recalculate_coordinates(t['hit_start'],
                                                    t['hit_end'],
                                                    self.overlap_value)
            x = set(range(x_start, x_end))
            res = x.intersection(y)
            if res: return True
        return False


    def get_best_hits_approximation(self, query_hits, sorting_class,sorting_type):
        '''
        this is just a lazy implementation when the amount of hits is too big and/or the hits are too small
        even with multiprocessing and cython we may run into computationally unfeasable calculations
        when this happens, we do generate a straightforward "best hit"
        Best hit will take the lowest evalue hit and add the next lowest evalue hit (without overlapping) until we reach a point where this cycle cant be repeated.
        This doesnt effectively calculate the "best hit", just a biased (since we start with lowest evalue as root) approximation
        Still, it is a pretty good approximation anyhow
        '''
        if sorting_class=='consensus': query_hits = self.sort_scaled_hits(query_hits,sorting_type=sorting_type)
        else: query_hits = self.sort_hits(query_hits, sorting_class,sorting_type=sorting_type)
        combo = []
        while query_hits:
            next_hit = query_hits.pop(0)
            if sorting_class == 'consensus':
                if not self.is_overlap_Consensus(combo, next_hit):
                    combo.append(next_hit)
            elif sorting_class == 'processor':
                if not self.is_overlap(combo, next_hit):
                    combo.append(next_hit)
        return combo

    def get_lowest_hit(self, query_hits, sorting_class,sorting_type):
        '''
        this will take the hit with the lowest evalue, regardless if there are multiple domains
        '''
        if not query_hits: return []
        if sorting_class=='consensus': query_hits = self.sort_scaled_hits(query_hits,sorting_type=sorting_type)
        else: query_hits = self.sort_hits(query_hits, sorting_class,sorting_type=sorting_type)
        lowest_hit = query_hits.pop(0)
        combo = [lowest_hit]
        return combo

    def cython_to_query_hits(self, cython_hits, conversion_dict):
        res = []
        for hit in cython_hits:
            res.append(conversion_dict[hit[0]])
        return res

    def get_min_max_dfs(self, cython_possible_combos, conversion_dict, sorting_class):
        min_val = None
        max_val = None
        for len_combo in cython_possible_combos:
            for cython_combo in cython_possible_combos[len_combo]:
                combo = self.cython_to_query_hits(cython_combo, conversion_dict)
                for hit in combo:
                    if sorting_class == 'processor':
                        hit_info = hit
                    elif sorting_class == 'consensus':
                        hmm_file, hmm_hit, hit_info = hit
                    # we dont consider 0 for the scaling, 0 will always be scaled to max/1
                    if hit_info[self.sorting_type]:
                        current_val = log10(hit_info[self.sorting_type])
                        if min_val is None: min_val = current_val
                        if max_val is None: max_val = current_val
                        if current_val > max_val: max_val = current_val
                        if current_val < min_val: min_val = current_val
        # lower is best
        if self.sorting_type == 'evalue':
            return max_val, min_val
        # higher is best
        elif self.sorting_type == 'bitscore':
            return min_val, max_val

    def get_combo_score(self, combo, query_length, min_val, max_val):

        average_value = 0
        average_hit_coverage = 0
        hit_ranges = []
        for hit in combo:
            hit_start, hit_end, ref_start, ref_end = hit['hit_start'], hit['hit_end'], \
                                                     hit['ref_start'], hit['ref_end']
            ref_len = hit['ref_len']
            hit_evalue = hit['evalue']
            hit_bitscore = hit['bitscore']
            hit_coverage = (hit_end - hit_start + 1) / query_length
            ref_coverage = (ref_end - ref_start + 1) / ref_len
            average_hit_coverage += hit_coverage
            # lower is better
            if self.sorting_type == 'evalue':
                if hit_evalue:
                    scaled_value = min_max_scale(log10(hit_evalue), min_val, max_val) + 0.01
                else:
                    scaled_value = 2
            # higher is better
            elif self.sorting_type == 'bitscore':
                if hit_bitscore:
                    scaled_value = min_max_scale(log10(hit_bitscore), min_val, max_val) + 0.01
                else:
                    scaled_value = 0
            average_value += scaled_value
            hit_ranges.append([hit_start, hit_end])

        # sequence hit quality
        average_value = average_value / len(combo)
        # average coverage of all hits
        average_hit_coverage /= len(combo)
        # coverage of combination
        combination_coverage = get_combination_ranges(hit_ranges)
        combination_coverage /= query_length
        '''
        average_hit_coverage 0-1
        combination_coverage 0-1
        average_value 0-1 and exceptionally above when e-value is 0
        
        the user is not meant to choose a formula, this was written for internal quality testing
        for 1e-3 formula 1 with bitscore and dfs algorithm produced the best results
        '''
        if self.best_combo_formula==1: #default
            combo_score = average_value
        elif self.best_combo_formula==2:
            combo_score = average_value * combination_coverage * average_hit_coverage
        elif self.best_combo_formula==3:
            combo_score = (average_value * combination_coverage * average_hit_coverage)**(1/3)
        elif self.best_combo_formula==4:
            combo_score = average_value * combination_coverage
        elif self.best_combo_formula==5:
            combo_score = (average_value * combination_coverage)**(1/2)
        elif self.best_combo_formula==6:
            combo_score = (average_value + average_hit_coverage + combination_coverage)/3
        elif self.best_combo_formula==7:
            combo_score = (2*average_value + average_hit_coverage + combination_coverage)/4
        elif self.best_combo_formula==8:
            combo_score = (2*average_value + combination_coverage)/3
        elif self.best_combo_formula==9:
            combo_score = (3*average_value + average_hit_coverage + combination_coverage)/5
        elif self.best_combo_formula==10:
            combo_score = (3*average_value + combination_coverage)/4
        elif self.best_combo_formula == 11:
            combo_score = average_value * (combination_coverage) ** 2 * (average_hit_coverage) ** 2
        elif self.best_combo_formula == 12:
            combo_score = average_value * (combination_coverage) ** 2

        return combo_score

    def get_best_combo_bit_score(self, good_combos):
        best_bitscore = None
        best_combo = None
        for combo in good_combos:
            combo_bitscore = 0
            for hit in combo:
                combo_bitscore += hit['bitscore']
            combo_bitscore /= len(combo)
            if best_bitscore is None:
                best_bitscore = combo_bitscore
                best_combo = combo
            else:
                if combo_bitscore > best_bitscore:
                    best_combo = combo
        return best_combo

    def get_best_hits(self, query_hits, query_length):
        '''
        we take into consideration the coverage of our hits and their evalue
        steps:
        1- get all possible combinations of hits
        2- check which possible combinations don't overlap
        3- get best combination (best evalue and coverage)
        '''
        ordered_query_hits = self.sort_hits(query_hits, sorting_class='processor',sorting_type=self.sorting_type)
        cython_hits, conversion_dict = self.query_hits_to_cython(ordered_query_hits)
        cython_possible_combos = get_non_overlapping_hits(cython_hits, time_limit=self.time_limit)
        # sometimes this calculation is not feasible (when we have too many small hits and the query sequence is too long, the calculation of conbinations would take forever - limit of 5mins)
        if not cython_possible_combos: return None
        best_hit_score = 0
        best_combo = None
        min_val, max_val = self.get_min_max_dfs(cython_possible_combos, conversion_dict, sorting_class='processor')
        good_combos = []
        for len_combo in cython_possible_combos:
            if len_combo:
                for cython_combo in cython_possible_combos[len_combo]:
                    combo = self.cython_to_query_hits(cython_combo, conversion_dict)
                    combo_score = self.get_combo_score(combo, query_length, min_val, max_val)
                    if combo_score > 1.5: good_combos.append(combo)
                    if not best_combo:
                        best_combo = combo
                        best_hit_score = combo_score
                    elif combo_score > best_hit_score:
                        best_hit_score = combo_score
                        best_combo = combo
        # some sequences will have extremely high scores against multiple profiles (e-value of 0), when this is the case, and the combo score is also high, we use bitscore to select the best
        # this can only be done here because we are comparing the same db
        if self.sorting_type == 'evalue':
            if good_combos:
                best_combo = self.get_best_combo_bit_score(good_combos)
        return best_combo

    def query_hits_to_cython(self, query_hits):
        conversion_dict = {}
        res = set()
        for hit_i in range(len(query_hits)):
            hit_start,hit_end=recalculate_coordinates(query_hits[hit_i]['hit_start'],
                                                      query_hits[hit_i]['hit_end'],
                                                      self.overlap_value)
            res.add(tuple([
                hit_i,
                hit_start,
                hit_end,
                query_hits[hit_i]['hit_name']
            ]))
            conversion_dict[hit_i] = query_hits[hit_i]
        return res, conversion_dict

    def sort_hits(self, query_hits, sorting_class,sorting_type):
        #we do a 2 step sorting to break ties (very frequent in evalue, less so in bitscore)
        if not query_hits: return query_hits
        #first we sort by the preferred sorting type
        if sorting_type=='bitscore': to_reverse=True
        else: to_reverse=False
        if sorting_class == 'processor':
            sorted_hits= sorted(query_hits, key=lambda k: k[sorting_type],reverse=to_reverse)
        elif sorting_class == 'consensus':
            sorted_hits= sorted(query_hits, key=lambda k: k[2][sorting_type],reverse=to_reverse)
        res=[]
        #then we separate by sorting value
        sorted_hits_groups=[]
        c=0
        for i in sorted_hits:
            if sorting_class == 'processor':    hit_value = i[sorting_type]
            elif sorting_class == 'consensus':  hit_value = i[2][sorting_type]
            if not sorted_hits_groups:
                sorted_hits_groups.append([])
                current=hit_value
            if hit_value==current: sorted_hits_groups[c].append(i)
            else:
                sorted_hits_groups.append([i])
                c+=1
                current=hit_value
        #then we sort each subgroup by the secondary sorting type
        sec_sorting_type= 'bitscore' if sorting_type=='evalue' else 'evalue'
        if sec_sorting_type=='bitscore': to_reverse=True
        else: to_reverse=False
        for sg in sorted_hits_groups:
            if sorting_class == 'processor':
                temp= sorted(sg, key=lambda k: k[sec_sorting_type],reverse=to_reverse)
            elif sorting_class == 'consensus':
                temp= sorted(sg, key=lambda k: k[2][sec_sorting_type],reverse=to_reverse)
            res.extend(temp)
        return res

    def calculate_evalue_threshold(self, query_len,diamond=False):
        #See :     https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3416-y
        #if the user sets evalue threshold, we use that
        if self.evalue_threshold and self.evalue_threshold!='dynamic': return self.evalue_threshold
        #if the user sets the evalue threshold to dynamic OR if the user doesnt set evalue threshold and we are analysing diamond's output
        elif self.evalue_threshold == 'dynamic' or diamond:
            #this somewhat translates to the findings in the paper
            res = float('1e-' + str(ceil(query_len / 10)))
            #for very small sequences, we set a hard threshold
            if res > self.default_evalue_threshold: res = self.default_evalue_threshold
            return res
        #if the user doesnt set the evalue threshold and we are analysing hmmer's output
        else:
            return self.default_evalue_threshold

    def recalculate_evalue(self, original_evalue, to_divide, to_multiply):
        '''
        keep in mind these calculations will result in minor rounding errors!
        For example in 3 different executions for the same query:
        When splitting the fasta into chunks:
         -one time the fasta had 20 sequences and an i-evalue of 5.98e-52
         - another time the fasta had 24 sequences and an i-evalue of 5.7e-52
         - when the fasta was not split we got 5.8 e-52
        Ideally we would not split the fasta into chunks; this would avoid rounding errors but would not allow parallelization
        '''
        if self.force_evalue:
            return original_evalue
        else:
            evalue = float(original_evalue)
            evalue /= to_divide
            evalue *= to_multiply
            if evalue>=self.minimum_evalue_threshold:
                evalue=self.minimum_evalue_threshold
            return evalue

    def read_searchout(self, output_path, count_seqs_chunk, count_seqs_original_file,count_residues_original_file, get_output=True):
        if '.dmndout' in output_path:
            return self.read_dmndout(output_path=output_path, count_seqs_chunk=count_seqs_chunk,
                                count_residues_original_file=count_residues_original_file,get_output=get_output)
        elif '.domtblout' in output_path:
            return self.read_domtblout(output_path=output_path, count_seqs_chunk=count_seqs_chunk,
                                count_seqs_original_file=count_seqs_original_file,get_output=get_output)

    def read_dmndout(self, output_path, count_seqs_chunk, count_residues_original_file, get_output=True):
        res = {}
        hit_counter = set()
        # reading the output file
        with open(output_path, 'r', encoding='utf8') as file:
            line = file.readline().strip('\n')
            while line:
                #try:
                if True:
                    line = line.split()
                    query_name = line[0]
                    query_len = int(line[1])
                    ref_name = line[2]
                    ref_len = int(line[3])
                    query_start = int(line[4])
                    query_end = int(line[5])
                    ref_start = int(line[6])
                    ref_end = int(line[7])
                    evalue = float(line[8])
                    bitscore = float(line[9])
                    #diamond's evalue is not scaled, so we need to multiply by residues count
                    evalue = self.recalculate_evalue(evalue, self.diamond_db_size, count_residues_original_file)
                    #hmmer's coordinates on the seq are always the same, regardless of direction
                    direction='+'
                    if query_start<query_end:
                        corrected_query_start=query_start
                        corrected_query_end=query_end
                    else:
                        direction='-'
                        corrected_query_start=query_end
                        corrected_query_end=query_start



                    if get_output:
                        if query_name not in res: res[query_name] = {'query_len': query_len, 'hits': []}
                    # when processing mg data, the evalue of hits for chunks should be higher than what it really is (because chunk has less seqs = more significant e-value)
                    if evalue <= self.calculate_evalue_threshold(int(query_len),diamond=True):
                        hit_dict = {'hit_name': ref_name,
                                    'hit_accession': '',
                                    'ref_len': ref_len,
                                    'evalue': evalue,
                                    'bitscore': bitscore,
                                    'direction': direction,
                                    'ref_start': ref_start if ref_start < ref_end else ref_end,
                                    'ref_end': ref_end if ref_end > ref_start else ref_start,
                                    'hit_start': corrected_query_start,
                                    'hit_end': corrected_query_end,
                                    }
                        if get_output:
                            res[query_name]['hits'].append(hit_dict)
                        hit_counter.add(query_name)
                #except:
                #    print(f'Could not read line: {line} in file', output_path, flush=True)
                line = file.readline().strip('\n')
        return res, hit_counter

    def read_domtblout(self, output_path, count_seqs_chunk, count_seqs_original_file, get_output=True):
        res = {}
        hit_counter = set()
        # reading the output file
        with open(output_path, 'r', encoding='utf8') as file:
            line = file.readline().strip('\n')
            while line:
                if line[0] != '#':
                    #try:
                    if True:
                        # after 22 it's just the free text of the annotation description
                        line = line.split()[0:22]
                        if len(line) == 22:
                            query_name = line[0]
                            query_len = int(line[2])
                            hmm_name = line[3]
                            hmm_accession = line[4]
                            hmm_len = int(line[5])
                            # for these two we accept multiple hits
                            if self.domain_algorithm in ['dfs', 'heuristic']:
                                evalue = float(line[12])
                                bitscore = float(line[13])
                            # here we want the full sequence best
                            else:
                                evalue = float(line[6])
                                bitscore = float(line[7])
                            evalue = self.recalculate_evalue(evalue, count_seqs_chunk, count_seqs_original_file)
                            hmm_coord_from = int(line[15])
                            hmm_coord_to = int(line[16])
                            #ali_coord_from = int(line[17])
                            #ali_coord_to = int(line[18])
                            # we will use envelop coords as per HMMER's manual recommendation
                            env_coord_from = int(line[19])
                            env_coord_to = int(line[20])
                            #hmmer's coordinates on the seq are always the same, regardless of direction
                            direction='+'
                            if env_coord_from<env_coord_to:
                                corrected_env_coord_from=env_coord_from
                                corrected_env_coord_to=env_coord_to
                            else:
                                direction='-'
                                corrected_env_coord_from=env_coord_to
                                corrected_env_coord_to=env_coord_from
                            if get_output:
                                if query_name not in res: res[query_name] = {'query_len': query_len, 'hits': []}
                            # when processing mg data, the evalue of hits for chunks should be higher than what it really is (because chunk has less seqs = more significant e-value)
                            if evalue <= self.calculate_evalue_threshold(int(query_len)):
                                hit_dict = {'hit_name': hmm_name,
                                            'hit_accession': hmm_accession,
                                            'ref_len': hmm_len,
                                            'evalue': evalue,
                                            'bitscore': bitscore,
                                            'direction': direction,
                                            'ref_start': hmm_coord_from if hmm_coord_from < hmm_coord_to else hmm_coord_to,
                                            'ref_end': hmm_coord_to if hmm_coord_to > hmm_coord_from else hmm_coord_from,
                                            'hit_start': corrected_env_coord_from,
                                            'hit_end': corrected_env_coord_to,
                                            }
                                if get_output:
                                    res[query_name]['hits'].append(hit_dict)
                                hit_counter.add(query_name)
                    #except:
                    #    print(f'Could not read line: {line} in file', output_path, flush=True)
                line = file.readline().strip('\n')
        return res, hit_counter

    def process_searchout(self, output_path, count_seqs_chunk, count_seqs_original_file,count_residues_original_file, stdout_path=None):
        '''
        this will read the domtblout file and export:
        query -> all hits with their respective evalue and hit start and end
        after we have this , we invoke get_best_hits to get our result

        we use i_evalue as per the hmmer manual

        this function can read raw and concatenated domtblout files
        '''
        if isinstance(stdout_path, str):
            stdout_file = open(stdout_path, 'a+')
        else:
            stdout_file = stdout_path
        print('------------------------------------------', flush=True, file=stdout_file)
        print('Processing hits for:\n' + output_path, flush=True, file=stdout_file)
        queries_searchout, hit_counter = self.read_searchout(output_path=output_path, count_seqs_chunk=count_seqs_chunk,
                                                             count_seqs_original_file=count_seqs_original_file,
                                                             count_residues_original_file=count_residues_original_file)
        hit_counter = len(hit_counter)
        # processing the hits and getting the best hit/non-overlapping hits
        print(f'Found {hit_counter} hits in:\n{output_path}\nWill now get best hits!', flush=True,file=stdout_file)
        approximated_hits = []
        hmm = get_path_level(output_path, remove_extension=True)
        res_annotation = {}
        for query in queries_searchout:
            if query not in res_annotation: res_annotation[query] = {}
            list_hits = queries_searchout[query].pop('hits')
            if self.domain_algorithm == 'heuristic':
                best_hit = self.get_best_hits_approximation(list(list_hits), sorting_class='processor',sorting_type=self.sorting_type)
            elif self.domain_algorithm == 'bpo':
                best_hit = self.get_lowest_hit(list(list_hits), sorting_class='processor',sorting_type=self.sorting_type)
            else:
                try:
                    best_hit = self.get_best_hits(list(list_hits), queries_searchout[query]['query_len'])
                except (TimeoutError, RecursionError) as e:
                    #heuristic with bitscore produces bad results, so we force evalue sorting
                    best_hit = self.get_best_hits_approximation(list(list_hits), sorting_class='processor',sorting_type='evalue')
                    approximated_hits.append(query)
            queries_searchout[query]['best_hit'] = best_hit
            res_annotation[query][hmm] = queries_searchout[query]
        if approximated_hits:
            approximated_hits = ' ; '.join(approximated_hits)
            print(f'Some hits in the output file:\n{output_path} were approximated. They were:\n{approximated_hits}',flush=True, file=stdout_file)
        print(f'Finished processing hits for:\n{output_path}', flush=True, file=stdout_file)
        print('------------------------------------------', flush=True, file=stdout_file)
        if isinstance(stdout_path, str): stdout_file.close()
        return res_annotation

    def save_processed_hits(self, annotation_output, output_dir, domtblout):
        annotation_output_file = output_dir + domtblout.replace('domtblout', 'pro')
        try:
            os.remove(annotation_output_file)
        except:
            pass
        with open(annotation_output_file, 'w') as output_file:
            first_line = ['Query',
                          'Ref_file',
                          'Ref_hit',
                          'Ref_hit_accession',
                          'evalue',
                          'bitscore',
                          'Direction',
                          'Query_length',
                          'Query_hit_start',
                          'Query_hit_end',
                          'Ref_hit_start',
                          'Ref_hit_end',
                          'Ref_length',
                          ]
            first_line = '\t'.join(first_line) + '\n'
            output_file.write(first_line)
            for query in annotation_output:
                for hmm_file in annotation_output[query]:
                    # some hmm hits won't pass the i-evalue threshold, so we dont consider them
                    if '_chunk_' in hmm_file:
                        str_hmm_file = hmm_file.split('_chunk_')[0]
                    else:
                        str_hmm_file = hmm_file
                    if annotation_output[query][hmm_file]['best_hit']:
                        for hit in annotation_output[query][hmm_file]['best_hit']:
                            line = [str(query),
                                    str(str_hmm_file),
                                    str(hit['hit_name']),
                                    str(hit['hit_accession']),
                                    str(hit['evalue']),
                                    str(hit['bitscore']),
                                    str(hit['direction']),
                                    str(annotation_output[query][hmm_file]['query_len']),
                                    str(hit['hit_start']),
                                    str(hit['hit_end']),
                                    str(hit['ref_start']),
                                    str(hit['ref_end']),
                                    str(hit['ref_len']),
                                    ]
                            line = '\t'.join(line)
                            output_file.write(line + '\n')

    def merge_gff_output(self,output_folder,output_file,chunks_path):
        all_seq_regions=[]
        gff_version=None
        for chunk_output in chunks_path:
            target_chunk_output = chunk_output + output_file
            target_chunk_output=target_chunk_output.replace('.tsv','.gff')
            with open(target_chunk_output, 'r') as chunk_file:
                chunk_line = chunk_file.readline()
                while chunk_line:
                    if '##gff-version' in chunk_line:
                        gff_version=chunk_line
                    if '##sequence-region' in chunk_line:
                        all_seq_regions.append(chunk_line)
                    chunk_line = chunk_file.readline()
        if not gff_version:
            kill_switch(InvalidGFFVersion,flush=True, file=self.redirect_verbose)
        output_gff_file=f'{output_folder}{output_file}'.replace('.tsv','.gff')
        with open(output_gff_file,'w+') as file:
            file.write(gff_version)
            for seq_reg in all_seq_regions:
                file.write(seq_reg)
            for chunk_output in chunks_path:
                target_chunk_output = chunk_output + output_file
                target_chunk_output = target_chunk_output.replace('.tsv', '.gff')
                with open(target_chunk_output, 'r') as chunk_file:
                    chunk_line = chunk_file.readline()
                    while chunk_line:
                        if not chunk_line.startswith('##'):
                            file.write(chunk_line)
                        chunk_line = chunk_file.readline()




    def merge_target_output(self, output_file, output_folder, chunks_path, stdout_file, same_output=True):
        print(f'Merging chunks to {output_folder}{output_file}',flush=True, file=stdout_file)
        if self.output_gff and ('consensus' in output_file or 'integrated' in output_file):
            self.merge_gff_output(output_folder,output_file,chunks_path)
        header = False
        with open(output_folder + output_file, 'w+') as file:
            for chunk_output in chunks_path:
                if same_output:
                    target_chunk_output = chunk_output + output_file
                else:
                    target_chunk_output = chunk_output
                with open(target_chunk_output, 'r') as chunk_file:
                    chunk_line = chunk_file.readline()
                    while chunk_line:
                        if 'Query' in chunk_line:
                            if not header:
                                file.write(chunk_line)
                                header = True
                        else:
                            file.write(chunk_line)
                        chunk_line = chunk_file.readline()




if __name__ == '__main__':
    hmm_pro = Homology_processor()
