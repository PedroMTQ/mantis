try:
    from source.MANTIS_Assembler import *
except:
    from MANTIS_Assembler import *

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
                    raise CythonNotCompiled
        # when cython's version is not compatible with the current python version, we need to recompile it
        else:
            compile_cython()
            try:
                from source.cython_src.get_non_overlapping_hits import get_non_overlapping_hits
            except:
                try:
                    from cython_src.get_non_overlapping_hits import get_non_overlapping_hits
                except:
                    raise CythonNotCompiled


# This class will process the domtblout output and get the best hit for our queries, it is inherited by the MANTIS

class MANTIS_Processor():

    def create_chunk_hmmer_dirs(self, chunk_dir):
        for hmmer_dir in ['output_hmmer', 'domtblout']:  # ,'tblout']:
            if not os.path.exists(chunk_dir + hmmer_dir):
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

    def taxon_annotation_finished(self, target_hmm, output_folder, chunks_n):
        c = 0
        domtblout_folder = add_slash(add_slash(output_folder) + 'domtblout')
        domtblout_files = os.listdir(domtblout_folder)
        for domtblout in domtblout_files:
            if target_hmm in domtblout:
                if '_finished' in domtblout: c += 1
        if c == chunks_n:
            for domtblout in domtblout_files:
                if target_hmm in domtblout:
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

    def group_output_chunks(self, domtblout_path, all_domtblout_with_chunks, chunk_suffix):
        # grouping the domtblout into the corresponding original 'hmm'
        res = {}
        for domtblout in all_domtblout_with_chunks:
            if 'chunk' in domtblout:
                hmm_name = domtblout.split('_chunk_')[0]
                hmm_domtblout = hmm_name + chunk_suffix
                if hmm_domtblout not in res: res[hmm_domtblout] = []
                res[hmm_domtblout].append(domtblout_path + domtblout)
        return res

    def merge_output_chunks(self, domtblout_path, all_domtblout_with_chunks, chunk_suffix, stdout_file=None):
        grouped_domtblout = self.group_output_chunks(domtblout_path, all_domtblout_with_chunks, chunk_suffix)
        for domtblout in grouped_domtblout:
            concat_files(output_file=domtblout_path + domtblout, list_file_paths=grouped_domtblout[domtblout],
                         stdout_file=stdout_file)
            for to_delete in grouped_domtblout[domtblout]:
                os.remove(to_delete)

    def get_domtblout_line(self, line):
        if line[0] != '#':
            res = line.strip('\n').split()[0:22]
            if '---' not in line:
                if len(res) == 22:
                    return res[0]
            else:
                print('BAD LINE', line)


    def split_hits(self, domtblout_path, domtblout, chunks_domtblout, current_chunk_dir, worker_count):
        # split the hits into different chunk allows for memory saving during hit processing
        split_domtblout_path = domtblout_path.replace('split_hits', 'domtblout')
        list_of_files = [f'{split_domtblout_path}{domtblout}_chunk_{i}' for i in range(worker_count)]
        Path(f'{current_chunk_dir}domtblout').mkdir(parents=True, exist_ok=True)
        hit_to_file = {}
        file_yielder = yield_file(list_of_files)
        for chunk_domtblout in chunks_domtblout:
            chunk_domtblout_path = domtblout_path + chunk_domtblout
            with open(chunk_domtblout_path) as file:
                line = file.readline()
                while line:
                    query_name = self.get_domtblout_line(line)
                    if query_name:
                        if query_name not in hit_to_file:
                            file_name = next(file_yielder)
                            hit_to_file[query_name] = file_name
                        out_file_path = hit_to_file[query_name]
                        with open(out_file_path, 'a+') as opened_file: opened_file.write(line)
                    line = file.readline()
            if not self.keep_files:
                os.remove(chunk_domtblout_path)

    ########To calculate new evalue

    def get_temp_fasta_length(self, chunk_dir, domtblout, db):
        hmm = domtblout.split('_')[0]
        with open(f'{chunk_dir}missing_annotation.{db}.length') as file:
            line = file.readline()
            while line:
                hmm_key, hmm_len = line.split('\t')
                if hmm == hmm_key: return int(hmm_len)
                line = file.readline()

    def remove_temp_fasta_length(self, chunk_dir, db):
        if file_exists(f'{chunk_dir}missing_annotation.{db}.length'):
            os.remove(f'{chunk_dir}missing_annotation.{db}.length')

    ########Processing protein fasta
    def remove_temp_fasta(self, temp_fasta_path, db):
        if f'missing_annotations.{db}.tmp' in temp_fasta_path:
            os.remove(temp_fasta_path)

    def remove_annotated_queries(self, missing_queries, annotated_queries):
        for p in annotated_queries:
            if p in missing_queries:
                del missing_queries[p]

    def process_protein_fasta_line(self, res, query, fasta_line, start_recording):
        # for the first > line
        if '>' in fasta_line and not start_recording:
            fasta_line = fasta_line.replace('\'', '')
            fasta_line = fasta_line.replace('>', '')
            fasta_line = fasta_line.replace('\"', '')
            new_query = fasta_line.split()[0]
            start_recording = True
            res[new_query] = ''
            return res, new_query, start_recording
        # for posterior > lines
        elif '>' in fasta_line and start_recording:
            fasta_line = fasta_line.replace('\'', '')
            fasta_line = fasta_line.replace('>', '')
            fasta_line = fasta_line.replace('\"', '')
            start_recording = True
            new_query = fasta_line.split()[0]
            res[new_query] = ''
            return res, new_query, start_recording
        # to get the sequences
        elif start_recording:
            fasta_line = fasta_line.replace('\"', '')
            fasta_line = fasta_line.replace('\'', '')
            res[query] += fasta_line.strip()
            return res, query, start_recording

    def read_protein_fasta(self, protein_fasta_path):
        res = {}
        with open(protein_fasta_path, 'r') as file:
            line = file.readline()
            if line[0] != '>':  raise FileNotFoundError
            start_recording = False
            query = None
            while line:
                line = line.strip('\n')
                if line:
                    res, query, start_recording = self.process_protein_fasta_line(res=res,
                                                                                  query=query,
                                                                                  fasta_line=line,
                                                                                  start_recording=start_recording)
                line = file.readline()
        return res

    def generate_temp_fasta(self, missing_queries, output_folder, db):
        temp_path = f'{output_folder}missing_annotations.{db}.tmp'
        if os.path.exists(temp_path):
            self.remove_temp_fasta(temp_path, db)
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
        y_start,y_end=recalculate_coordinates(current_query['env_coord_from'],
                                              current_query['env_coord_to'],
                                              self.overlap_value)
        y=set(range(y_start, y_end))

        for t in temp_queries:
            if t['hmm_name'] == current_query['hmm_name']:  return True
            x_start,x_end = recalculate_coordinates(t['env_coord_from'],
                                                    t['env_coord_to'],
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
        query_hits = self.sort_hits(query_hits, sorting_class,sorting_type=sorting_type)
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
        query_hits = self.sort_hits(query_hits, sorting_class,sorting_type=sorting_type)
        lowest_hit = query_hits.pop(0)
        combo = [lowest_hit]
        return combo

    def cython_to_query_hits(self, cython_hits, conversion_dict):
        res = []
        for hit in cython_hits:
            res.append(conversion_dict[hit[0]])
        return res

    def get_min_max(self, cython_possible_combos, conversion_dict, sorting_class):
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
            hit_start, hit_end, hmm_start, hmm_end = hit['env_coord_from'], hit['env_coord_to'], \
                                                     hit['hmm_coord_from'], hit['hmm_coord_to']
            hmm_len = hit['hmm_len']
            hit_evalue = hit['evalue']
            hit_bitscore = hit['bitscore']
            hit_coverage = (hit_end - hit_start + 1) / query_length
            hmm_coverage = (hmm_end - hmm_start + 1) / hmm_len
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
        if self.best_combo_formula==1: #best
            combo_score = average_value
        elif self.best_combo_formula==2:    #default
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
        min_val, max_val = self.get_min_max(cython_possible_combos, conversion_dict, sorting_class='processor')
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
            hit_start,hit_end=recalculate_coordinates(query_hits[hit_i]['env_coord_from'],
                                                      query_hits[hit_i]['env_coord_to'],
                                                      self.overlap_value)
            res.add(tuple([
                hit_i,
                hit_start,
                hit_end,
                query_hits[hit_i]['hmm_name']
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

    def calculate_evalue_threshold(self, query_len):
        '''
        See :
        https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-3416-y

        #These values were for diamond, but we assume they would also work with HMMER
        "the best choice of e-value threshold was dependent on read length.
         For read lengths 100 bp and higher, the best choice of e-value was not the default (Fig. 2c-d).
        For 100-150 bp, a threshold of 1-e10 yielded the highest accuracy (Additional file 1, Figure S3).
        For 200-250 bp read lengths, 1e-25 was most accurate (Additional file 1, Figure S3).
        In general, the optimal e-value threshold increases with increasing read length."
        '''
        if self.evalue_threshold == 'dynamic':
            if query_len <= 150:
                return 1e-10
            elif query_len >= 250:
                return 1e-25
            else:
                return float('1e-' + str(ceil(query_len / 10)))
        else:
            return self.evalue_threshold

    def recalculate_evalue(self, i_evalue, count_seqs_chunk, count_seqs_original_file):
        '''
        keep in mind these calculations will result in minor rounding errors!
        For example in 3 different executions for the same query:
        When splitting the fasta into chunks:
         -one time the fasta had 20 sequences and an i-evalue of 5.98e-52
         - another time the fasta had 24 sequences and an i-evalue of 5.7e-52
         - when the fasta was not split we got 5.8 e-52
        Ideally we would not split the fasta into chunks; this would avoid rounding errors but would not allow parallelization
        '''
        evalue = float(i_evalue)
        evalue /= count_seqs_chunk
        evalue *= count_seqs_original_file
        return evalue

    def read_domtblout(self, output_path, count_seqs_chunk, count_seqs_original_file, get_output=True):
        res = {}
        hit_counter = set()
        # reading the output file
        with open(output_path, 'r', encoding='utf8') as file:
            line = file.readline().strip('\n')
            while line:
                if line[0] != '#':
                    try:
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
                            ali_coord_from = int(line[17])
                            ali_coord_to = int(line[18])
                            # we will use envelop coords as per HMMER's manual recommendation
                            env_coord_from = int(line[19])
                            env_coord_to = int(line[20])
                            #hmmer's coordinates on the seq are always the same, regardless of direction
                            direction='Forward'
                            if env_coord_from<env_coord_to:
                                corrected_env_coord_from=env_coord_from
                                corrected_env_coord_to=env_coord_to
                            else:
                                direction='Backward'
                                corrected_env_coord_from=env_coord_to
                                corrected_env_coord_to=env_coord_from

                            #hit_range = corrected_env_coord_to - corrected_env_coord_from
                            #hit_overlap = ceil(self.overlap_value * hit_range / 2)
                            #final_env_coord_from = corrected_env_coord_from# + hit_overlap
                            #final_env_coord_to = corrected_env_coord_to# - hit_overlap


                            if get_output:
                                if query_name not in res: res[query_name] = {'query_len': query_len, 'hits': []}
                            # when processing mg data, the evalue of hits for chunks should be higher than what it really is (because chunk has less seqs = more significant e-value)
                            if evalue <= self.calculate_evalue_threshold(int(query_len)):
                                hit_dict = {'hmm_name': hmm_name,
                                            'hmm_accession': hmm_accession,
                                            'hmm_len': hmm_len,
                                            'evalue': evalue,
                                            'bitscore': bitscore,
                                            'direction': direction,
                                            'hmm_coord_from': hmm_coord_from if hmm_coord_from < hmm_coord_to else hmm_coord_to,
                                            'hmm_coord_to': hmm_coord_to if hmm_coord_to > hmm_coord_from else hmm_coord_from,
                                            'ali_coord_from': ali_coord_from if ali_coord_from < ali_coord_to else ali_coord_to,
                                            'ali_coord_to': ali_coord_to if ali_coord_to > ali_coord_from else ali_coord_from,
                                            'env_coord_from': corrected_env_coord_from,
                                            'env_coord_to': corrected_env_coord_to,
                                            }
                                if get_output:
                                    res[query_name]['hits'].append(hit_dict)
                                hit_counter.add(query_name)
                    except:
                        print(f'Could not read line: {line} in file', output_path, flush=True)
                line = file.readline().strip('\n')
        return res, hit_counter

    def process_domtblout(self, output_path, count_seqs_chunk, count_seqs_original_file, stdout_path=None):
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
        queries_domtblout, hit_counter = self.read_domtblout(output_path, count_seqs_chunk, count_seqs_original_file)
        hit_counter = len(hit_counter)
        # processing the hits and getting the best hit/non-overlapping hits
        print(f'Found {hit_counter} hits in:\n{output_path}\nWill now get best hits!', flush=True,file=stdout_file)
        approximated_hits = []
        hmm = get_path_level(output_path, remove_extension=True)
        res_annotation = {}
        for query in queries_domtblout:
            if query not in res_annotation: res_annotation[query] = {}
            list_hits = queries_domtblout[query].pop('hits')
            if self.domain_algorithm == 'heuristic':
                best_hit = self.get_best_hits_approximation(list(list_hits), sorting_class='processor',sorting_type=self.sorting_type)
            elif self.domain_algorithm == 'bpo':
                best_hit = self.get_lowest_hit(list(list_hits), sorting_class='processor',sorting_type=self.sorting_type)
            else:
                try:
                    best_hit = self.get_best_hits(list(list_hits), queries_domtblout[query]['query_len'])
                except (TimeoutError, RecursionError) as e:
                    #heuristic with bitscore produces bad results, so we force evalue sorting
                    best_hit = self.get_best_hits_approximation(list(list_hits), sorting_class='processor',sorting_type='evalue')
                    approximated_hits.append(query)
            queries_domtblout[query]['best_hit'] = best_hit
            res_annotation[query][hmm] = queries_domtblout[query]
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
                          'HMM_file',
                          'HMM_hit',
                          'HMM_hit_accession',
                          'evalue',
                          'bitscore',
                          'Direction',
                          'Query_length',
                          'Query_hit_start',
                          'Query_hit_end',
                          'HMM_hit_start',
                          'HMM_hit_end',
                          'HMM_length',
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
                                    str(hit['hmm_name']),
                                    str(hit['hmm_accession']),
                                    # all calculations with evalues are done, so just making it more readable
                                    # str(float(format(hit['i_evalue'],'.6g'))),
                                    str(hit['evalue']),
                                    str(hit['bitscore']),
                                    str(hit['direction']),
                                    str(annotation_output[query][hmm_file]['query_len']),
                                    str(hit['env_coord_from']),
                                    str(hit['env_coord_to']),
                                    str(hit['hmm_coord_from']),
                                    str(hit['hmm_coord_to']),
                                    str(hit['hmm_len']),
                                    ]
                            line = '\t'.join(line)
                            output_file.write(line + '\n')

    def merge_target_output(self, output_file, output_folder, chunks_path, stdout_file, same_output=True):
        header = False
        print(f'Merging chunks to {output_folder}{output_file}',flush=True, file=stdout_file)
        #print('Merging output', output_folder + output_file, 'from chunks:', [get_path_level(i) for i in chunks_path],flush=True, file=stdout_file)
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
    hmm_pro = MANTIS_Processor()
    query_hits = [{'hmm_name': 'hit1', 'hmm_accession': 'NF038151.1', 'hmm_len': 833, 'evalue': 6.1e-13,
                   'bitscore': 126.8, 'direction': 'Forward', 'hmm_coord_from': 234, 'hmm_coord_to': 458,
                   'ali_coord_from': 597, 'ali_coord_to': 820, 'env_coord_from': 580, 'env_coord_to': 870},
                  {'hmm_name': 'hit2', 'hmm_accession': 'NF033483.0', 'hmm_len': 563, 'evalue': 6.9e-41,
                   'bitscore': 126.8, 'direction': 'Forward', 'hmm_coord_from': 13, 'hmm_coord_to': 206,
                   'ali_coord_from': 589, 'ali_coord_to': 788, 'env_coord_from': 580, 'env_coord_to': 795},
                  {'hmm_name': 'hit3', 'hmm_accession': 'TIGR00606.1', 'hmm_len': 1311, 'evalue': 3e-05,
                   'bitscore': 7.8, 'direction': 'Forward', 'hmm_coord_from': 90, 'hmm_coord_to': 171,
                   'ali_coord_from': 607, 'ali_coord_to': 689, 'env_coord_from': 603, 'env_coord_to': 702},
                  {'hmm_name': 'hit4', 'hmm_accession': 'NF038150.1', 'hmm_len': 851, 'evalue': 6.7e-23,
                   'bitscore': 126.8, 'direction': 'Forward', 'hmm_coord_from': 203, 'hmm_coord_to': 454,
                   'ali_coord_from': 597, 'ali_coord_to': 843, 'env_coord_from': 582, 'env_coord_to': 880},
                  {'hmm_name': 'hit5', 'hmm_accession': 'TIGR03903.1', 'hmm_len': 1267, 'evalue': 5.7e-23,
                   'bitscore': 66.6, 'direction': 'Forward', 'hmm_coord_from': 3, 'hmm_coord_to': 190,
                   'ali_coord_from': 607, 'ali_coord_to': 791, 'env_coord_from': 605, 'env_coord_to': 799},
                  {'hmm_name': 'hit6', 'hmm_accession': 'NF033442.0', 'hmm_len': 1391, 'evalue': 5.7e-23,
                   'bitscore': 90, 'direction': 'Forward', 'hmm_coord_from': 258, 'hmm_coord_to': 426,
                   'ali_coord_from': 629, 'ali_coord_to': 788, 'env_coord_from': 626, 'env_coord_to': 794}]
    hmm_pro.sort_hits(query_hits,sorting_type='bitscore',sorting_class='processor')
