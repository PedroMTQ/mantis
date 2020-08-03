try:
    from source.cython_src.get_non_overlapping_hits import get_non_overlapping_hits
    from source.MANTIS_Assembler import *
except:
    from cython_src.get_non_overlapping_hits import get_non_overlapping_hits
    from MANTIS_Assembler import *

#This class will process the domtblout output and get the best hit for our queries, it is inherited by the MANTIS

class MANTIS_Processor():


    def create_chunk_hmmer_dirs(self,chunk_dir):
        for hmmer_dir in ['output_hmmer','domtblout']:#,'tblout']:
           if not os.path.exists(chunk_dir+hmmer_dir):
                Path(chunk_dir+hmmer_dir).mkdir(parents=True, exist_ok=True)

    ######SPLITTING TARGET FASTA INTO CHUNKS######
    def chunk_dict_generator(self,protein_seqs,seq_chunks):
        res={}
        for seq in seq_chunks:
            res[seq]=protein_seqs[seq]
        return res



    ######CHECKING HMMER JOBS######

    def get_taxon_chunks(self,taxon_id,domtblout_path):
        all_domtblout_with_chunks = os.listdir(domtblout_path)
        res=[]
        for domtblout in all_domtblout_with_chunks:
            if 'NOGT'+str(taxon_id) in domtblout:
                res.append(domtblout)
        return res

    def taxon_annotation_finished(self,target_hmm,output_folder,chunks_n):
        c=0
        domtblout_folder=add_slash(add_slash(output_folder)+'domtblout')
        domtblout_files=os.listdir(domtblout_folder)
        for domtblout in domtblout_files:
            if target_hmm in domtblout:
                if '_finished' in domtblout: c+=1
        if c==chunks_n:
            for domtblout in domtblout_files:
                if target_hmm in domtblout:
                    if '_finished' in domtblout:
                        os.rename(domtblout_folder+domtblout,domtblout_folder+domtblout.strip('_finished'))
            return True
        else: return False


    def save_temp_fasta_length(self, chunk_dir, hmm, count_seqs):
        with open(chunk_dir + 'missing_annotation.length', 'a+') as file:
            file.write(hmm + '\t' + str(count_seqs) + '\n')


    ########To merge domtblout

    def group_output_chunks(self, domtblout_path, all_domtblout_with_chunks,chunk_suffix):
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
        grouped_domtblout = self.group_output_chunks(domtblout_path, all_domtblout_with_chunks,chunk_suffix)
        for domtblout in grouped_domtblout:
            concat_files(output_file=domtblout_path + domtblout, list_file_paths=grouped_domtblout[domtblout],
                         stdout_file=stdout_file)
            for to_delete in grouped_domtblout[domtblout]:
                os.remove(to_delete)

    def split_hits(self,domtblout_path,worker_count):
        #split the hits into different chunk allows for memory saving during hit processing
        list_of_files=[domtblout_path+'_chunk_'+str(i) for i in range(worker_count)]
        hit_to_file={}
        file_yielder = yield_file(list_of_files)
        with open(domtblout_path) as file:
            original_line = file.readline()
            line = str(original_line).strip('\n')
            while line:
                if line[0] != '#':
                    line = line.split()[0:22]
                    if len(line) == 22:
                        query_name = line[0]
                        if query_name not in hit_to_file:
                            file_name = next(file_yielder)
                            hit_to_file[query_name]=file_name
                        opened_file= open(hit_to_file[query_name],'a+')
                        opened_file.write(original_line)

                original_line=file.readline()
                line = str(original_line).strip('\n')
        os.remove(domtblout_path)

    ########To calculate new evalue

    def get_temp_fasta_length(self,chunk_dir,domtblout):
        hmm=domtblout.split('_')[0]
        with open(chunk_dir+'missing_annotation.length') as file:
            line=file.readline()
            while line:
                hmm_key,hmm_len=line.split('\t')
                if hmm ==hmm_key: return int(hmm_len)
                line=file.readline()

    def remove_temp_fasta_length(self,chunk_dir):
        if self.file_exists(chunk_dir+'missing_annotation.length'):
            os.remove(chunk_dir+'missing_annotation.length')


    ########Processing protein fasta
    def remove_temp_fasta(self, temp_fasta_path):
        if 'missing_annotations.fa.tmp' in temp_fasta_path:
            os.remove(temp_fasta_path)

    def remove_annotated_queries(self, missing_queries,annotated_queries):
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
        with open(protein_fasta_path,'r') as file:
            line=file.readline()
            if line[0]!='>':
                raise FileNotFoundError
            start_recording = False
            query = None
            while line:
                line=line.strip('\n')
                if line:
                    res, query, start_recording = self.process_protein_fasta_line(res=res,
                                                                                  query=query,
                                                                                  fasta_line=line,
                                                                                  start_recording=start_recording)
                line=file.readline()
        return res


    def generate_temp_fasta(self, missing_queries,output_folder):
        temp_path = output_folder + 'missing_annotations.fa.tmp'
        if os.path.exists(temp_path):
            self.remove_temp_fasta(temp_path)
        with open(temp_path,'w+') as file:
            for mq in missing_queries:
                chunks = [missing_queries[mq][x:x + 60] for x in range(0, len(missing_queries[mq]), 60)]
                chunk_str = '\n'.join(i for i in chunks)
                file.write('>'+mq+'\n'+chunk_str+'\n')
        return temp_path




    ########Processing HMMER hits

    def is_overlap(self, temp_queries,  current_query):
        if not temp_queries or not current_query: return False
        y = range(current_query['env_coord_from'], current_query['env_coord_to'] + 1)
        for t in temp_queries:
            x = range(t['env_coord_from'], t['env_coord_to'] + 1)
            xs = set(x)
            res = xs.intersection(y)
            if res: return True
            if t['hmm_name']==current_query['hmm_name']:  return True
        return False

    def get_best_hits_approximation(self, query_hits):
        '''
        this is just a lazy implementation when the amount of hits is too big and/or the hits are too small
        even with multiprocessing and cython we may run into computationally unfeasable calculations
        when this happens, we do generate a straightforward "best hit"
        Best hit will take the lowest evalue hit and add the next lowest evalue hit (without overlapping) until we reach a point where this cycle cant be repeated.
        This doesnt effectively calculate the "best hit", just a biased (since we start with lowest evalue as root) approximation
        Still, it is a pretty good approximation anyhow
        '''
        query_hits = self.sort_by_evalue(query_hits)
        combo=[]
        while query_hits:
            next_hit= query_hits.pop(0)
            if not self.is_overlap(combo,next_hit):
                combo.append(next_hit)
        return combo

    def get_lowest_hit(self,query_hits):
        '''
        this will take the hit with the lowest evalue, regardless if there are multiple domains
        '''
        if not query_hits: return []
        query_hits = self.sort_by_evalue(query_hits)
        lowest_hit= query_hits.pop(0)
        combo=[lowest_hit]
        return combo

    def get_min_max(self,cython_possible_combos,conversion_dict):
        #highest evalue
        minX_evalue=-999
        #lowest evalue
        maxX_evalue=999
        for len_combo in cython_possible_combos:
            for cython_combo in cython_possible_combos[len_combo]:
                combo=self.cython_to_query_hits(cython_combo,conversion_dict)
                for hit in combo:
                    #we dont consider 0 for the scaling, 0 will always be scaled to max/1
                    if hit['i_evalue']:
                        evalue = log10(hit['i_evalue'])
                        if evalue>minX_evalue: minX_evalue=evalue
                        if evalue<maxX_evalue:maxX_evalue=evalue
        return minX_evalue,maxX_evalue

    def min_max_scale_normal(self, X, minX, maxX):
        if minX == maxX: return 1
        return (X - minX) / (maxX - minX)

    def get_combo_score(self,combo,query_length,minX_evalue,maxX_evalue):
        total_evalue=0
        total_coverage=0
        for hit in combo:
            coverage = hit['env_coord_to'] - hit['env_coord_from']+1
            if hit['i_evalue']: evalue = self.min_max_scale_normal(log10(hit['i_evalue']),minX_evalue,maxX_evalue)+0.01
            else:               evalue = 2
            total_coverage+=coverage
            total_evalue+=evalue
        hit_score= total_evalue*(total_coverage/query_length)
        hit_score /= len(combo)
        return hit_score



    def get_best_hits(self, query_hits,query_length,time_limit):
        '''
        we take into consideration the coverage of our hits and their evalue
        steps:
        1- get all possible combinations of hits
        2- check which possible combinations don't overlap
        3- get best combination (best evalue and coverage)
        '''
        ordered_query_hits = self.sort_by_evalue(query_hits)
        cython_hits,conversion_dict=self.query_hits_to_cython(ordered_query_hits)
        cython_possible_combos=get_non_overlapping_hits(cython_hits,time_limit=time_limit)
        #sometimes this calculation is not feasible (when we have too many small hits and the query sequence is too long, the calculation of conbinations would take forever - limit of 5mins)
        if not cython_possible_combos: return None
        best_hit_score=0
        best_combo = None
        minX_evalue,maxX_evalue= self.get_min_max(cython_possible_combos,conversion_dict)
        for len_combo in cython_possible_combos:
            if len_combo:
                for cython_combo in cython_possible_combos[len_combo]:
                    combo=self.cython_to_query_hits(cython_combo,conversion_dict)
                    hit_score=self.get_combo_score(combo,query_length,minX_evalue,maxX_evalue)
                    if not best_combo:
                        best_combo = combo
                        best_hit_score = hit_score
                    elif hit_score > best_hit_score :
                        best_hit_score = hit_score
                        best_combo = combo
        return best_combo

    def query_hits_to_cython(self,query_hits):
        conversion_dict={}
        res=set()
        for hit_i in range(len(query_hits)):
            res.add(tuple([
                    hit_i,
                    query_hits[hit_i]['env_coord_from'],
                    query_hits[hit_i]['env_coord_to'],
                    query_hits[hit_i]['hmm_name']
                    ]))
            conversion_dict[hit_i]=query_hits[hit_i]
        return res,conversion_dict

    def cython_to_query_hits(self,cython_hits,conversion_dict):
        res=[]
        for hit in cython_hits:
            res.append(conversion_dict[hit[0]])
        return res

    def sort_by_evalue(self,query_hits):
        return sorted(query_hits, key=lambda k: k['i_evalue'])

    def calculate_evalue_threshold(self,query_len):
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
        if self.evalue_threshold =='dynamic':
            if query_len<=150:      return 1e-10
            elif query_len>=250:    return 1e-25
            else:                   return float('1e-'+str(ceil(query_len/10)))
        else: return self.evalue_threshold


    def recalculate_evalue(self,i_evalue,count_seqs_chunk,count_seqs_original_file):
        '''
        keep in mind these calculations will result in minor rounding errors!
        For example in 3 different executions for the same query:
        When splitting the fasta into chunks:
         -one time the fasta had 20 sequences and an i-evalue of 5.98e-52
         - another time the fasta had 24 sequences and an i-evalue of 5.7e-52
         - when the fasta was not split we got 5.8 e-52
        Ideally we would not split the fasta into chunks; this would avoid rounding errors but would not allow parallelization
        '''
        e_value=float(i_evalue)
        e_value/=count_seqs_chunk
        e_value*=count_seqs_original_file
        return e_value

    def read_domtblout(self,output_path,count_seqs_chunk,count_seqs_original_file):
        res = {}
        hit_counter = 0
        # reading the output file
        with open(output_path, 'r', encoding='utf8') as file:
            line = file.readline().strip('\n')
            while line:
                if line[0] != '#':
                    # after 22 it's just the free text of the annotation description
                    line = line.split()[0:22]
                    if len(line) == 22:
                        query_name = line[0]
                        query_len = int(line[2])
                        hmm_name = line[3]
                        hmm_accession = line[4]
                        hmm_len = line[5]
                        i_evalue = float(line[12])
                        e_value = self.recalculate_evalue(i_evalue, count_seqs_chunk, count_seqs_original_file)
                        hmm_coord_from = int(line[15])
                        hmm_coord_to = int(line[16])
                        ali_coord_from = int(line[17])
                        ali_coord_to = int(line[18])
                        # we will use envevelop coords as per HMMER's manual recommendation
                        env_coord_from = int(line[19])
                        env_coord_to = int(line[20])
                        # correcting coordinates for 5'-3' and 3'-5'
                        corrected_env_coord_from = env_coord_from if env_coord_from < env_coord_to else env_coord_to
                        corrected_env_coord_to = env_coord_to if env_coord_to > env_coord_from else env_coord_from
                        # set the coordinates according to the allowed overlap value
                        if self.overlap_value > 0.3: self.overlap_value = 0.3
                        final_env_coord_from = ceil(corrected_env_coord_from + corrected_env_coord_from * self.overlap_value / 2)
                        final_env_coord_to = ceil(corrected_env_coord_to - corrected_env_coord_to * self.overlap_value / 2)
                        if query_name not in res: res[query_name] = {'query_len': query_len, 'hits': []}
                        # when processing mg data, the evalue of hits for chunks should be higher than what it really is (because chunk has less seqs = more significant e-value)
                        if e_value <= self.calculate_evalue_threshold(int(query_len)):
                            hit_dict = {'hmm_name': hmm_name,
                                        'hmm_accession': hmm_accession,
                                        'hmm_len': hmm_len,
                                        'i_evalue': e_value,
                                        'hmm_coord_from': hmm_coord_from,
                                        'hmm_coord_to': hmm_coord_to,
                                        'ali_coord_from': ali_coord_from if ali_coord_from < ali_coord_to else ali_coord_to,
                                        'ali_coord_to': ali_coord_to if ali_coord_to > ali_coord_from else ali_coord_from,
                                        'env_coord_from': final_env_coord_from,
                                        'env_coord_to': final_env_coord_to,
                                        }
                            res[query_name]['hits'].append(hit_dict)
                            hit_counter += 1
                line = file.readline().strip('\n')
        return res,hit_counter

    def process_domtblout(self, output_path,count_seqs_chunk,count_seqs_original_file,stdout_path=None):
        '''
        this will read the domtblout file and export:
        query -> all hits with their respective evalue and hit start and end
        after we have this , we invoke get_best_hits to get our result

        we use i_evalue as per the hmmer manual

        this function can read raw and concatenated domtblout files
        '''
        if isinstance(stdout_path,str):  stdout_file=open(stdout_path ,'a+')
        else: stdout_file=stdout_path
        print('------------------------------------------', flush=True, file=stdout_file)
        print('Processing hits for:\n' + output_path, flush=True, file=stdout_file)
        queries_domtblout,hit_counter=self.read_domtblout(output_path,count_seqs_chunk,count_seqs_original_file)
        # processing the hits and getting the best hit/non-overlapping hits
        print('Found '+str(hit_counter)+' hits in:\n'+output_path+'\nWill now get best hits!',flush=True,file=stdout_file)
        approximated_hits=[]
        hmm=get_path_level(output_path,remove_extension=True)
        res_annotation = {}
        for query in queries_domtblout:
            if query not in res_annotation: res_annotation[query] = {}
            list_hits=queries_domtblout[query].pop('hits')
            #self.domain_algorithm defines which algorithm to use, exhaustive for get_best_hits, heuristic for get_best_hits_approximation, and lowest for get_lowest_hit
            if self.domain_algorithm=='heuristic':
                best_hit = self.get_best_hits_approximation(list(list_hits))
            elif self.domain_algorithm=='lowest':
                best_hit = self.get_lowest_hit(list(list_hits))
            else:
                try:
                    best_hit = self.get_best_hits(list(list_hits),queries_domtblout[query]['query_len'],time_limit=60)
                except (TimeoutError,RecursionError):
                    approximated_hits.append(query)
                    best_hit = self.get_best_hits_approximation(list(list_hits))
            queries_domtblout[query]['best_hit']=best_hit
            res_annotation[query][hmm] = queries_domtblout[query]
        if approximated_hits:
            approximated_hits=' ; '.join(approximated_hits)
            print('Some hits in the output file:\n'+output_path+' were approximated. They were:\n'+approximated_hits,flush=True, file=stdout_file)
        print('Finished processing hits for:\n'+ output_path,flush=True, file=stdout_file)
        print('------------------------------------------', flush=True, file=stdout_file)
        if isinstance(stdout_path,str): stdout_file.close()
        return res_annotation



    def save_processed_hits(self, annotation_output,output_dir,domtblout):
        annotation_output_file=output_dir+domtblout.replace('domtblout','pro')
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
                          'Query_length',
                          'Query_hit_start',
                          'Query_hit_end',
                          'HMM_hit_start',
                          'HMM_hit_end',
                          ]
            first_line = '\t'.join(first_line)+'\n'
            output_file.write(first_line)
            for query in annotation_output:
                for hmm_file in annotation_output[query]:
                    #some hmm hits won't pass the i-evalue threshold, so we dont consider them
                    if '_chunk_' in hmm_file: str_hmm_file=hmm_file.split('_chunk_')[0]
                    else: str_hmm_file=hmm_file
                    if annotation_output[query][hmm_file]['best_hit']:
                        for hit in annotation_output[query][hmm_file]['best_hit']:
                            line = [str(query),
                                    str(str_hmm_file),
                                    str(hit['hmm_name']),
                                    str(hit['hmm_accession']),
                                    #all calculations with evalues are done, so just making it more readable
                                    #str(float(format(hit['i_evalue'],'.6g'))),
                                    str(hit['i_evalue']),
                                    str(annotation_output[query][hmm_file]['query_len']),
                                    str(hit['env_coord_from']),
                                    str(hit['env_coord_to']),
                                    str(hit['hmm_coord_from']),
                                    str(hit['hmm_coord_to']),
                                    ]
                            line = '\t'.join(line)
                            output_file.write(line + '\n')

    def merge_target_output(self,output_file,output_folder,chunks_path,stdout_file,same_output=True):
        header=False
        print('Merging output', output_folder+output_file,'from chunks:',[get_path_level(i) for i in chunks_path],flush=True,file=stdout_file)
        with open(output_folder+output_file,'w+') as file:
            for chunk_output in chunks_path:
                if same_output:  target_chunk_output=chunk_output+output_file
                else:            target_chunk_output= chunk_output
                with open(target_chunk_output, 'r') as chunk_file:
                    chunk_line = chunk_file.readline()
                    while chunk_line:
                        if 'Query' in chunk_line:
                            if not header:
                                file.write(chunk_line)
                                header=True
                        else:
                            file.write(chunk_line)
                        chunk_line = chunk_file.readline()



if __name__ == '__main__':
    tests=[
            [
            {'i_evalue':1e-9,'env_coord_from':0,'env_coord_to':10},
            {'i_evalue':1e-10,'env_coord_from':11,'env_coord_to':20},
            {'i_evalue':1e-15,'env_coord_from':21,'env_coord_to':30},
            ],
            [
            {'i_evalue':1e-9,'env_coord_from':0,'env_coord_to':10},
            {'i_evalue':1e-10,'env_coord_from':11,'env_coord_to':20},
            {'i_evalue':1e-50,'env_coord_from':21,'env_coord_to':30},
            ],
            [
            {'i_evalue':1e-9,'env_coord_from':0,'env_coord_to':10},
            {'i_evalue':1e-10,'env_coord_from':11,'env_coord_to':20},
            {'i_evalue':1e-50,'env_coord_from':21,'env_coord_to':51},
            ]    ,
            [
            {'i_evalue':1e-9,'env_coord_from':0,'env_coord_to':10},
            {'i_evalue':1e-10,'env_coord_from':11,'env_coord_to':20},
            {'i_evalue':1e-50,'env_coord_from':21,'env_coord_to':51},
            ],
            [
            {'i_evalue':1e-10,'env_coord_from':11,'env_coord_to':20},
            {'i_evalue':1e-50,'env_coord_from':21,'env_coord_to':51},
            ]
        ]
    maxX=log10(1e-50)
    minX=log10(1e-9)
    hmm_pro=MANTIS_Processor()
    f = '/home/pedroq/Python_projects/DRAX/Data/Annotation_Output/metag_test.faa'
    o='/home/pedroq/Python_projects/DRAX/Data/Annotation_Output/'
    #p = hmm_pro.split_metagenome(f,o)
    hmm_pro.original_target_path='/home/pedroq/Desktop/mantis_mg/custom_mg_loop/fasta_chunks/p0_c0/p0_c0.faa'
    hmm_pro.redirect_verbose=None
    hmm_pro.overlap_value=0.1
    hmm_pro.evalue_threshold=1e-9
    hmm_pro.acceptable_range=0.05
    hmm_pro.process_domtblout('/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis_free/hmm/test_hmm/sample/out/query3/fasta_chunks/p1_c0/domtblout/Resfams-full.domtblout',1000,2000)
    #hmm_pro.read_protein_fasta(hmm_pro.original_target_path)
    #hmm_pro.target_path='/home/pedroq/Desktop/mantis_mg/custom_mg/fasta_chunks/p1_c0/p1_c0.faa'
    #hmm_pro.evalue_threshold=1e-6
    #a=hmm_pro.hit_is_significant(1e-7)
    #print(a)
