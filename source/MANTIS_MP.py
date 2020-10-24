try:
    from source.MANTIS_Assembler import *
    from source.MANTIS_Processor import MANTIS_Processor
    from source.MANTIS_Interpreter import MANTIS_Interpreter
    from source.MANTIS_Consensus import MANTIS_Consensus
except:
    from MANTIS_Assembler import *
    from MANTIS_Processor import MANTIS_Processor
    from MANTIS_Interpreter import MANTIS_Interpreter
    from MANTIS_Consensus import MANTIS_Consensus



class MANTIS_MP(MANTIS_Assembler,MANTIS_Processor,MANTIS_Interpreter,MANTIS_Consensus):

    def prepare_queue_split_sample(self,protein_seqs,seq_chunks,chunk_dir):
        c=0
        for chunk in seq_chunks:
            self.queue.append([self.chunk_dict_generator(protein_seqs,chunk),c,chunk_dir])
            c+=1
        return c

    def worker_split_sample(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            chunk_seqs, chunk_number,chunk_dir = record
            self.generate_fasta_chunks(chunk_seqs, chunk_number,chunk_dir)


    def generate_fasta_chunks(self,protein_seqs,chunk_number,chunk_dir):
        process_number = str(current_process()._name).split('-')[-1]
        chunk_name = 'p' + str(process_number) + '_c' + str(chunk_number)
        current_chunk_dir = add_slash(chunk_dir + chunk_name)
        chunk_path = current_chunk_dir + chunk_name + '.faa'
        Path(current_chunk_dir).mkdir(parents=True, exist_ok=True)
        with open(chunk_path, 'w+') as file:
            for seq_id in protein_seqs:
                chunks = [protein_seqs[seq_id][x:x + 60] for x in range(0, len(protein_seqs[seq_id]), 60)]
                chunk_str = '\n'.join(i for i in chunks)
                file.write('>' + seq_id + '\n' + chunk_str + '\n')


    def set_chunks_to_annotate(self):
        for file_path,output_path,organism_lineage,count_seqs_original_file in self.fastas_to_annotate:
            chunk_dir=output_path+'fasta_chunks'+splitter
            all_chunks = os.listdir(chunk_dir)
            for chunk in all_chunks:
                current_chunk_dir =add_slash(chunk_dir+chunk)
                chunk_name=chunk
                chunk_path=current_chunk_dir +chunk_name+'.faa'
                count_seqs_chunk=get_seqs_count(chunk_path)
                self.chunks_to_annotate.append([chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path])
                if output_path not in self.chunks_to_fasta: self.chunks_to_fasta[output_path]=[]
                self.chunks_to_fasta[output_path].append(current_chunk_dir)

    def split_sample(self,minimum_worker_load=20000,time_limit=300,load_balancing_limit=200000):
        '''
        minimum_worker_load we assume a minumum amount of workload for the workers, since generating results in process spawning overhead
        #if we have a chunk size of 5k, each worker will produce minimum_worker_load/5000
        it will use 5 processes here to split faster

        On generating chunks:
        HMMER performance is determined by several factors, one of which is the length of the query sequence.
        A naive approach but efficient approach to generating chunks would be to split chunks by a moving window <chunk_generator>
        A better approach would be to load balance the chunks by sequence length <chunk_generator_load_balanced>. This will effectively create more even chunks. This is slower and more ram consuming though since we are pop each sequence at a time and also saving all seq keys in memory

        Just a comparison on the chunk size (bytes) of 1307 chunks between the non balanced and balanced generator (metagenomic sample with 1306215 sequences and 325.1 MB)
                        non_balanced        balanced
        sum	            209394320	        209421151
        average	        160332.557427259	160230.413925019
        min	            118979	            159526
        max	            197024	            163953
        stdev	        11069.918226454	    602.515771601041
        There's obviously some variation, even in the balanced generator, but it's way lower (18 times more deviation)
        Load balancing affects the file size by very little too (around 0.013%)

        After some testing I've found that load balancing doesn't affect the total speed when dealing with metagenomic samples.
        With so many chunks, the theoretical speed gained by balancing the chunks doesn't come into play since we never get idle processes.
        This load balancing will only be useful for smaller sample sizes
        Ordering around 200k sequences length should take no longer than 10 seconds

        A inconvenient outcome of using the <chunk_generator_load_balanced> is the fact that it uses more RAM (since it's not a generator) and requires more time. This is hardly perceptible with smaller samples though.

        '''
        print_cyan('Splitting samples into chunks!',flush=True,file=self.redirect_verbose)
        worker_count=1
        n_chunks=0
        for file_path,output_path,organism_lineage,count_seqs_original_file in self.fastas_to_annotate:
            protein_seqs=self.read_protein_fasta(file_path)
            chunk_dir=output_path+'fasta_chunks'+splitter
            if not os.path.exists(chunk_dir):
                Path(chunk_dir).mkdir(parents=True, exist_ok=True)
            current_worker_count= estimate_number_workers_split_sample(minimum_worker_load,len(protein_seqs))
            chunk_size=estimate_chunk_size(total_n_seqs=len(protein_seqs),
                                           annotation_workers=self.estimate_number_workers_annotation(split_sample=True),
                                           chunk_size=self.chunk_size,
                                           )
            if current_worker_count> worker_count: worker_count=current_worker_count
            if len(protein_seqs)<load_balancing_limit:
                proteins_seqs_keys_len={i:len(protein_seqs[i]) for i in protein_seqs}
                list_ordered = sorted(proteins_seqs_keys_len, key=proteins_seqs_keys_len.__getitem__)
                seq_chunks= chunk_generator_load_balanced(list_ordered, chunk_size,time_limit=time_limit)
            else:
                proteins_seqs_keys=list(protein_seqs.keys())
                seq_chunks= chunk_generator(proteins_seqs_keys, chunk_size)
            current_chunks= self.prepare_queue_split_sample(protein_seqs,seq_chunks,chunk_dir)
            n_chunks+=current_chunks
            stdout_file= open(output_path + 'Mantis.out','a+')
            print( 'The current sample: '+file_path+' will be split into ' + str(current_chunks) + ' chunks (up to '+str(chunk_size)+' sequences each), which will be stored at:\n'+chunk_dir,flush=True, file=stdout_file)
            stdout_file.close()
        if worker_count<environment_cores*worker_per_core:
            if len(self.fastas_to_annotate)<environment_cores*worker_per_core:
                worker_count=len(self.fastas_to_annotate)
            else:
                worker_count=environment_cores*worker_per_core
        print_cyan('Samples will be split into '+str(n_chunks)+' chunks with '+str(worker_count)+' workers',flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_split_sample,worker_count)




    ####To run HMMER


    def compile_annotation_job(self, hmm_path, target_path,output_folder, output_initials=''):
        hmm = get_path_level(hmm_path)
        hmm = hmm.split('.')[0]
        # what is more efficient? hmmsearch or hmmscan? hmmsearch: https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/
        command = 'hmmsearch '
        # summarized output
        #command += ' --tblout ' + output_folder + 'tblout'+splitter+output_initials +hmm+'.tblout'
        # output with coordinates of hits
        domtblout_path=output_folder + 'domtblout' + splitter + output_initials + hmm + '.domtblout'
        command += ' --domtblout ' + domtblout_path
        # hmm requires a master thread, so we deduct it from the number of threads
        command += ' --cpu ' + str(self.hmmer_threads)
        #since we split the original fasta into chunks, hmmer might remove some hits ( I correct for this further down the line, but not in hmmer)
        #even when using the default evalue threshold, there isn't much of a loss
        if self.evalue_threshold=='dynamic':
            command += ' -E ' + str(1e-8)
        elif self.evalue_threshold:
            command += ' -E ' + str(self.evalue_threshold/10)
        else:
            #hmmers default evalue threshold will be 1e-6, Mantis will be more strict further down the line - this is just to allow splitting up the file into chunks
            command += ' -E ' + str(1e-3)
        command += ' --notextw '
        command += ' ' + hmm_path
        command += ' ' + target_path
        console_stdout = output_folder + 'output_hmmer' +splitter+ output_initials+hmm+'.out'
        return command, domtblout_path, console_stdout


    #####For the general HMMS

    def calculate_total_hmms_annotation(self):
        #some lineage taxons (since we might fully annotate a fasta before the top taxon level) might be unnecessary but this should provide a good estimate
        n_hmms = 0
        hmm_list = self.compile_hmms_list()
        #this will take into account hmms that are split into chunks
        for hmm in hmm_list:
            n_hmms += len(compile_hmm_chunks_path(hmm))
        if self.mantis_paths['NOGG'][0:2] != 'NA': n_hmms+=len(compile_hmm_chunks_path(self.mantis_paths['NOGG']))
        if self.mantis_paths['NOGT'][0:2] != 'NA':
            tax_hmms=0
            for file_path,output_path,organism_lineage,count_seqs_original_file in self.fastas_to_annotate:
                organism_lineage_temp = list(organism_lineage)
                if organism_lineage_temp:
                    current_taxon = organism_lineage_temp.pop(-1)
                    hmm_path = self.get_lineage_hmm_path(current_taxon)
                    # to skip taxons without an hmm
                    while not hmm_path and organism_lineage_temp:
                        current_taxon = organism_lineage_temp.pop(-1)
                        hmm_path = self.get_lineage_hmm_path(current_taxon)
                    if hmm_path:
                        chunks_path = compile_hmm_chunks_path(hmm_path)
                        tax_hmms+=len(chunks_path)
            n_hmms+=int(tax_hmms/len(self.fastas_to_annotate))
        self.total_hmms_annotation = n_hmms
        return n_hmms

    def estimate_number_workers_annotation(self,split_sample=False):
        if not hasattr(self,'total_hmms_annotation'):
            n_hmms=self.calculate_total_hmms_annotation()
        else:
            n_hmms = self.total_hmms_annotation
        return estimate_number_workers_annotation(n_chunks=len(self.chunks_to_annotate),
                                                    n_hmms=n_hmms,
                                                    default_workers=self.default_workers,
                                                    user_cores=self.user_cores,
                                                    split_sample=split_sample,
                                                    )

    def run_hmmer(self):
        worker_count=self.estimate_number_workers_annotation()
        self.prepare_queue_hmmer()
        print_cyan('Running HMMER with '+str(worker_count)+' workers. HMMER will run around '+str(len(self.chunks_to_annotate)*self.total_hmms_annotation)+' times (lineage annotation may change this number)',flush=True,file=self.redirect_verbose)
        self.processes_handler(self.worker_hmmer,worker_count)




    def run_hmmer_annotation(self,hmmer_command,hmmer_stdout_path,stdout_path,master_pid,output_file=None):
        hmmer_stdout_file = open(hmmer_stdout_path, 'w+')
        stdout_file = open(stdout_path, 'a+')
        print('Running HMMER command:\n', hmmer_command, flush=True, file=stdout_file)
        start_time = time()
        run_command(hmmer_command, stdout_file=hmmer_stdout_file,master_pid=master_pid,wanted_child='hmmsearch',user_memory=self.user_memory)
        print('Finished running HMMER (' + str(round(time() - start_time,3)) + ' seconds):\n', hmmer_command, flush=True,file=stdout_file)
        hmmer_stdout_file.close()
        stdout_file.close()
        if output_file:
            move_file(output_file,output_file+'_finished')



    def worker_hmmer(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            if record[0]=='General':
                _,hmmer_command,hmmer_stdout_path,stdout_path = record
                self.run_hmmer_annotation(hmmer_command=hmmer_command,hmmer_stdout_path=hmmer_stdout_path,stdout_path=stdout_path,master_pid=master_pid)
            #for taxon runs
            elif record[0]=='NOGT':
                _,hmmer_command, hmmer_stdout_path, stdout_path, output_file = record
                self.run_hmmer_annotation(hmmer_command=hmmer_command, hmmer_stdout_path=hmmer_stdout_path, stdout_path=stdout_path, master_pid=master_pid,output_file=output_file)
            elif record[0]=='NOGG_checkpoint':
                _,current_chunk_dir,fasta_path,count_seqs_original_file,chunks_n = record
                if self.taxon_annotation_finished('NOGG',current_chunk_dir,chunks_n):
                    protein_sequences = self.read_protein_fasta(fasta_path)
                    count_seqs_chunk = len(protein_sequences)
                    self.remove_temp_fasta(fasta_path)
                else:
                   self.queue.insert(0, record)
            elif record[0]=='NOGT_checkpoint':
                _,current_chunk_dir,fasta_path,taxon_id,organism_lineage,count_seqs_original_file,chunks_n,stdout_path = record
                if self.taxon_annotation_finished('NOGT'+str(taxon_id),current_chunk_dir,chunks_n):
                    protein_sequences = self.read_protein_fasta(fasta_path)
                    count_seqs_chunk = len(protein_sequences)
                    domtblout_path = current_chunk_dir + 'domtblout' + splitter
                    taxon_domtblouts = self.get_taxon_chunks(taxon_id,domtblout_path)
                    stdout_file= open(stdout_path,'a+')
                    #while merging here is not optimal, we are only using NOGT for small samples (taxonomically classified) so it shouldnt consume too much memory anyhow
                    self.merge_output_chunks(domtblout_path,taxon_domtblouts,chunk_suffix='.domtblout',stdout_file=stdout_file)
                    stdout_file.close()
                    annotated_queries = self.process_domtblout(output_path=domtblout_path+'NOGT'+str(taxon_id)+'_merged.domtblout',count_seqs_chunk=count_seqs_chunk,count_seqs_original_file=count_seqs_original_file,stdout_path=stdout_path)
                    # in each iteration we only annotate the missing sequences
                    self.remove_annotated_queries(protein_sequences,annotated_queries)
                    if protein_sequences:
                        fasta_path = self.generate_temp_fasta(protein_sequences,current_chunk_dir)
                        self.add_to_queue_lineage_annotation(fasta_path,current_chunk_dir,list(organism_lineage),count_seqs_original_file,stdout_path)
                #if annotations havent finished, we add the checker back into the queue
                else:
                    self.queue.insert(0, record)


    def prepare_queue_hmmer(self):
        hmms_list=self.compile_hmms_list()
        chunked_hmms_list=[]
        for hmm_path in hmms_list:
            chunked_hmms_list.extend(compile_hmm_chunks_path(hmm_path))
        chunked_hmms_list,chunked_hmms_list_size=self.order_by_size_descending(chunked_hmms_list)
        # this will build the hmmer processes to run as well as give the domtblout we want
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            self.create_chunk_hmmer_dirs(current_chunk_dir)
            self.add_to_queue_lineage_annotation(chunk_path,current_chunk_dir, list(organism_lineage),count_seqs_original_file, output_path + 'Mantis.out')
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            for hmm_path in chunked_hmms_list:
                # full hmmer command to be run with subprocess
                command, output_file, console_stdout = self.compile_annotation_job(hmm_path,target_path=chunk_path,output_folder=current_chunk_dir)
                # adding our hmmer command to be consumed by the hmmer processes later on
                self.queue.append(['General',command,console_stdout,output_path + 'Mantis.out'])

    ####For the lineage HMMs


    def add_to_queue_lineage_annotation(self,fasta_path,current_chunk_dir,organism_lineage,count_seqs_original_file,stdout_path):
        protein_sequences = self.read_protein_fasta(fasta_path)
        count_seqs_chunk = len(protein_sequences)
        if count_seqs_chunk:
            if organism_lineage:
                if self.mantis_paths['NOGT'][0:2]!='NA':
                    current_taxon = organism_lineage.pop(-1)
                    hmm_path = self.get_lineage_hmm_path(current_taxon)
                    #to skip taxons without an hmm
                    while not hmm_path and organism_lineage:
                        current_taxon = organism_lineage.pop(-1)
                        hmm_path = self.get_lineage_hmm_path(current_taxon)
                    if hmm_path:
                        chunks_path = compile_hmm_chunks_path(hmm_path)
                        for chunk_hmm in chunks_path:
                            command, output_file,console_stdout = self.compile_annotation_job(chunk_hmm,target_path=fasta_path,output_initials='NOGT',output_folder=current_chunk_dir)
                            self.queue.insert(0,['NOGT',command, console_stdout, stdout_path,output_file])
                        #will be used for checking whether chunks have been annotated
                        self.queue.insert(len(chunks_path),['NOGT_checkpoint',current_chunk_dir,fasta_path,current_taxon,organism_lineage,count_seqs_original_file,len(chunks_path),stdout_path])
                        self.save_temp_fasta_length(current_chunk_dir, 'NOGT'+str(current_taxon) , count_seqs_chunk)

                    else:
                        # if there are still missing annotations from the lineage annotation or there's not taxonomic classification we query against the whole nog database
                        if self.mantis_paths['NOGG'][0:2] != 'NA':
                            hmm_path = get_hmm_in_folder(self.mantis_paths['NOGG'])
                            chunks_path = compile_hmm_chunks_path(hmm_path)
                            for chunk_hmm in chunks_path:
                                command, output_file, console_stdout = self.compile_annotation_job(chunk_hmm, target_path=fasta_path, output_folder=current_chunk_dir)
                                self.queue.insert(0, ['NOGT', command, console_stdout, stdout_path, output_file])
                            self.queue.insert(len(chunks_path),['NOGG_checkpoint', current_chunk_dir, fasta_path, count_seqs_original_file, len(chunks_path)])
                            self.save_temp_fasta_length(current_chunk_dir, 'NOGG', count_seqs_chunk)
            else:
                # if there are still missing annotations from the lineage annotation or there's not taxonomic classification we query against the whole nog database
                if self.mantis_paths['NOGG'][0:2] != 'NA':
                    hmm_path = get_hmm_in_folder(self.mantis_paths['NOGG'])
                    chunks_path = compile_hmm_chunks_path(hmm_path)
                    for chunk_hmm in chunks_path:
                        command, output_file, console_stdout = self.compile_annotation_job(chunk_hmm, target_path=fasta_path, output_folder=current_chunk_dir)
                        self.queue.insert(0,['NOGT',command, console_stdout, stdout_path,output_file])
                    self.queue.insert(len(chunks_path),['NOGG_checkpoint', current_chunk_dir, fasta_path, count_seqs_original_file, len(chunks_path)])
                    self.save_temp_fasta_length(current_chunk_dir, 'NOGG' , count_seqs_chunk)

    ####Merging hmmer output






    def estimate_domtblouts_per_chunk(self):
        n_hmms = len(self.compile_hmms_list())
        if self.mantis_paths['NOGT'][0:2] != 'NA': n_hmms+=1
        if self.mantis_paths['NOGG'][0:2] != 'NA': n_hmms+=1
        return n_hmms

    def process_output(self):
        domtblout_per_chunks=self.estimate_domtblouts_per_chunk()
        worker_count = estimate_number_workers_process_output(n_chunks=len(self.chunks_to_annotate),domtblout_per_chunks=domtblout_per_chunks)
        print_cyan('Processing output with ' + str(worker_count) + ' workers.', flush=True, file=self.redirect_verbose)

        #this needs to be merged at this point so that we can properly get the best hits
        #Since an HMM might be split into chunks we merge all the results from the same HMM to a single file
        #this means that before we could have hmm_chunk_1 - hit_1, hmm_chunk_2- hit_1. Now all this info is in one file
        self.prepare_queue_merge_domtblout()
        self.processes_handler(self.worker_merge_domtblout, worker_count)

        #However when we have a lot of hits for that HMM file the hit processing can be quite memory heavy, so instead we now split hits into chunks
        #This process is quite light since it only stores the file the hit should be stored at, all the hit information is read and discarded from memory
        #this also allows for parallelization
        self.prepare_queue_split_hits(worker_count)
        self.processes_handler(self.worker_split_hits, worker_count)


        self.prepare_queue_process_output()
        self.processes_handler(self.worker_process_output, worker_count)


        self.prepare_queue_merge_output()
        self.processes_handler(self.worker_merge_output, worker_count)


        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            self.remove_temp_fasta_length(current_chunk_dir)


    def prepare_queue_merge_domtblout(self):
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            chunks_output_path= current_chunk_dir+'domtblout'+splitter
            all_output_with_chunks = os.listdir(chunks_output_path)
            self.queue.append([chunks_output_path,all_output_with_chunks,output_path+'Mantis.out'])

    def worker_merge_domtblout(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            chunks_output_path,all_output_with_chunks,stdout_path = record
            stdout_file=open(stdout_path,'a+')
            self.merge_output_chunks(chunks_output_path,all_output_with_chunks,chunk_suffix='.domtblout',stdout_file=stdout_file)
            stdout_file.close()

    def prepare_queue_merge_processed_output(self):
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            chunks_output_path= current_chunk_dir+'processed_output'+splitter
            all_output_with_chunks = os.listdir(chunks_output_path)
            self.queue.append([chunks_output_path,all_output_with_chunks,output_path+'Mantis.out'])

    def worker_merge_processed_output(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            chunks_output_path,all_output_with_chunks,stdout_path = record
            stdout_file=open(stdout_path,'a+')
            self.merge_output_chunks(chunks_output_path,all_output_with_chunks,chunk_suffix='.pro',stdout_file=stdout_file)
            stdout_file.close()

    def prepare_queue_split_hits(self,worker_count):
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            domtblout_path= current_chunk_dir+'domtblout'+splitter
            all_domtblout = os.listdir(domtblout_path)
            for domtblout in all_domtblout:
                self.queue.append([domtblout_path+domtblout,current_chunk_dir,worker_count,output_path + 'Mantis.out'])

    def worker_split_hits(self,queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            domtblout_path,current_chunk_dir,worker_count,stdout_path = record
            self.split_hits(domtblout_path,worker_count)




    def prepare_queue_process_output(self):
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            domtblout_path= current_chunk_dir+'domtblout'+splitter
            all_domtblout = os.listdir(domtblout_path)
            if not os.path.exists(add_slash(add_slash(current_chunk_dir)+'processed_output')):
                Path(add_slash(add_slash(current_chunk_dir)+'processed_output')).mkdir(parents=True, exist_ok=True)
            for domtblout in all_domtblout:
                if 'NOGT' in domtblout or 'NOGG' in domtblout:
                    count_seqs_chunk_domtblout =self.get_temp_fasta_length(current_chunk_dir,domtblout)
                else:
                    count_seqs_chunk_domtblout= int(count_seqs_chunk)
                self.queue.append([domtblout_path+domtblout,current_chunk_dir,count_seqs_chunk_domtblout,count_seqs_original_file,output_path + 'Mantis.out'])

    def worker_process_output(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            domtblout_path,current_chunk_dir,count_seqs_chunk,count_seqs_original_file,stdout_path = record
            processed_hits= self.process_domtblout(output_path=domtblout_path,count_seqs_chunk=count_seqs_chunk,count_seqs_original_file=count_seqs_original_file,stdout_path=stdout_path)
            self.save_processed_hits(processed_hits,add_slash(add_slash(current_chunk_dir)+'processed_output'),domtblout=get_path_level(domtblout_path))

    def prepare_queue_merge_output(self):
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            self.queue.append([current_chunk_dir,output_path])

    def worker_merge_output(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            current_chunk_dir,output_path = record
            chunks_path= add_slash(current_chunk_dir+'processed_output')
            chunks_to_merge=[chunks_path+i for i in os.listdir(chunks_path)]
            stdout_file=open(output_path+'Mantis.out','a+')
            self.merge_target_output('output_annotation.tsv',current_chunk_dir,chunks_to_merge,stdout_file,same_output=False)
            stdout_file.close()






    ###Interpreting output

    def interpret_output(self):
        worker_count = estimate_number_workers_process_output(n_chunks=len(self.chunks_to_annotate))
        self.prepare_queue_interpret_output()
        print_cyan('Interpreting output with '+str(worker_count)+' workers.', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_interpret_output, worker_count)


    def prepare_queue_interpret_output(self):
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            output_annotation_tsv= current_chunk_dir+'output_annotation.tsv'
            self.queue.append([output_annotation_tsv,current_chunk_dir])


    def worker_interpret_output(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            output_annotation_tsv,current_chunk_dir = record
            interpreted_annotation_tsv= current_chunk_dir+'integrated_annotation.tsv'
            self.generate_interpreted_output(output_annotation_tsv,interpreted_annotation_tsv)



    ###Generate consensus output

    def get_consensus_output(self):
        MANTIS_Consensus.__init__(self)
        worker_count = estimate_number_workers_process_output(n_chunks=len(self.chunks_to_annotate))
        self.prepare_queue_generate_consensus()
        print_cyan('Generating consensus output with '+str(worker_count)+' workers.', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_consensus_output, worker_count)


    def prepare_queue_generate_consensus(self):
        for chunk_name,chunk_path,current_chunk_dir,organism_lineage,count_seqs_chunk,count_seqs_original_file,output_path in self.chunks_to_annotate:
            interepreted_annotation_tsv= current_chunk_dir+'integrated_annotation.tsv'
            stdout_file_path=output_path + 'Mantis.out'
            self.queue.append([interepreted_annotation_tsv,current_chunk_dir,stdout_file_path])

    def worker_consensus_output(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            interpreted_annotation_tsv,current_chunk_dir,stdout_file_path = record
            consensus_annotation_tsv= current_chunk_dir+'consensus_annotation.tsv'
            self.generate_consensus_output(interpreted_annotation_tsv,consensus_annotation_tsv,stdout_file_path)

    #Merging Mantis output

    def merge_mantis_output(self):
        worker_count = estimate_number_workers_process_output(n_chunks=len(self.chunks_to_annotate))
        self.prepare_queue_merge_mantis_output()
        print_cyan('Merging output with '+str(worker_count)+' workers.', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_merge_mantis_output, worker_count)


    def prepare_queue_merge_mantis_output(self):
        for output_path in self.chunks_to_fasta:
            self.queue.append([output_path,self.chunks_to_fasta[output_path]])

    def worker_merge_mantis_output(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            output_path,chunks_path = record
            chunks_output=[i+'output_annotation.tsv' for i in chunks_path]
            self.merge_chunks_outputs(output_path,chunks_path)

    def merge_chunks_outputs(self,output_path,chunks_path):
        stdout_file = open(output_path + 'Mantis.out', 'a+')
        self.merge_target_output(output_file='output_annotation.tsv',output_folder=output_path,chunks_path=chunks_path,stdout_file=stdout_file)
        self.merge_target_output(output_file='integrated_annotation.tsv',output_folder=output_path,chunks_path=chunks_path,stdout_file=stdout_file)
        if not self.skip_consensus:
            self.merge_target_output(output_file='consensus_annotation.tsv',output_folder=output_path,chunks_path=chunks_path,stdout_file=stdout_file)
        print('------------------------------------------', flush=True, file=stdout_file)
        print_cyan('This sample has been sucessfully annotated!', flush=True, file=stdout_file)
        stdout_file.close()



