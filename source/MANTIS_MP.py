try:
    from source.MANTIS_Assembler import *
    from source.MANTIS_Processor import MANTIS_Processor
    from source.MANTIS_Metadata import MANTIS_Metadata
    from source.MANTIS_Consensus import MANTIS_Consensus
except:
    from MANTIS_Assembler import *
    from MANTIS_Processor import MANTIS_Processor
    from MANTIS_Metadata import MANTIS_Metadata
    from MANTIS_Consensus import MANTIS_Consensus


class MANTIS_MP(MANTIS_Assembler, MANTIS_Processor, MANTIS_Metadata, MANTIS_Consensus):

    def prepare_queue_split_sample(self, protein_seqs, seq_chunks, chunk_dir):
        c = 0
        for chunk in seq_chunks:
            self.queue.append([self.chunk_dict_generator(protein_seqs, chunk), c, chunk_dir])
            c += 1
        return c

    def worker_split_sample(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            chunk_seqs, chunk_number, chunk_dir = record
            self.generate_fasta_chunks(chunk_seqs, chunk_number, chunk_dir)

    def generate_fasta_chunks(self, protein_seqs, chunk_number, chunk_dir):
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
        for file_path, output_path, organism_lineage, count_seqs_original_file in self.fastas_to_annotate:
            chunk_dir = add_slash(output_path + 'fasta_chunks')
            all_chunks = os.listdir(chunk_dir)
            for chunk in all_chunks:
                current_chunk_dir = add_slash(chunk_dir + chunk)
                chunk_name = chunk
                chunk_path = current_chunk_dir + chunk_name + '.faa'
                count_seqs_chunk = get_seqs_count(chunk_path)
                self.chunks_to_annotate.append(
                    [chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk,
                     count_seqs_original_file, output_path])
                if output_path not in self.chunks_to_fasta: self.chunks_to_fasta[output_path] = []
                self.chunks_to_fasta[output_path].append(current_chunk_dir)

    def split_sample(self, minimum_worker_load=20000, time_limit=300, load_balancing_limit=200000):
        '''
        minimum_worker_load we assume a minumum amount of workload for the workers, since generating results in process spawning overhead
        #if we have a chunk size of 5k, each worker will produce minimum_worker_load/5000
        it will use 5 processes here to split faster

        On generating chunks:
        HMMER performance is determined by several factors, one of which is the length of the query sequence.
        A naive approach but efficient approach to generating chunks would be to split chunks by a moving window <chunk_generator>
        A better approach would be to load balance the chunks by sequence length <chunk_generator_load_balanced>. This will effectively create more even chunks. This is slower and more ram consuming though since we pop each sequence at a time and also saving all seq keys in memory

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
        print_cyan('Splitting samples into chunks!', flush=True, file=self.redirect_verbose)
        worker_count = 1
        n_chunks = 0
        for file_path, output_path, organism_lineage, count_seqs_original_file in self.fastas_to_annotate:
            protein_seqs = self.read_protein_fasta(file_path)
            chunk_dir = add_slash(output_path + 'fasta_chunks')
            if not os.path.exists(chunk_dir):
                Path(chunk_dir).mkdir(parents=True, exist_ok=True)
            current_worker_count = estimate_number_workers_split_sample(minimum_worker_load, len(protein_seqs))
            chunk_size = estimate_chunk_size(total_n_seqs=len(protein_seqs),
                                             annotation_workers=self.estimate_number_workers_annotation(split_sample=True),
                                             chunk_size=self.chunk_size,
                                             )
            if current_worker_count > worker_count: worker_count = current_worker_count
            if len(protein_seqs) < load_balancing_limit:
                proteins_seqs_keys_len = {i: len(protein_seqs[i]) for i in protein_seqs}
                list_ordered = sorted(proteins_seqs_keys_len, key=proteins_seqs_keys_len.__getitem__)
                seq_chunks = chunk_generator_load_balanced(list_ordered, chunk_size, time_limit=time_limit)
            else:
                proteins_seqs_keys = list(protein_seqs.keys())
                seq_chunks = chunk_generator(proteins_seqs_keys, chunk_size)
            current_chunks = self.prepare_queue_split_sample(protein_seqs, seq_chunks, chunk_dir)
            n_chunks += current_chunks
            stdout_file = open(output_path + 'Mantis.out', 'a+')
            print('The current sample: ' + file_path + ' will be split into ' + str(
                current_chunks) + ' chunks (up to ' + str(
                chunk_size) + ' sequences each), which will be stored at:\n' + chunk_dir, flush=True, file=stdout_file)
            stdout_file.close()
        if worker_count < ENVIRONMENT_CORES * WORKER_PER_CORE:
            if len(self.fastas_to_annotate) < ENVIRONMENT_CORES * WORKER_PER_CORE:
                worker_count = len(self.fastas_to_annotate)
            else:
                worker_count = ENVIRONMENT_CORES * WORKER_PER_CORE
        print_cyan('Samples will be split into ' + str(n_chunks) + ' chunks with ' + str(worker_count) + ' workers',
                   flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_split_sample, worker_count)

    ####To run HMMER

    def compile_annotation_job(self, hmm_path, target_path, output_folder, output_initials=''):
        hmm = get_path_level(hmm_path)
        hmm = hmm.split('.')[0]
        # what is more efficient? hmmsearch or hmmscan? hmmsearch: https://cryptogenomicon.org/2011/05/27/hmmscan-vs-hmmsearch-speed-the-numerology/
        command = 'hmmsearch '
        console_stdout = output_folder + 'output_hmmer' + SPLITTER + output_initials + hmm + '.out'
        if self.keep_files:
            command += '-o '+console_stdout
        else:
            command += '-o /dev/null'
        # summarized output
        # command += ' --tblout ' + output_folder + 'tblout'+SPLITTER+output_initials +hmm+'.tblout'
        # output with coordinates of hits
        domtblout_path = output_folder + 'domtblout' + SPLITTER + output_initials + hmm + '.domtblout'
        command += ' --domtblout ' + domtblout_path
        # hmm requires a master thread, so we deduct it from the number of threads
        command += ' --cpu ' + str(self.hmmer_threads)
        # since we split the original fasta into chunks, hmmer might remove some hits ( I correct for this further down the line, but not in hmmer)
        # even when using the default evalue threshold, there isn't much of a loss
        # we use domE because we accept multiple hits with these two algorithms
        if self.domain_algorithm in ['dfs', 'heuristic']:
            threshold_type = '--domE'
        # whereas bpo we only accept one
        else:
            threshold_type = '-E'
        if self.evalue_threshold == 'dynamic':
            command += ' ' + threshold_type + ' ' + str(1e-6)
        else:
            command += ' ' + threshold_type + ' ' + str(self.evalue_threshold * 10)
        command += ' --notextw '
        command += ' ' + hmm_path
        command += ' ' + target_path
        return command, domtblout_path

    #####For the general HMMS

    def calculate_total_hmms_annotation(self):
        # some lineage taxons (since we might fully annotate a fasta before the top taxon level) might be unnecessary but this should provide a good estimate
        n_hmms = 0
        hmm_list = self.compile_hmms_list()
        # this will take into account hmms that are split into chunks
        for hmm in hmm_list:
            n_hmms += len(compile_hmm_chunks_path(hmm))
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            n_hmms += len(compile_hmm_chunks_path(add_slash(self.mantis_paths['NCBI'] + 'NCBIG')))
        if self.mantis_paths['NOG'][0:2] != 'NA':
            n_hmms += len(compile_hmm_chunks_path(add_slash(self.mantis_paths['NOG'] + 'NOGG')))
        if self.mantis_paths['NOG'][0:2] != 'NA':
            tax_hmms = 0
            for file_path, output_path, organism_lineage, count_seqs_original_file in self.fastas_to_annotate:
                organism_lineage_temp = list(organism_lineage)
                if organism_lineage_temp:
                    current_taxon = organism_lineage_temp.pop(-1)
                    hmm_path = self.get_lineage_hmm_path(current_taxon, db='NOG')
                    # to skip taxons without an hmm
                    while not hmm_path and organism_lineage_temp:
                        current_taxon = organism_lineage_temp.pop(-1)
                        hmm_path = self.get_lineage_hmm_path(current_taxon, db='NOG')
                    if hmm_path:
                        chunks_path = compile_hmm_chunks_path(hmm_path)
                        tax_hmms += len(chunks_path)
            n_hmms += int(tax_hmms / len(self.fastas_to_annotate))
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            tax_hmms = 0
            for file_path, output_path, organism_lineage, count_seqs_original_file in self.fastas_to_annotate:
                organism_lineage_temp = list(organism_lineage)
                if organism_lineage_temp:
                    current_taxon = organism_lineage_temp.pop(-1)
                    hmm_path = self.get_lineage_hmm_path(current_taxon, db='NCBI')
                    # to skip taxons without an hmm
                    while not hmm_path and organism_lineage_temp:
                        current_taxon = organism_lineage_temp.pop(-1)
                        hmm_path = self.get_lineage_hmm_path(current_taxon, db='NCBI')
                    if hmm_path:
                        chunks_path = compile_hmm_chunks_path(hmm_path)
                        tax_hmms += len(chunks_path)
            n_hmms += int(tax_hmms / len(self.fastas_to_annotate))
        self.total_hmms_annotation = n_hmms
        return n_hmms

    def estimate_number_workers_annotation(self, split_sample=False):
        if not hasattr(self, 'total_hmms_annotation'):
            n_hmms = self.calculate_total_hmms_annotation()
        else:
            n_hmms = self.total_hmms_annotation
        return estimate_number_workers_annotation(n_chunks=len(self.chunks_to_annotate),
                                                  n_hmms=n_hmms,
                                                  default_workers=self.default_workers,
                                                  user_cores=self.user_cores,
                                                  split_sample=split_sample,
                                                  )

    def run_hmmer(self):
        worker_count = self.estimate_number_workers_annotation()
        self.prepare_queue_hmmer()
        print_cyan('Running HMMER with ' + str(worker_count) + ' workers.', flush=True, file=self.redirect_verbose)
        if worker_count < 1:
            print('Invalid number of workers in run_hmmer')
            raise BadNumberWorkers
        self.processes_handler(self.worker_hmmer, worker_count)

    def run_hmmer_annotation(self, hmmer_command,stdout_path, master_pid, output_file=None):
        stdout_file = open(stdout_path, 'a+')
        print('Running HMMER command:\n', hmmer_command, flush=True, file=stdout_file)
        start_time = time()
        if self.skip_managed_memory:
            run_command(hmmer_command, master_pid=None, wanted_child='hmmsearch',user_memory=self.user_memory)
        else:
            run_command(hmmer_command, master_pid=master_pid, wanted_child='hmmsearch',user_memory=self.user_memory)
        print('Finished running HMMER (' + str(round(time() - start_time, 3)) + ' seconds):\n', hmmer_command,
              flush=True, file=stdout_file)
        stdout_file.close()
        if output_file:
            move_file(output_file, output_file + '_finished')

    def worker_hmmer(self, queue, master_pid):
        '''
        this is a bit tricky to understand but it works the following way:
        There is a queue which is stacked with jobs (i.e. hmmer commands to execute), each process/worker takes one job and runs it
        when the worker finishes the job it just takes another job from the queue
        if the worker received a record that is None (which is a sentinel, see processes_handler in mantis_assembler.py) then the while cycle breaks and the worker is done


        Now, onto the different types of records/jobs
        -if a record is None, the worker stops
        -if a record is not None the worker has to work (^.^):
            -if the record_type is General we just run the usual hmmer command (e.g. annotate against a Pfam chunk)
            -if the record_type is NOGT or NCBIT we just run the usual hmmer command (e.g. annotate against a Pfam chunk)
            -if the record_type is NOGT/G_checkpoint or NCBIT/G_checkpoint then we dont execute hmmer, we do something else:
                 you might have noticed that when we run with TSHMMs we dont run the whole sample against the whole TSHMM, instead,
                 in each iteration we need to check which sequences still need to be annotated, so thats why we have the functions
                 that reads and generates fastas, as well as the read_domtblout - these are basically processing results during the
                 queue (unlike the others where results are processed in the end).
                 Another important thing about these checkpoints is that they also check whether a certain sample has finished annotated
                 with the taxon, you might be wondering, "is this just a normal execution?", and the answer is no.
                 Unfortunately some taxons have too many HMMs to run them all at once, so we had them split into chunks, which then forces
                 us to make sure all the HMM chunks have finished (so we can generate the fasta for the next taxon). And so with these checkpoints,
                 if the hmmer execution finished running we proceed to the next taxon, otherwise we add another checkpoint to the queue (see add_to_queue_lineage_annotation)....
                 these checkpoints will keep being added and processed by the worker until all the other workers have finished annotating all the required HMM chunks!


        '''
        while True:
            record = queue.pop(0)
            if record is None: break
            record_type = record[0]
            if record_type == 'General':
                _, hmmer_command, stdout_path = record
                self.run_hmmer_annotation(hmmer_command=hmmer_command,
                                          stdout_path=stdout_path, master_pid=master_pid)
            # for taxon runs
            elif record_type in ['NOGT', 'NCBIT']:
                _, hmmer_command, stdout_path, output_file = record
                self.run_hmmer_annotation(hmmer_command=hmmer_command,
                                          stdout_path=stdout_path, master_pid=master_pid, output_file=output_file)
            elif record_type in ['NOGG_checkpoint', 'NCBIG_checkpoint']:
                _, current_chunk_dir, fasta_path, count_seqs_original_file, chunks_n = record
                target_hmm = record_type.split('_')[0]
                target_db = target_hmm[0:-1]
                if self.taxon_annotation_finished(target_hmm, current_chunk_dir, chunks_n):
                    protein_sequences = self.read_protein_fasta(fasta_path)
                    count_seqs_chunk = len(protein_sequences)
                    self.remove_temp_fasta(fasta_path, target_db)
                else:
                    self.queue.insert(0, record)

            elif record_type in ['NOGT_checkpoint', 'NCBIT_checkpoint']:
                _, current_chunk_dir, fasta_path, taxon_id, organism_lineage, original_lineage, count_seqs_original_file, chunks_n, stdout_path = record
                target_hmm = record_type.split('_')[0]
                target_db = target_hmm[0:-1]

                if self.taxon_annotation_finished(target_hmm + str(taxon_id), current_chunk_dir, chunks_n):
                    protein_sequences = self.read_protein_fasta(fasta_path)
                    count_seqs_chunk = len(protein_sequences)
                    domtblout_path = add_slash(current_chunk_dir + 'domtblout')
                    taxon_domtblouts = self.get_taxon_chunks(taxon_id, domtblout_path, target_hmm)
                    annotated_queries = set()
                    for dp in taxon_domtblouts:
                        _, hit_counter = self.read_domtblout(output_path=domtblout_path + dp,
                                                             count_seqs_chunk=count_seqs_chunk,
                                                             count_seqs_original_file=count_seqs_original_file,
                                                             get_output=False)
                        annotated_queries.update(hit_counter)
                    # in each iteration we only annotate the missing sequences
                    self.remove_annotated_queries(protein_sequences, annotated_queries)
                    if protein_sequences:
                        fasta_path = self.generate_temp_fasta(protein_sequences, current_chunk_dir, target_db)
                        self.add_to_queue_lineage_annotation(fasta_path=fasta_path, current_chunk_dir=current_chunk_dir,
                                                             organism_lineage=list(organism_lineage),
                                                             original_lineage=list(original_lineage),
                                                             count_seqs_original_file=count_seqs_original_file,
                                                             stdout_path=stdout_path,
                                                             dbs_to_use=[target_db])
                # if annotations havent finished, we add the checker back into the queue
                else:
                    self.queue.insert(0, record)

    def prepare_queue_hmmer(self):
        hmms_list = self.compile_hmms_list()
        chunked_hmms_list = []
        for hmm_path in hmms_list:
            chunked_hmms_list.extend(compile_hmm_chunks_path(hmm_path))
        chunked_hmms_list, chunked_hmms_list_size = self.order_by_size_descending(chunked_hmms_list)
        # this will build the hmmer processes to run as well as give the domtblout we want
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            self.create_chunk_hmmer_dirs(current_chunk_dir)
            self.add_to_queue_lineage_annotation(fasta_path=chunk_path, current_chunk_dir=current_chunk_dir,
                                                 organism_lineage=list(organism_lineage),
                                                 original_lineage=list(organism_lineage),
                                                 count_seqs_original_file=count_seqs_original_file,
                                                 stdout_path=output_path + 'Mantis.out')
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            for hmm_path in chunked_hmms_list:
                # full hmmer command to be run with subprocess
                command, output_file = self.compile_annotation_job(hmm_path, target_path=chunk_path,
                                                                                   output_folder=current_chunk_dir)
                # adding our hmmer command to be consumed by the hmmer processes later on
                self.queue.append(['General', command, output_path + 'Mantis.out'])

    ####For the lineage HMMs

    def add_to_queue_lineage_annotation(self, fasta_path, current_chunk_dir, organism_lineage, original_lineage,
                                        count_seqs_original_file, stdout_path, dbs_to_use=['NCBI', 'NOG']):
        protein_sequences = self.read_protein_fasta(fasta_path)
        count_seqs_chunk = len(protein_sequences)
        # taxa specific dbs include nog and ncbi
        for db in dbs_to_use:
            db_tax = db + 'T'
            db_general = db + 'G'
            path_tax, path_general, hmm_db_general = db, db, db_general
            current_lineage = list(organism_lineage)
            if count_seqs_chunk:
                # we use this one with NCBI and NOG if there's a lineage
                if current_lineage and original_lineage:
                    if self.mantis_paths[path_tax][0:2] != 'NA':
                        current_taxon = current_lineage.pop(-1)
                        hmm_path = self.get_lineage_hmm_path(current_taxon, db=path_tax)
                        # to skip taxons without an hmm
                        while not hmm_path and current_lineage:
                            current_taxon = current_lineage.pop(-1)
                            hmm_path = self.get_lineage_hmm_path(current_taxon, db=path_tax)
                        if hmm_path:
                            chunks_path = compile_hmm_chunks_path(hmm_path)
                            for chunk_hmm in chunks_path:
                                command, output_file = self.compile_annotation_job(chunk_hmm,
                                                                                   target_path=fasta_path,
                                                                                   output_initials=db_tax,
                                                                                   output_folder=current_chunk_dir)
                                self.queue.insert(0, [db_tax, command, stdout_path, output_file])
                            # will be used for checking whether chunks have been annotated
                            self.queue.insert(len(chunks_path),
                                              [db_tax + '_checkpoint', current_chunk_dir, fasta_path, current_taxon,
                                               current_lineage, original_lineage, count_seqs_original_file,
                                               len(chunks_path), stdout_path])
                            self.save_temp_fasta_length(current_chunk_dir, db_tax + str(current_taxon),
                                                        count_seqs_chunk, db)
                        # we only use the general NCBI, NOGG is too big and would slow down annotation considerably
                        elif not hmm_path and db == 'NCBI':
                            # if there are still missing annotations from the lineage annotation or there's no taxonomic classification we query against the whole TSHMM database
                            if self.mantis_paths[path_general][0:2] != 'NA':
                                hmm_path = get_hmm_in_folder(self.mantis_paths[path_general] + hmm_db_general)
                                chunks_path = compile_hmm_chunks_path(hmm_path)
                                for chunk_hmm in chunks_path:
                                    command, output_file = self.compile_annotation_job(chunk_hmm,
                                                                                                       target_path=fasta_path,
                                                                                                       output_folder=current_chunk_dir)
                                    self.queue.insert(0, [db_tax, command, stdout_path, output_file])
                                self.queue.insert(len(chunks_path),
                                                  [db_general + '_checkpoint', current_chunk_dir, fasta_path,
                                                   count_seqs_original_file, len(chunks_path)])
                                self.save_temp_fasta_length(current_chunk_dir, db_general, count_seqs_chunk, db)
                # we use this for any situation for NCBI, for NOG only when there's no original_lineage
                else:
                    if (db == 'NOG' and not original_lineage) or (db == 'NCBI'):
                        # if there are still missing annotations from the lineage annotation or there's not taxonomic classification we query against the whole TSHMM database
                        if self.mantis_paths[path_general][0:2] != 'NA':
                            hmm_path = get_hmm_in_folder(self.mantis_paths[path_general] + hmm_db_general)
                            chunks_path = compile_hmm_chunks_path(hmm_path)
                            for chunk_hmm in chunks_path:
                                command, output_file = self.compile_annotation_job(chunk_hmm,
                                                                                                   target_path=fasta_path,
                                                                                                   output_folder=current_chunk_dir)
                                self.queue.insert(0, [db_tax, command, stdout_path, output_file])
                            self.queue.insert(len(chunks_path),
                                              [db_general + '_checkpoint', current_chunk_dir, fasta_path,
                                               count_seqs_original_file, len(chunks_path)])
                            self.save_temp_fasta_length(current_chunk_dir, db_general, count_seqs_chunk, db)

    ####Merging hmmer output

    def estimate_domtblouts_per_chunk(self):
        n_hmms = len(self.compile_hmms_list())
        if self.mantis_paths['NOG'][0:2] != 'NA': n_hmms += 2
        if self.mantis_paths['NCBI'][0:2] != 'NA': n_hmms += 2
        return n_hmms

    def process_output(self):
        domtblout_per_chunks = self.estimate_domtblouts_per_chunk()
        worker_count = estimate_number_workers_process_output(n_chunks=len(self.chunks_to_annotate),
                                                              domtblout_per_chunks=domtblout_per_chunks)
        print_cyan('Processing output with ' + str(worker_count) + ' workers.', flush=True, file=self.redirect_verbose)
        if worker_count < 1:
            print('Invalid number of workers in process_output')
            raise BadNumberWorkers

        # when we have a lot of hits in the HMM file the hit processing can be quite memory heavy, so instead we now split hits into chunks
        # This process is quite light since it only stores the file the hit should be stored at, all the hit information is read and discarded from memory
        # this also allows for parallelization
        self.prepare_queue_split_hits(worker_count)
        self.processes_handler(self.worker_split_hits, worker_count)

        self.prepare_queue_process_output()
        self.processes_handler(self.worker_process_output, worker_count)


        self.prepare_queue_merge_output()
        self.processes_handler(self.worker_merge_output, worker_count)

        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            self.remove_temp_fasta_length(current_chunk_dir, 'NOG')
            self.remove_temp_fasta_length(current_chunk_dir, 'NCBI')

    def prepare_queue_merge_processed_output(self):
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            chunks_output_path = add_slash(current_chunk_dir + 'processed_output')
            all_output_with_chunks = os.listdir(chunks_output_path)
            self.queue.append([chunks_output_path, all_output_with_chunks, output_path + 'Mantis.out'])

    def worker_merge_processed_output(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            chunks_output_path, all_output_with_chunks, stdout_path = record
            stdout_file = open(stdout_path, 'a+')
            self.merge_output_chunks(chunks_output_path, all_output_with_chunks, chunk_suffix='.pro',stdout_file=stdout_file)
            stdout_file.close()




    def prepare_queue_split_hits(self, worker_count):
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            domtblout_path = add_slash(current_chunk_dir + 'domtblout')
            raw_domtblout_path = add_slash(current_chunk_dir + 'split_hits')
            move_file(domtblout_path, raw_domtblout_path)
            all_domtblout = os.listdir(raw_domtblout_path)
            same_db_chunks = {}
            for domtblout in all_domtblout:
                chunk_n = re.search('_chunk_\d+', domtblout)
                if chunk_n: temp = domtblout.replace(chunk_n.group(), '')
                else:       temp = domtblout
                if temp not in same_db_chunks: same_db_chunks[temp] = set()
                same_db_chunks[temp].add(domtblout)
            for domtblout in same_db_chunks:
                self.queue.append(
                    [raw_domtblout_path, domtblout, same_db_chunks[domtblout], current_chunk_dir, worker_count,output_path + 'Mantis.out'])

    def worker_split_hits(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            domtblout_path, domtblout, chunks_domtblout, current_chunk_dir, worker_count, stdout_path = record
            self.split_hits(domtblout_path, domtblout, chunks_domtblout, current_chunk_dir, worker_count)


    def prepare_queue_process_output(self):
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            domtblout_path = add_slash(current_chunk_dir + 'domtblout')
            all_domtblout = os.listdir(domtblout_path)
            if not os.path.exists(add_slash(add_slash(current_chunk_dir) + 'processed_output')):
                Path(add_slash(add_slash(current_chunk_dir) + 'processed_output')).mkdir(parents=True, exist_ok=True)
            for domtblout in all_domtblout:
                if 'NOG' in domtblout:
                    count_seqs_chunk_domtblout = self.get_temp_fasta_length(current_chunk_dir, domtblout, 'NOG')
                elif 'NCBIT' in domtblout or 'NCBIG' in domtblout:
                    count_seqs_chunk_domtblout = self.get_temp_fasta_length(current_chunk_dir, domtblout, 'NCBI')
                else:
                    count_seqs_chunk_domtblout = int(count_seqs_chunk)
                self.queue.append([domtblout_path + domtblout, current_chunk_dir, count_seqs_chunk_domtblout,
                                   count_seqs_original_file, output_path + 'Mantis.out'])

    def worker_process_output(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            domtblout_path, current_chunk_dir, count_seqs_chunk, count_seqs_original_file, stdout_path = record
            processed_hits = self.process_domtblout(output_path=domtblout_path, count_seqs_chunk=count_seqs_chunk,
                                                    count_seqs_original_file=count_seqs_original_file,
                                                    stdout_path=stdout_path)
            self.save_processed_hits(processed_hits, add_slash(add_slash(current_chunk_dir) + 'processed_output'),
                                     domtblout=get_path_level(domtblout_path))
            if not self.keep_files:
                os.remove(domtblout_path)

    def prepare_queue_merge_output(self):
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            self.queue.append([current_chunk_dir, output_path])

    def worker_merge_output(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            current_chunk_dir, output_path = record
            chunks_path = add_slash(current_chunk_dir + 'processed_output')
            chunks_to_merge = [chunks_path + i for i in os.listdir(chunks_path)]
            stdout_file = open(output_path + 'Mantis.out', 'a+')
            self.merge_target_output('output_annotation.tsv', current_chunk_dir, chunks_to_merge, stdout_file,
                                     same_output=False)
            stdout_file.close()

    ###Interpreting output

    def interpret_output(self):
        worker_count = estimate_number_workers_process_output(n_chunks=len(self.chunks_to_annotate))
        self.prepare_queue_interpret_output()
        print_cyan('Interpreting output with ' + str(worker_count) + ' workers.', flush=True,
                   file=self.redirect_verbose)
        if worker_count < 1:
            print('Invalid number of workers in interpret_output')
            raise BadNumberWorkers

        self.processes_handler(self.worker_interpret_output, worker_count)

    def prepare_queue_interpret_output(self):
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            output_annotation_tsv = current_chunk_dir + 'output_annotation.tsv'
            self.queue.append([output_annotation_tsv, current_chunk_dir])

    def worker_interpret_output(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            output_annotation_tsv, current_chunk_dir = record
            interpreted_annotation_tsv = current_chunk_dir + 'integrated_annotation.tsv'
            self.generate_interpreted_output(output_annotation_tsv, interpreted_annotation_tsv)

    ###Generate consensus output

    def get_consensus_output(self):
        MANTIS_Consensus.__init__(self)
        worker_count = estimate_number_workers_process_output(n_chunks=len(self.chunks_to_annotate))
        self.prepare_queue_generate_consensus()
        print_cyan('Generating consensus output with ' + str(worker_count) + ' workers.', flush=True,
                   file=self.redirect_verbose)
        if worker_count < 1:
            print('Invalid number of workers in get_consensus_output')
            raise BadNumberWorkers

        self.processes_handler(self.worker_consensus_output, worker_count)

    def prepare_queue_generate_consensus(self):
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            interepreted_annotation_tsv = current_chunk_dir + 'integrated_annotation.tsv'
            stdout_file_path = output_path + 'Mantis.out'
            self.queue.append([interepreted_annotation_tsv, current_chunk_dir, stdout_file_path])

    def worker_consensus_output(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            interpreted_annotation_tsv, current_chunk_dir, stdout_file_path = record
            consensus_annotation_tsv = current_chunk_dir + 'consensus_annotation.tsv'
            self.generate_consensus_output(interpreted_annotation_tsv, consensus_annotation_tsv, stdout_file_path)

    # Merging Mantis output

    def merge_mantis_output(self):
        worker_count = estimate_number_workers_process_output(n_chunks=len(self.chunks_to_annotate))
        self.prepare_queue_merge_mantis_output()
        print_cyan('Merging output with ' + str(worker_count) + ' workers.', flush=True, file=self.redirect_verbose)
        if worker_count < 1:
            print('Invalid number of workers in merge_mantis_output')
            raise BadNumberWorkers

        self.processes_handler(self.worker_merge_mantis_output, worker_count)

    def prepare_queue_merge_mantis_output(self):
        for output_path in self.chunks_to_fasta:
            self.queue.append([output_path, self.chunks_to_fasta[output_path]])

    def worker_merge_mantis_output(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            output_path, chunks_path = record
            chunks_output = [i + 'output_annotation.tsv' for i in chunks_path]
            self.merge_chunks_outputs(output_path, chunks_path)

    def merge_chunks_outputs(self, output_path, chunks_path):
        stdout_file = open(output_path + 'Mantis.out', 'a+')
        self.merge_target_output(output_file='output_annotation.tsv', output_folder=output_path,
                                 chunks_path=chunks_path, stdout_file=stdout_file)
        self.merge_target_output(output_file='integrated_annotation.tsv', output_folder=output_path,
                                 chunks_path=chunks_path, stdout_file=stdout_file)
        if not self.skip_consensus:
            self.merge_target_output(output_file='consensus_annotation.tsv', output_folder=output_path,
                                     chunks_path=chunks_path, stdout_file=stdout_file)
        print('------------------------------------------', flush=True, file=stdout_file)
        print_cyan('This sample has been sucessfully annotated!', flush=True, file=stdout_file)
        stdout_file.close()
