try:
    from source.MANTIS_MP import *
except:
    from MANTIS_MP import *


def run_mantis(target_path,
               output_folder,
               mantis_config=None,
               evalue_threshold=None,
               overlap_value=None,
               minimum_consensus_overlap=None,
               organism_details=None,
               domain_algorithm=None,
               best_combo_formula=None,
               sorting_type=None,
               keep_files=False,
               skip_consensus=False,
               skip_managed_memory=False,
               no_consensus_expansion=False,
               no_unifunc=False,
               kegg_matrix=False,
               verbose=True,
               default_workers=None,
               chunk_size=None,
               time_limit=None,
               hmmer_threads=None,
               cores=None,
               memory=None,
               ):
    if evalue_threshold:
        if evalue_threshold != 'dynamic':   evalue_threshold = float(evalue_threshold)
    if overlap_value:                       overlap_value = float(overlap_value)
    if minimum_consensus_overlap:           minimum_consensus_overlap = float(minimum_consensus_overlap)
    if best_combo_formula:                  best_combo_formula = int(best_combo_formula)
    if default_workers:                     default_workers = int(default_workers)
    if chunk_size:                          chunk_size = int(chunk_size)
    if time_limit:                          time_limit = int(time_limit)
    if hmmer_threads:                       hmmer_threads = int(hmmer_threads)
    if cores:                               cores = int(cores)
    if memory:                              memory = int(memory)
    mantis = MANTIS(
        target_path=target_path,
        output_folder=output_folder,
        mantis_config=mantis_config,
        evalue_threshold=evalue_threshold,
        overlap_value=overlap_value,
        minimum_consensus_overlap=minimum_consensus_overlap,
        organism_details=organism_details,
        domain_algorithm=domain_algorithm,
        best_combo_formula=best_combo_formula,
        sorting_type=sorting_type,
        keep_files=keep_files,
        skip_consensus=skip_consensus,
        skip_managed_memory=skip_managed_memory,
        no_consensus_expansion=no_consensus_expansion,
        no_unifunc=no_unifunc,
        kegg_matrix=kegg_matrix,
        verbose=verbose,
        default_workers=default_workers,
        chunk_size=chunk_size,
        time_limit=time_limit,
        hmmer_threads=hmmer_threads,
        user_cores=cores,
        user_memory=memory,
    )
    mantis.run_mantis()


def run_mantis_test(target_path,
                    output_folder,
                    mantis_config,
                    ):
    mantis = MANTIS(
        target_path=target_path,
        output_folder=output_folder,
        mantis_config=mantis_config,
        keep_files=True)
    mantis.run_mantis_test()

'''
The MANTIS class is just a wrapper:
    -   Exceptions.py contains custom exceptions raised by Mantis
    -   MANTIS_Assembler.py contains methods responsible for the basic initialization of Mantis (i.e. reading config file, getting taxa lineage, and checking installation)
    -   MANTIS_Consensus.py contains methods responsible for generating the consensus_annotation.tsv
    -   MANTIS_DB.py contains methods responsible for setting up the reference databases
    -   MANTIS_Metadata.py contains methods responsible for integrating the metadata/generating the integrated_annotation.tsv
    -   MANTIS_MP.py contains most methods responsible for annotation, it wraps around other classes and uses their methods in a multi processing environment
    -   MANTIS_NLP.py is a wrapper for UniFunc, it has no methods
    -   MANTIS_Processor.py contains methods responsible for the processing of hits (e.g. hit processing algorithms)
    -   utils.py contains functions common to most classes
    -   __main__.py contains the front end for user interaction
    
This is how the classes are connected
UniFunc----------->MANTIS_NLP------>MANTIS_DB-->MANTIS_Assembler--v
                                ^-->MANTIS_Consensus-------------------->MANTIS_MP--->MANTIS--->__main__
MANTIS_Metadata---------------------------------------------^  ^
MANTIS_Processor-----------------------------------------------^

utils and Exceptions are only general funtions and therefore don't respect this hierarchy, they are imported as needed

'''
class MANTIS(MANTIS_MP):
    def __init__(self,
                 target_path=None,
                 output_folder=None,
                 mantis_config=None,
                 evalue_threshold=None,
                 overlap_value=None,
                 minimum_consensus_overlap=None,
                 domain_algorithm=None,
                 sorting_type=None,
                 best_combo_formula=None,
                 organism_details={},
                 redirect_verbose=None,
                 keep_files=False,
                 skip_consensus=False,
                 skip_managed_memory=False,
                 no_consensus_expansion=False,
                 no_unifunc=False,
                 kegg_matrix=False,
                 verbose=True,
                 default_workers=None,
                 chunk_size=None,
                 time_limit=None,
                 hmmer_threads=None,
                 user_cores=None,
                 user_memory=None,
                 ):
        self.output_folder = add_slash(output_folder)
        self.redirect_verbose = redirect_verbose
        print('------------------------------------------', flush=True, file=self.redirect_verbose)
        print_cyan('Setting up Mantis!', flush=True, file=self.redirect_verbose)
        print('------------------------------------------', flush=True, file=self.redirect_verbose)
        self.target_path = target_path
        self.mantis_config = mantis_config

        #Prediction parameters
        if evalue_threshold:            self.evalue_threshold = evalue_threshold
        else:                           self.evalue_threshold = 1e-3

        if overlap_value:               self.overlap_value = overlap_value
        else:                           self.overlap_value = 0.1

        if minimum_consensus_overlap:   self.minimum_consensus_overlap=minimum_consensus_overlap
        else:                           self.minimum_consensus_overlap = 0.7

        if domain_algorithm:            self.domain_algorithm = domain_algorithm
        else:                           self.domain_algorithm = 'dfs'

        if  best_combo_formula:         self.best_combo_formula = best_combo_formula
        else:                           self.best_combo_formula = 1
        if hmmer_threads:               self.hmmer_threads = hmmer_threads
        # 1 should be ideal if we are already using the maximum amount of cores with Mantis
        else:                           self.hmmer_threads = 1
        #the user can force the sorting type
        if sorting_type:                self.sorting_type = sorting_type
        else:
        #but we recommend using bitscore for dfs, evalue for bpo or heuristic
            if self.domain_algorithm !='dfs':   self.sorting_type='evalue'
            else:                               self.sorting_type='bitscore'
        self.organism_details = organism_details
        #Execution parameters
        self.skip_consensus = skip_consensus
        self.skip_managed_memory = skip_managed_memory
        self.no_consensus_expansion = no_consensus_expansion
        self.no_unifunc = no_unifunc
        self.kegg_matrix = kegg_matrix
        self.default_workers = default_workers
        self.user_cores = user_cores
        self.user_memory = user_memory
        # chunk size is highly relevant in the execution time
        self.chunk_size = chunk_size
        if time_limit:
            self.time_limit = time_limit
        else:
            self.time_limit = 60
        print_cyan('Reading config file and setting up paths', flush=True, file=self.redirect_verbose)
        MANTIS_Assembler.__init__(self, verbose=verbose, redirect_verbose=redirect_verbose,
                                  mantis_config=mantis_config,keep_files=keep_files)
        datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        if self.target_path:
            print_cyan(f'This MANTIS process started running at {datetime_str}',flush=True, file=self.redirect_verbose)
        self.chunks_to_annotate = []
        self.chunks_to_fasta = {}
        self.fastas_to_annotate = []
        self.print_available_hardware()

    def print_available_hardware(self):
        if self.user_cores:
            print(f'Cores allocated: {self.user_cores}')
        else:
            print(f'Cores allocated: {ENVIRONMENT_CORES}')
        if self.user_memory:
            print(f'Memory allocated: {self.user_memory}')
        else:
            print(f'Memory allocated: {round(AVAILABLE_RAM, 2)}')
        print(f'Workers per core: {WORKER_PER_CORE}')

    def __str__(self):
        output_list = [
            'Output folder:\t\t\t' + str(self.output_folder) + '\n' if self.output_folder else '',
            'Mantis config:\t\t\t' + str(self.mantis_config) + '\n' if self.mantis_config else '',
            'Target path:\t\t\t' + str(self.target_path) + '\n' if self.target_path else '',
            'E-value threshold:\t\t' + str(self.evalue_threshold) + '\n' if self.evalue_threshold else '',
            'Overlap value:\t\t\t' + str(self.overlap_value) + '\n' if self.overlap_value else '',
            'Default workers:\t\t' + str(self.default_workers) + '\n' if self.default_workers else '',
            'User cores:\t\t\t\t' + str(self.user_cores) + '\n' if self.user_cores else '',
            'HMMER threads:\t\t\t' + str(self.hmmer_threads) + '\n' if self.hmmer_threads else '',
            'Chunk size:\t\t\t' + str(self.chunk_size) + '\n' if self.chunk_size else '',
            'Algorithm:\t\t\t' + str(self.domain_algorithm) + '\n' if self.domain_algorithm else '',
            'Formula:\t\t\t' + str(self.best_combo_formula) + '\n' if self.best_combo_formula else '',
            'Sorting type:\t\t\t' + str(self.sorting_type) + '\n' if self.sorting_type else '',
            'Skip consensus:\t\t' + str(self.skip_consensus) + '\n' if self.skip_consensus else '',
            'Skip memory management:\t\t' + str(self.skip_managed_memory) + '\n' if self.skip_managed_memory else '',
            'Skip consensus expansion:\t' + str(self.no_consensus_expansion) + '\n' if self.no_consensus_expansion else '',
            'Skip text similarity analysis:\t' + str(self.no_unifunc) + '\n' if self.no_unifunc else '',
            'Generate KEGG modules matrix:\t' + str(self.kegg_matrix) + '\n' if self.kegg_matrix else '',
            '------------------------------------------']
        return 'User configuration:' + '\n' + '------------------------------------------' + '\n' + ''.join(output_list)

    def generate_fastas_to_annotate(self):
        if '.' not in self.target_path:
            print('Your file does not have an extension, so Mantis can\'t detect the file format',flush=True, file=self.redirect_verbose)
            raise InvalidTargetFile
        else:
            if os.path.isdir(self.target_path):
                self.annotate_directory()
            elif self.target_path.endswith('.tsv'):
                self.annotate_multiple_samples()
            elif self.target_path.split('.')[-2] in ['tar']:
                self.annotate_compressed_sample()
            elif self.target_path.endswith('.gz') or self.target_path.endswith('.zip'):
                self.annotate_compressed_sample()
            elif is_fasta(self.target_path):
                self.annotate_one_sample()
            else:
                print('Your file does not appear to be a fasta. If you want to annotate multiple samples, make sure your file has the <.tsv> extension.',flush=True, file=self.redirect_verbose)
                raise InvalidTargetFile
        if not self.fastas_to_annotate: raise InvalidTargetFile
        for file_path, output_path, organism_details, count_seqs_original_file in self.fastas_to_annotate:
            Path(output_path).mkdir(parents=True, exist_ok=True)

    def annotate_multiple_samples(self):
        try:
            with open(self.target_path) as file:
                line = file.readline()
                if SPLITTER not in line:
                    line = file.readline()
                while line:
                    line = line.strip('\n').split()
                    if len(line) >= 2:
                        query_name = line[0]
                        line_path = line[1]
                        if len(line) > 2:
                            organism_details = ' '.join(line[2:])
                            if organism_details=='None': organism_details=''
                        else:
                            organism_details = ''
                        count_seqs_original_file = get_seqs_count(line_path)
                        self.fastas_to_annotate.append([line_path, add_slash(self.output_folder + query_name), organism_details,count_seqs_original_file])
                    line = file.readline()
        except:
            print('If you want to annotate multiple samples, make sure your file is correctly formatted. Please see the examples in the <tests> folder.',
                flush=True, file=self.redirect_verbose)
            raise InvalidTargetFile

    def annotate_directory(self):
        try:
            list_dir = os.listdir(self.target_path)
            for file in list_dir:
                if 'faa' in file.split('.')[-1]:
                    query_name = '.'.join(file.split('.')[0:-1])
                    query_path = self.target_path + file
                    count_seqs_original_file = get_seqs_count(query_path)
                    self.fastas_to_annotate.append(
                        [query_path, add_slash(self.output_folder + query_name), None, count_seqs_original_file])
        except:
            raise InvalidTargetFile

    def annotate_compressed_sample(self):
        try:
            uncompressed_path = self.output_folder + 'uncompressed_samples/'
            Path(uncompressed_path).mkdir(parents=True, exist_ok=True)
            uncompressing_function = uncompress_archive(source_filepath=self.target_path,
                                                        extract_path=uncompressed_path)
            list_dir = os.listdir(uncompressed_path)
            for file in list_dir:
                if os.path.isdir(uncompressed_path + file):
                    sub_list_dir = os.listdir(uncompressed_path + file)
                    for sub_file in sub_list_dir:
                        if 'faa' in sub_file.split('.')[-1]:
                            query_name = '.'.join(sub_file.split('.')[0:-1])
                            query_path = add_slash(uncompressed_path + file) + sub_file
                            count_seqs_original_file = get_seqs_count(query_path)
                            self.fastas_to_annotate.append(
                                [query_path, add_slash(self.output_folder + query_name), None,
                                 count_seqs_original_file])
                if 'faa' in file.split('.')[-1]:
                    query_name = '.'.join(file.split('.')[0:-1])
                    query_path = uncompressed_path + file
                    count_seqs_original_file = get_seqs_count(query_path)
                    self.fastas_to_annotate.append(
                        [query_path, add_slash(self.output_folder + query_name), None, count_seqs_original_file])
        except:
            raise InvalidTargetFile

    def annotate_one_sample(self):
        count_seqs_original_file = get_seqs_count(self.target_path)
        self.fastas_to_annotate.append(
            [self.target_path, self.output_folder, self.organism_details, count_seqs_original_file])

    def setup_organism_lineage(self, organism_details, stdout_file):
        # when running with DRAX
        if 'organism_instance' in organism_details:
            print_cyan('Setting up organism lineage from provided organism instance', flush=True, file=stdout_file)
            organism_lineage = self.get_organism_lineage(organism_details['organism_instance'].get_detail('ncbi_taxon'),
                                                         stdout_file=stdout_file)
        elif 'organism_lineage' in organism_details:
            print_cyan('Setting up organism lineage from provided organism lineage', flush=True, file=stdout_file)
            organism_lineage = organism_details['organism_lineage']
        # when running standalone
        elif 'ncbi_taxon' in organism_details:
            print_cyan('Setting up organism lineage from provided NCBI taxon id', flush=True, file=stdout_file)
            organism_lineage = self.get_organism_lineage(organism_details['ncbi_taxon'], stdout_file=stdout_file)
        elif 'synonyms' in organism_details:
            print_cyan('Setting up organism lineage from provided taxon synonym', flush=True, file=stdout_file)
            ncbi_taxon_id = self.get_taxa_ncbi(organism_details['synonyms'])
            organism_lineage = self.get_organism_lineage(ncbi_taxon_id, stdout_file=stdout_file)
        else:
            print_cyan('No data provided for organism lineage!', flush=True, file=stdout_file)
            organism_lineage = []
        return organism_lineage

    def generate_sample_lineage(self):
        for i in range(len(self.fastas_to_annotate)):
            file_path, output_path, organism_details, count_seqs_original_file = self.fastas_to_annotate[i]
            if not organism_details:
                organism_details_dict = {}
            elif re.match('\d+', organism_details):
                organism_details_dict = {'ncbi_taxon': organism_details}
            else:
                organism_details_dict = {'synonyms': organism_details}
            stdout_file = open(output_path + 'Mantis.out', 'a+')
            organism_lineage = self.setup_organism_lineage(organism_details_dict, stdout_file)
            self.fastas_to_annotate[i][2] = organism_lineage
            print('------------------------------------------', flush=True, file=stdout_file)
            if organism_lineage:
                print_cyan('Target file:\n' + file_path + '\n has the following taxonomy lineage: ' + ' > '.join(
                    organism_lineage), flush=True, file=stdout_file)
            else:
                print_cyan('Target file:\n' + file_path + '\n has no organism lineage!', flush=True, file=stdout_file)
            print('------------------------------------------', flush=True, file=stdout_file)

            stdout_file.close()

    def remove_non_essential_files(self):
        for file_path, output_path, organism_details, count_seqs_original_file in self.fastas_to_annotate:
            if os.path.exists(output_path + 'fasta_chunks/'):
                shutil.rmtree(output_path + 'fasta_chunks/')

    def remove_uncompressed_files(self):
        if os.path.exists(self.output_folder + 'uncompressed_samples/'):
            shutil.rmtree(self.output_folder + 'uncompressed_samples/')

    @timeit_class
    def run_mantis_test(self):
        self.mantis_paths['custom']=MANTIS_FOLDER + 'tests/test_hmm/'
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        self.generate_fastas_to_annotate()
        self.generate_sample_lineage()
        self.split_sample()

        self.set_chunks_to_annotate()

        worker_count = 1
        for hmm_path in [MANTIS_FOLDER + 'tests/test_hmm/test1/test1.hmm', MANTIS_FOLDER + 'tests/test_hmm/test2/test2.hmm']:
            for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
                self.create_chunk_hmmer_dirs(current_chunk_dir)
                command, output_file = self.compile_annotation_job(hmm_path, target_path=chunk_path,output_folder=current_chunk_dir)
                self.queue.append(['General', command, output_path + 'Mantis.out'])

        self.processes_handler(self.worker_hmmer, worker_count)

        self.prepare_queue_split_hits(worker_count)
        self.processes_handler(self.worker_split_hits, worker_count)

        self.prepare_queue_process_output()
        self.processes_handler(self.worker_process_output, worker_count)

        self.prepare_queue_merge_output()
        self.processes_handler(self.worker_merge_output, worker_count)

        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            output_annotation_tsv = current_chunk_dir + 'output_annotation.tsv'
            self.queue.append([output_annotation_tsv, current_chunk_dir])
        self.processes_handler(self.worker_interpret_output, worker_count)
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, output_path in self.chunks_to_annotate:
            interepreted_annotation_tsv = current_chunk_dir + 'integrated_annotation.tsv'
            stdout_file_path = output_path + 'Mantis.out'
            self.queue.append([interepreted_annotation_tsv, current_chunk_dir, stdout_file_path])
        MANTIS_Consensus.__init__(self)
        self.processes_handler(self.worker_consensus_output, worker_count)

        for output_path in self.chunks_to_fasta:
            self.queue.append([output_path, self.chunks_to_fasta[output_path]])
        self.processes_handler(self.worker_merge_mantis_output, worker_count)
        # if os.path.exists(self.output_folder):
        #    shutil.rmtree(self.output_folder)

        print('Mantis ran sucessfully!', flush=True)

    @timeit_class
    def run_mantis(self):
        self.check_installation(verbose=False)
        if not self.passed_check:
            raise InstallationCheckNotPassed

        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        print_cyan('Reading input file', flush=True, file=self.redirect_verbose)
        self.generate_fastas_to_annotate()
        print_cyan('Determining lineage', flush=True, file=self.redirect_verbose)
        self.generate_sample_lineage()
        self.split_sample()
        self.set_chunks_to_annotate()
        start_time = time()
        self.run_hmmer()
        print(f'HMMER took {int(time() - start_time)} seconds to run', flush=True, file=self.redirect_verbose)
        processing_start_time = time()
        start_time = time()
        self.process_output()
        print(f'Output processing took {int(time() - start_time)} seconds to run', flush=True,
              file=self.redirect_verbose)
        start_time = time()

        self.interpret_output()
        print(f'Output interpretation took {int(time() - start_time)} seconds to run', flush=True,
              file=self.redirect_verbose)
        if not self.skip_consensus:
            start_time = time()
            self.get_consensus_output()
            print(f'Consensus generation took {int(time() - start_time)} seconds to run', flush=True,file=self.redirect_verbose)
        start_time = time()
        self.merge_mantis_output()
        print(f'Output merge took {int(time() - start_time)} seconds to run', flush=True, file=self.redirect_verbose)
        start_time = time()
        self.remove_uncompressed_files()
        print(f'In total, Mantis took {int(time() - processing_start_time)} seconds to process HMMER\'s output',flush=True, file=self.redirect_verbose)
        if not self.keep_files:
            print_cyan('Removing temporary files', flush=True, file=self.redirect_verbose)
            self.remove_non_essential_files()
        if self.kegg_matrix and not self.skip_consensus:
            print('Generating KEGG modules matrix!',flush=True, file=self.redirect_verbose)
            self.generate_matrix()
        print(f'Mantis has finished running! It took {str(int(time() - self.start_time))} seconds to run.',flush=True, file=self.redirect_verbose)


if __name__ == '__main__':
    m = MANTIS()
    essential_genes = m.get_essential_genes_list()
    print(essential_genes)
