try:
    from mantis.multiprocessing import *
except:
    from multiprocessing import *


def print_citation_mantis():
    paper_doi = 'https://doi.org/10.1093/gigascience/giab042'
    separator = '##########################################################################################################################'
    res = f'{separator}\n# Thank you for using Mantis, please make sure you cite the respective paper {paper_doi} #\n{separator}'
    print(res)


def print_version(user, project):
    import requests
    response = requests.get(f"https://api.github.com/repos/{user}/{project}/releases/latest")
    json = response.json()
    if 'name' in json:
        print(f'{project}\'s latest release is:', json['name'])
    else:
        print('No release available')


def run_mantis(input_path,
               output_folder,
               mantis_config=None,
               evalue_threshold=None,
               overlap_value=None,
               minimum_consensus_overlap=None,
               organism_details=None,
               genetic_code=None,
               domain_algorithm=None,
               best_combo_formula=None,
               sorting_type=None,
               keep_files=False,
               skip_consensus=False,
               skip_managed_memory=False,
               force_evalue=False,
               no_consensus_expansion=False,
               no_taxonomy=False,
               no_unifunc=False,
               kegg_matrix=False,
               verbose_kegg_matrix=False,
               output_gff=False,
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
    if genetic_code:                        genetic_code = int(genetic_code)
    mantis = MANTIS(
        input_path=input_path,
        output_folder=output_folder,
        mantis_config=mantis_config,
        evalue_threshold=evalue_threshold,
        overlap_value=overlap_value,
        minimum_consensus_overlap=minimum_consensus_overlap,
        organism_details=organism_details,
        genetic_code=genetic_code,
        domain_algorithm=domain_algorithm,
        best_combo_formula=best_combo_formula,
        sorting_type=sorting_type,
        keep_files=keep_files,
        skip_consensus=skip_consensus,
        skip_managed_memory=skip_managed_memory,
        force_evalue=force_evalue,
        no_consensus_expansion=no_consensus_expansion,
        no_taxonomy=no_taxonomy,
        no_unifunc=no_unifunc,
        kegg_matrix=kegg_matrix,
        verbose_kegg_matrix=verbose_kegg_matrix,
        output_gff=output_gff,
        verbose=verbose,
        default_workers=default_workers,
        chunk_size=chunk_size,
        time_limit=time_limit,
        hmmer_threads=hmmer_threads,
        user_cores=cores,
        user_memory=memory,
    )
    mantis.run_mantis()


def run_mantis_test(input_path,
                    output_folder,
                    mantis_config,
                    ):
    mantis = MANTIS(
        input_path=input_path,
        output_folder=output_folder,
        mantis_config=mantis_config,
        keep_files=True)
    mantis.run_mantis_test()


class MANTIS(Multiprocessing):
    def __init__(self,
                 input_path=None,
                 output_folder=None,
                 mantis_config=None,
                 evalue_threshold=None,
                 overlap_value=None,
                 minimum_consensus_overlap=None,
                 domain_algorithm=None,
                 sorting_type=None,
                 best_combo_formula=None,
                 organism_details={},
                 genetic_code=None,
                 redirect_verbose=None,
                 keep_files=False,
                 skip_consensus=False,
                 skip_managed_memory=False,
                 force_evalue=False,
                 no_consensus_expansion=False,
                 no_taxonomy=False,
                 no_unifunc=False,
                 kegg_matrix=False,
                 verbose_kegg_matrix=False,
                 output_gff=False,
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
        self.input_path = input_path
        self.mantis_config = mantis_config

        # Prediction parameters
        self.evalue_threshold = evalue_threshold
        self.default_evalue_threshold = 1e-3  # 1e-6 might be better a better default
        self.minimum_evalue_threshold = 1e-2
        self.force_evalue = force_evalue

        if overlap_value:
            self.overlap_value = overlap_value
        else:
            self.overlap_value = 0.1

        if minimum_consensus_overlap:
            self.minimum_consensus_overlap = minimum_consensus_overlap
        else:
            self.minimum_consensus_overlap = 0.7

        if domain_algorithm:
            self.domain_algorithm = domain_algorithm
        else:
            self.domain_algorithm = 'heuristic'

        if best_combo_formula:
            self.best_combo_formula = best_combo_formula
        else:
            self.best_combo_formula = 1
        if hmmer_threads:
            self.hmmer_threads = hmmer_threads
        # 1 should be ideal if we are already using the maximum amount of cores with Mantis
        else:
            self.hmmer_threads = 1
        # the user can force the sorting type
        if sorting_type:
            self.sorting_type = sorting_type
        else:
            # but we recommend using bitscore for dfs, evalue for bpo or heuristic
            if self.domain_algorithm == 'dfs':
                self.sorting_type = 'bitscore'
            else:
                self.sorting_type = 'evalue'
        self.organism_details = organism_details
        self.genetic_code = genetic_code
        # Execution parameters
        self.skip_consensus = skip_consensus
        self.skip_managed_memory = skip_managed_memory
        self.no_consensus_expansion = no_consensus_expansion
        self.no_unifunc = no_unifunc
        self.kegg_matrix = kegg_matrix
        self.verbose_kegg_matrix = verbose_kegg_matrix
        self.output_gff = output_gff
        if self.verbose_kegg_matrix: self.kegg_matrix = True
        self.default_workers = default_workers
        self.user_memory = user_memory
        # chunk size is highly relevant in the execution time
        self.chunk_size = chunk_size
        if time_limit:
            self.time_limit = time_limit
        else:
            self.time_limit = 60
        # diamond db size for scaling. we increase the db size to avoid overly good e-values, i.e., 0 where sample scaling by multiplication wouldn't change anything
        self.diamond_db_size = 1e6
        print_cyan('Reading config file and setting up paths', flush=True, file=self.redirect_verbose)
        Assembler.__init__(self, verbose=verbose, redirect_verbose=redirect_verbose, mantis_config=mantis_config,
                           keep_files=keep_files, user_cores=user_cores, no_taxonomy=no_taxonomy)
        datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        if self.input_path:
            print_cyan(f'This MANTIS process started running at {datetime_str}', flush=True, file=self.redirect_verbose)
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
        if self.kegg_matrix and not self.verbose_kegg_matrix:
            kegg_matrix_str = 'Generate KEGG modules matrix:\t' + str(self.kegg_matrix) + '\n'
        elif self.kegg_matrix and self.verbose_kegg_matrix:
            kegg_matrix_str = 'Generate KEGG modules matrix in verbose mode:\t' + str(self.verbose_kegg_matrix) + '\n'
        else:
            kegg_matrix_str = ''

        output_list = [
            'Output folder:\t\t\t' + str(self.output_folder) + '\n' if self.output_folder else '',
            'Mantis config:\t\t\t' + str(self.mantis_config) + '\n' if self.mantis_config else '',
            'Target path:\t\t\t' + str(self.input_path) + '\n' if self.input_path else '',
            'E-value threshold:\t\t' + str(self.evalue_threshold) + '\n' if self.evalue_threshold else '',
            'E-value threshold:\t\t' + str(self.default_evalue_threshold) + '\n' if not self.evalue_threshold else '',
            'Forcing e-value:\t\t' + str(self.force_evalue) + '\n' if self.force_evalue else '',
            'Overlap value:\t\t\t' + str(self.overlap_value) + '\n' if self.overlap_value else '',
            'Default workers:\t\t' + str(self.default_workers) + '\n' if self.default_workers else '',
            'User cores:\t\t\t' + str(self.user_cores) + '\n' if self.user_cores else '',
            'HMMER threads:\t\t\t' + str(self.hmmer_threads) + '\n' if self.hmmer_threads else '',
            'Chunk size:\t\t\t' + str(self.chunk_size) + '\n' if self.chunk_size else '',
            'Algorithm:\t\t\t' + str(self.domain_algorithm) + '\n' if self.domain_algorithm else '',
            # 'Formula:\t\t\t' + str(self.best_combo_formula) + '\n' if self.best_combo_formula else '',
            # 'Sorting type:\t\t\t' + str(self.sorting_type) + '\n' if self.sorting_type else '',
            'Outputting GFF:\t\t\t' + str(self.output_gff) + '\n' if self.output_gff else '',
            'Skip consensus:\t\t' + str(self.skip_consensus) + '\n' if self.skip_consensus else '',
            'Skip memory management:\t\t' + str(self.skip_managed_memory) + '\n' if self.skip_managed_memory else '',
            'Skip consensus expansion:\t' + str(
                self.no_consensus_expansion) + '\n' if self.no_consensus_expansion else '',
            'Use taxonomy:\t' + str(self.use_taxonomy) + '\n' if not self.use_taxonomy else '',
            'Skip text similarity analysis:\t' + str(self.no_unifunc) + '\n' if self.no_unifunc else '',
            kegg_matrix_str,
            '------------------------------------------']
        return 'User configuration:' + '\n' + '------------------------------------------' + '\n' + ''.join(output_list)

    def generate_fastas_to_annotate(self):
        if '.' not in self.input_path:
            kill_switch(InvalidTargetFile,
                        'Your file does not have an extension, so Mantis can\'t detect the file format.\nPlease provide a valid target file',
                        flush=True, file=self.redirect_verbose)
        else:
            if os.path.isdir(self.input_path):
                self.annotate_directory()
            elif self.input_path.endswith('.tsv'):
                self.annotate_multiple_samples()
            elif self.input_path.split('.')[-2] in ['tar']:
                self.annotate_compressed_sample()
            elif self.input_path.endswith('.gz') or self.input_path.endswith('.zip'):
                self.annotate_compressed_sample()
            elif is_fasta(self.input_path):
                self.annotate_one_sample()
            else:
                kill_switch(InvalidTargetFile,
                            'Your file does not appear to be a fasta. If you want to annotate multiple samples, make sure your file has the <.tsv> extension.',
                            flush=True, file=self.redirect_verbose)
        if not self.fastas_to_annotate: kill_switch(NoValidFiles)
        for file_path, output_path, organism_details, genetic_code, count_seqs_original_file, count_residues_original_file in self.fastas_to_annotate:
            Path(output_path).mkdir(parents=True, exist_ok=True)

    def annotate_multiple_samples(self):
        try:
            with open(self.input_path) as file:
                line = file.readline()
                if SPLITTER not in line:
                    line = file.readline()
                while line:
                    line = line.strip('\n').split('\t')
                    if len(line) >= 2:
                        query_name = line[0]
                        line_path = line[1]
                        if len(line) == 4:
                            genetic_code = ' '.join(line[-1])
                            organism_details = ' '.join(line[2:-1])
                            if organism_details == 'None': organism_details = ''
                        elif len(line) == 3:
                            organism_details = ' '.join(line[2:])
                            if organism_details == 'None': organism_details = ''
                            genetic_code = None
                        else:
                            organism_details = None
                            genetic_code = None
                        count_seqs_original_file = get_seqs_count(line_path)
                        count_residues_original_file = count_residues(line_path)
                        if os.path.exists(line_path):
                            self.fastas_to_annotate.append([line_path, add_slash(self.output_folder + query_name),
                                                            organism_details, genetic_code,
                                                            count_seqs_original_file, count_residues_original_file])
                        else:
                            kill_switch(TargetFileNotFound,
                                        flush=True, file=self.redirect_verbose)
                    line = file.readline()
        except:
            kill_switch(InvalidTargetFile,
                        'If you want to annotate multiple samples, make sure your file is correctly formatted. Please see the examples in the <tests> folder.',
                        flush=True, file=self.redirect_verbose)

    def annotate_directory(self):
        try:
            list_dir = os.listdir(self.input_path)
            for file in list_dir:
                if 'faa' in file.split('.')[-1]:
                    query_name = '.'.join(file.split('.')[0:-1])
                    query_path = self.input_path + file
                    count_seqs_original_file = get_seqs_count(query_path)
                    count_residues_original_file = count_residues(query_path)
                    self.fastas_to_annotate.append([query_path, add_slash(self.output_folder + query_name), None, None,
                                                    count_seqs_original_file, count_residues_original_file])
        except:
            kill_switch(InvalidTargetFile, 'Something went wrong when annotating the provided directory!', flush=True,
                        file=self.redirect_verbose)

    def annotate_compressed_sample(self):
        try:
            uncompressed_path = self.output_folder + 'uncompressed_samples/'
            Path(uncompressed_path).mkdir(parents=True, exist_ok=True)
            uncompressing_function = uncompress_archive(source_filepath=self.input_path,
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
                            count_residues_original_file = count_residues(query_path)
                            self.fastas_to_annotate.append(
                                [query_path, add_slash(self.output_folder + query_name), None, None,
                                 count_seqs_original_file, count_residues_original_file])
                if 'faa' in file.split('.')[-1]:
                    query_name = '.'.join(file.split('.')[0:-1])
                    query_path = uncompressed_path + file
                    count_seqs_original_file = get_seqs_count(query_path)
                    count_residues_original_file = count_residues(query_path)
                    self.fastas_to_annotate.append(
                        [query_path, add_slash(self.output_folder + query_name), None, None,
                         count_seqs_original_file, count_residues_original_file])
        except:
            kill_switch(InvalidTargetFile, 'Something went wrong when annotating the provided compressed file!',
                        flush=True, file=self.redirect_verbose)

    def annotate_one_sample(self):
        count_seqs_original_file = get_seqs_count(self.input_path)
        count_residues_original_file = count_residues(self.input_path)
        self.fastas_to_annotate.append(
            [self.input_path, self.output_folder, self.organism_details, self.genetic_code,
             count_seqs_original_file, count_residues_original_file])

    def setup_organism_lineage(self, organism_details, stdout_file):
        if not self.use_taxonomy: return []
        if not organism_details:
            print_cyan('No data provided for organism lineage!', flush=True, file=stdout_file)
            return []
        if re.match('\d+', organism_details):
            print_cyan('Setting up organism lineage from provided NCBI taxon id', flush=True, file=stdout_file)
            organism_lineage = self.fetch_ncbi_lineage(organism_details)
        else:
            print_cyan('Setting up organism lineage from provided taxon synonym or GTDB lineage', flush=True,
                       file=stdout_file)
            organism_details_dict = {'synonyms': organism_details}
            ncbi_taxon_id = self.get_taxa_ncbi(organism_details)
            organism_lineage = self.fetch_ncbi_lineage(ncbi_taxon_id)
        return organism_lineage

    def generate_translated_sample(self):
        ncbi_resources = add_slash(self.mantis_paths['resources'] + 'NCBI')
        translation_tables = parse_translation_tables(ncbi_resources + 'gc.prt')
        for i in range(len(self.fastas_to_annotate)):
            file_path, output_path, organism_details, genetic_code, count_seqs_original_file, count_residues_original_file = \
            self.fastas_to_annotate[i]
            sample_type = check_sample_type(file_path)
            if sample_type == 'dna' or sample_type == 'rna':
                if not genetic_code:
                    genetic_code = 11
                translated_fasta_path = f'{output_path}translated_gc_{genetic_code}.faa'
                try:
                    write_translated_fasta(original_fasta_path=file_path, translated_fasta_path=translated_fasta_path,
                                           translation_table=translation_tables[genetic_code], sample_type=sample_type)
                    self.fastas_to_annotate[i][0] = translated_fasta_path
                    self.fastas_to_annotate[i][5] = count_residues(translated_fasta_path)
                except Exception as e:
                    kill_switch(InvalidTranslation, file_path)

    def generate_sample_lineage(self):
        self.start_taxonomy_connection()
        for i in range(len(self.fastas_to_annotate)):
            file_path, output_path, organism_details, genetic_code, count_seqs_original_file, count_residues_original_file = self.fastas_to_annotate[i]
            stdout_file = open(output_path + 'Mantis.out', 'a+')
            organism_lineage = self.setup_organism_lineage(organism_details, stdout_file)
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
        for file_path, output_path, organism_details, genetic_code, count_seqs_original_file, count_residues_original_file in self.fastas_to_annotate:
            if file_exists(output_path + 'fasta_chunks/'):
                shutil.rmtree(output_path + 'fasta_chunks/')

    def remove_uncompressed_files(self):
        if file_exists(self.output_folder + 'uncompressed_samples/'):
            shutil.rmtree(self.output_folder + 'uncompressed_samples/')

    @timeit_class
    def run_mantis_test(self):
        self.mantis_paths['custom'] = MANTIS_FOLDER + 'tests/test_hmm/'
        self.output_gff = True
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        self.generate_fastas_to_annotate()
        self.generate_translated_sample()
        self.generate_sample_lineage()

        self.split_sample()

        self.set_chunks_to_annotate()

        worker_count = 1
        for hmm_path in [MANTIS_FOLDER + 'tests/test_hmm/test1/test1.hmm',
                         MANTIS_FOLDER + 'tests/test_hmm/test2/test2.hmm',
                         MANTIS_FOLDER + 'tests/test_hmm/test3/test3.hmm',
                         ]:
            for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, count_residues_original_file, output_path in self.chunks_to_annotate:
                self.create_chunk_output_dirs(current_chunk_dir)
                command, output_file = self.compile_annotation_job(hmm_path,
                                                                   target_path=chunk_path,
                                                                   output_folder=current_chunk_dir,
                                                                   count_seqs_chunk=count_seqs_chunk,
                                                                   count_seqs_original_file=count_seqs_original_file)
                self.queue.append(['General', command, output_path + 'Mantis.out'])

        self.processes_handler(self.worker_annotation, worker_count)

        self.prepare_queue_split_hits(worker_count)
        self.processes_handler(self.worker_split_hits, worker_count)

        self.prepare_queue_process_output()
        self.processes_handler(self.worker_process_output, worker_count)

        self.prepare_queue_merge_output()
        self.processes_handler(self.worker_merge_output, worker_count)

        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, count_residues_original_file, output_path in self.chunks_to_annotate:
            output_annotation_tsv = current_chunk_dir + 'output_annotation.tsv'
            self.queue.append([output_annotation_tsv, current_chunk_dir])
        self.processes_handler(self.worker_interpret_output, worker_count)
        for chunk_name, chunk_path, current_chunk_dir, organism_lineage, count_seqs_chunk, count_seqs_original_file, count_residues_original_file, output_path in self.chunks_to_annotate:
            interepreted_annotation_tsv = current_chunk_dir + 'integrated_annotation.tsv'
            stdout_file_path = output_path + 'Mantis.out'
            self.queue.append([interepreted_annotation_tsv, current_chunk_dir, stdout_file_path])
        Consensus.__init__(self)
        self.processes_handler(self.worker_consensus_output, worker_count)

        for output_path in self.chunks_to_fasta:
            self.queue.append([output_path, self.chunks_to_fasta[output_path]])
        self.processes_handler(self.worker_merge_mantis_output, worker_count)
        self.generate_matrix()

        print('Mantis ran sucessfully!', flush=True)

    @timeit_class
    def run_mantis(self):
        self.check_installation(verbose=True)
        if not self.passed_check:
            kill_switch(InstallationCheckNotPassed)

        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        print_cyan('Reading input file', flush=True, file=self.redirect_verbose)
        self.generate_fastas_to_annotate()
        print_cyan('Determining lineage', flush=True, file=self.redirect_verbose)
        self.generate_translated_sample()
        self.generate_sample_lineage()
        if self.use_taxonomy:
            try:
                self.close_taxonomy_connection()
            except:
                pass
        self.split_sample()
        self.set_chunks_to_annotate()
        start_time = time()
        self.run_homology_search()
        print(f'Homology search took {int(time() - start_time)} seconds to run', flush=True, file=self.redirect_verbose)
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
            print(f'Consensus generation took {int(time() - start_time)} seconds to run', flush=True,
                  file=self.redirect_verbose)
        start_time = time()
        self.merge_mantis_output()
        print(f'Output merge took {int(time() - start_time)} seconds to run', flush=True, file=self.redirect_verbose)
        start_time = time()
        self.remove_uncompressed_files()
        print(f'In total, Mantis took {int(time() - processing_start_time)} seconds to process homology search results',
              flush=True, file=self.redirect_verbose)
        if not self.keep_files:
            print_cyan('Removing temporary files', flush=True, file=self.redirect_verbose)
            self.remove_non_essential_files()
        if self.kegg_matrix and not self.skip_consensus:
            print('Generating KEGG modules matrix!', flush=True, file=self.redirect_verbose)
            self.generate_matrix()
        print(f'Mantis has finished running! It took {str(int(time() - self.start_time))} seconds to run.', flush=True,
              file=self.redirect_verbose)


if __name__ == '__main__':
    m = MANTIS()
