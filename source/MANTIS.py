try:
    from source.MANTIS_MP import *
except:
    from MANTIS_MP import *

'''
The way this class works is:
1- If we have taxonomically classified the sample we go to step 2 else go to step 6
2- We get tax id for organism we are annotating
3- We get the lineage of the organism
4- We get the most specific taxa id of the organism
5- We find the NOG hmm corresponding to the current taxa id (per_tax_level). If its not found we go one level up (less specific) until we find an hmm.
    If an hmm is found hmmscan with it, else got to step 6. If the hmmscan has a hit below the threshold, return hit, else go to step 6
6- Get all global hmms available (NOG_hmm and others)
7- Iterate over the global hmms
8- Retrieve all hits below the threshold
9- If hit seq start--end overlap get the hit with the best evalue. Repeat until there are no overlaps in our query
'''


def run_mantis(  target_path,
                 output_folder,
                 mantis_config=None,
                 evalue_threshold=None,
                 acceptable_range=None,
                 overlap_value=None,
                 organism_details=None,
                 domain_algorithm=None,
                 keep_files=False,
                 skip_consensus=False,
                 target_hmm=None,
                 verbose=True,
                 default_workers=None,
                 chunk_size=None,
                 hmmer_threads=None,
                 cores=None,
                 memory=None,
                 ):
    if not acceptable_range:        acceptable_range=0.05
    if not overlap_value:           overlap_value=0.1
    if not domain_algorithm:        domain_algorithm='dfs'
    if evalue_threshold:    evalue_threshold=float(evalue_threshold)
    if overlap_value:       overlap_value=float(overlap_value)
    if default_workers:     default_workers=int(default_workers)
    if chunk_size:          chunk_size=int(chunk_size)
    if hmmer_threads:       hmmer_threads=int(hmmer_threads)
    if cores:               cores=int(cores)
    if memory:              memory=int(memory)
    mantis = MANTIS(
                        target_path=target_path,
                        output_folder=output_folder,
                        mantis_config=mantis_config,
                        evalue_threshold=evalue_threshold,
                        acceptable_range=acceptable_range,
                        overlap_value=overlap_value,
                        organism_details=organism_details,
                        domain_algorithm=domain_algorithm,
                        keep_files=keep_files,
                        skip_consensus=skip_consensus,
                        target_hmm=target_hmm,
                        verbose=verbose,
                        default_workers=default_workers,
                        chunk_size=chunk_size,
                        hmmer_threads=hmmer_threads,
                        user_cores=cores,
                        user_memory=memory,
    )
    mantis.run_mantis()

class MANTIS(MANTIS_MP):
    def __init__(self,
                 target_path=None,
                 output_folder=None,
                 mantis_config=None,
                 evalue_threshold=None,
                 acceptable_range=0.05,
                 overlap_value=0.1,
                 organism_details={},
                 domain_algorithm='dfs',
                 redirect_verbose=None,
                 keep_files=False,
                 skip_consensus=False,
                 target_hmm=None,
                 verbose=True,
                 default_workers=None,
                 chunk_size=None,
                 hmmer_threads=None,
                 user_cores=None,
                 user_memory=None,
                 ):
        self.output_folder = add_slash(output_folder)
        self.redirect_verbose=redirect_verbose
        print('------------------------------------------', flush=True, file=self.redirect_verbose)
        print_cyan('Setting up Mantis!', flush=True, file=self.redirect_verbose)
        print('------------------------------------------', flush=True, file=self.redirect_verbose)
        print_cyan('Reading config file and setting up paths', flush=True, file=self.redirect_verbose)
        self.target_path = target_path
        self.mantis_config = mantis_config
        if evalue_threshold:    self.evalue_threshold = evalue_threshold
        else:                   self.evalue_threshold = 1e-3
        self.acceptable_range = acceptable_range
        self.overlap_value = overlap_value
        self.organism_details = organism_details
        self.domain_algorithm = domain_algorithm
        self.keep_files = keep_files
        self.skip_consensus = skip_consensus
        self.default_workers=default_workers
        self.user_cores=user_cores
        self.user_memory=user_memory
        #chunk size is highly relevant in the execution time
        self.chunk_size=chunk_size
        MANTIS_Assembler.__init__(self,verbose=verbose,hmmer_threads=hmmer_threads,redirect_verbose=redirect_verbose,mantis_config=mantis_config)
        self.target_hmm = target_hmm
        datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        if self.target_path:
            print_cyan('This MANTIS process started running at '+ datetime_str, flush=True, file=self.redirect_verbose)
        self.chunks_to_annotate=[]
        self.chunks_to_fasta={}
        self.fastas_to_annotate=[]
        self.print_available_hardware()



    def print_available_hardware(self):
        if self.user_cores:
            print('Cores allocated:', self.user_cores)
        else:
            print('Cores allocated:', environment_cores)
        if self.user_memory:
            print('Memory allocated:', self.user_memory)
        else:
            print('Memory allocated:', round(available_ram,2))
        print('Workers per core:',worker_per_core)


    def __str__(self):
        output_list =[
        'Output folder:\t\t\t' + str(self.output_folder) + '\n' if self.output_folder else '',
        'Mantis config:\t\t\t'+str(self.mantis_config)+'\n' if self.mantis_config else '',
        'Target path:\t\t\t'+str(self.target_path)+'\n' if self.target_path else '',
        'E-value threshold:\t\t'+str(self.evalue_threshold)+'\n' if self.evalue_threshold else '',
        'Acceptable range:\t\t'+str(self.acceptable_range)+'\n' if self.acceptable_range else '',
        'Overlap value:\t\t\t'+str(self.overlap_value)+'\n' if self.overlap_value else '',
        'Target HMM:\t\t\t'+str(self.target_hmm)+'\n' if self.target_hmm else '',
        'Default workers:\t\t'+str(self.default_workers)+'\n' if self.default_workers else '',
        'User cores:\t\t\t\t'+str(self.user_cores)+'\n' if self.user_cores else '',
        'Chunk size:\t\t\t'+str(self.chunk_size)+'\n' if self.chunk_size else '',
        'Domain algorithm:\t\t'+str(self.domain_algorithm)+'\n' if self.domain_algorithm else '',
        '------------------------------------------']
        return 'User configuration:'+'\n'+'------------------------------------------'+'\n'+''.join(output_list)

    def generate_fastas_to_annotate(self):
        if self.target_path.split('.')[-1] in ['tsv']:
            self.annotate_multiple_samples()
        elif self.target_path[-1]==splitter:
            self.annotate_directory()
        elif self.target_path.split('.')[-2] in ['tar']:
            self.annotate_compressed_sample()
        elif self.target_path.split('.')[-1] in ['gz','zip']:
            self.annotate_compressed_sample()
        elif is_fasta(self.target_path):
            self.annotate_one_sample()
        else:
            print('Your file does not appear to be a fasta. If you want to annotate multiple samples, make sure your file has the <.tsv> extension.', flush=True, file=self.redirect_verbose)
            raise InvalidTargetFile
        for file_path, output_path, organism_details, count_seqs_original_file in self.fastas_to_annotate:
            Path(output_path).mkdir(parents=True, exist_ok=True)

    def annotate_multiple_samples(self):
        try:
            with open(self.target_path) as file:
                line = file.readline()
                line = file.readline()
                while line:
                    line = line.strip('\n').split()
                    if len(line) >= 2:
                        query_name = line[0]
                        line_path = line[1]
                        if len(line)>2: organism_details = ' '.join(line[2:])
                        else:           organism_details=''
                        count_seqs_original_file = get_seqs_count(line_path)
                        self.fastas_to_annotate.append([line_path,add_slash(self.output_folder+query_name),organism_details,count_seqs_original_file])
                    line = file.readline()
        except:
            print('If you want to annotate multiple samples, make sure your file is correctly formatted. Please see the examples in the <tests> folder.', flush=True, file=self.redirect_verbose)
            raise InvalidTargetFile

    def annotate_directory(self):
        try:
            list_dir=os.listdir(self.target_path)
            for file in list_dir:
                if 'faa' in file.split('.')[-1]:
                    query_name = '.'.join(file.split('.')[0:-1])
                    query_path = self.target_path+file
                    count_seqs_original_file = get_seqs_count(query_path)
                    self.fastas_to_annotate.append([query_path,add_slash(self.output_folder+query_name),None,count_seqs_original_file])
        except:
            raise InvalidTargetFile

    def annotate_compressed_sample(self):
        try:
            uncompressed_path=self.output_folder+'uncompressed_samples/'
            Path(uncompressed_path).mkdir(parents=True, exist_ok=True)
            uncompressing_function=uncompress_archive(source_filepath=self.target_path, extract_path=uncompressed_path)
            list_dir=os.listdir(uncompressed_path)
            for file in list_dir:
                if os.path.isdir(uncompressed_path+file):
                    sub_list_dir = os.listdir(uncompressed_path+file)
                    for sub_file in sub_list_dir:
                        if 'faa' in sub_file.split('.')[-1]:
                            query_name = '.'.join(sub_file.split('.')[0:-1])
                            query_path = add_slash(uncompressed_path+file)+sub_file
                            count_seqs_original_file = get_seqs_count(query_path)
                            self.fastas_to_annotate.append([query_path,add_slash(self.output_folder+query_name),None,count_seqs_original_file])
                if 'faa' in file.split('.')[-1]:
                    query_name = '.'.join(file.split('.')[0:-1])
                    query_path = uncompressed_path+file
                    count_seqs_original_file = get_seqs_count(query_path)
                    self.fastas_to_annotate.append([query_path,add_slash(self.output_folder+query_name),None,count_seqs_original_file])
        except:
            raise InvalidTargetFile

    def annotate_one_sample(self):
        count_seqs_original_file = get_seqs_count(self.target_path)
        self.fastas_to_annotate.append([self.target_path, self.output_folder, self.organism_details,count_seqs_original_file])

    def get_organism_lineage(self, taxon_id,stdout_file=None):
        lineage_file_path = self.mantis_paths['ncbi'] + 'taxidlineage.dmp'
        try:
            lineage_file = open(lineage_file_path,'r')
        except:
            print_cyan('Lineage dump is not present! If you\'d like to run taxonomic lineage annotation, please run < setup_databases >',flush=True,file=stdout_file)
            return []
        line = lineage_file.readline().strip('\n').replace('|', '')
        while line:
            line = line.split()
            if str(taxon_id) == str(line[0]):
                lineage_file.close()
                return line[1:]
            line = lineage_file.readline().strip('\n').replace('|', '')
        lineage_file.close()
        return []

    def check_internet_connection(self):
        try:
            response = requests.get("http://www.google.com")
            return True
        except requests.ConnectionError:
            print("Could not connect to internet!\nIf you would like to run offline make sure you introduce organism NCBI IDs instead of synonyms!")
            return False

    def get_taxa_ncbi(self, organism_name):
        if not self.check_internet_connection(): return
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=' + organism_name
        webpage = None
        c = 0
        while not webpage and c <= 10:
            req = requests.get(url)
            try:
                webpage = req.text
            except:
                c += 1
        taxa_id = re.search('<Id>\d+</Id>', webpage)
        if taxa_id: return re.search('\d+', taxa_id.group()).group()



    def setup_organism_lineage(self, organism_details,stdout_file):
        # when running with DRAX
        if 'organism_instance' in organism_details:
            print_cyan('Setting up organism lineage from provided organism instance', flush=True,file=stdout_file)
            organism_lineage = self.get_organism_lineage(organism_details['organism_instance'].get_detail('ncbi_taxon'),stdout_file=stdout_file)
        elif 'organism_lineage' in organism_details:
            print_cyan('Setting up organism lineage from provided organism lineage', flush=True,file=stdout_file)
            organism_lineage = organism_details['organism_lineage']
        # when running standalone
        elif 'ncbi_taxon' in organism_details:
            print_cyan('Setting up organism lineage from provided NCBI taxon id', flush=True,file=stdout_file)
            organism_lineage = self.get_organism_lineage(organism_details['ncbi_taxon'],stdout_file=stdout_file)
        elif 'synonyms' in organism_details:
            print_cyan('Setting up organism lineage from provided taxon synonym', flush=True,file=stdout_file)
            ncbi_taxon_id = self.get_taxa_ncbi(organism_details['synonyms'])
            organism_lineage = self.get_organism_lineage(ncbi_taxon_id,stdout_file=stdout_file)
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
            stdout_file=open(output_path+'Mantis.out','a+')
            organism_lineage=self.setup_organism_lineage(organism_details_dict,stdout_file)
            self.fastas_to_annotate[i][2]=organism_lineage
            print('------------------------------------------', flush=True, file=stdout_file)
            if organism_lineage:
                print_cyan('Target file:\n'+file_path+'\n has the following taxonomy lineage: ' + ' > '.join(organism_lineage), flush=True,file=stdout_file)
            else:
                print_cyan('Target file:\n'+file_path+'\n has no organism lineage!', flush=True, file=stdout_file)
            print('------------------------------------------', flush=True,file=stdout_file)

            stdout_file.close()

    def remove_non_essential_files(self):
        for file_path,output_path,organism_details,count_seqs_original_file in self.fastas_to_annotate:
            if os.path.exists(output_path + 'fasta_chunks/'):
                shutil.rmtree(output_path + 'fasta_chunks/')

    def remove_uncompressed_files(self):
        if os.path.exists(self.output_folder + 'uncompressed_samples/'):
            shutil.rmtree(self.output_folder + 'uncompressed_samples/')

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
        start_time=time()
        self.run_hmmer()
        print('HMMER took',int(time()-start_time),'seconds to run',flush=True,file=self.redirect_verbose)
        processing_start_time=time()
        start_time=time()
        self.process_output()
        print('Output processing took',int(time()-start_time),'seconds to run',flush=True,file=self.redirect_verbose)
        start_time=time()
        self.interpret_output()
        print('Output interpretation took',int(time()-start_time),'seconds to run',flush=True,file=self.redirect_verbose)
        if not self.skip_consensus:
            start_time = time()
            self.get_consensus_output()
            print('Consensus generation took',int(time()-start_time),'seconds to run',flush=True,file=self.redirect_verbose)
        start_time=time()
        self.merge_mantis_output()
        print('Output merge took',int(time()-start_time),'seconds to run',flush=True,file=self.redirect_verbose)
        start_time=time()
        self.remove_uncompressed_files()
        print('In total, Mantis took',int(time()-processing_start_time),'seconds to process HMMER\'s output',flush=True,file=self.redirect_verbose)
        if not self.keep_files:
            print_cyan('Removing temporary files', flush=True, file=self.redirect_verbose)
            self.remove_non_essential_files()
        stdout_file=open(self.output_folder+'Mantis.out','a+')
        print('Mantis has finished running! It took '+str(int(time() - self.start_time))+' seconds to run.', flush=True, file=stdout_file)
        stdout_file.close()

if __name__ == '__main__':
    ts='/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis/tests/test_Here/metag_test.faa'
    #ts='/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis/tests/test_Here/partial_mg_sample.faa'
    to='/home/pedroq/Desktop/test_hmm/test_general_annot'
    tsv_file='/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis/tests/test_sample.faa'
    tsv_file='/home/pedroq/Desktop/test_hmm/test2.faa'
    tsv_file='/home/pedroq/Desktop/test_hmm/test1.faa'
    #l=12158
    l=1306215
    mm=MANTIS(target_path=tsv_file,output_folder=to,mantis_config='/home/pedroq/Desktop/test_hmm/t.config',organism_details='escherichia coli')#,chunk_size=1000)
    print(mm.check_internet_connection())
    #mm.keep_files=True
    #mm.run_mantis()