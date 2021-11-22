try:
    from source.Exceptions import *
    from source.utils import *
    from source.MANTIS_NLP import MANTIS_NLP
except:
    from Exceptions import *
    from utils import *
    from MANTIS_NLP import MANTIS_NLP

if not unifunc_downloaded():
    download_unifunc()
if not diamond_downloaded():
    download_diamond()


class MANTIS_DB(MANTIS_NLP):

#####################   Main function
    @timeit_class
    def setup_databases(self, force_download=False):
        if not cython_compiled():
            compile_cython()
        if not unifunc_downloaded():
            download_unifunc()
        if not diamond_downloaded():
            download_diamond()
        self.output_folder = f'{MANTIS_FOLDER}setup_databases/'
        self.mantis_out = f'{self.output_folder}Mantis.out'
        if file_exists(self.mantis_out):
            os.remove(self.mantis_out)
        Path(self.mantis_paths['default']).mkdir(parents=True, exist_ok=True)
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        dbs_list = [
            'ncbi_res',
            'gtdb_res',
            'pfam' if self.mantis_paths['pfam'][0:2] != 'NA' else None,
            'kofam' if self.mantis_paths['kofam'][0:2] != 'NA' else None,
            'tcdb' if self.mantis_paths['tcdb'][0:2] != 'NA' else None,
            'NCBI' if self.mantis_paths['NCBI'][0:2] != 'NA' else None,
        ]
        # DOWNLOADING
        self.prepare_queue_setup_databases(dbs_list, force_download)
        # for unzipping tax specific hmms
        self.prepare_queue_setup_databases_tax(force_download)
        worker_count = estimate_number_workers_setup_database(len(self.queue), user_cores=self.user_cores)
        if worker_count: print(f'Database will be setup with {worker_count} workers!', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_setup_databases, worker_count)
        print_cyan('Finished downloading all data!', flush=True, file=self.redirect_verbose)
        # METADATA
        self.prepare_queue_extract_metadata()
        worker_count = estimate_number_workers_setup_database(len(self.queue), minimum_jobs_per_worker=2, user_cores=self.user_cores)
        if worker_count:print(f'Metadata will be extracted with {worker_count} workers!', flush=True,file=self.redirect_verbose)
        self.processes_handler(self.worker_extract_NOG_metadata, worker_count)
        print('NOGG will now be compiled',flush=True,file=self.redirect_verbose)
        self.compile_NOGG(force_download)
        # SPLITTING
        if self.hmm_chunk_size:
            print_cyan('Will now split data into chunks!', flush=True, file=self.redirect_verbose)
            self.prepare_queue_split_hmms()
            worker_count = estimate_number_workers_setup_database(len(self.queue), user_cores=self.user_cores)
            print(f'Database will be split with {worker_count} workers!', flush=True,
                  file=self.redirect_verbose)
            self.processes_handler(self.worker_split_hmms, worker_count)
        self.prepare_queue_press_custom_hmms()
        worker_count = estimate_number_workers_setup_database(len(self.queue), user_cores=self.user_cores)
        print(f'HMMs will be pressed with {worker_count} workers!', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_press_custom_hmms, worker_count)
        print('Preparing NLP Resources!', flush=True, file=self.redirect_verbose)
        MANTIS_NLP.__init__(self)
        print_cyan('Finished setting up databases!', flush=True, file=self.redirect_verbose)

    #this is a standalone execution, this is just for the developers to update the NOG tars and then upload them to Github
    @timeit_class
    def extract_nog_metadata(self, metadata_path):
        # metadata_path where we are extracting metadata to. This is an auxilliary function for compillation of NOG's metadata
        self.mantis_paths['NOG'] = add_slash(metadata_path + 'NOG')
        self.output_folder = metadata_path + 'extract_nog_metadata/'
        self.mantis_out = self.output_folder + 'Mantis.out'
        self.mantis_paths['resources'] = add_slash(metadata_path + 'Resources')
        self.mantis_paths['default'] = add_slash(metadata_path)

        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        Path(self.mantis_paths['resources']).mkdir(parents=True, exist_ok=True)

        if file_exists(self.mantis_out):
            os.remove(self.mantis_out)
        stdout_file = open(self.mantis_out, 'a+')
        taxon_ids = self.get_ref_taxon_ids('NOG')
        self.download_and_unzip_eggnogdb()

        for t in taxon_ids:
            Path(self.mantis_paths['NOG'] + t).mkdir(parents=True, exist_ok=True)
            url = f'http://eggnog5.embl.de/download/latest/per_tax_level/{t}/{t}_annotations.tsv.gz'
            download_file(url, output_folder=add_slash(self.mantis_paths['NOG'] + t), stdout_file=stdout_file)
            uncompress_archive(source_filepath=add_slash(self.mantis_paths['NOG'] + t) + f'{t}_annotations.tsv.gz',
                               stdout_file=stdout_file, remove_source=True)
        self.prepare_queue_extract_metadata()
        worker_count = estimate_number_workers_setup_database(len(self.queue), minimum_jobs_per_worker=2, user_cores=self.user_cores)
        print(f'Metadata will be extracted with {worker_count} workers!', flush=True,file=self.redirect_verbose)
        print(f'The following taxons will be extracted:{taxon_ids}', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_extract_NOG_metadata, worker_count)

#####################   Filling queue with jobs

    def prepare_queue_setup_databases(self, dbs_list, force_download):
        for database in dbs_list:
            if database:
                self.queue.append([database, force_download, self.mantis_out])

    def prepare_queue_setup_databases_tax(self, force_download):
        stdout_file = open(self.mantis_out, 'a+')
        if self.mantis_paths['NOG'][0:2] != 'NA':
            list_taxon_ids = self.get_taxon_for_queue_NOGT(stdout_file=stdout_file)
            for taxon_id in list_taxon_ids:
                self.queue.append(['NOG', force_download, taxon_id, self.mantis_out])


    def prepare_queue_split_hmms(self):
        print('Checking which HMMs need to be split, this may take a while...', flush=True, file=self.redirect_verbose)
        taxon_nog_hmms = self.get_taxon_for_queue_NOG_split_hmms()
        taxon_ncbi_hmms = self.get_taxon_for_queue_ncbi_split_hmms()
        general_hmms = self.get_hmm_for_queue_split_hmms()
        print('Will split: ', [get_path_level(i) for i in taxon_nog_hmms + taxon_ncbi_hmms + general_hmms], flush=True,file=self.redirect_verbose)
        for hmm_path in taxon_nog_hmms + taxon_ncbi_hmms + general_hmms:
            self.queue.append([hmm_path, self.mantis_out])

    def prepare_queue_press_custom_hmms(self):
        print('Checking which custom hmms need to be pressed', flush=True, file=self.redirect_verbose)
        hmms_list=[]
        for hmm_path in self.get_custom_refs_paths(folder=False):
            if hmm_path.endswith('.hmm'):
                hmm_folder=add_slash(SPLITTER.join(hmm_path.split(SPLITTER)[:-1]))
                hmm_name=hmm_path.split(SPLITTER)[-1]
                if 'chunks' not in os.listdir(hmm_folder):
                    if hmm_folder[0:2] != 'NA':
                        if not self.check_missing_chunk_files(hmm_name,hmm_folder):
                            hmms_list.append(hmm_path)
        if hmms_list:
            print(f'Will hmmpress: {hmms_list}', flush=True,file=self.redirect_verbose)
        for hmm_path in hmms_list:
            self.queue.append([hmm_path, self.mantis_out])

    def prepare_queue_extract_metadata(self):
        if self.mantis_paths['NOG'][0:2] != 'NA':
            print_cyan('Will now extract metadata from NOG!', flush=True, file=self.redirect_verbose)
            stdout_file = open(self.mantis_out, 'a+')
            print('Checking which NOGs we need to extract metadata from', flush=True, file=stdout_file)
            self.unpack_NOG_sql(stdout_file=stdout_file)
            to_download=False
            for taxon_id in os.listdir(self.mantis_paths['NOG']):
                if os.path.isdir(self.mantis_paths['NOG'] + taxon_id) and taxon_id!='NOGG':
                    target_annotation_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_annotations.tsv'
                    target_sql_file = add_slash(self.mantis_paths['NOG'] + taxon_id) +  f'{taxon_id}_sql_annotations.tsv'
                    if file_exists(target_sql_file):
                        if len(self.get_hmms_annotation_file(target_sql_file, hmm_col=0)) != len(self.get_hmms_annotation_file(target_annotation_file, hmm_col=1)):
                            to_download=True
                            self.queue.append([target_sql_file, target_annotation_file, taxon_id, self.mantis_out])
                        else:
                            print(f'Skipping metadata extraction for NOGT {taxon_id}', flush=True, file=stdout_file)
                    else:
                        to_download = True
                        self.queue.append([target_sql_file, target_annotation_file, taxon_id, self.mantis_out])
            if to_download: self.download_and_unzip_eggnogdb()
            stdout_file.close()

#####################   Executing queue jobs

    def worker_setup_databases(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            if len(record) == 3:
                database, force_download, stdout_path = record
                taxon_id = None
            else:
                database, force_download, taxon_id, stdout_path = record
            self.download_database(database, force_download, taxon_id, stdout_path)
            if not self.check_reference_exists(database, taxon_id):
                stdout_file = open(stdout_path, 'a+')
                print(f'Setup failed on {database} {taxon_id}', flush=True, file=stdout_file)
                stdout_file.close()
                queue.insert(0, record)

    def worker_split_hmms(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            hmm_path, stdout_path = record
            self.split_hmm_into_chunks(hmm_path, stdout_path)

    def worker_press_custom_hmms(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            hmm_path, stdout_path = record
            self.press_custom_hmms(hmm_path, stdout_path)

    def worker_extract_NOG_metadata(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            target_sql_file, target_annotation_file, taxon_id, stdout_path = record
            self.get_metadata_hmms(target_annotation_file=target_annotation_file, target_sql_file=target_sql_file,
                                   taxon_id=taxon_id, stdout_path=stdout_path)

#####################   Main download function

    def download_database(self, database, force_download=False, taxon_id=None, stdout_path=None):
        stdout_file = open(stdout_path, 'a+')
        if database == 'ncbi_res':      self.download_ncbi_resources(force_download=force_download, stdout_file=stdout_file)
        elif database == 'gtdb_res':    self.download_gtdb_resources(force_download=force_download, stdout_file=stdout_file)
        elif database == 'pfam':        self.download_pfam(force_download=force_download, stdout_file=stdout_file)
        elif database == 'kofam':       self.download_kofam(force_download=force_download, stdout_file=stdout_file)
        elif database == 'tcdb':        self.download_tcdb(force_download=force_download, stdout_file=stdout_file)
        elif database == 'NCBI':        self.download_NCBI(force_download=force_download, stdout_file=stdout_file)
        elif database == 'NOG':         self.download_NOGT(taxon_id=taxon_id, stdout_file=stdout_file)
        if taxon_id:
            print(f'Finished downloading {database} with taxon {taxon_id}', flush=True, file=stdout_file)
        else:
            print(f'Finished downloading {database}', flush=True, file=stdout_file)
        stdout_file.close()

    def download_ncbi_resources(self, force_download=False, stdout_file=None):
        ncbi_resources=add_slash(self.mantis_paths['resources']+'NCBI')
        Path(ncbi_resources).mkdir(parents=True, exist_ok=True)
        if file_exists(ncbi_resources + 'taxidlineage.dmp', force_download) and file_exists(ncbi_resources + 'gc.prt', force_download):
            print('NCBI taxonomic lineage file already exists! Skipping...', flush=True, file=stdout_file)
            return
        try:
            os.remove(ncbi_resources + 'taxidlineage.dmp')
        except:
            pass
        taxonomy_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz'
        translation_tables_url='ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt'
        with open(ncbi_resources + 'readme.md', 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'These files were downloaded on {datetime_str} from:\n{taxonomy_url}\n{translation_tables_url}\nThey are used to traceback organism lineage for taxonomically resolved annotations and to translate CDS')
        print_cyan('Downloading and unzipping NCBI resources', flush=True, file=stdout_file)
        for url in [taxonomy_url,translation_tables_url]:
            download_file(url, output_folder=ncbi_resources, stdout_file=stdout_file)
        uncompress_archive(source_filepath=ncbi_resources + 'new_taxdump.tar.gz',
                           extract_path=ncbi_resources, stdout_file=stdout_file, remove_source=True)
        # deleting the extra files that come with the tar dump
        files_dir = os.listdir(ncbi_resources)
        for file in files_dir:
            if file not in ['taxidlineage.dmp','gc.prt','readme.md']:
                os.remove(ncbi_resources + file)

    def download_gtdb_resources(self, force_download=False, stdout_file=None):
        gtdb_resources=add_slash(self.mantis_paths['resources']+'GTDB')
        self.launch_gtdb_connector(resources_folder=gtdb_resources)

    def get_common_links_metadata(self, string, res):
        ec = find_ecs(string)
        if ec:
            if 'enzyme_ec' not in res: res['enzyme_ec'] = set()
            res['enzyme_ec'].update(ec)
        tc = find_tcdb(string)
        if tc:
            if 'tcdb' not in res: res['tcdb'] = set()
            res['tcdb'].update(tc)
        tigr = find_tigrfam(string)
        if tigr:
            if 'tigrfam' not in res: res['tigrfam'] = set()
            res['tigrfam'].update(tigr)
        ko = find_ko(string)
        if ko:
            if 'kegg_ko' not in res: res['kegg_ko'] = set()
            res['kegg_ko'].update(ko)
        pfam = find_pfam(string)
        if pfam:
            if 'pfam' not in res: res['pfam'] = set()
            res['pfam'].update(pfam)
        cog = find_cog(string)
        if cog:
            if 'cog' not in res: res['cog'] = set()
            res['cog'].update(cog)
        go = find_go(string)
        if go:
            if 'go' not in res: res['go'] = set()
            res['go'].update(go)

    def write_metadata(self,metadata,metadata_file):
        with open(metadata_file,'w+') as file:
            for seq in metadata:
                link_line = f'{seq}\t|'
                for link_type in metadata[seq]:
                    for inner_link in metadata[seq][link_type]:
                        if inner_link:
                            link_line += f'\t{link_type}:{inner_link}'
                link_line+='\n'
                file.write(link_line)



#####################   PFAM

    def read_pfam2go(self):
        res = {}
        with open(self.mantis_paths['pfam'] + 'pfam2go') as pfam2go_file:
            line = pfam2go_file.readline()
            while line:
                line = line.strip('\n')
                if '!' not in line[0]:
                    line = line.split('>')
                    line = [i.strip() for i in line]
                    pfam_id = line[0].split()[0].replace('Pfam:', '')
                    go_annots = line[1].split(';')
                    go_description = [i.replace('GO:', '').strip() for i in go_annots if not re.search('GO:\d{3,}', i)]
                    go_ids = [i.replace('GO:', '').strip() for i in go_annots if re.search('GO:\d{3,}', i)]
                    if pfam_id not in res:
                        res[pfam_id] = {'go': set(go_ids), 'description': set(go_description)}
                    else:
                        res[pfam_id]['go'].update(go_ids)
                        res[pfam_id]['description'].update(go_description)
                line = pfam2go_file.readline()
        return res

    def is_good_description(self, hmm, row_description):
        temp = [i.lower() for i in row_description.split()]
        if hmm.lower() in temp and 'protein' in temp and len(temp) == 2:
            return False
        if re.search(' [uU]nknown [Ff]unction', row_description): return False
        return True

    def build_pfam_line(self, hmm, metadata):
        link_line = ''
        pfam_ids = set(metadata['pfam'])
        metadata['pfam'] = set()
        for p_id in pfam_ids:
            metadata['pfam'].add(p_id.split('.')[0])
        for link_type in metadata:
            for inner_link in metadata[link_type]:
                link_line += f'\t{link_type}:{inner_link}'
        return hmm + f'\t|{link_line}\n'

    def compile_pfam_metadata(self):
        pfam2go = self.read_pfam2go()
        with open(self.mantis_paths['pfam'] + 'metadata.tsv', 'w+') as pfam_metadata_file:
            with open(self.mantis_paths['pfam'] + 'Pfam-A.hmm.dat') as pfam_dat_file:
                stop = False
                line = pfam_dat_file.readline()
                while line:
                    line = line.strip('\n').split('   ')
                    if len(line) == 2:
                        row_header, row_description = line
                        if row_header == '#=GF ID':
                            stop = True
                            hmm = str(row_description)
                        elif row_header == '#=GF AC':
                            pfam_accession = str(row_description)
                            pfam_accession = pfam_accession.split('.')[0].strip()
                        elif row_header == '#=GF DE' and stop:
                            if pfam_accession in pfam2go:
                                current_metadata = pfam2go[pfam_accession]
                            else:
                                current_metadata = {'description': set()}
                            current_metadata['pfam'] = set()
                            current_metadata['pfam'].add(pfam_accession)
                            if self.is_good_description(hmm, row_description):
                                current_metadata['description'].add(row_description)
                            self.get_common_links_metadata(row_description, current_metadata)
                            stop = False
                            metadata_line = self.build_pfam_line(hmm, current_metadata)
                            pfam_metadata_file.write(metadata_line)
                    line = pfam_dat_file.readline()
        remove_file(self.mantis_paths['pfam'] + 'Pfam-A.hmm.dat')
        remove_file(self.mantis_paths['pfam'] + 'pfam2go')

    def download_pfam(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['pfam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('pfam', force_download) and \
                file_exists(self.mantis_paths['pfam'] + 'metadata.tsv', force_download):
            print('Pfam hmm already exists! Skipping...', flush=True, file=stdout_file)
            return
        pfam_hmm = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
        pfam_metadata = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
        pfam2go = 'http://current.geneontology.org/ontology/external2go/pfam2go'
        with open(self.mantis_paths['pfam'] + 'readme.md', 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'This hmm was downloaded on {datetime_str} from:\n{pfam_hmm}\nMetadata was downloaded from:\n{pfam_metadata}\nPfam2GO was downloaded from:\n{pfam2go}')
        print_cyan('Downloading and unzipping Pfam hmms ', flush=True, file=stdout_file)
        to_download = []
        to_unzip = []
        if not file_exists(self.mantis_paths['pfam'] + 'Pfam-A.hmm', force_download):
            to_download.append(pfam_hmm)
            to_unzip.append('Pfam-A.hmm.gz')
        if not file_exists(self.mantis_paths['pfam'] + 'metadata.tsv', force_download):
            to_unzip.append('Pfam-A.hmm.dat.gz')
            to_download.append(pfam_metadata)
            to_download.append(pfam2go)
        for url in to_download:
            download_file(url, output_folder=self.mantis_paths['pfam'], stdout_file=stdout_file)
        for file_to_unzip in to_unzip:
            uncompress_archive(source_filepath=self.mantis_paths['pfam'] + file_to_unzip, stdout_file=stdout_file,
                               remove_source=True)
        if not self.check_reference_exists('pfam'):
            run_command('hmmpress ' + self.mantis_paths['pfam'] + 'Pfam-A.hmm', stdout_file=stdout_file)
        self.compile_pfam_metadata()

#####################   KOFAM

    def compile_kofam_metadata(self):
        metadata_to_write={}
        self.get_link_kofam_ko_list(metadata_to_write)
        self.get_link_kofam_ko_to_binary(metadata_to_write, target_file='ko2cog.xl')
        self.get_link_kofam_ko_to_binary(metadata_to_write, target_file='ko2go.xl')
        self.get_link_kofam_ko_to_binary(metadata_to_write, target_file='ko2tc.xl')
        self.get_link_kofam_ko_to_binary(metadata_to_write, target_file='ko2cazy.xl')
        self.write_metadata(metadata_to_write, self.mantis_paths['kofam'] + 'metadata.tsv')
        remove_file(self.mantis_paths['kofam'] + 'ko_list')
        remove_file(self.mantis_paths['kofam'] + 'ko2cazy.xl')
        remove_file(self.mantis_paths['kofam'] + 'ko2cog.xl')
        remove_file(self.mantis_paths['kofam'] + 'ko2go.xl')
        remove_file(self.mantis_paths['kofam'] + 'ko2tc.xl')

    def get_link_kofam_ko_list(self, res):
        file_path = self.mantis_paths['kofam'] + 'ko_list'
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                ko, description = line[0], line[-1]
                if ko not in res: res[ko]={}
                if '[EC:' in description:
                    description, temp_links = description.split('[EC:')
                else:
                    temp_links = description
                self.get_common_links_metadata(temp_links, res[ko])
                if 'kegg_ko' not in res[ko]: res[ko]['kegg_ko']=set()
                res[ko]['kegg_ko'].add(ko)
                if 'description' not in res[ko]:
                    res[ko]['description']=set()
                res[ko]['description'].add(description)
                line = file.readline()

    def get_link_kofam_ko_to_binary(self, res, target_file):
        file_path = self.mantis_paths['kofam'] + target_file
        if 'ko2tc' in target_file:
            target_link = 'tcdb'
        else:
            target_link = target_file.replace('ko2', '').replace('.xl', '')
        with open(file_path) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                ko, link = line
                link = link.strip('[]').split(':')[1].split()
                if ko not in res: res[ko]={}
                if target_link not in res[ko]: res[ko][target_link]=set()
                res[ko][target_link].update(link)
                line = file.readline()



    def download_kofam(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['kofam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('kofam', force_download) and \
                file_exists(self.mantis_paths['kofam'] + 'metadata.tsv', force_download):
            print('KOfam HMM already exists! Skipping...', flush=True, file=stdout_file)
            return
        kofam_hmm = 'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz'
        ko_list = 'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz'
        ko_to_cog = 'https://www.kegg.jp/kegg/files/ko2cog.xl'
        ko_to_go = 'https://www.kegg.jp/kegg/files/ko2go.xl'
        ko_to_tc = 'https://www.kegg.jp/kegg/files/ko2tc.xl'
        ko_to_cazy = 'https://www.kegg.jp/kegg/files/ko2cazy.xl'
        with open(self.mantis_paths['kofam'] + 'readme.md', 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'This hmm was downloaded on {datetime_str} from:\n{kofam_hmm}\nMetadata was downloaded from:\n{ko_list}\n{ko_to_cog}\n{ko_to_go}\n{ko_to_tc}\n{ko_to_cazy}')
        print_cyan('Downloading and unzipping KOfam hmms ', flush=True, file=stdout_file)
        for url in [kofam_hmm, ko_list, ko_to_cog, ko_to_go, ko_to_tc, ko_to_cazy]:
            download_file(url, output_folder=self.mantis_paths['kofam'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'profiles.tar.gz',
                           extract_path=self.mantis_paths['kofam'], stdout_file=stdout_file, remove_source=True)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'ko_list.gz', stdout_file=stdout_file,
                           remove_source=True)
        merge_profiles(self.mantis_paths['kofam'] + 'profiles/', self.mantis_paths['kofam'] + 'kofam_merged.hmm',stdout_file=stdout_file)
        run_command('hmmpress ' + self.mantis_paths['kofam'] + 'kofam_merged.hmm', stdout_file=stdout_file)
        self.compile_kofam_metadata()

#####################   NCBI


    def download_NCBI(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['NCBI']).mkdir(parents=True, exist_ok=True)
        ncbi_hmm = 'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz'
        metadata = 'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv'

        # we cant verify a priori which foulders we should have, so you need to delete the folder to restart

        if self.check_reference_exists('NCBI', force_download) and \
                file_exists(self.mantis_paths['NCBI'] + 'readme.md', force_download):
            print('NCBI hmm folder already exists! Skipping...', flush=True, file=stdout_file)
            return

        with open(self.mantis_paths['NCBI'] + 'readme.md', 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'This hmm was downloaded on {datetime_str} from:\n{ncbi_hmm}\nMetadata was downloaded from:\n{metadata}')

        print_cyan('Downloading and unzipping NCBI hmms ', flush=True, file=stdout_file)
        for url in [ncbi_hmm, metadata]:
            download_file(url, output_folder=self.mantis_paths['NCBI'], stdout_file=stdout_file)
        move_file(self.mantis_paths['NCBI'] + 'hmm_PGAP.HMM.tgz', self.mantis_paths['NCBI'] + 'hmm_PGAP.HMM.tar.gz')
        uncompress_archive(source_filepath=self.mantis_paths['NCBI'] + 'hmm_PGAP.HMM.tar.gz',
                           extract_path=self.mantis_paths['NCBI'], stdout_file=stdout_file, remove_source=True)
        self.compile_hmms_NCBI(stdout_file=stdout_file)
        remove_file(self.mantis_paths['NCBI'] + 'hmm_PGAP.tsv')
        if file_exists(self.mantis_paths['NCBI'] + 'hmm_PGAP'):
            shutil.rmtree(self.mantis_paths['NCBI'] + 'hmm_PGAP')

    # sorting profiles by the taxonomic id from ncbi metadata file
    def compile_hmms_NCBI(self,stdout_file=None):
        sorted_metadata = self.sort_hmms_NCBI()
        self.write_metadata_ncbi(sorted_metadata)
        self.assign_hmm_profiles(sorted_metadata, stdout_file=stdout_file)

    def assign_hmm_profiles(self, sorted_metadata, stdout_file=None):
        for taxa in sorted_metadata:
            for hmm, hmm_label, description, enzyme_ec,go_terms,common_links in sorted_metadata[taxa]:
                try:
                    copy_file(add_slash(self.mantis_paths['NCBI'] + 'hmm_PGAP') + hmm + '.HMM',
                              add_slash(add_slash(self.mantis_paths['NCBI'] + taxa) + 'to_merge') + hmm + '.hmm')
                except:
                    pass
            if os.listdir(add_slash(add_slash(self.mantis_paths['NCBI'] + taxa) + 'to_merge')):
                merge_profiles(add_slash(add_slash(self.mantis_paths['NCBI'] + taxa) + 'to_merge'),
                               add_slash(self.mantis_paths['NCBI'] + taxa) + taxa + '_merged.hmm')
                run_command('hmmpress ' + add_slash(self.mantis_paths['NCBI'] + taxa) + taxa + '_merged.hmm',stdout_file=stdout_file)
            else:
                shutil.rmtree(self.mantis_paths['NCBI'] + taxa)




    def write_metadata_ncbi(self, sorted_metadata):
        for taxa in sorted_metadata:
            Path(self.mantis_paths['NCBI'] + taxa).mkdir(parents=True, exist_ok=True)
            Path(add_slash(self.mantis_paths['NCBI'] + taxa) + 'to_merge').mkdir(parents=True, exist_ok=True)
            with open(add_slash(self.mantis_paths['NCBI'] + taxa) + 'metadata.tsv', 'w+') as metadata_file:
                for hmm, hmm_label, description, enzyme_ec,go_terms,common_links in sorted_metadata[taxa]:
                    line = [hmm_label, '|', f'description:{description}']

                    for db in common_links:
                        for db_id in common_links[db]:
                            if f'{db}:{db_id}' not in line:
                                line.append(f'{db}:{db_id}')
                    for ec in enzyme_ec:
                        if f'enzyme_ec:{ec}' not in line:
                            line.append(f'enzyme_ec:{ec}')
                    for go_term in go_terms:
                        if f'go:{go_term}' not in line:
                            line.append(f'go:{go_term}')
                    #ncbi also contains tigrfam hmms
                    if hmm_label.startswith('TIGR'):
                        if f'tigrfam:{hmm_label}' not in line:
                            line.append(f'tigrfam:{hmm_label}')

                    line = '\t'.join(line) + '\n'
                    metadata_file.write(line)

    def get_ncbi_domains(self):
        # this is from NCBI's top level taxonomy page
        return [
            '2157',  # Archaea
            '2',  # Bacteria
            '2759',  # Eukaryota
            '10239',  # Viruses
            '28384',  # Others
            '12908',  # Unclassified
        ]

    def sort_hmms_NCBI(self):
        general_taxon_ids = self.get_ncbi_domains()
        res = {'NCBIG': []}
        metadata = self.mantis_paths['NCBI'] + 'hmm_PGAP.tsv'
        with open(metadata) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.split('\t')
                hmm, hmm_label, description, enzyme_ec, go_terms, taxa_id = line[0], line[2], line[10], line[12], line[13], line[15]
                common_links={}
                self.get_common_links_metadata(description,common_links)
                for db in common_links:
                    for db_id in common_links[db]:
                        description=description.replace(db_id,'').strip()
                description = description.replace('(Provisional)', '')
                description = description.strip()
                enzyme_ec = [i for i in enzyme_ec.split(',') if i]
                go_terms = [i.replace('GO:','') for i in go_terms.split(',') if i]
                line = file.readline()
                if taxa_id:
                    if taxa_id not in res: res[taxa_id] = []
                    res[taxa_id].append([hmm, hmm_label, description, enzyme_ec,go_terms,common_links])
                    if taxa_id in general_taxon_ids:
                        res['NCBIG'].append([hmm, hmm_label, description, enzyme_ec,go_terms,common_links])
                else:
                    res['NCBIG'].append([hmm, hmm_label, description, enzyme_ec,go_terms,common_links])
        return res

#####################   TCDB


    def parse_tsv_tcdb(self,file_path,key_col,val_col,res_key,res,val_clean_function=None):
        with open(file_path) as file:
            line=file.readline()
            while line:
                line=line.strip('\n')
                line=line.split('\t')
                if line[key_col] not in res: res[line[key_col]]={}
                if res_key not in res[line[key_col]]:res[line[key_col]][res_key]=set()
                if val_clean_function:
                    line[val_col]=val_clean_function(line[val_col])
                res[line[key_col]][res_key].add(line[val_col])
                line = file.readline()

    def read_tcdb_headers(self):
        file_path=self.mantis_paths['tcdb'] + 'tcdb'
        res={}
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                if '>' in line:
                    line=line.split('|')
                    uniprot_accession=line[2]
                    tc_id=line[3].split()[0]
                    description=line[3].replace(tc_id,'').strip()
                    if re.search('([Hh]ypothetical|[Uu]ncharacterized|[Uu]ndetermined)',description): description=''
                    description=description.split('[')[0]
                    if re.search('[A-Z]+=',description):
                        description=description.split(re.search('[A-Z]+=',description).group())[0]
                    if ' - ' in description:
                        description=description.split(' - ')[0]
                    description=description.strip()
                    if description:
                        res[uniprot_accession]={'tcdb':tc_id,'description': {description.strip()}}
                    else:
                        res[uniprot_accession]={'tcdb':tc_id}
                line=file.readline()
        return res

    def remove_bad_entries(self,all_seqs):
        #some proteins are not in uniprot, but are in some other dbs.... Im not sure if they should be removed, but we will consider uniprot as the central repo, so we will remove them
        chunks_post=chunk_generator(all_seqs,500)
        all_found=set()
        for seqs_chunk in chunks_post:
            c = 0
            url='https://www.uniprot.org/uploadlists/'
            #post_obj={'query':'P13368 P20806 ZP_02735776','format':'tab','from':'ACC','to':'ACC'}
            post_obj={'query':' '.join(seqs_chunk),'format':'tab','from':'ACC','to':'ACC'}
            headers = {'Content_Type': 'form-data'}
            response=requests.post(url, data=post_obj,headers=headers)
            while response.status_code!=200:
                response = requests.post(url, data=post_obj, headers=headers)
                if c==10:
                    kill_switch(ConnectionError,url)
                c += 1

            response_text=response.text
            for seq in seqs_chunk:
                if seq in response_text:
                    all_found.add(seq)
        return all_found

    def yield_tcdb_seqs(self,tcdb_fasta,seqs_found):
        res = {}
        tcdb_seqs=read_protein_fasta_generator(tcdb_fasta)
        for seq_id,seq in tcdb_seqs:
            uniprot_accession = seq_id.split('|')[2]
            if uniprot_accession in seqs_found:
                yield uniprot_accession,seq

    def generate_tcdb_fasta(self,tcdb_fasta,seqs_found):
        write_fasta_generator(self.yield_tcdb_seqs(tcdb_fasta,seqs_found),self.mantis_paths['tcdb']+'tcdb.faa')



    def compile_tcdb_metadata(self):
        #acession will be the key
        all_seqs=self.read_tcdb_headers()
        seqs_found=list(all_seqs.keys())

        seqs_found=self.remove_bad_entries(list(all_seqs.keys()))
        self.generate_tcdb_fasta(self.mantis_paths['tcdb']+'tcdb',seqs_found)
        #here tc ids will be keys
        metadata={}
        self.parse_tsv_tcdb(self.mantis_paths['tcdb']+'go.py',1,0,'go',metadata,val_clean_function=lambda a : a.replace('GO:',''))
        self.parse_tsv_tcdb(self.mantis_paths['tcdb']+'pfam.py',1,0,'pfam',metadata)
        #we now add the tc specific metadata to each acession
        metadata_to_write={}
        for seq in seqs_found:
            tc_id=all_seqs[seq]['tcdb']
            if tc_id in metadata:
                metadata_to_write[seq] = metadata[tc_id]
            else: metadata_to_write[seq]={}
            if 'description' in all_seqs[seq]: metadata_to_write[seq]['description']=all_seqs[seq]['description']
            metadata_to_write[seq]['tcdb']={tc_id}


        self.write_metadata(metadata_to_write,self.mantis_paths['tcdb']+'metadata.tsv')

    def download_tcdb(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['tcdb']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('tcdb', force_download) and \
                file_exists(self.mantis_paths['tcdb'] + 'metadata.tsv', force_download):
            print('TCDB sequences already exists! Skipping...', flush=True, file=stdout_file)
            return
        tcdb_seqs = 'http://www.tcdb.org/public/tcdb'
        tcdb_go='https://www.tcdb.org/cgi-bin/projectv/public/go.py'
        tcdb_pfam='https://www.tcdb.org/cgi-bin/projectv/public/pfam.py'

        with open(self.mantis_paths['tcdb'] + 'readme.md', 'w+') as  file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'These sequences were downloaded on {datetime_str} from:\n{tcdb_seqs}\n'
                           f'Metadata was downloaded from:\n{tcdb_go}\n{tcdb_pfam}\n')
        print_cyan('Downloading and unzipping TCDB sequences', flush=True, file=stdout_file)
        for link in [tcdb_seqs,
                     tcdb_go,
                     tcdb_pfam,
                     ]:
            download_file(link, output_folder=self.mantis_paths['tcdb'], stdout_file=stdout_file)
        self.compile_tcdb_metadata()

        if not file_exists(self.mantis_paths['tcdb']+'tcdb.dmnd'):
            run_command(f'{DIAMOND_PATH} makedb --in ' + self.mantis_paths['tcdb'] + 'tcdb.faa -d '+self.mantis_paths['tcdb']+'tcdb', stdout_file=stdout_file)

#####################   NOG


    def check_completeness_NOGG(self, nogg_file, list_file_paths, force_download):
        if force_download:
            return False
        if not file_exists(nogg_file):
            return False
        nogg_hmms = set()
        nogt_hmms = set()
        with open(nogg_file) as file:
            line = file.readline()
            while line:
                hmm_name = line.split('\t')[0]
                nogg_hmms.add(hmm_name)
                line = file.readline()
        for f in list_file_paths:
            with open(f) as file:
                line = file.readline()
                while line:
                    hmm_name = line.split('\t')[0]
                    if hmm_name not in nogg_hmms:
                        return False
                    nogt_hmms.add(hmm_name)
                    line = file.readline()
        if not nogg_hmms == nogt_hmms:
            return False
        return True

    def compile_NOGG(self, force_download=False):
        # this is not feasible for all HMMs, so we just select the most general taxa aka domains
        if self.mantis_paths['NOG'][0:2] != 'NA':
            stdout_file = open(self.mantis_out, 'a+')
            target_annotation_file = add_slash(self.mantis_paths['NOG'] + 'NOGG') + 'NOGG_sql_annotations.tsv'
            target_merged_hmm = add_slash(self.mantis_paths['NOG'] + 'NOGG') + 'NOGG_merged.hmm'
            all_sql = set()
            all_hmm = set()
            for taxon_id in self.get_ncbi_domains():
                if taxon_id != 'NOGG':
                    if os.path.isdir(self.mantis_paths['NOG'] + taxon_id):
                        target_hmm_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_merged.hmm'
                        target_sql_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_sql_annotations.tsv'
                        all_sql.add(target_sql_file)
                        all_hmm.add(target_hmm_file)

            if not self.check_completeness_NOGG(target_annotation_file, all_sql,force_download) or\
                    not self.check_reference_exists('NOGG',force_download):
                if file_exists(self.mantis_paths['NOG'] + 'NOGG'): shutil.rmtree(self.mantis_paths['NOG'] + 'NOGG')
                Path(self.mantis_paths['NOG'] + 'NOGG').mkdir(parents=True, exist_ok=True)
            else:
                print('NOGG already compiled, skipping...', flush=True, file=stdout_file)
                return
            print_cyan('Compiling global NOG hmms ', flush=True, file=stdout_file)
            # now to merge all tshmms into one
            merge_redundant_profiles(output_file=target_merged_hmm, list_file_paths=all_hmm, stdout_file=stdout_file)
            merge_redundant_sql_annotations(output_file=target_annotation_file, list_file_paths=all_sql,stdout_file=stdout_file)
            run_command('hmmpress ' + target_merged_hmm, stdout_file=stdout_file)
            stdout_file.close()

    def download_NOGT(self, taxon_id, stdout_file=None):
        folder_path = add_slash(self.mantis_paths['NOG'] + taxon_id)
        Path(folder_path).mkdir(parents=True, exist_ok=True)
        target_annotation_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_annotations.tsv'
        target_merged_hmm = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_merged.hmm'
        eggnog_downloads_page = f'http://eggnog5.embl.de/download/latest/per_tax_level/{taxon_id}/'
        for file in ['_annotations.tsv.gz', '_hmms.tar.gz']:
            url = f'{eggnog_downloads_page}{taxon_id}{file}'
            download_file(url, output_folder=folder_path, stdout_file=stdout_file)
        if file_exists(f'{folder_path}profiles'): shutil.rmtree(f'{folder_path}profiles')
        uncompress_archive(source_filepath= f'{folder_path}{taxon_id}_hmms.tar.gz', extract_path= f'{folder_path}profiles',stdout_file=stdout_file, remove_source=True)
        uncompress_archive(source_filepath= f'{folder_path}{taxon_id}_annotations.tsv.gz', stdout_file=stdout_file,remove_source=True)
        # removing gz extension
        hmm_files = os.listdir(folder_path + 'profiles')
        for hmm_profile in hmm_files:
            if '.hmm' in hmm_profile: move_file(hmm_profile, hmm_profile.strip('.gz'))
        merge_profiles(f'{folder_path}profiles/{taxon_id}', f'{folder_path}{taxon_id}_merged.hmm',
                       stdout_file=stdout_file)
        if file_exists(f'{folder_path}profiles'):      shutil.rmtree(f'{folder_path}profiles')
        run_command(f'hmmpress {target_merged_hmm}', stdout_file=stdout_file)

    def download_and_unzip_eggnogdb(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['default']).mkdir(parents=True, exist_ok=True)
        if file_exists(self.mantis_paths['default'] + 'eggnog.db', force_download):
            print('eggnog.db already exists! Skipping...', flush=True, file=stdout_file)
            return
        url = 'http://eggnogdb.embl.de/download/emapperdb-'+self.get_latest_version_eggnog()+'/eggnog.db.gz'
        download_file(url, output_folder=self.mantis_paths['default'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['default'] + 'eggnog.db.gz', stdout_file=stdout_file,
                           remove_source=True)

    def unpack_NOG_sql(self, stdout_file=None):
        '''
        splitting into different folders for github pushing (without lfs restrictions)
        x=("NOGT1" "NOGT2" "NOGT3" "NOGT4" "NOGT5")
        c=0
        for f in *
        do
            mv "$f" "../${x[c]}/"
            c=$(( (c+1)%5 ))
        done

        split into 6 subfolders(5 NOGT + NOGG) which should be compressed with:
        export GZIP=-9
        then in the NOGT folder do:
        for i in *; do tar cvzf $i.tar.gz $i; done
        '''
        passed = True
        if self.mantis_paths['NOG'][0:2] == 'NA': return
        if self.mantis_nogt_tax:
            for tax_id in self.mantis_nogt_tax:
                if not file_exists(add_slash(self.mantis_paths['NOG'] + str(tax_id)) + f'{tax_id}_sql_annotations.tsv'): passed = False
            if passed:
                print('All chosen NOGT were already extracted! Skipping...', flush=True, file=stdout_file)
                return
        resources_path = self.mantis_paths['resources']
        NOG_sql_resources_folder = add_slash(resources_path + 'NOG_sql')
        NOGT_sql_folder = add_slash(self.mantis_paths['NOG'] + 'NOG_sql')
        Path(NOGT_sql_folder).mkdir(parents=True, exist_ok=True)
        Path(NOG_sql_resources_folder).mkdir(parents=True, exist_ok=True)


        if not file_exists(NOG_sql_resources_folder): return
        for compressed_sql in os.listdir(NOG_sql_resources_folder):
            if 'NOG' in compressed_sql and '.tar.gz' in compressed_sql:
                uncompress_archive(source_filepath=NOG_sql_resources_folder + compressed_sql,
                                   extract_path=NOGT_sql_folder, stdout_file=stdout_file, remove_source=False)
        for sql_folder in os.listdir(NOGT_sql_folder):
            # sql_folder - e.g. NOGT1
            if re.search('NOGT\d+', sql_folder) and 'tar.gz' not in sql_folder:
                for tax_sql in os.listdir(NOGT_sql_folder + sql_folder):
                    tax_id = tax_sql.split('_')[0]
                    passed = False
                    if self.mantis_nogt_tax:
                        if str(tax_id) in self.mantis_nogt_tax:
                            passed = True
                    else:
                        passed = True
                    if passed:
                        Path(add_slash(self.mantis_paths['NOG'] + tax_id)).mkdir(parents=True, exist_ok=True)
                        move_file(add_slash(NOGT_sql_folder + sql_folder) + tax_sql,
                                  add_slash(self.mantis_paths['NOG'] + tax_id) + tax_sql)
        # removing empty folders
        if file_exists(NOGT_sql_folder):
            shutil.rmtree(NOGT_sql_folder)

    ### Support functions for setting up queue to download NOGT hmms

    def get_latest_version_eggnog(self):
        url = 'http://eggnog5.embl.de/download/'
        webpage = None
        c = 0
        while not webpage and c <= 10:
            try:
                req = requests.get(url)
                webpage = req.text
            except:
                c += 1
        versions_id = re.findall('href="emapperdb-.*"', webpage)
        res=[]
        for v in versions_id:
            version_id=v.replace('href="emapperdb-','')
            version_id=version_id.replace('/"','')
            res.append(version_id)
        if res:
            return max(res)
        else:
            #hard coded just in case the user doesnt have access to internet, should be updated
            current_version='5.0.2'
            print(f'eggNOG database version retrieval failed, so returning hardcoded eggnog database version {current_version}, please change <current_version> if a newer version is available')
            return current_version

    def get_taxon_ids_eggNOG(self):
        #this will get all available taxon ids from the eggnog webpage
        url = 'http://eggnog5.embl.de/download/latest/per_tax_level/'
        webpage = None
        c = 0
        while not webpage and c <= 10:
            try:
                req = requests.get(url)
                webpage = req.text
            except:
                c += 1
        if isinstance(webpage,str):
            taxons_search = re.findall('href="\d+/"', webpage)
            taxons = [re.search('\d+', i).group() for i in taxons_search]
        else: taxons=[]
        #if no connection is established, we just return the local taxon ids (will be empty if no connection is available during setup)
        if not taxons: return self.get_local_ref_taxon_ids('NOG')
        return taxons

    def get_taxon_for_queue_NOGT(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['NOG']).mkdir(parents=True, exist_ok=True)
        # we will download tax specific hmms, which can be used for more specific hmmscans
        if not file_exists(self.mantis_paths['NOG'] + 'readme.md'):
            with open(self.mantis_paths['NOG'] + 'readme.md', 'w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write(f'These hmms were downloaded on {datetime_str} from:\nhttp://eggnog5.embl.de/download/latest/per_tax_level/')
        taxon_ids = self.get_ref_taxon_ids('NOG')
        res = []
        for taxon_id in taxon_ids:
            if taxon_id != '1':
                target_annotation_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + taxon_id + '_annotations.tsv'
                if self.check_reference_exists('NOG', taxon_id, force_download) and file_exists(target_annotation_file,force_download):
                    hmm_path = get_ref_in_folder(self.mantis_paths['NOG'] + taxon_id)
                    profile_count = get_hmm_profile_count(hmm_path)
                    annotations_count = len(self.get_hmms_annotation_file(target_annotation_file, 1))
                    # should be the same but some HMMs dont have an annotation (probably an error from NOG)
                    if profile_count in range(annotations_count - 10, annotations_count + 10 + 1):
                        print(f'Tax {taxon_id} hmm and annotations.tsv files already exist, and they were merged correctly. Skipping setup...',flush=True, file=stdout_file)
                    else:
                        print(f'Tax {taxon_id} hmm already exists but was not merged correctly!', flush=True,file=stdout_file)
                        res.append(taxon_id)
                else:
                    res.append(taxon_id)
        print('Will compile data for the following NOGT:\n' + ','.join(res), flush=True, file=stdout_file)
        return res

    ### Support functions for setting up queue for split hmms

    def get_taxon_for_queue_NOG_split_hmms(self):
        res = []
        stdout_file = open(self.mantis_out, 'a+')
        if self.mantis_paths['NOG'][0:2] != 'NA':
            for taxon_id in self.get_ref_taxon_ids('NOG'):
                if taxon_id != '1' and os.path.isdir(self.mantis_paths['NOG'] + taxon_id):
                    hmm_path = get_ref_in_folder(self.mantis_paths['NOG'] + taxon_id)
                    print('Checking NOG for splitting:', hmm_path, flush=True, file=stdout_file)
                    if 'chunks' not in os.listdir(self.mantis_paths['NOG'] + taxon_id):
                        profile_count = get_hmm_profile_count(hmm_path)
                        if profile_count > self.hmm_chunk_size:
                            res.append(hmm_path)
                        else:
                            print(f'NOG already split: {hmm_path} {taxon_id}', flush=True, file=stdout_file)
                    else:
                        print(f'NOG already split: {hmm_path} {taxon_id}', flush=True, file=stdout_file)
        stdout_file.close()
        return res

    def get_taxon_for_queue_ncbi_split_hmms(self):
        res = []
        stdout_file = open(self.mantis_out, 'a+')
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            for taxon_id in self.get_ref_taxon_ids('NCBI'):
                if os.path.isdir(self.mantis_paths['NCBI'] + taxon_id):
                    hmm_path = get_ref_in_folder(self.mantis_paths['NCBI'] + taxon_id)
                    print(f'Checking NCBI for splitting: {hmm_path} {taxon_id}', flush=True, file=stdout_file)
                    if 'chunks' not in os.listdir(self.mantis_paths['NCBI'] + taxon_id):
                        profile_count = get_hmm_profile_count(hmm_path)
                        if profile_count > self.hmm_chunk_size:
                            res.append(hmm_path)
                        else:
                            print(f'NCBI already split: {hmm_path} {taxon_id}', flush=True, file=stdout_file)
                    else:
                        print(f'NCBI already split: {hmm_path} {taxon_id}', flush=True, file=stdout_file)

        stdout_file.close()
        return res


    def get_hmm_for_queue_split_hmms(self):
        hmms_list = self.compile_refs_list()
        res = []
        for hmm in hmms_list:
            hmm_folder = get_folder(hmm)
            if 'chunks' not in os.listdir(hmm_folder):
                hmm_profile_count = get_hmm_profile_count(hmm)
                if hmm_profile_count > self.hmm_chunk_size:
                    res.append(hmm)
        return res

    ### Support functions for splitting hmms

    def split_hmm_into_chunks(self, hmm_path, stdout_path=None):
        # we dont really load balance the profiles based on profile length because we would need to keep too much data in memory
        # so we just balance out the number of profiles
        stdout_file = open(stdout_path, 'a+')
        print(f'Counting profiles in {hmm_path}', flush=True, file=stdout_file)
        profile_count = get_hmm_profile_count(hmm_path, stdout=stdout_file)
        print(f'{hmm_path} has {profile_count} profiles', flush=True, file=stdout_file)
        hmm_folder = get_folder(hmm_path)
        hmm_chunks_folder = f'{hmm_folder}chunks/'
        if file_exists(hmm_chunks_folder):
            shutil.rmtree(hmm_chunks_folder)
        Path(hmm_chunks_folder).mkdir(parents=True, exist_ok=True)
        print('Load balancing chunks', flush=True, file=stdout_file)
        load_balanced_chunks = chunk_generator_load_balanced([i for i in range(profile_count)],
                                                             get_hmm_chunk_size(total_profiles=profile_count,
                                                                                current_chunk_size=self.hmm_chunk_size,
                                                                                max_chunks=100), time_limit=None)
        print(f'Splitting {hmm_path} into {len(load_balanced_chunks)} chunks.', flush=True,file=stdout_file)
        chunks_name = get_path_level(hmm_path, remove_extension=True) + '_chunk_'
        with open(hmm_path) as hmm_file:
            for chunk_i in range(len(load_balanced_chunks)):
                chunk_file_path = f'{hmm_chunks_folder}{chunks_name}{chunk_i}.hmm'
                with open(chunk_file_path, 'w+') as chunk_file:
                    for profile_i in range(len(load_balanced_chunks[chunk_i])):
                        profile = ''.join(read_profile(hmm_file))
                        chunk_file.write(profile)
        chunks_dir = os.listdir(hmm_chunks_folder)
        for chunk in chunks_dir:
            run_command(f'hmmpress {hmm_chunks_folder}{chunk}', stdout_file=stdout_file)
        stdout_file.close()


    def press_custom_hmms(self, hmm_path, stdout_path=None):
        stdout_file = open(stdout_path, 'a+')
        hmm_folder = add_slash(SPLITTER.join(hmm_path.split(SPLITTER)[:-1]))
        old_files=['h3f', 'h3i', 'h3m', 'h3p']
        for inner_file in os.listdir(hmm_folder):
            inner_file_ending=inner_file.split('.')[-1]
            if inner_file_ending in old_files:
                os.remove(hmm_folder+inner_file)
        run_command(f'hmmpress {hmm_path}', stdout_file=stdout_file)
        stdout_file.close()


    ### Support functions for extracting metadata

    def get_hmms_annotation_file(self, annotations_file, hmm_col):
        res = set()
        if not file_exists(annotations_file): return res
        with open(annotations_file, 'r') as file:
            line = file.readline()
            while line:
                line = line.split('\t')
                res.add(line[hmm_col])
                line = file.readline()
        return res

    def get_free_text_NOG(self, hmm, hmm_annotations_path, description_col):
        with open(hmm_annotations_path, 'r') as file:
            line = file.readline()
            while line:
                line = line.split('\t')
                if line[1] == hmm:
                    description = line[description_col].strip('\n').replace('NA', '')
                    return description
                line = file.readline()
        return ''

    def pfam_id_to_acc(self):
        res = {}
        pfam_metadata = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
        download_file(pfam_metadata, output_folder=self.mantis_paths['resources'])
        uncompress_archive(source_filepath=self.mantis_paths['resources'] + 'Pfam-A.hmm.dat.gz',remove_source=True)
        with open(self.mantis_paths['resources'] + 'Pfam-A.hmm.dat') as pfam_dat_file:
            line = pfam_dat_file.readline()
            while line:
                line = line.strip('\n').split('   ')
                if len(line) == 2:
                    row_header, row_description = line
                    if row_header == '#=GF ID':
                        hmm = str(row_description)
                    elif row_header == '#=GF AC':
                        pfam_accession = str(row_description)
                        pfam_accession = pfam_accession.split('.')[0].strip()
                        res[hmm] = pfam_accession
                line = pfam_dat_file.readline()
        remove_file(self.mantis_paths['resources'] + 'Pfam-A.hmm.dat')
        return res

    def get_metadata_hmms_thread(self, hmm, taxon_id, hmm_annotations_path, eggnog_db):
        '''
        some OG HMMs contain a lot of protein sequences... when that happens we might get a ton of IDs.
        the current eggnog mapper annotation will be different since they now use diamond which instead of matching against
        a whole hmm (multiple sequences), they match against a single sequence. So with mantis the annotation will contain all the annotations from the OG sequences

        notes on the sql command:
        >SELECT< will retrieve all the metadata (gos,kegg_ec,kegg_ko,kegg_pathway...)
        >LEFT JOIN< each db (i.e. gene_ontology,kegg,bigg) to get all seqs (from each db) where the db.name matches with the eggnog.name
        this doesnt guarantee that each db will have an entry, but if 1 or more do, then db1.name = to db2.name
        >WHERE< will search only for the seqs that belong to the group for the hmm name in the TSHMM
        final >AND< is just to restrict results to only retrieve sequences with metadata

        '''

        codes_to_exclude = ['IEA','ND']
        #this will convert protein names to pfam ids (which is typically what is used with mantis)
        pfam_id_to_acc=self.pfam_id_to_acc()
        taxon_query = '@' + str(taxon_id)
        connection = sqlite3.connect(eggnog_db)
        sql_command = 'SELECT eggnog.groups,pfam.pfam,gene_ontology.gos,kegg.ec, kegg.ko, kegg.pathway, kegg.module, kegg.reaction, kegg.rclass, kegg.brite, kegg.tc, kegg.cazy, bigg.reaction FROM eggnog ' \
                      ' LEFT JOIN gene_ontology on gene_ontology.name = eggnog.name ' \
                      ' LEFT JOIN kegg on kegg.name = eggnog.name ' \
                      ' LEFT JOIN bigg on bigg.name = eggnog.name ' \
                      ' LEFT JOIN pfam on pfam.name = eggnog.name ' \
                      ' WHERE eggnog.groups LIKE "%' + hmm + taxon_query + '%"' \
                       'AND (gene_ontology.name IS NOT NULL OR kegg.name IS NOT NULL OR bigg.name IS NOT NULL OR pfam.name IS NOT NULL);'
        cursor = connection.cursor()
        cursor.execute(sql_command)
        rows = cursor.fetchall()
        res = {'go': set(), 'enzyme_ec': set(), 'kegg_ko': set(), 'kegg_pathway': set(), 'kegg_module': set(),
               'kegg_reaction': set(), 'kegg_rclass': set(), 'kegg_brite': set(), 'kegg_cazy': set(),
               'bigg_reaction': set(), 'description': set(), 'pfam':set()}
        for row in rows:
            hmms, pfam_ids,gene_ontology_gos, kegg_ec, kegg_ko, kegg_pathway, kegg_module, kegg_reaction, kegg_rclass, kegg_brite, kegg_tc, kegg_cazy, bigg_reaction = row
            if kegg_ko: kegg_ko = kegg_ko.replace('ko:', '')
            for query_hmm_taxon in hmms.split(','):
                query_hmm, query_taxon = query_hmm_taxon.split('@')
                if int(taxon_id) == int(query_taxon):
                    passed_test = True
                else:
                    passed_test = False
                if passed_test:
                    if gene_ontology_gos:
                        gene_ontology_gos_copy = str(gene_ontology_gos)
                        gene_ontology_gos_copy = gene_ontology_gos_copy.split(',')
                        for go_group in gene_ontology_gos_copy:
                            if '|' not in go_group: print(sql_command, '\n', go_group)
                            _, go, evidence_code = go_group.split('|')
                            if evidence_code not in codes_to_exclude:
                                res['go'].add(go.split(':')[1])
                    if pfam_ids:
                        temp_pfam_ids=pfam_ids.split(',')
                        for tpi in temp_pfam_ids:
                            if tpi in pfam_id_to_acc:
                                res['pfam'].add(pfam_id_to_acc[tpi])
                    if kegg_ec: res['enzyme_ec'].update(kegg_ec.split(','))
                    if kegg_ko: res['kegg_ko'].update(kegg_ko.split(','))
                    if kegg_pathway: res['kegg_pathway'].update(kegg_pathway.split(','))
                    if kegg_module: res['kegg_module'].update(kegg_module.split(','))
                    if kegg_reaction: res['kegg_reaction'].update(kegg_reaction.split(','))
                    if kegg_rclass: res['kegg_rclass'].update(kegg_rclass.split(','))
                    if kegg_brite: res['kegg_brite'].update(kegg_brite.split(','))
                    if kegg_cazy: res['kegg_cazy'].update(kegg_cazy.split(','))
                    if bigg_reaction: res['bigg_reaction'].update(bigg_reaction.split(','))
        description = self.get_free_text_NOG(hmm=hmm, hmm_annotations_path=hmm_annotations_path, description_col=3)
        if description:   res['description'].add(description)
        return [hmm, res]

    def launch_pool(self, n_processes, function_to_process, function_args):
        pool = IO_Pool(n_processes)
        res = pool.starmap(function_to_process, function_args)
        pool.close()
        pool.join()
        return res

    def get_metadata_hmms(self, target_annotation_file, target_sql_file, taxon_id=None, stdout_path=None):
        max_threads = 100
        chunk_size = 500
        eggnog_db_path = self.mantis_paths['default'] + "eggnog.db"
        hmm_list = self.get_hmms_annotation_file(target_annotation_file, hmm_col=1)
        stdout_file = open(stdout_path, 'a+')

        if file_exists(target_sql_file):
            already_compiled_hmms_list = self.get_hmms_annotation_file(target_sql_file, hmm_col=0)
            if len(already_compiled_hmms_list) == len(hmm_list):
                print('Target sql annotation file already finished, skipping!', target_sql_file, flush=True,
                      file=stdout_file)
                return
            else:
                for hmm in already_compiled_hmms_list:
                    if hmm in hmm_list:
                        hmm_list.remove(hmm)
                if not hmm_list:
                    print('Target sql annotation file already complete, skipping!', target_sql_file, flush=True,
                          file=stdout_file)
                    return
                else:
                    print('Target sql annotation file already exists but is incomplete, continuing...', target_sql_file,
                          flush=True, file=stdout_file)

        print('Starting metadata extraction into', target_sql_file, flush=True, file=stdout_file)
        # http://geneontology.org/docs/guide-go-evidence-codes/
        while not file_exists(self.mantis_paths['default'] + "eggnog.db"): sleep(5)
        print(f'Exporting metadata for NOG {taxon_id}', flush=True, file=stdout_file)
        hmm_annotations_path = self.mantis_paths['NOG'] + f'{taxon_id}/{taxon_id}_annotations.tsv'
        pool_args = []
        for hmm in hmm_list:
            pool_args.append([hmm, taxon_id, hmm_annotations_path, eggnog_db_path])
        if len(pool_args) < max_threads:
            n_threads = len(pool_args)
        else:
            n_threads = max_threads
        n_threads = 5
        # to avoid keeping everything in memory, we chunk it
        chunks = chunk_generator(pool_args, chunk_size)
        for chunk in chunks:
            hmm_metadata = self.launch_pool(n_threads, self.get_metadata_hmms_thread, chunk)
            with open(target_sql_file, "a+") as file:
                for hmm_links in hmm_metadata:
                    link_line = ''
                    hmm, metadata = hmm_links
                    for link_type in metadata:
                        for inner_link in metadata[link_type]:
                            link_line += '\t' + link_type + ':' + inner_link
                    file.write(hmm + '\t|' + link_line + '\n')
        print(f'Finished exporting metadata for NOGT {taxon_id}', flush=True, file=stdout_file)
        stdout_file.close()

    ###### NOG diamond (not implemented yet)

    def download_and_unzip_eggnogdmnd(self, force_download=False, stdout_file=None):
        folder_path = add_slash(self.mantis_paths['NOG'])
        if file_exists(folder_path + 'eggnog_proteins.dmnd', force_download):
            print('eggnog_proteins.dmnd already exists! Skipping...', flush=True, file=stdout_file)
            return
        url = 'http://eggnogdb.embl.de/download/emapperdb-'+self.get_latest_version_eggnog()+'/eggnog_proteins.dmnd.gz'
        download_file(url, output_folder=folder_path, stdout_file=stdout_file)
        uncompress_archive(source_filepath= f'{folder_path}eggnog_proteins.dmnd.gz', extract_path= folder_path,stdout_file=stdout_file, remove_source=False)



    def download_NOGT_diamond(self, stdout_file=None):
        folder_path = add_slash(self.mantis_paths['NOG'])
        Path(folder_path).mkdir(parents=True, exist_ok=True)
        self.download_and_unzip_eggnogdb()
        self.download_and_unzip_eggnogdmnd()
        diamond_db=f'{folder_path}eggnog_proteins'
        extract_seqs_command= f'{DIAMOND_PATH} getseq -d {diamond_db} > {folder_path}eggnog_seqs.faa'
        eggnog_proteins_path=f'{folder_path}eggnog_seqs.faa'
        print('Extracting seqs from Diamond database')
        #run_command(extract_seqs_command,join_command=True,shell=True)
        #self.get_metadata_diamond(eggnog_proteins_path)
        print('Making taxa specific dmnds')
        for taxon_id in os.listdir(self.mantis_paths['NOG']):
            taxon_folder = self.mantis_paths['NOG'] + taxon_id + SPLITTER
            if re.search('\d+',taxon_id) and os.path.isdir(taxon_folder):
                if taxon_id!='1':
                    taxon_folder = self.mantis_paths['NOG'] + taxon_id + SPLITTER
                    taxon_dmnd=f'{taxon_folder}{taxon_id}'
                    taxon_faa=f'{taxon_folder}{taxon_id}.faa'
                    if not file_exists(self.mantis_paths['tcdb'] + taxon_dmnd):
                        print(f'Creating dmnd with {taxon_faa}')
                        run_command(f'{DIAMOND_PATH} makedb --in {taxon_faa} -d {taxon_dmnd}', stdout_file=stdout_file)


    def export_metadata_and_seqs_nog_diamond(self,all_seqs,sequence,taxons,metadata):
        for taxon_id in taxons:
            if taxon_id != '1':
                taxon_folder=self.mantis_paths['NOG'] + taxon_id + SPLITTER
                if not file_exists(taxon_folder): Path(taxon_folder).mkdir(parents=True, exist_ok=True)
                if sequence in all_seqs:
                    with open(f'{taxon_folder}{taxon_id}.faa','a+') as seq_file:
                        seq_file.write(f'>{sequence}\n{all_seqs[sequence]}\n')
                    with open(f'{taxon_folder}metadata.tsv','a+') as metadata_file:
                        link_line=f'{sequence}\t|\t'
                        for link_type in metadata:
                            for inner_link in metadata[link_type]:
                                link_line += '\t' + link_type + ':' + inner_link
                        metadata_file.write(f'{link_line}\n')

    def get_metadata_diamond(self,eggnog_proteins_path,stdout_path=None):
        print('Reading eggNOG protein fasta')
        eggnog_proteins=read_protein_fasta(eggnog_proteins_path)
        eggnog_db_path = self.mantis_paths['default'] + "eggnog.db"
        target_sql_file=self.mantis_paths['NOG']+'metadata.tsv'
        #stdout_file = open(stdout_path, 'a+')
        print('Generating metadata tsv')
        self.generate_metadata_diamond(eggnog_db_path,target_sql_file)
        print('Generating taxon specific metadata tsvs')
        with open(target_sql_file) as file:
            line=file.readline()
            while line:
                current_metadata={}
                line = line.strip('\n')
                line = line.split('\t')
                current_seq = line[0]
                current_taxons=set()
                annotations = line[1:]
                for link in annotations:
                    if link:
                        temp_link = link.split(':')
                        link_type = temp_link[0]
                        if link_type == 'kegg_cazy': link_type = 'cazy'
                        if link_type == 'kegg_ec': link_type = 'enzyme_ec'
                        link_text = ':'.join(temp_link[1:])
                        link_text = link_text.strip()
                        if link_type=='taxon': current_taxons.add(link_text)
                        else:
                            if link_type not in current_metadata: current_metadata[link_type] = set()
                            if link_type == 'description' and link_text == 'NA':   link_text = ''
                            if link_text and link_type == 'description': self.get_common_links(link_text, res=current_metadata)
                            current_metadata[link_type].add(link_text)
                self.export_metadata_and_seqs_nog_diamond(all_seqs=eggnog_proteins,sequence=current_seq,taxons=current_taxons,metadata=current_metadata)
                line = file.readline()





    #we generate a metadata.tsv with the functional annotation of each sequence
    def generate_metadata_diamond(self,eggnog_db,target_sql_file):
        codes_to_exclude = ['IEA','ND']
        #this will convert protein names to pfam ids (which is typically what is used with mantis)
        pfam_id_to_acc=self.pfam_id_to_acc()
        print('Querying SQL')
        connection = sqlite3.connect(eggnog_db)
        sql_command = 'SELECT eggnog.name, eggnog.groups,pfam.pfam,gene_ontology.gos,kegg.ec, kegg.ko, kegg.pathway, kegg.module, kegg.reaction, kegg.rclass, kegg.brite, kegg.tc, kegg.cazy, bigg.reaction FROM eggnog ' \
                      ' LEFT JOIN gene_ontology on gene_ontology.name = eggnog.name ' \
                      ' LEFT JOIN kegg on kegg.name = eggnog.name ' \
                      ' LEFT JOIN bigg on bigg.name = eggnog.name ' \
                      ' LEFT JOIN pfam on pfam.name = eggnog.name '
        print(f'SQL command:\n{sql_command}')
        cursor = connection.cursor()
        cursor.execute(sql_command)
        rows = cursor.fetchall()
        print('Reading SQL results')
        with open(target_sql_file,'w+') as file:
            for row in rows:
                res = {'go': set(), 'enzyme_ec': set(), 'kegg_ko': set(), 'kegg_pathway': set(), 'kegg_module': set(),
                       'kegg_reaction': set(), 'kegg_rclass': set(), 'kegg_brite': set(), 'kegg_cazy': set(),
                       'bigg_reaction': set(), 'description': set(), 'pfam': set()}
                seq_name,hmms, pfam_ids,gene_ontology_gos, kegg_ec, kegg_ko, kegg_pathway, kegg_module, kegg_reaction, kegg_rclass, kegg_brite, kegg_tc, kegg_cazy, bigg_reaction = row
                if kegg_ko: kegg_ko = kegg_ko.replace('ko:', '')
                taxon_ids=set()
                for query_hmm_taxon in hmms.split(','):
                    query_hmm, query_taxon = query_hmm_taxon.split('@')
                    taxon_ids.add(query_taxon)
                if gene_ontology_gos:
                    gene_ontology_gos_copy = str(gene_ontology_gos)
                    gene_ontology_gos_copy = gene_ontology_gos_copy.split(',')
                    for go_group in gene_ontology_gos_copy:
                        if '|' not in go_group: print(sql_command, '\n', go_group)
                        _, go, evidence_code = go_group.split('|')
                        if evidence_code not in codes_to_exclude:
                            res['go'].add(go.split(':')[1])
                if pfam_ids:
                    temp_pfam_ids=pfam_ids.split(',')
                    for tpi in temp_pfam_ids:
                        if tpi in pfam_id_to_acc:
                            res['pfam'].add(pfam_id_to_acc[tpi])
                if kegg_ec: res['enzyme_ec'].update(kegg_ec.split(','))
                if kegg_ko: res['kegg_ko'].update(kegg_ko.split(','))
                if kegg_pathway: res['kegg_pathway'].update(kegg_pathway.split(','))
                if kegg_module: res['kegg_module'].update(kegg_module.split(','))
                if kegg_reaction: res['kegg_reaction'].update(kegg_reaction.split(','))
                if kegg_rclass: res['kegg_rclass'].update(kegg_rclass.split(','))
                if kegg_brite: res['kegg_brite'].update(kegg_brite.split(','))
                if kegg_cazy: res['kegg_cazy'].update(kegg_cazy.split(','))
                if bigg_reaction: res['bigg_reaction'].update(bigg_reaction.split(','))
                link_line = f'{seq_name}'
                for taxon_id in taxon_ids:
                    link_line+=f'\ttaxon:{taxon_id}'
                for link_type in res:
                    for inner_link in res[link_type]:
                        link_line += '\t' + link_type + ':' + inner_link
                file.write(f'{link_line}\n')

