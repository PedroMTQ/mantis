try:
    from source.Exceptions import *
    from source.utils import *
    from source.MANTIS_NLP import MANTIS_NLP
except:
    from Exceptions import *
    from utils import *
    from MANTIS_NLP import MANTIS_NLP


class MANTIS_DB(MANTIS_NLP):

    ###Main function
    @timeit_class
    def setup_databases(self, force_download=False):
        if not cython_compiled():
            compile_cython()
        if not unifunc_downloaded():
            download_unifunc()
        self.output_folder = f'{MANTIS_FOLDER}setup_databases/'
        self.mantis_out = f'{self.output_folder}Mantis.out'
        if file_exists(self.mantis_out):
            os.remove(self.mantis_out)
        Path(self.mantis_paths['default']).mkdir(parents=True, exist_ok=True)
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        dbs_list = [
            'ncbi_tax' if self.mantis_paths['ncbi_tax'][0:2] != 'NA' else None,
            'pfam' if self.mantis_paths['pfam'][0:2] != 'NA' else None,
            'kofam' if self.mantis_paths['kofam'][0:2] != 'NA' else None,
            'tigrfam' if self.mantis_paths['tigrfam'][0:2] != 'NA' else None,
            'NCBI' if self.mantis_paths['NCBI'][0:2] != 'NA' else None,
        ]
        # DOWNLOADING
        self.prepare_queue_setup_databases(dbs_list, force_download)
        # for unzipping tax specific hmms
        self.prepare_queue_setup_databases_tax(force_download)
        worker_count = estimate_number_workers_setup_database(len(self.queue))
        if worker_count: print(f'Database will be setup with {worker_count} workers!', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_setup_databases, worker_count)
        print_cyan('Finished downloading all data!', flush=True, file=self.redirect_verbose)
        # METADATA
        print_cyan('Will now extract metadata from NOG!', flush=True, file=self.redirect_verbose)
        self.prepare_queue_extract_metadata()
        worker_count = estimate_number_workers_setup_database(len(self.queue), minimum_jobs_per_worker=2)
        if worker_count:print(f'Metadata will be extracted with {worker_count} workers!', flush=True,file=self.redirect_verbose)
        self.processes_handler(self.worker_extract_NOG_metadata, worker_count)
        print('NOGG will now be compiled',flush=True,file=self.redirect_verbose)
        self.compile_NOGG(force_download)
        # SPLITTING
        if self.hmm_chunk_size:
            print_cyan('Will now split data into chunks!', flush=True, file=self.redirect_verbose)
            self.prepare_queue_split_hmms()
            worker_count = estimate_number_workers_setup_database(len(self.queue))
            print(f'Database will be split with {worker_count} workers!', flush=True,
                  file=self.redirect_verbose)
            self.processes_handler(self.worker_split_hmms, worker_count)
        print('Preparing NLP Resources!', flush=True, file=self.redirect_verbose)
        MANTIS_NLP.__init__(self)
        print_cyan('Finished setting up databases!', flush=True, file=self.redirect_verbose)

    #this is a standalone execution, this is just for the developers to update the NOG tars and then upload them on Github
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

        eggnog_downloads_page = 'http://eggnog5.embl.de/download/latest/per_tax_level/'
        taxon_ids = self.get_taxon_ids_NOGT(eggnog_downloads_page)

        url = 'http://eggnogdb.embl.de/download/emapperdb-5.0.1/eggnog.db.gz'
        download_file(url, output_folder=self.mantis_paths['default'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['default'] + 'eggnog.db.gz', stdout_file=stdout_file,
                           remove_source=False)

        for t in taxon_ids:
            Path(self.mantis_paths['NOG'] + t).mkdir(parents=True, exist_ok=True)
            url = f'http://eggnog5.embl.de/download/latest/per_tax_level/{t}/{t}_annotations.tsv.gz'
            download_file(url, output_folder=add_slash(self.mantis_paths['NOG'] + t), stdout_file=stdout_file)
            uncompress_archive(source_filepath=add_slash(self.mantis_paths['NOG'] + t) + f'{t}_annotations.tsv.gz',
                               stdout_file=stdout_file, remove_source=True)

        self.prepare_queue_extract_metadata()
        worker_count = estimate_number_workers_setup_database(len(self.queue), minimum_jobs_per_worker=2)
        print(f'Metadata will be extracted with {worker_count} workers!', flush=True,
              file=self.redirect_verbose)
        print(f'The following taxons will be extracted:{taxon_ids}', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_extract_NOG_metadata, worker_count)

    ###Filling queue with jobs

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
        print('Checking which hmms need to be split, this may take a while...', flush=True, file=self.redirect_verbose)
        taxon_nog_hmms = self.get_taxon_for_queue_NOG_split_hmms()
        taxon_ncbi_hmms = self.get_taxon_for_queue_ncbi_split_hmms()
        general_hmms = self.get_hmm_for_queue_split_hmms()
        print('Will split: ', [get_path_level(i) for i in taxon_nog_hmms + taxon_ncbi_hmms + general_hmms], flush=True,
              file=self.redirect_verbose)
        for hmm_path in taxon_nog_hmms + taxon_ncbi_hmms + general_hmms:
            self.queue.append([hmm_path, self.mantis_out])

    def prepare_queue_extract_metadata(self):
        stdout_file = open(self.mantis_out, 'a+')
        print('Checking which NOGs we need to extract metadata from', flush=True, file=stdout_file)
        if self.mantis_paths['NOG'][0:2] != 'NA':
            self.unpack_NOG_sql(stdout_file=stdout_file)
        to_download=False
        if self.mantis_paths['NOG'][0:2] != 'NA':
            for taxon_id in os.listdir(self.mantis_paths['NOG']):
                if os.path.isdir(self.mantis_paths['NOG'] + taxon_id) and taxon_id!='NOGG':
                    target_annotation_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_annotations.tsv'
                    target_sql_file = add_slash(self.mantis_paths['NOG'] + taxon_id) +  f'{taxon_id}_sql_annotations.tsv'
                    if file_exists(target_sql_file):
                        if len(self.get_hmms_annotation_file(target_sql_file, hmm_col=0)) != len(
                                self.get_hmms_annotation_file(target_annotation_file, hmm_col=1)):
                            to_download=True
                            self.queue.append([target_sql_file, target_annotation_file, taxon_id, self.mantis_out])
                        else:
                            print(f'Skipping metadata extraction for NOGT {taxon_id}', flush=True, file=stdout_file)
                    else:
                        to_download = True
                        self.queue.append([target_sql_file, target_annotation_file, taxon_id, self.mantis_out])
        if to_download: self.download_and_unzip_eggnogdb()
        stdout_file.close()

    ###Executing queue jobs

    def worker_setup_databases(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            if len(record) == 3:
                database, force_download, stdout_path = record
                taxon_id = None
            else:
                database, force_download, taxon_id, stdout_path = record
            self.run_download_and_unzip(database, force_download, taxon_id, stdout_path)
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

    def worker_extract_NOG_metadata(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            target_sql_file, target_annotation_file, taxon_id, stdout_path = record
            self.get_metadata_hmms(target_annotation_file=target_annotation_file, target_sql_file=target_sql_file,
                                   taxon_id=taxon_id, stdout_path=stdout_path)

    ###Main download function

    def run_download_and_unzip(self, database, force_download=False, taxon_id=None, stdout_path=None):
        stdout_file = open(stdout_path, 'a+')
        if database == 'ncbi_tax':      self.download_and_unzip_tax_lineage_dmp(force_download=force_download, stdout_file=stdout_file)
        elif database == 'pfam':        self.download_and_unzip_pfam_hmm(force_download=force_download, stdout_file=stdout_file)
        elif database == 'kofam':       self.download_and_unzip_kofam_hmm(force_download=force_download, stdout_file=stdout_file)
        elif database == 'tigrfam':     self.download_and_unzip_tigrfam_hmm(force_download=force_download, stdout_file=stdout_file)
        elif database == 'NCBI':        self.download_and_unzip_NCBI(force_download=force_download, stdout_file=stdout_file)
        elif database == 'NOG':         self.download_and_unzip_NOGT(taxon_id=taxon_id, stdout_file=stdout_file)
        if taxon_id:
            print(f'Finished downloading {database} with taxon {taxon_id}', flush=True, file=stdout_file)
        else:
            print(f'Finished downloading {database}', flush=True, file=stdout_file)
        stdout_file.close()

    def download_and_unzip_tax_lineage_dmp(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['ncbi_tax']).mkdir(parents=True, exist_ok=True)
        if file_exists(self.mantis_paths['ncbi_tax'] + 'taxidlineage.dmp', force_download):
            print('NCBI taxonomic lineage file already exists! Skipping...', flush=True, file=stdout_file)
            return
        try:
            os.remove(self.mantis_paths['ncbi_tax'] + 'taxidlineage.dmp')
        except:
            pass
        url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz'
        if not file_exists(self.mantis_paths['ncbi_tax'] + 'readme.md'):
            with open(self.mantis_paths['ncbi_tax'] + 'readme.md', 'w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write(f'This file was downloaded on {datetime_str} from:\n{url}\nIt is used to traceback organism lineage for taxonomically resolved annotations')
        print_cyan('Downloading and unzipping ncbi taxon lineage dump', flush=True, file=stdout_file)
        download_file(url, output_folder=self.mantis_paths['ncbi_tax'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['ncbi_tax'] + 'new_taxdump.tar.gz',
                           extract_path=self.mantis_paths['ncbi_tax'], stdout_file=stdout_file, remove_source=True)
        # deleting the extra files that come with the tar dump
        files_dir = os.listdir(self.mantis_paths['ncbi_tax'])
        for file in files_dir:
            if file != 'taxidlineage.dmp':
                os.remove(self.mantis_paths['ncbi_tax'] + file)


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
                    res[pfam_id] = {'go': set(go_ids), 'description': set(go_description)}
                line = pfam2go_file.readline()
        return res

    def is_good_description(self, hmm, row_description):
        temp = [i.lower() for i in row_description.split()]
        if hmm.lower() in temp and 'protein' in temp and len(temp) == 2:
            return False
        if re.search(' [uU]nknown [Ff]unction', row_description): return False
        return True

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

    def build_pfam_line(self, hmm, metadata):
        link_line = ''
        pfam_ids = set(metadata['pfam'])
        metadata['pfam'] = set()
        for p_id in pfam_ids:
            metadata['pfam'].add(p_id.split('.')[0])
        for link_type in metadata:
            # print('link',link,flush=True)
            for inner_link in metadata[link_type]:
                link_line += f'\t{link_type}:{inner_link}'
        return hmm + f'\t|{link_line}\n'

    def compile_pfam_metadata(self):
        pfam2go = self.read_pfam2go()
        with open(self.mantis_paths['pfam'] + 'pfam_metadata.tsv', 'w+') as pfam_metadata_file:
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
                            if hmm in pfam2go:
                                current_metadata = pfam2go[hmm]
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

    def download_and_unzip_pfam_hmm(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['pfam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists(self.mantis_paths['pfam'] + 'Pfam-A.hmm', force_download) and \
                file_exists(self.mantis_paths['pfam'] + 'Pfam-A.hmm.dat', force_download) and \
                file_exists(self.mantis_paths['pfam'] + 'pfam_metadata.tsv', force_download) and \
                file_exists(self.mantis_paths['pfam'] + 'pfam2go', force_download):
            print('Pfam hmm already exists! Skipping...', flush=True, file=stdout_file)
            return
        pfam_hmm = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
        pfam_metadata = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
        pfam2go = 'http://current.geneontology.org/ontology/external2go/pfam2go'
        if not file_exists(self.mantis_paths['pfam'] + 'readme.md'):
            with open(self.mantis_paths['pfam'] + 'readme.md', 'w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write(f'This hmm was downloaded on {datetime_str} from:\n{pfam_hmm}\nMetadata was downloaded from:\n{pfam_metadata}\nPfam2GO was downloaded from:\n{pfam2go}')
        print_cyan('Downloading and unzipping Pfam hmms ', flush=True, file=stdout_file)
        to_download = []
        to_unzip = []
        if not file_exists(self.mantis_paths['pfam'] + 'Pfam-A.hmm.dat', force_download):
            to_download.append(pfam_hmm)
            to_unzip.append('Pfam-A.hmm.gz')
        if not file_exists(self.mantis_paths['pfam'] + 'pfam_metadata.tsv', force_download):
            to_unzip.append('Pfam-A.hmm.dat.gz')
            to_download.append(pfam_metadata)
        if not file_exists(self.mantis_paths['pfam'] + 'pfam2go', force_download): to_download.append(pfam2go)
        for url in to_download:
            download_file(url, output_folder=self.mantis_paths['pfam'], stdout_file=stdout_file)
        for file_to_unzip in to_unzip:
            uncompress_archive(source_filepath=self.mantis_paths['pfam'] + file_to_unzip, stdout_file=stdout_file,
                               remove_source=True)
        if not self.check_reference_exists('pfam'):
            run_command('hmmpress ' + self.mantis_paths['pfam'] + 'Pfam-A.hmm', stdout_file=stdout_file)
        self.compile_pfam_metadata()

    def download_and_unzip_kofam_hmm(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['kofam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('kofam', force_download) and \
                file_exists(self.mantis_paths['kofam'] + 'ko_list', force_download) and \
                file_exists(self.mantis_paths['kofam'] + 'ko2cog.xl', force_download) and \
                file_exists(self.mantis_paths['kofam'] + 'ko2go.xl', force_download) and \
                file_exists(self.mantis_paths['kofam'] + 'ko2tc.xl', force_download) and \
                file_exists(self.mantis_paths['kofam'] + 'ko2cazy.xl', force_download) and \
                file_exists(self.mantis_paths['kofam'] + 'ko_to_path', force_download) and \
                file_exists(self.mantis_paths['kofam'] + 'map_description', force_download):
            print('KOfam hmm already exists! Skipping...', flush=True, file=stdout_file)
            return
        kofam_hmm = 'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz'
        ko_list = 'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz'
        ko_to_cog = 'https://www.kegg.jp/kegg/files/ko2cog.xl'
        ko_to_go = 'https://www.kegg.jp/kegg/files/ko2go.xl'
        ko_to_tc = 'https://www.kegg.jp/kegg/files/ko2tc.xl'
        ko_to_cazy = 'https://www.kegg.jp/kegg/files/ko2cazy.xl'
        ko_to_path = 'http://rest.kegg.jp/link/pathway/ko'
        map_description = 'https://www.genome.jp/kegg/pathway.html'
        if not file_exists(self.mantis_paths['kofam'] + 'readme.md'):
            with open(self.mantis_paths['kofam'] + 'readme.md', 'w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write(f'This hmm was downloaded on {datetime_str} from:\n{kofam_hmm}\nMetadata was downloaded from:\n{ko_list}\n{ko_to_cog}\n{ko_to_go}\n{ko_to_tc}\n{ko_to_cazy}')
        print_cyan('Downloading and unzipping KOfam hmms ', flush=True, file=stdout_file)
        for url in [kofam_hmm, ko_list, ko_to_cog, ko_to_go, ko_to_tc, ko_to_cazy, ko_to_path, map_description]:
            download_file(url, output_folder=self.mantis_paths['kofam'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'profiles.tar.gz',
                           extract_path=self.mantis_paths['kofam'], stdout_file=stdout_file, remove_source=True)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'ko_list.gz', stdout_file=stdout_file,
                           remove_source=True)
        move_file(self.mantis_paths['kofam'] + 'ko', self.mantis_paths['kofam'] + 'ko_to_path')
        move_file(self.mantis_paths['kofam'] + 'pathway.html', self.mantis_paths['kofam'] + 'map_description')
        merge_profiles(self.mantis_paths['kofam'] + 'profiles/', self.mantis_paths['kofam'] + 'kofam_merged.hmm',
                       stdout_file=stdout_file)
        run_command('hmmpress ' + self.mantis_paths['kofam'] + 'kofam_merged.hmm', stdout_file=stdout_file)

    def download_and_unzip_NCBI(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['NCBI']).mkdir(parents=True, exist_ok=True)
        ncbi_hmm = 'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz'
        metadata = 'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv'
        # we cant verify a priori which foulders we should have, so you need to delete the folder to restart

        if file_exists(self.mantis_paths['NCBI'] + 'readme.md', force_download):
            print('NCBI hmm folder already exists! Skipping...', flush=True, file=stdout_file)
            return

        if not file_exists(self.mantis_paths['NCBI'] + 'readme.md'):
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
        if os.path.exists(self.mantis_paths['NCBI'] + 'hmm_PGAP'):
            shutil.rmtree(self.mantis_paths['NCBI'] + 'hmm_PGAP')

    # sorting profiles by the taxonomic id from ncbi metadata file
    def compile_hmms_NCBI(self, stdout_file=None):
        sorted_metadata = self.sort_hmms_NCBI()
        self.write_metadata(sorted_metadata)
        self.assign_hmm_profiles(sorted_metadata, stdout_file=stdout_file)

    def assign_hmm_profiles(self, sorted_metadata, stdout_file=None):
        for taxa in sorted_metadata:
            for hmm, hmm_label, description, enzyme_ec,common_links in sorted_metadata[taxa]:
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

    def write_metadata(self, sorted_metadata):
        for taxa in sorted_metadata:
            Path(self.mantis_paths['NCBI'] + taxa).mkdir(parents=True, exist_ok=True)
            Path(add_slash(self.mantis_paths['NCBI'] + taxa) + 'to_merge').mkdir(parents=True, exist_ok=True)
            with open(add_slash(self.mantis_paths['NCBI'] + taxa) + 'metadata.tsv', 'w+') as metadata_file:
                for hmm, hmm_label, description, enzyme_ec,common_links in sorted_metadata[taxa]:
                    line = [hmm_label, '|', f'description:{description}']
                    for db in common_links:
                        for db_id in common_links[db]:
                            line.append(f'{db}:{db_id}')
                    for ec in enzyme_ec:
                        line.append(f'enzyme_ec:{ec}')
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
                hmm, hmm_label, description, enzyme_ec, taxa_id = line[0], line[2], line[10], line[12], line[14]
                common_links={}
                self.get_common_links_metadata(description,common_links)
                for db in common_links:
                    for db_id in common_links[db]:
                        description=description.replace(db_id,'').strip()
                description = description.replace('(Provisional)', '')
                description = description.strip()
                enzyme_ec = [i for i in enzyme_ec.split(',') if i]
                line = file.readline()
                if taxa_id:
                    if taxa_id not in res: res[taxa_id] = []
                    res[taxa_id].append([hmm, hmm_label, description, enzyme_ec,common_links])
                    if taxa_id in general_taxon_ids:
                        res['NCBIG'].append([hmm, hmm_label, description, enzyme_ec,common_links])
                else:
                    res['NCBIG'].append([hmm, hmm_label, description, enzyme_ec,common_links])
        return res

    def download_and_unzip_tigrfam_hmm(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['tigrfam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('tigrfam', force_download) and \
                file_exists(self.mantis_paths['tigrfam'] + 'gpl.html', force_download) and \
                file_exists(self.mantis_paths['tigrfam'] + 'COPYRIGHT', force_download) and \
                file_exists(self.mantis_paths['tigrfam'] + 'TIGRFAMS_GO_LINK', force_download) and \
                file_exists(self.mantis_paths['tigrfam'] + 'TIGRFAMS_ROLE_LINK', force_download) and \
                file_exists(self.mantis_paths['tigrfam'] + 'TIGR_ROLE_NAMES', force_download):
            print('TIGRfam hmm already exists! Skipping...', flush=True, file=stdout_file)
            return
        tigrfam_hmm = 'ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz'
        tigrfam_go_link = 'ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMS_GO_LINK'
        tigrfam_role_link = 'ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMS_ROLE_LINK'
        tigrfam_role_names = 'ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGR_ROLE_NAMES'
        copyright = 'ftp://ftp.jcvi.org/pub/data/TIGRFAMs/COPYRIGHT'
        license = 'http://www.gnu.org/licenses/gpl.html'
        if not file_exists(self.mantis_paths['tigrfam'] + 'readme.md'):
            with open(self.mantis_paths['tigrfam'] + 'readme.md', 'w+') as  file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write(f'This hmm was downloaded on {datetime_str} from:\ntigrfam_hmm\nMetadata was downloaded from:\n{tigrfam_go_link}\n{tigrfam_role_link}\n{tigrfam_role_names}\n{copyright}\nLicense was downloaded from:\n{license}')
        print_cyan('Downloading and unzipping TIGRfam hmms', flush=True, file=stdout_file)
        for link in [tigrfam_hmm,
                     tigrfam_go_link,
                     tigrfam_role_link,
                     tigrfam_role_names,
                     copyright,
                     license,
                     ]:
            download_file(link, output_folder=self.mantis_paths['tigrfam'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['tigrfam'] + 'TIGRFAMs_15.0_HMM.tar.gz',
                           extract_path=self.mantis_paths['tigrfam'] + 'profiles', stdout_file=stdout_file,
                           remove_source=True)
        merge_profiles(self.mantis_paths['tigrfam'] + 'profiles/', self.mantis_paths['tigrfam'] + 'tigrfam_merged.hmm',
                       stdout_file=stdout_file)
        run_command('hmmpress ' + self.mantis_paths['tigrfam'] + 'tigrfam_merged.hmm ', stdout_file=stdout_file)



    def check_completeness_NOGG(self, nogg_file, list_file_paths, force_download):
        if force_download:
            return False
        if not os.path.exists(nogg_file):
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
                if os.path.exists(self.mantis_paths['NOG'] + 'NOGG'): shutil.rmtree(self.mantis_paths['NOG'] + 'NOGG')
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

    def download_and_unzip_NOGT(self, taxon_id, stdout_file=None):
        folder_path = add_slash(self.mantis_paths['NOG'] + taxon_id)
        Path(folder_path).mkdir(parents=True, exist_ok=True)
        target_annotation_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_annotations.tsv'
        target_merged_hmm = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_merged.hmm'
        eggnog_downloads_page = f'http://eggnog5.embl.de/download/latest/per_tax_level/{taxon_id}/'
        for file in ['_annotations.tsv.gz', '_hmms.tar']:
            url = f'{eggnog_downloads_page}{taxon_id}{file}'
            download_file(url, output_folder=folder_path, stdout_file=stdout_file)
        if os.path.exists(f'{folder_path}profiles'): shutil.rmtree(f'{folder_path}profiles')
        uncompress_archive(source_filepath= f'{folder_path}{taxon_id}_hmms.tar', extract_path= f'{folder_path}profiles',stdout_file=stdout_file, remove_source=True)
        uncompress_archive(source_filepath= f'{folder_path}{taxon_id}_annotations.tsv.gz', stdout_file=stdout_file,remove_source=True)
        # removing gz extension
        hmm_files = os.listdir(folder_path + 'profiles')
        for hmm_profile in hmm_files:
            if '.hmm' in hmm_profile: move_file(hmm_profile, hmm_profile.strip('.gz'))
        merge_profiles(f'{folder_path}profiles/{taxon_id}', f'{folder_path}{taxon_id}_merged.hmm',
                       stdout_file=stdout_file)
        if os.path.exists(f'{folder_path}profiles'):      shutil.rmtree(f'{folder_path}profiles')
        run_command(f'hmmpress {target_merged_hmm}', stdout_file=stdout_file)

    def download_and_unzip_eggnogdb(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['default']).mkdir(parents=True, exist_ok=True)
        if file_exists(self.mantis_paths['default'] + 'eggnog.db', force_download):
            print('eggnog.db already exists! Skipping...', flush=True, file=stdout_file)
            return
        url = 'http://eggnogdb.embl.de/download/emapperdb-5.0.1/eggnog.db.gz'
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
                if not os.path.exists(add_slash(self.mantis_paths['NOG'] + str(tax_id)) + f'{tax_id}_sql_annotations.tsv'): passed = False
            if passed:
                print('All chosen NOGT were already extracted! Skipping...', flush=True, file=stdout_file)
                return
        resources_path = self.mantis_paths['resources']
        NOG_sql_resources_folder = add_slash(resources_path + 'NOG_sql')
        NOGT_sql_folder = add_slash(self.mantis_paths['NOG'] + 'NOG_sql')
        Path(NOGT_sql_folder).mkdir(parents=True, exist_ok=True)
        Path(NOG_sql_resources_folder).mkdir(parents=True, exist_ok=True)


        if not os.path.exists(NOG_sql_resources_folder): return
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
        if os.path.exists(NOGT_sql_folder):
            shutil.rmtree(NOGT_sql_folder)

    ### Support functions for setting up queue to download NOGT hmms

    def get_taxon_ids_NOGT(self):
        url = 'http://eggnog5.embl.de/download/latest/per_tax_level/'
        webpage = None
        c = 0
        while not webpage and c <= 10:
            req = requests.get(url)
            try:
                webpage = req.text
            except:
                c += 1
        taxons_search = re.findall('href="\d+/"', webpage)
        taxons = [re.search('\d+', i).group() for i in taxons_search]
        return taxons

    def get_taxon_for_queue_NOGT(self, force_download=False, stdout_file=None):
        Path(self.mantis_paths['NOG']).mkdir(parents=True, exist_ok=True)
        # we will download tax specific hmms, which can be used for more specific hmmscans
        if not file_exists(self.mantis_paths['NOG'] + 'readme.md'):
            with open(self.mantis_paths['NOG'] + 'readme.md', 'w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write(f'These hmms were downloaded on {datetime_str} from:\nhttp://eggnog5.embl.de/download/latest/per_tax_level/')
        available_taxon_ids = self.get_taxon_ids_NOGT()
        if self.mantis_nogt_tax:
            taxon_ids = self.mantis_nogt_tax
        else:
            taxon_ids=available_taxon_ids
        res = []
        for taxon_id in taxon_ids:
            if taxon_id != '1' and taxon_id in available_taxon_ids:
                target_annotation_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + taxon_id + '_annotations.tsv'
                if self.check_reference_exists('NOG', taxon_id, force_download) and file_exists(target_annotation_file,force_download):
                    hmm_path = get_hmm_in_folder(self.mantis_paths['NOG'] + taxon_id)
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
            taxon_ids = os.listdir(self.mantis_paths['NOG'])
            for taxon_id in taxon_ids:
                if taxon_id != '1' and os.path.isdir(self.mantis_paths['NOG'] + taxon_id):
                    hmm_path = get_hmm_in_folder(self.mantis_paths['NOG'] + taxon_id)
                    print('Checking NOG for splitting:', hmm_path, taxon_id, flush=True, file=stdout_file)
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
            taxon_ids = os.listdir(self.mantis_paths['NCBI'])
            for taxon_id in taxon_ids:
                if os.path.isdir(self.mantis_paths['NCBI'] + taxon_id):
                    hmm_path = get_hmm_in_folder(self.mantis_paths['NCBI'] + taxon_id)
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
        hmms_list = self.compile_hmms_list()
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
        if os.path.exists(hmm_chunks_folder):
            shutil.rmtree(hmm_chunks_folder)
        Path(hmm_chunks_folder).mkdir(parents=True, exist_ok=True)
        print('Load balancing chunks', flush=True, file=stdout_file)
        load_balanced_chunks = chunk_generator_load_balanced([i for i in range(profile_count)],
                                                             get_hmm_chunk_size(total_profiles=profile_count,
                                                                                current_chunk_size=self.hmm_chunk_size,
                                                                                max_chunks=100), time_limit=None)
        print(f'Splitting {hmm_path} into {len(load_balanced_chunks)} chunks.', flush=True,
              file=stdout_file)
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
        with open(self.mantis_paths['pfam'] + 'Pfam-A.hmm.dat') as pfam_dat_file:
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
                        res['enzyme_ec'].update(kegg_ec.split(','))
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

        if os.path.exists(target_sql_file):
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
        while not os.path.exists(self.mantis_paths['default'] + "eggnog.db"): sleep(5)
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
