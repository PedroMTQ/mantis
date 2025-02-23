import os
import shutil
from pathlib import Path

from mantis.src.metadata.utils import get_common_links_metadata
from mantis.src.setup.ncbi import SetupNcbi
from mantis.src.setup.pfam import SetupPfam
from mantis.src.setup.resources_ncbi import SetupResourcesNcbi
from mantis.src.setup.resources_taxonomy import SetupResourcesTaxonomy
from mantis.src.setup.tcdb import SetupTcdb
from mantis.src.utils import estimate_number_workers_setup_database, print_cyan, remove_file
from mantis.src.utils.logger import logger


class DatabaseGenerator():

    #####################   Main function
    @timeit_class
    def setup_databases(self):
        self.output_folder = f'{MANTIS_FOLDER}setup_databases/'
        self.mantis_out = f'{self.output_folder}Mantis.out'
        if os.path.exists(self.mantis_out):
            os.remove(self.mantis_out)
        Path(self.mantis_paths['default']).mkdir(parents=True, exist_ok=True)
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        dbs_list = [
            'ncbi_res',
            'pfam' if self.mantis_paths['pfam'][0:2] != 'NA' else None,
            'kofam' if self.mantis_paths['kofam'][0:2] != 'NA' else None,
            'tcdb' if self.mantis_paths['tcdb'][0:2] != 'NA' else None,
            'NCBI' if self.mantis_paths['NCBI'][0:2] != 'NA' else None,
        ]
        if self.use_taxonomy:
            dbs_list.insert(0, 'taxonomy')
        # DOWNLOADING
        self.prepare_queue_setup_databases(dbs_list)
        # for unzipping tax specific hmms
        passed_tax_check = self.prepare_queue_setup_databases_tax()
        worker_count = estimate_number_workers_setup_database(len(self.queue), user_cores=self.user_cores)
        if worker_count:
            print(f'Database will be setup with {worker_count} workers!',
                  flush=True,
                  file=self.redirect_verbose)

        self.processes_handler(self.worker_setup_databases, worker_count)
        print_cyan('Finished downloading all data!', flush=True, file=self.redirect_verbose)

        if not passed_tax_check:
            # METADATA for NOG
            if self.nog_db == 'hmm':
                # for NOG HMM we just recompile incomplete ones
                self.prepare_queue_extract_metadata_NOG_HMM()
                worker_count = estimate_number_workers_setup_database(len(self.queue),
                                                                      minimum_jobs_per_worker=2,
                                                                      user_cores=self.user_cores)
                if worker_count: print(f'NOG metadata will be extracted with {worker_count} workers!',
                                       flush=True,
                                       file=self.redirect_verbose)
                self.processes_handler(self.worker_extract_NOG_metadata_HMM, worker_count)
                print('NOGG will now be compiled', flush=True, file=self.redirect_verbose)
                self.compile_NOGG_HMM()
                remove_file(self.mantis_paths['NOG'] + 'Pfam-A.hmm.dat')
                for file in ['eggnog.db', 'Pfam-A.hmm.dat']:
                    file_path = self.mantis_paths['NOG'] + file
                    remove_file(file_path)
            else:
                # for NOG diamond we recompile everything
                seqs_taxons = self.compile_NOG_DMND()
                self.prepare_queue_extract_fastas_DMND(seqs_taxons)
                worker_count = estimate_number_workers_setup_database(len(self.queue),
                                                                      minimum_jobs_per_worker=2,
                                                                      user_cores=self.user_cores)
                if worker_count: print(f'NOG diamond fastas will be extracted with {worker_count} workers!',
                                       flush=True,
                                       file=self.redirect_verbose)
                self.processes_handler(self.worker_extract_fastas_DMND, worker_count)
                print('Creating NOGT diamond databases', flush=True, file=self.redirect_verbose)
                self.compile_NOGT_DMND()
                self.compile_NOGG_DMND()
                for file in ['eggnog_proteins.dmnd', 'eggnog_proteins.dmnd.gz', 'eggnog.db', 'eggnog_seqs.faa',
                             'metadata.tsv', 'Pfam-A.hmm.dat']:
                    file_path = self.mantis_paths['NOG'] + file
                    remove_file(file_path)
            # we also remove the 1 folder since it doesn't actually exist in NCBI, this taxa is just a general taxon which we already added to NOGG anyway
            if os.path.exists(self.mantis_paths['NOG'] + '1/'):
                shutil.rmtree(self.mantis_paths['NOG'] + '1/')

        # SPLITTING
        if self.hmm_chunk_size:
            print_cyan('Will now split data into chunks!', flush=True, file=self.redirect_verbose)
            self.prepare_queue_split_hmms()
            worker_count = estimate_number_workers_setup_database(len(self.queue), user_cores=self.user_cores)
            print(f'Database will be split with {worker_count} workers!', flush=True, file=self.redirect_verbose)
            self.processes_handler(self.worker_split_hmms, worker_count)
        self.prepare_queue_press_custom_hmms()
        worker_count = estimate_number_workers_setup_database(len(self.queue), user_cores=self.user_cores)
        print(f'HMMs will be pressed with {worker_count} workers!', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_press_custom_hmms, worker_count)
        print('Preparing NLP Resources!', flush=True, file=self.redirect_verbose)
        UniFunc_wrapper.__init__(self)
        print_cyan('Finished setting up databases!', flush=True, file=self.redirect_verbose)

    #####################   Filling queue with jobs

    def prepare_queue_setup_databases(self, dbs_list):
        for database in dbs_list:
            if database:
                self.queue.append([database, self.mantis_out])

    def prepare_queue_setup_databases_tax(self):
        if not self.use_taxonomy:
            return True
        stdout_file = open(self.mantis_out, 'a+')
        passed_tax_check = True
        if self.mantis_paths['NOG'][0:2] != 'NA':
            Path(self.mantis_paths['NOG']).mkdir(parents=True, exist_ok=True)
            list_taxon_ids = self.get_taxon_for_queue_NOGT()
            if not os.path.exists(self.mantis_paths['NOG']): passed_tax_check = False
            for taxon_id in list_taxon_ids:
                if taxon_id != '1':
                    if not self.check_reference_exists('NOGT', taxon_id=taxon_id):
                        passed_tax_check = False
                        break
            if self.nog_db == 'hmm':
                if not passed_tax_check:
                    for taxon_id in list_taxon_ids:
                        self.queue.append(['NOG_HMM', taxon_id, self.mantis_out])
            else:
                if not passed_tax_check:
                    self.queue.append(['NOG_DMND', self.mantis_out])

        if not passed_tax_check:
            if os.path.exists(self.mantis_paths['NOG']):
                shutil.rmtree(self.mantis_paths['NOG'])
            Path(self.mantis_paths['NOG']).mkdir(parents=True, exist_ok=True)
            with open(self.mantis_paths['NOG'] + 'readme.md', 'w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write(f'This data was downloaded on {datetime_str}')
            self.download_and_unzip_eggnogdb()
            self.download_pfam_id_to_acc()

        return passed_tax_check

    def prepare_queue_split_hmms(self):
        print('Checking which HMMs need to be split, this may take a while...', flush=True, file=self.redirect_verbose)
        taxon_nog_hmms = self.get_taxon_for_queue_NOG_split_hmms()
        taxon_ncbi_hmms = self.get_taxon_for_queue_ncbi_split_hmms()
        general_hmms = self.get_hmm_for_queue_split_hmms()
        print('Will split: ', [get_path_level(i) for i in taxon_nog_hmms + taxon_ncbi_hmms + general_hmms],
              flush=True,
              file=self.redirect_verbose)
        for hmm_path in taxon_nog_hmms + taxon_ncbi_hmms + general_hmms:
            self.queue.append([hmm_path, self.mantis_out])

    def prepare_queue_press_custom_hmms(self):
        print('Checking which custom hmms need to be pressed', flush=True, file=self.redirect_verbose)
        hmms_list = []
        for hmm_path in self.get_custom_reference_paths(folder=False):
            if hmm_path.endswith('.hmm'):
                hmm_folder = add_slash(SPLITTER.join(hmm_path.split(SPLITTER)[:-1]))
                hmm_name = hmm_path.split(SPLITTER)[-1]
                if 'chunks' not in os.listdir(hmm_folder):
                    if hmm_folder[0:2] != 'NA':
                        if not self.check_missing_chunk_files(hmm_name, hmm_folder):
                            hmms_list.append(hmm_path)
        if hmms_list:
            print(f'Will hmmpress: {hmms_list}', flush=True, file=self.redirect_verbose)
        for hmm_path in hmms_list:
            self.queue.append([hmm_path, self.mantis_out])

    def prepare_queue_extract_metadata_NOG_HMM(self):
        if self.mantis_paths['NOG'][0:2] != 'NA':
            print_cyan('Will now extract metadata from NOG!', flush=True, file=self.redirect_verbose)
            stdout_file = open(self.mantis_out, 'a+')
            for taxon_id in os.listdir(self.mantis_paths['NOG']):
                if os.path.isdir(self.mantis_paths['NOG'] + taxon_id) and taxon_id != 'NOGG':
                    target_annotation_file = add_slash(
                        self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_annotations.tsv'
                    target_sql_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + 'metadata.tsv'
                    if os.path.exists(target_sql_file):
                        if len(self.get_hmms_annotation_file(target_sql_file, hmm_col=0)) != len(
                                self.get_hmms_annotation_file(target_annotation_file, hmm_col=1)):
                            self.queue.append([target_sql_file, target_annotation_file, taxon_id, self.mantis_out])
                        else:
                            print(f'Skipping metadata extraction for NOGT {taxon_id}', flush=True, file=stdout_file)
                    else:
                        self.queue.append([taxon_id, self.mantis_out])
            stdout_file.close()

    #####################   Executing queue jobs

    def worker_setup_databases(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None:
                break
            if len(record) == 2:
                database, stdout_path = record
                taxon_id = None
            else:
                database, taxon_id, stdout_path = record
            self.download_database(database, taxon_id, stdout_path)
            if not self.check_reference_exists(database, taxon_id):
                stdout_file = open(stdout_path, 'a+')
                if taxon_id:
                    print(f'Setup failed on {database} {taxon_id}', flush=True, file=stdout_file)
                else:
                    print(f'Setup failed on {database}', flush=True, file=stdout_file)
                stdout_file.close()
                queue.insert(0, record)

    def worker_split_hmms(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None:
                break
            hmm_path, stdout_path = record
            self.split_hmm_into_chunks(hmm_path, stdout_path)

    def worker_press_custom_hmms(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            hmm_path, stdout_path = record
            self.press_custom_hmms(hmm_path, stdout_path)

    def worker_extract_NOG_metadata_HMM(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            taxon_id, stdout_path = record
            self.get_metadata_hmms(taxon_id=taxon_id, stdout_path=stdout_path)

    #####################   Main download function

    def download_database(self, database, taxon_id=None, stdout_path=None):
        stdout_file = open(stdout_path, 'a+')
        if database == 'taxonomy':
            SetupResourcesTaxonomy(config_file=self.config_file).run()
        elif database == 'ncbi_res':
            SetupResourcesNcbi(config_file=self.config_file).run()
        elif database == 'pfam':
            SetupPfam(config_file=self.config_file).run()
        elif database == 'kofam':
            self.download_kofam(stdout_file=stdout_file)
        elif database == 'tcdb':
            SetupTcdb(config_file=self.config_file).run()
        elif database == 'NCBI':
            SetupNcbi(config_file=self.config_file).run()
        elif database == 'NOG_HMM':
            self.download_NOGT(taxon_id=taxon_id, stdout_file=stdout_file)
        elif database == 'NOG_DMND':
            self.download_NOG_DMND(stdout_file=stdout_file)
        if taxon_id:
            logger.info(f'Finished downloading {database} with taxon {taxon_id}')

    #####################   KOFAM

    def compile_kofam_metadata(self):
        metadata_to_write = {}
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
                if ko not in res: res[ko] = {}
                if '[EC:' in description:
                    description, temp_links = description.split('[EC:')
                else:
                    temp_links = description
                get_common_links_metadata(input_string=temp_links,
                                          metadata_dict=res[ko])
                if 'kegg_ko' not in res[ko]: res[ko]['kegg_ko'] = set()
                res[ko]['kegg_ko'].add(ko)
                if 'description' not in res[ko]:
                    res[ko]['description'] = set()
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
                if ko not in res: res[ko] = {}
                if target_link not in res[ko]: res[ko][target_link] = set()
                res[ko][target_link].update(link)
                line = file.readline()

    def download_kofam(self, stdout_file=None):
        Path(self.mantis_paths['kofam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('kofam') and \
                os.path.exists(self.mantis_paths['kofam'] + 'metadata.tsv'):
            print('KOfam HMM already exists! Skipping...', flush=True, file=stdout_file)
            return
        kofam_hmm = 'https://www.genome.jp/ftp/db/kofam/profiles.tar.gz'
        ko_list = 'https://www.genome.jp/ftp/db/kofam/ko_list.gz'
        ko_to_cog = 'https://www.kegg.jp/kegg/files/ko2cog.xl'
        ko_to_go = 'https://www.kegg.jp/kegg/files/ko2go.xl'
        ko_to_tc = 'https://www.kegg.jp/kegg/files/ko2tc.xl'
        ko_to_cazy = 'https://www.kegg.jp/kegg/files/ko2cazy.xl'
        with open(self.mantis_paths['kofam'] + 'readme.md', 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(
                f'This hmm was downloaded on {datetime_str} from:\n{kofam_hmm}\nMetadata was downloaded from:\n{ko_list}\n{ko_to_cog}\n{ko_to_go}\n{ko_to_tc}\n{ko_to_cazy}')
        print_cyan('Downloading and unzipping KOfam hmms ', flush=True, file=stdout_file)
        for url in [kofam_hmm, ko_list, ko_to_cog, ko_to_go, ko_to_tc, ko_to_cazy]:
            download_file(url, output_folder=self.mantis_paths['kofam'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'profiles.tar.gz',
                           extract_path=self.mantis_paths['kofam'], stdout_file=stdout_file, remove_source=True)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'ko_list.gz', stdout_file=stdout_file,
                           remove_source=True)
        merge_profiles(self.mantis_paths['kofam'] + 'profiles/', self.mantis_paths['kofam'] + 'kofam_merged.hmm',
                       stdout_file=stdout_file)
        run_command('hmmpress ' + self.mantis_paths['kofam'] + 'kofam_merged.hmm', stdout_file=stdout_file)
        self.compile_kofam_metadata()

    #####################   NOG

    def check_completeness_NOGG(self, nogg_file, list_file_paths):
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

    def compile_NOGG_HMM(self):
        # this is not feasible for all HMMs, so we just select the most general taxa (i.e., domains)
        if self.mantis_paths['NOG'][0:2] != 'NA':
            stdout_file = open(self.mantis_out, 'a+')
            target_annotation_file = add_slash(self.mantis_paths['NOG'] + 'NOGG') + 'metadata.tsv'
            target_merged_hmm = add_slash(self.mantis_paths['NOG'] + 'NOGG') + 'NOGG_merged.hmm'
            all_sql = set()
            all_hmm = set()
            for taxon_id in self.get_ncbi_domains():
                if taxon_id != 'NOGG':
                    if os.path.isdir(self.mantis_paths['NOG'] + taxon_id):
                        target_hmm_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_merged.hmm'
                        target_sql_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + 'metadata.tsv'
                        all_sql.add(target_sql_file)
                        all_hmm.add(target_hmm_file)

            if not self.check_completeness_NOGG(target_annotation_file, all_sql) or \
                    not self.check_reference_exists('NOGG'):
                if os.path.exists(self.mantis_paths['NOG'] + 'NOGG'): shutil.rmtree(self.mantis_paths['NOG'] + 'NOGG')
                Path(self.mantis_paths['NOG'] + 'NOGG').mkdir(parents=True, exist_ok=True)
            else:
                print('NOGG already compiled, skipping...', flush=True, file=stdout_file)
                return
            print_cyan('Compiling global NOG hmms ', flush=True, file=stdout_file)
            # now to merge all tshmms into one
            merge_redundant_profiles(output_file=target_merged_hmm, list_file_paths=all_hmm, stdout_file=stdout_file)
            merge_redundant_sql_annotations(output_file=target_annotation_file, list_file_paths=all_sql,
                                            stdout_file=stdout_file)
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
        if os.path.exists(f'{folder_path}profiles'): shutil.rmtree(f'{folder_path}profiles')
        uncompress_archive(source_filepath=f'{folder_path}{taxon_id}_hmms.tar.gz',
                           extract_path=f'{folder_path}profiles', stdout_file=stdout_file, remove_source=True)
        uncompress_archive(source_filepath=f'{folder_path}{taxon_id}_annotations.tsv.gz',
                           stdout_file=stdout_file,
                           remove_source=True)
        # removing gz extension
        hmm_files = os.listdir(folder_path + 'profiles')
        for hmm_profile in hmm_files:
            if '.hmm' in hmm_profile: move_file(hmm_profile, hmm_profile.strip('.gz'))
        merge_profiles(f'{folder_path}profiles/{taxon_id}', f'{folder_path}{taxon_id}_merged.hmm',
                       stdout_file=stdout_file)
        if os.path.exists(f'{folder_path}profiles'):      shutil.rmtree(f'{folder_path}profiles')
        run_command(f'hmmpress {target_merged_hmm}', stdout_file=stdout_file)

    def download_and_unzip_eggnogdb(self, stdout_file=None):
        Path(self.mantis_paths['default']).mkdir(parents=True, exist_ok=True)
        if os.path.exists(self.mantis_paths['NOG'] + 'eggnog.db'):
            print('eggnog.db already exists! Skipping...', flush=True, file=stdout_file)
            return
        else:
            if os.path.exists(self.mantis_paths['default'] + 'eggnog.NOG'):
                remove_file(self.mantis_paths['NOG'] + 'eggnog.db')
        url = 'http://eggnogdb.embl.de/download/emapperdb-' + self.get_latest_version_eggnog() + '/eggnog.db.gz'
        download_file(url, output_folder=self.mantis_paths['NOG'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['NOG'] + 'eggnog.db.gz',
                           stdout_file=stdout_file,
                           remove_source=True)

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
        res = []
        for v in versions_id:
            version_id = v.replace('href="emapperdb-', '')
            version_id = version_id.replace('/"', '')
            res.append(version_id)
        if res:
            return max(res)
        else:
            # hard coded just in case the user doesnt have access to internet, should be updated
            current_version = '5.0.2'
            print(
                f'eggNOG database version retrieval failed, so returning hardcoded eggnog database version {current_version}, please change <current_version> if a newer version is available')
            return current_version

    def get_taxon_ids_eggNOG(self):
        # this will get all available taxon ids from the eggnog webpage
        url = 'http://eggnog5.embl.de/download/latest/per_tax_level/'
        webpage = None
        c = 0
        while not webpage and c <= 10:
            try:
                req = requests.get(url)
                webpage = req.text
            except:
                c += 1
        if isinstance(webpage, str):
            taxons_search = re.findall(r'href="\d+/"', webpage)
            taxons = [re.search(r'\d+', i).group() for i in taxons_search]
        else:
            taxons = []
        # if no connection is established, we just return the local taxon ids (will be empty if no connection is available during setup)
        if not taxons: return self.get_local_ref_taxon_ids('NOG')
        return taxons

    def get_taxon_for_queue_NOGT(self, stdout_file=None):
        # we will download tax specific hmms, which can be used for more specific hmmscans
        taxon_ids = self.get_ref_taxon_ids('NOG')
        res = []
        for taxon_id in taxon_ids:
            target_annotation_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + taxon_id + '_annotations.tsv'
            if self.check_reference_exists('NOG', taxon_id) and os.path.exists(target_annotation_file):
                hmm_path = get_ref_in_folder(self.mantis_paths['NOG'] + taxon_id)
                profile_count = get_hmm_profile_count(hmm_path)
                annotations_count = len(self.get_hmms_annotation_file(target_annotation_file, 1))
                # should be the same but some HMMs dont have an annotation (probably an error from NOG)
                if profile_count in range(annotations_count - 10, annotations_count + 10 + 1):
                    print(
                        f'Tax {taxon_id} hmm and annotations.tsv files already exist, and they were merged correctly. Skipping setup...',
                        flush=True, file=stdout_file)
                else:
                    print(f'Tax {taxon_id} hmm already exists but was not merged correctly!', flush=True,
                          file=stdout_file)
                    res.append(taxon_id)
            else:
                res.append(taxon_id)
        print('Will check data for the following NOGT:\n' + ','.join(res), flush=True, file=stdout_file)
        return res

    ### Support functions for setting up queue for split hmms

    def get_taxon_for_queue_NOG_split_hmms(self):
        res = []
        if self.mantis_paths['NOG'][0:2] != 'NA' and self.nog_db == 'hmm':
            stdout_file = open(self.mantis_out, 'a+')

            for taxon_id in self.get_ref_taxon_ids('NOG'):
                if os.path.isdir(self.mantis_paths['NOG'] + taxon_id):
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
        if os.path.exists(hmm_chunks_folder):
            shutil.rmtree(hmm_chunks_folder)
        Path(hmm_chunks_folder).mkdir(parents=True, exist_ok=True)
        print('Load balancing chunks', flush=True, file=stdout_file)
        load_balanced_chunks = chunk_generator_load_balanced([i for i in range(profile_count)],
                                                             get_hmm_chunk_size(total_profiles=profile_count,
                                                                                current_chunk_size=self.hmm_chunk_size,
                                                                                max_chunks=100), time_limit=None)
        print(f'Splitting {hmm_path} into {len(load_balanced_chunks)} chunks.', flush=True, file=stdout_file)
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
        old_files = ['h3f', 'h3i', 'h3m', 'h3p']
        for inner_file in os.listdir(hmm_folder):
            inner_file_ending = inner_file.split('.')[-1]
            if inner_file_ending in old_files:
                os.remove(hmm_folder + inner_file)
        run_command(f'hmmpress {hmm_path}', stdout_file=stdout_file)
        stdout_file.close()

    ### Support functions for extracting metadata

    def get_hmms_annotation_file(self, annotations_file, hmm_col):
        res = set()
        if not os.path.exists(annotations_file): return res
        with open(annotations_file, 'r') as file:
            line = file.readline()
            while line:
                line = line.split('\t')
                res.add(line[hmm_col])
                line = file.readline()
        return res

    def download_pfam_id_to_acc(self):
        if not os.path.exists(self.mantis_paths['NOG'] + 'Pfam-A.hmm.dat'):
            pfam_metadata = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
            download_file(pfam_metadata, output_folder=self.mantis_paths['NOG'])
            uncompress_archive(source_filepath=self.mantis_paths['NOG'] + 'Pfam-A.hmm.dat.gz', remove_source=True)

    def pfam_id_to_acc(self):
        res = {}
        with open(self.mantis_paths['NOG'] + 'Pfam-A.hmm.dat') as pfam_dat_file:
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

    def clean_up_sql_results_description(self, sql_row, ids_res):
        og, level, description = sql_row
        description = description.strip()
        if description:
            if og not in ids_res: ids_res[og] = {'description': set(), 'cog': set(), 'arcog': set()}
            ids_res[og]['description'].add(description)
            get_common_links_metadata(input_string=description,
                                      metadata_dict=ids_res[og])

    def clean_up_sql_results_ids_hmm(self, sql_row, taxon_id, pfam_id_to_acc):
        codes_to_exclude = ['IEA', 'ND']
        res = {
            'go': set(),
            'pfam': set(),
            'kegg_ko': set(),
            'cog': set(),
            'arcog': set(),
            'enzyme_ec': set(),
            'kegg_pathway': set(),
            'kegg_module': set(),
            'kegg_reaction': set(),
            'kegg_rclass': set(),
            'kegg_brite': set(),
            'cazy': set(),
            'bigg_reaction': set(),
            'description': set(),
            'eggnog': set(),
            'tcdb': set(),
        }
        eggnog_hmms = set()
        ogs, gene_ontology_gos, pfam_ids, kegg_ko, kegg_cog, kegg_ec, kegg_brite, kegg_rclass, kegg_tc, kegg_cazy, kegg_pathway, kegg_module, kegg_reaction, kegg_go = sql_row
        for query_hmm_taxon in ogs.split(','):
            query_hmm, query_taxon = query_hmm_taxon.split('@')
            if int(taxon_id) == int(query_taxon):
                eggnog_hmms.add(query_hmm)
        if eggnog_hmms:
            if gene_ontology_gos:
                # http://geneontology.org/docs/guide-go-evidence-codes/
                gene_ontology_gos_copy = str(gene_ontology_gos)
                gene_ontology_gos_copy = gene_ontology_gos_copy.split(',')
                for go_group in gene_ontology_gos_copy:
                    if '|' not in go_group: print(ids_command, '\n', go_group)
                    _, go, evidence_code = go_group.split('|')
                    if evidence_code not in codes_to_exclude:
                        res['go'].add(go.split(':')[1])
            if pfam_ids:
                temp_pfam_ids = pfam_ids.split(',')
                for tpi in temp_pfam_ids:
                    if tpi in pfam_id_to_acc:
                        res['pfam'].add(pfam_id_to_acc[tpi])
            if kegg_ko:
                kegg_ko = kegg_ko.replace('ko:', '')
                res['kegg_ko'].update(kegg_ko.split(','))
            if kegg_cog: res['cog'].update(kegg_cog.split(','))
            if kegg_ec: res['enzyme_ec'].update(kegg_ec.split(','))
            if kegg_brite: res['kegg_brite'].update(kegg_brite.split(','))
            if kegg_rclass: res['kegg_rclass'].update(kegg_rclass.split(','))
            if kegg_tc: res['tcdb'].update(kegg_tc.split(','))
            if kegg_cazy: res['cazy'].update(kegg_cazy.split(','))
            if kegg_pathway: res['kegg_pathway'].update(kegg_pathway.split(','))
            if kegg_module: res['kegg_module'].update(kegg_module.split(','))
            if kegg_reaction: res['kegg_reaction'].update(kegg_reaction.split(','))
            if kegg_go: res['go'].update(kegg_go.split(','))
        return eggnog_hmms, res

    def yield_sql_rows(self, cursor, step_size=1000):
        while True:
            results = cursor.fetchmany(step_size)
            if not results:
                break
            for result in results:
                yield result

    def fetch_eggnog_metadata_hmm(self, cursor, taxon_id, ids_command, description_command, pfam_id_to_acc):
        res = {}
        cursor.execute(ids_command)
        rows = self.yield_sql_rows(cursor)
        for row in rows:
            eggnog_hmms, processed_row = self.clean_up_sql_results_ids_hmm(row, taxon_id, pfam_id_to_acc)
            for hmm in eggnog_hmms:
                if hmm not in res:
                    res[hmm] = dict(processed_row)
                    res[hmm]['eggnog'].add(hmm)
                    if find_cog(hmm): res[hmm]['cog'].add(hmm)
                    if find_arcog(hmm): res[hmm]['arcog'].add(hmm)
                for link_type in res[hmm]:
                    res[hmm][link_type].update(processed_row[link_type])
        cursor.execute(description_command)
        rows = self.yield_sql_rows(cursor)
        for row in rows:
            self.clean_up_sql_results_description(row, res)
        return res

    def get_metadata_hmms(self, taxon_id, stdout_path=None):
        stdout_file = open(stdout_path, 'a+')

        print(f'Exporting metadata for NOG {taxon_id}', flush=True, file=stdout_file)
        metadata_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + 'metadata.tsv'
        eggnog_db_path = self.mantis_paths['NOG'] + "eggnog.db"
        connection = sqlite3.connect(eggnog_db_path)
        cursor = connection.cursor()
        pfam_id_to_acc = self.pfam_id_to_acc()
        ids_command = f'SELECT ogs, gos, pfam, kegg_ko,kegg_cog,kegg_ec,kegg_brite,kegg_rclass,kegg_tc,kegg_cazy,kegg_pathway,kegg_module,kegg_reaction,kegg_go FROM prots WHERE ogs LIKE "%@{taxon_id}%";'
        description_command = f'SELECT og,level,description FROM og WHERE level="{taxon_id}"'
        taxon_metadata = self.fetch_eggnog_metadata_hmm(cursor, taxon_id, ids_command, description_command,
                                                        pfam_id_to_acc)

        with open(metadata_file, "a+") as file:
            for hmm in taxon_metadata:
                link_line = ''
                for link_type in taxon_metadata[hmm]:
                    for inner_link in taxon_metadata[hmm][link_type]:
                        inner_link = inner_link.strip()
                        if inner_link:
                            link_line += '\t' + link_type + ':' + inner_link
                file.write(hmm + '\t|' + link_line + '\n')
        print(f'Finished exporting metadata for NOGT {taxon_id}', flush=True, file=stdout_file)
        stdout_file.close()

        ###### NOG diamond (not implemented yet)

    def download_NOG_DMND(self, stdout_file=None):
        folder_path = add_slash(self.mantis_paths['NOG'])
        if os.path.exists(folder_path + 'eggnog_proteins.dmnd'):
            print('eggnog_proteins.dmnd already exists! Skipping...', flush=True, file=stdout_file)
            return
        url = 'http://eggnogdb.embl.de/download/emapperdb-' + self.get_latest_version_eggnog() + '/eggnog_proteins.dmnd.gz'
        download_file(url, output_folder=folder_path, stdout_file=stdout_file)
        uncompress_archive(source_filepath=f'{folder_path}eggnog_proteins.dmnd.gz', extract_path=folder_path,
                           stdout_file=stdout_file, remove_source=False)

    def merge_fasta_files(self, target_merged_faa, all_faa):
        already_added = set()
        with open(target_merged_faa, 'w+') as file:
            for faa_file in all_faa:
                for query, seq in read_protein_fasta_generator(faa_file):
                    if query not in already_added:
                        line = f'>{query}\n{seq}\n'
                        file.write(line)
                        already_added.add(query)

    def merge_metadata_files(self, target_metadata_file, all_metadata):
        already_added = set()
        with open(target_merged_faa, 'w+') as outfile:
            for metadata_file in all_metadata:
                with open(metadata_file) as infile:
                    for line in infile:
                        query = line.split('\t')[0]
                        if query not in already_added:
                            outfile.write(line)
                            already_added.add(query)

    def compile_NOGT_DMND(self):
        if self.mantis_paths['NOG'][0:2] != 'NA':
            for taxon in self.get_ref_taxon_ids('NOG'):
                taxon_folder = self.mantis_paths['NOG'] + taxon + SPLITTER
                taxon_fasta = f'{taxon_folder}{taxon}_merged.faa'
                taxon_dmnd = f'{taxon_folder}{taxon}'
                if not os.path.exists(f'{taxon_dmnd}.dmnd') and os.path.exists(taxon_fasta):
                    run_command(f'diamond makedb --in {taxon_fasta} -d {taxon_dmnd}')

    def compile_NOGG_DMND(self):
        if self.mantis_paths['NOG'][0:2] != 'NA':
            stdout_file = open(self.mantis_out, 'a+')
            nogg_folder_path = add_slash(self.mantis_paths['NOG'] + 'NOGG')
            target_metadata_file = f'{nogg_folder_path}metadata.tsv'
            target_merged_faa = f'{nogg_folder_path}NOGG_merged.faa'
            all_metadata = set()
            all_faa = set()
            for taxon_id in self.get_ncbi_domains():
                if taxon_id != 'NOGG':
                    if os.path.isdir(self.mantis_paths['NOG'] + taxon_id):
                        taxon_faa_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_merged.faa'
                        all_faa.add(taxon_faa_file)
                        taxon_metadata_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + 'metadata.tsv'
                        all_metadata.add(taxon_metadata_file)

            if not self.check_reference_exists('NOGG'):
                if os.path.exists(nogg_folder_path): shutil.rmtree(nogg_folder_path)
                Path(nogg_folder_path).mkdir(parents=True, exist_ok=True)
            else:
                print('NOGG already compiled, skipping...', flush=True, file=stdout_file)
                return
            print_cyan('Compiling global NOG diamond database ', flush=True, file=stdout_file)
            self.merge_fasta_files(target_merged_faa, all_faa)
            concat_files(target_metadata_file, all_metadata)

            nogg_dmnd = f'{nogg_folder_path}NOGG_merged'
            if not os.path.exists(f'{nogg_dmnd}.dmnd'):
                run_command(f'diamond makedb --in {target_merged_faa} -d {nogg_dmnd}', stdout_file=stdout_file)

            stdout_file.close()

    def compile_NOG_DMND(self):
        if self.mantis_paths['NOG'][0:2] != 'NA':
            folder_path = add_slash(self.mantis_paths['NOG'])
            diamond_db = f'{folder_path}eggnog_proteins'
            eggnog_proteins_path = f'{folder_path}eggnog_seqs.faa'
            extract_seqs_command = f'diamond getseq -d {diamond_db} > {eggnog_proteins_path}'
            if not os.path.exists(eggnog_proteins_path):
                print('Extracting sequences from NOG Diamond database', flush=True, file=self.redirect_verbose)
                run_command(extract_seqs_command, join_command=True, shell=True)
            return self.create_fastas_NOG_DMND()

    def create_fastas_NOG_DMND(self, stdout_path=None):
        eggnog_db_path = self.mantis_paths['NOG'] + "eggnog.db"
        all_metadata = self.mantis_paths['NOG'] + 'metadata.tsv'
        stdout_file = open(self.mantis_out, 'a+')
        print('Generating metadata tsv', flush=True, file=stdout_file)
        seqs_taxons = self.generate_metadata_diamond(eggnog_db_path, all_metadata)
        print('Generating taxon specific metadata tsvs', flush=True, file=stdout_file)
        stdout_file.close()
        return seqs_taxons

    def clean_up_sql_results_ids_dmnd(self, sql_row, pfam_id_to_acc):
        codes_to_exclude = ['IEA', 'ND']
        res = {
            'go': set(),
            'pfam': set(),
            'kegg_ko': set(),
            'cog': set(),
            'arcog': set(),
            'enzyme_ec': set(),
            'kegg_pathway': set(),
            'kegg_module': set(),
            'kegg_reaction': set(),
            'kegg_rclass': set(),
            'kegg_brite': set(),
            'cazy': set(),
            'bigg_reaction': set(),
            'description': set(),
            'eggnog': set(),
            'tcdb': set(),
        }
        seq_name, ogs, gene_ontology_gos, pfam_ids, kegg_ko, kegg_cog, kegg_ec, kegg_brite, kegg_rclass, kegg_tc, kegg_cazy, kegg_pathway, kegg_module, kegg_reaction, kegg_go = sql_row
        if seq_name == '#member': return None, None, None
        eggnog_hmms = set()
        cogs = set()
        arcogs = set()
        taxons = set()
        if ogs:
            for query_hmm_taxon in ogs.split(','):
                query_hmm, query_taxon = query_hmm_taxon.split('@')
                eggnog_hmms.add(query_hmm)
                query_cogs = find_cog(query_hmm)
                query_arcogs = find_arcog(query_hmm)
                if query_cogs: cogs.update(query_cogs)
                if query_arcogs: arcogs.update(query_arcogs)
                taxons.add(query_taxon)
        if gene_ontology_gos:
            # http://geneontology.org/docs/guide-go-evidence-codes/
            gene_ontology_gos_copy = str(gene_ontology_gos)
            gene_ontology_gos_copy = gene_ontology_gos_copy.split(',')
            for go_group in gene_ontology_gos_copy:
                if '|' not in go_group: print(ids_command, '\n', go_group)
                _, go, evidence_code = go_group.split('|')
                if evidence_code not in codes_to_exclude:
                    res['go'].add(go.split(':')[1])
        if pfam_ids:
            temp_pfam_ids = pfam_ids.split(',')
            for tpi in temp_pfam_ids:
                if tpi in pfam_id_to_acc:
                    res['pfam'].add(pfam_id_to_acc[tpi])
        if kegg_ko:
            kegg_ko = kegg_ko.replace('ko:', '')
            res['kegg_ko'].update(kegg_ko.split(','))
        if kegg_cog: res['cog'].update(kegg_cog.split(','))
        if kegg_ec: res['enzyme_ec'].update(kegg_ec.split(','))
        if kegg_brite: res['kegg_brite'].update(kegg_brite.split(','))
        if kegg_rclass: res['kegg_rclass'].update(kegg_rclass.split(','))
        if kegg_tc: res['tcdb'].update(kegg_tc.split(','))
        if kegg_cazy: res['cazy'].update(kegg_cazy.split(','))
        if kegg_pathway: res['kegg_pathway'].update(kegg_pathway.split(','))
        if kegg_module: res['kegg_module'].update(kegg_module.split(','))
        if kegg_reaction: res['kegg_reaction'].update(kegg_reaction.split(','))
        if kegg_go: res['go'].update(kegg_go.split(','))
        if eggnog_hmms: res['eggnog'].update(eggnog_hmms)
        if cogs: res['cog'].update(cogs)
        if arcogs: res['arcog'].update(arcogs)

        return seq_name, taxons, res

    # we generate a metadata.tsv with the functional annotation of each sequence
    def generate_metadata_diamond(self, eggnog_db, target_sql_file):
        stdout_file = open(self.mantis_out, 'a+')
        codes_to_exclude = ['IEA', 'ND']
        res = {}
        # this will convert protein names to pfam ids (which is typically what is used with mantis)
        pfam_id_to_acc = self.pfam_id_to_acc()
        connection = sqlite3.connect(eggnog_db)
        sql_command = 'SELECT name,ogs, gos, pfam, kegg_ko,kegg_cog,kegg_ec,kegg_brite,kegg_rclass,kegg_tc,kegg_cazy,kegg_pathway,kegg_module,kegg_reaction,kegg_go FROM prots;'
        print(f'Querying SQL:\n{sql_command}', flush=True, file=stdout_file)
        cursor = connection.cursor()
        cursor.execute(sql_command)
        rows = self.yield_sql_rows(cursor)
        with open(target_sql_file, 'w+') as file:
            for row in rows:
                seq_name, taxons, row_info = self.clean_up_sql_results_ids_dmnd(row, pfam_id_to_acc)
                if seq_name:
                    for t in taxons:
                        if t not in res: res[t] = set()
                        res[t].add(seq_name)

                    link_line = [seq_name, '|']
                    for link_type in row_info:
                        for inner_link in row_info[link_type]:
                            inner_link = inner_link.strip()
                            if inner_link:
                                link_line.append(f'{link_type}:{inner_link}')
                    if len(link_line) > 2:
                        link_line = '\t'.join(link_line)
                        file.write(f'{link_line}\n')
        stdout_file.close()
        return res

    def prepare_queue_extract_fastas_DMND(self, seqs_taxons):
        if self.mantis_paths['NOG'][0:2] != 'NA':
            for taxon in self.get_ref_taxon_ids('NOG'):
                if taxon in seqs_taxons:
                    current_seqs = seqs_taxons[taxon]
                    taxon_folder = self.mantis_paths['NOG'] + taxon + SPLITTER
                    taxon_fasta = f'{taxon_folder}{taxon}_merged.faa'
                    taxon_dmnd = f'{taxon_folder}{taxon}_merged.dmnd'
                    taxon_metadata = f'{taxon_folder}metadata.tsv'
                    if not os.path.exists(taxon_fasta) or \
                            not os.path.exists(taxon_metadata) or \
                            not os.path.exists(taxon_dmnd):
                        self.queue.append([taxon, current_seqs, self.mantis_out])

    def worker_extract_fastas_DMND(self, queue, master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            folder_path = add_slash(self.mantis_paths['NOG'])
            all_proteins_path = f'{folder_path}eggnog_seqs.faa'
            all_metadata_path = f'{folder_path}metadata.tsv'
            taxon, taxon_seqs, stdout_path = record
            self.create_fasta_DMND(taxon, taxon_seqs, all_proteins_path, all_metadata_path, stdout_path)

    def create_fasta_DMND(self, taxon, taxon_seqs, all_proteins_path, all_metadata_path, stdout_path):

        stdout_file = open(stdout_path, 'a+')
        sequences_generator = yield_target_seqs(all_proteins_path, taxon_seqs)
        metadata_generator = yield_target_metadata(all_metadata_path, taxon_seqs)
        taxon_folder = self.mantis_paths['NOG'] + taxon + SPLITTER
        taxon_faa = f'{taxon_folder}{taxon}_merged.faa'
        taxon_metadata = f'{taxon_folder}metadata.tsv'

        Path(taxon_folder).mkdir(parents=True, exist_ok=True)

        with open(taxon_faa, 'w+') as seq_file:
            for seq in sequences_generator:
                seq_file.write(seq)

        with open(taxon_metadata, 'w+') as met_file:
            for met_seq in metadata_generator:
                met_file.write(met_seq)

        stdout_file.close()
