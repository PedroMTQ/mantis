import os

from mantis.src.metadata.metadata_sqlite_connector import MetadataSqliteConnector
from mantis.src.taxonomy.taxonomy_sqlite_connector import TaxonomySqliteConnector
from mantis.src.utils.logger import logger


class MantisInstaller(MetadataSqliteConnector, TaxonomySqliteConnector):
    def __init__(self, verbose=True,
                 redirect_verbose=None,
                 no_taxonomy=False,
                 mantis_config=None,
                 hmm_chunk_size=None,
                 keep_files=False,
                 user_cores=None):
        self.redirect_verbose = redirect_verbose
        self.keep_files = keep_files
        self.verbose = verbose
        if no_taxonomy:
            self.use_taxonomy = False
        else:
            self.use_taxonomy = True
        self.mantis_config = mantis_config
        self.user_cores = user_cores
        self.broken_merged_hmms = set()
        self.clean_merged_hmms = set()
        self.start_time = time()
        # to speed up hmm search we split very big hmms into smaller chunks - better job distribution
        if hmm_chunk_size == 0:
            self.hmm_chunk_size = None
        elif hmm_chunk_size is None:
            self.hmm_chunk_size = 5000
        else:
            self.hmm_chunk_size = hmm_chunk_size
        self.read_config_file()

        # I use manager instead of queue since I need to be able to add records to the end and start of the 'queue' (actually a list) which is not possible with the multiprocessing.Queue
        # we need manager.list because a normal list cant communicate during multiprocessing
        self.manager = Manager()
        self.queue = self.manager.list()




    def setup_paths_config_file(self):
        mantis_config = MantisConfig(config_file=self.config_file)
        self.mantis_paths = mantis_config.mantis_paths
        self.nog_db = mantis_config.nog_db
        # setting up which taxa we need to have references for
        Taxonomy_SQLITE_Connector.__init__(self, resources_folder=self.mantis_paths['resources'])
        if self.use_taxonomy:
            if mantis_config.nogt_line:
                if self.launch_taxonomy_connector():
                    self.set_nogt_line(mantis_config.nogt_line)




    def order_by_size_descending(self, refs_list):
        res = {}
        for ref in refs_list:
            res[ref] = os.stat(ref).st_size
        # mixing big and low size HMMs so that we try not to run out of memory, might lead to more idle time.
        sorted_res = sorted(res, key=res.get, reverse=True)
        resorted_res = []
        c = 1
        while sorted_res:
            if c == 1:
                resorted_res.append(sorted_res.pop(0))
                c = -1
            elif c == -1:
                resorted_res.append(sorted_res.pop(-1))
                c = 1
        return resorted_res

    def compile_refs_list(self, folder=False):
        # doesnt include NOG or NCBI
        refs_list = []
        default_list = [
            get_ref_in_folder(self.mantis_paths['pfam']) if not folder else self.mantis_paths['pfam'],
            get_ref_in_folder(self.mantis_paths['kofam']) if not folder else self.mantis_paths['kofam'],
            get_ref_in_folder(self.mantis_paths['tcdb']) if not folder else self.mantis_paths['tcdb'],
        ]
        for ref_path in self.get_custom_refs_paths(folder):
            if ref_path[0:2] != 'NA':
                refs_list.append(ref_path)
        for ref_path in default_list:
            if ref_path and ref_path[0:2] != 'NA':
                refs_list.append(ref_path)
        return refs_list

    #####SETTING UP DATABASE#####

    def get_path_default_ref(self, database, taxon_id=None):
        target_file = None
        if 'kofam' in database.lower():
            target_file = get_ref_in_folder(self.mantis_paths['kofam'])
        elif 'pfam' in database.lower():
            target_file = get_ref_in_folder(self.mantis_paths['pfam'])
        elif 'tcdb' in database.lower():
            target_file = get_ref_in_folder(self.mantis_paths['tcdb'])
        elif 'NOG'.lower() in database.lower():
            if not taxon_id: taxon_id = 'NOGG'
            target_file = get_ref_in_folder(self.mantis_paths['NOG'] + taxon_id)
        elif 'NCBI'.lower() in database.lower():
            if not taxon_id: taxon_id = 'NCBIG'
            target_file = get_ref_in_folder(self.mantis_paths['NCBI'] + taxon_id)
        return target_file


    #####LISTING HMMS DATABASE#####

    def check_installation_extras(self, res, verbose=True):
        ncbi_resources = add_slash(self.mantis_paths['resources'] + 'NCBI')
        essential_genes = f'{MANTIS_FOLDER}Resources{SPLITTER}essential_genes/essential_genes.txt'
        taxonomy_db = self.mantis_paths['resources'] + 'taxonomy.db'

        if verbose:
            yellow('Checking extra files', flush=True, file=self.redirect_verbose)

        if not os.path.exists(essential_genes):
            red('Essential genes list is missing, it should be in the github repo!')
            if verbose:
                red(f'Failed installation check on [files missing]: {essential_genes}',
                    flush=True,
                    file=self.redirect_verbose)
            res.append(self.mantis_paths['resources'] + 'essential_genes/')
        else:
            if verbose:
                green(f'Passed installation check on: {self.mantis_paths["resources"]}essential_genes',
                      flush=True, file=self.redirect_verbose)

        if not os.path.exists(ncbi_resources + 'gc.prt'):
            if verbose:
                red(f'Failed installation check on [files missing]: {ncbi_resources}gc.prt.dmp',
                    flush=True,
                    file=self.redirect_verbose)
            res.append(ncbi_resources)
        else:
            if verbose:
                green(f'Passed installation check on: {ncbi_resources}gc.prt.dmp',
                      flush=True,
                      file=self.redirect_verbose)

        if self.use_taxonomy:
            if not os.path.exists(taxonomy_db):
                if verbose:
                    red(f'Failed installation check on [files missing]: {taxonomy_db}',
                        flush=True,
                        file=self.redirect_verbose)
                res.append(taxonomy_db)
            else:
                if verbose:
                    green(f'Passed installation check on: {taxonomy_db}',
                          flush=True,
                          file=self.redirect_verbose)
        return res

    def check_chunks_dir(self, chunks_dir):
        all_chunks = []
        for hmm in os.listdir(chunks_dir):
            if hmm.endswith('.hmm'):
                all_chunks.append(hmm)
        for hmm in all_chunks:
            if not self.check_missing_chunk_files(hmm, chunks_dir):
                return False
        return True

    def check_missing_chunk_files(self, hmm, chunks_dir):
        missing_files = ['.h3f', '.h3i', '.h3m', '.h3p']
        res = 0
        for inner_file in os.listdir(chunks_dir):
            for mf in missing_files:
                if inner_file == f'{hmm}{mf}':
                    res += 1
        if res == len(missing_files):
            return True
        red(f'Failed installation check on [files missing]: {hmm} in chunks folder: {chunks_dir}',
            flush=True,
            file=self.redirect_verbose)
        return False

    def check_installation_folder(self, ref_folder_path, res, verbose=True, extra_requirements=[]):
        missing_files = set(extra_requirements)
        try:
            files_dir = os.listdir(ref_folder_path)
        except:
            if verbose:
                red(f'Failed installation check on [path unavailable]: {ref_folder_path}',
                    flush=True,
                    file=self.redirect_verbose)
            res.append(ref_folder_path)
            self.passed_check = False
            return
        ref_type = None
        for file in files_dir:
            if file.endswith('.dmnd'):
                ref_type = 'dmnd'
                missing_files.update(['.dmnd'])
            elif file.endswith('.hmm'):
                ref_type = 'hmm'
                missing_files.update(['.hmm', '.h3f', '.h3i', '.h3m', '.h3p'])

        if not ref_type:
            if verbose:
                red(f'Failed installation check on [invalid referecence type]: {ref_folder_path}',
                    flush=True,
                    file=self.redirect_verbose)
            res.append(ref_folder_path)
            self.passed_check = False
            return
        check = len(missing_files)

        if 'chunks' in files_dir:
            if not self.check_chunks_dir(f'{ref_folder_path}chunks'):
                self.passed_check = False
                return
            else:
                missing_files = set(extra_requirements)
                check = len(missing_files)
        for file in files_dir:
            if ref_type == 'hmm':
                if file.endswith('.hmm'):
                    if '.hmm' in missing_files:
                        check -= 1
                        missing_files.remove('.hmm')
                elif file.endswith('.h3f'):
                    if '.h3f' in missing_files:
                        check -= 1
                        missing_files.remove('.h3f')
                elif file.endswith('.h3i'):
                    if '.h3i' in missing_files:
                        check -= 1
                        missing_files.remove('.h3i')
                elif file.endswith('.h3m'):
                    if '.h3m' in missing_files:
                        check -= 1
                        missing_files.remove('.h3m')
                elif file.endswith('.h3p'):
                    if '.h3p' in missing_files:
                        check -= 1
                        missing_files.remove('.h3p')
            elif ref_type == 'dmnd':
                if file.endswith('.dmnd'):
                    check -= 1
                    missing_files.remove('.dmnd')
            if file in extra_requirements:
                check -= 1
                missing_files.remove(file)
        if check != 0:
            missing_files_str = '; '.join(missing_files)
            red(f'Failed installation check on [files missing]: {ref_folder_path}\n{missing_files_str}',
                flush=True,
                file=self.redirect_verbose)
            res.append(ref_folder_path)
        else:
            if verbose:
                green(f'Passed installation check on: {ref_folder_path}',
                      flush=True,
                      file=self.redirect_verbose)

    def compile_sql_metadata(self):
        all_files = set()
        for ref in self.compile_refs_list(folder=True):
            metadata_file = f'{ref}metadata.tsv'
            all_files.add(metadata_file)
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            ncbi_tax = self.get_taxon_refs('NCBI', folder=True)
            for ref in ncbi_tax:
                metadata_file = f'{ref}metadata.tsv'
                all_files.add(metadata_file)
        if self.mantis_paths['NOG'][0:2] != 'NA':
            nog_tax = self.get_taxon_refs('NOG', folder=True)
            for ref in nog_tax:
                metadata_file = f'{ref}metadata.tsv'
                all_files.add(metadata_file)
        for metadata_file in all_files:

            if not os.path.exists(metadata_file.replace('.tsv', '.db')):
                cursor = MetadataSqliteConnector(metadata_file)
                cursor.close_sql_connection()

    def check_sql_databases(self, ref_dbs):
        broken_refs = set()
        broken_ids = {}
        for db in ref_dbs:
            yellow(f'Checking {db}metadata.db', flush=True, file=self.redirect_verbose)
            cursor = MetadataSqliteConnector(f'{db}metadata.tsv')
            db_res = cursor.test_database()
            if db_res: broken_refs.add(db)
            if db_res: broken_ids[db] = db_res
            cursor.close_sql_connection()
        for db in broken_ids:
            red(f'Failed SQL check in {db} for the following IDs:\n{broken_ids[db]}', flush=True,
                file=self.redirect_verbose)
        if not broken_refs:
            green('------------------------------------------', flush=True, file=self.redirect_verbose)
            green('-------------SQL CHECK PASSED-------------', flush=True, file=self.redirect_verbose)
            green('------------------------------------------', flush=True, file=self.redirect_verbose)
        else:
            red('------------------------------------------', flush=True, file=self.redirect_verbose)
            red('-------------SQL CHECK FAILED-------------', flush=True, file=self.redirect_verbose)
            red('------------------------------------------', flush=True, file=self.redirect_verbose)
            raise SQLCheckNotPassed

    def check_installation(self, verbose=True, check_sql=False):
        # we use the verbose mode when running the check_installation directly
        self.compile_sql_metadata()
        self.passed_check = True
        ref_dbs = set()
        if not cython_compiled():
            self.passed_check = False
            if verbose:
                red('Cython needs to be compiled!', flush=True, file=self.redirect_verbose)
        else:
            if verbose:
                green('Cython correctly compiled!', flush=True, file=self.redirect_verbose)
        res = []
        res = self.check_installation_extras(res, verbose)

        if verbose:
            yellow('Checking references installation', flush=True, file=self.redirect_verbose)
        requirements = {
            self.mantis_paths['pfam']: ['metadata.tsv'],
            self.mantis_paths['tcdb']: ['metadata.tsv'],
            self.mantis_paths['kofam']: ['metadata.tsv'],
        }
        # per tax level FOR EGGNOG
        if self.mantis_paths['NOG'][0:2] != 'NA':
            tax_refs = self.get_taxon_refs(db='NOG', folder=True)
            if not tax_refs:
                if verbose:
                    red(F'Failed installation check on [path unavailable]: {self.mantis_paths["NOG"]}',
                        flush=True,
                        file=self.redirect_verbose)
                res.append(self.mantis_paths['NOG'])
            for tax_ref_folder in tax_refs:
                self.check_installation_folder(tax_ref_folder, res, verbose=False, extra_requirements=['metadata.tsv'])
                ref_dbs.add(tax_ref_folder)
            nogt_check = [i for i in res if self.mantis_paths['NOG'] in i]
            if not nogt_check:
                if verbose: green('Passed installation check on: ' + self.mantis_paths['NOG'], flush=True,
                                  file=self.redirect_verbose)

        # per tax level FOR NCBI
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            # checking those already present
            tax_refs = self.get_taxon_refs(db='NCBI', folder=True)
            if not tax_refs:
                if verbose:
                    red('Failed installation check on [path unavailable]: self.mantis_paths["NCBI"]',
                        flush=True,
                        file=self.redirect_verbose)
                res.append(self.mantis_paths['NCBI'])
            for tax_ref_folder in tax_refs:
                # we skip the taxon 1 since it has no hmms
                self.check_installation_folder(tax_ref_folder, res, verbose=False, extra_requirements=['metadata.tsv'])
                ref_dbs.add(tax_ref_folder)

            ncbi_check = [i for i in res if self.mantis_paths['NCBI'] in i]
            if not ncbi_check:
                if verbose: green('Passed installation check on: ' + self.mantis_paths['NCBI'], flush=True,
                                  file=self.redirect_verbose)
        for ref_folder in self.compile_refs_list(folder=True):
            if ref_folder in requirements:
                self.check_installation_folder(ref_folder, res, verbose, extra_requirements=requirements[ref_folder])
            else:
                self.check_installation_folder(ref_folder, res, verbose)
            ref_dbs.add(ref_folder)
        if res:
            self.passed_check = False
            fail_res = ''
            for i in res:
                fail_res += f'{i}\n'
            if verbose:
                red(f'Installation check failed on:\n{fail_res}', flush=True, file=self.redirect_verbose)
        if self.passed_check:
            logger.info('Installation check passed', flush=True, file=self.redirect_verbose)

        else:
            logger.error('Installation check failed', flush=True, file=self.redirect_verbose)
            raise InstallationCheckNotPassed
        if check_sql: self.check_sql_databases(ref_dbs)


    def get_taxon_ref_path(self, taxon_id, db):
        tax_refs = self.get_local_ref_taxon_ids(db=db)
        if taxon_id in tax_refs:
            if db == 'NOG' and self.nog_db == 'dmnd':
                return add_slash(self.mantis_paths[db] + taxon_id) + f'{taxon_id}.dmnd'
            else:
                return add_slash(self.mantis_paths[db] + taxon_id) + f'{taxon_id}_merged.hmm'
        else:
            return None

    def get_ref_taxon_ids(self, db):
        res = set()
        if not os.path.exists(self.mantis_paths[db]): return res
        if db == 'NOG':
            available_taxon_ids = self.get_taxon_ids_eggNOG()
            if self.mantis_nogt_tax:
                for tax_id in self.mantis_nogt_tax:
                    if tax_id in available_taxon_ids:
                        res.add(tax_id)
                return res
            else:
                return set(available_taxon_ids)
        else:
            for i in os.listdir(self.mantis_paths[db]):
                if re.search(r'\d+', i): res.add(i)
            return res

    def get_local_ref_taxon_ids(self, db):
        res = set()
        if os.path.exists(self.mantis_paths[db]):
            if db == 'NOG':
                if self.mantis_nogt_tax:
                    for i in self.mantis_nogt_tax:
                        res.add(i)
            for i in os.listdir(self.mantis_paths[db]):
                if re.search(r'\d+', i): res.add(i)
        return res

    def get_taxon_refs(self, db, folder=False):
        # locally available taxon ids
        local_taxon_ids = self.get_local_ref_taxon_ids(db)
        # all taxon ids
        taxon_ids = self.get_ref_taxon_ids(db)
        res = []
        for t in taxon_ids:
            if t in local_taxon_ids:
                if folder:
                    res.append(add_slash(self.mantis_paths[db] + t))
                else:
                    if self.nog_db == 'hmm':
                        res.append(add_slash(self.mantis_paths[db] + t) + f'{t}_merged.hmm')
                    else:
                        res.append(add_slash(self.mantis_paths[db] + t) + f'{t}_merged.dmnd')
        global_folder = add_slash(self.mantis_paths[db] + db + 'G')
        if folder:
            if os.path.exists(global_folder):
                res.append(global_folder)
        else:
            if self.nog_db == 'hmm':
                if os.path.exists(f'{global_folder}{db}G_merged.hmm'):
                    res.append(f'{global_folder}{db}G_merged.hmm')
            else:
                if os.path.exists(f'{global_folder}{db}G_merged.dmnd'):
                    res.append(f'{global_folder}{db}G_merged.dmnd')

        return res

    def processes_handler(self, target_worker_function, worker_count, add_sentinels=True):
        '''
        this will first generate one process per worker, then we add sentinels to the end of the list which will basically tell us when the queue is empty
        if we need to add new work (e.g. when doing taxa annotation) we just add the new work to the start of the list
        '''
        # os.getpid to add the master_pid
        processes = [Process(target=target_worker_function,
                             args=(self.queue, os.getpid(),)) for _ in range(worker_count)]
        # adding sentinel record since queue can be signaled as empty when its really not
        if add_sentinels:
            for _ in range(worker_count):   self.queue.append(None)
        for process in processes:
            process.start()
        # we could manage the processes memory here with a while cycle
        for process in processes:
            process.join()
            # exitcode 0 for sucessful exists
            if process.exitcode != 0:
                sleep(5)
                print('Ran into an issue, check the log for details. Exitting!')
                os._exit(1)

