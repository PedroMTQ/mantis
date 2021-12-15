try:
    from source.Database_generator import *
    from source.Taxonomy_SQLITE_Connector import Taxonomy_SQLITE_Connector
    from source.Metadata_SQLITE_Connector import Metadata_SQLITE_Connector

except:
    from Database_generator import *
    from Taxonomy_SQLITE_Connector import Taxonomy_SQLITE_Connector
    from Metadata_SQLITE_Connector import Metadata_SQLITE_Connector



def setup_databases(force_download=False, chunk_size=None,no_taxonomy=False,mantis_config=None,cores=None):
    print_cyan('Setting up databases')
    if force_download == 'None': force_download = None
    if chunk_size: chunk_size = int(chunk_size)
    if cores: cores=int(cores)
    mantis = Assembler(hmm_chunk_size=chunk_size, mantis_config=mantis_config,user_cores=cores,no_taxonomy=no_taxonomy)
    mantis.setup_databases(force_download)




def check_installation(mantis_config=None,no_taxonomy=False,check_sql=False):
    yellow('Checking installation')
    mantis = Assembler(mantis_config=mantis_config,no_taxonomy=no_taxonomy)
    mantis.check_installation(check_sql=check_sql)



class Assembler(Database_generator,Taxonomy_SQLITE_Connector):
    def __init__(self, verbose=True, redirect_verbose=None,no_taxonomy=False, mantis_config=None,
                 hmm_chunk_size=None,keep_files=False,user_cores=None):
        self.redirect_verbose = redirect_verbose
        self.keep_files = keep_files
        self.verbose = verbose
        if no_taxonomy: self.use_taxonomy = False
        else:           self.use_taxonomy=True
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
        Taxonomy_SQLITE_Connector.__init__(self,resources_folder=self.mantis_paths['resources'])

        #self.requirements_met()
        # I use manager instead of queue since I need to be able to add records to the end and start of the 'queue' (actually a list) which is not possible with the multiprocessing.Queue
        # we need manager.list because a normal list cant communicate during multiprocessing
        self.manager = Manager()
        self.queue = self.manager.list()


    def __str__(self):
        custom_refs = self.get_custom_refs_paths(folder=True)
        custom_refs_str = ''
        custom_res=''
        for cref in custom_refs:
            custom_refs_str += cref + '\n'
        if custom_refs_str:
            custom_res = f'# Custom references:\n{custom_refs_str}'

        res=[]
        if hasattr(self, 'output_folder'):
            res.append(f'Output folder:\nself.output_folder')
        res.append(f'Default references folder:\n{self.mantis_paths["default"]}')
        res.append(f'Custom references folder:\n{self.mantis_paths["custom"]}')
        if self.mantis_paths['NOG'][0:2] != 'NA':
            res.append(f'TAX NOG references folder:\n{self.mantis_paths["NOG"]}')
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            res.append(f'TAX NCBI references folder:\n{self.mantis_paths["NCBI"]}')
        if self.mantis_paths['pfam'][0:2] != 'NA':
            res.append(f'Pfam reference folder:\n{self.mantis_paths["pfam"]}')
        if self.mantis_paths['kofam'][0:2] != 'NA':
            res.append(f'KOfam reference folder:\n{self.mantis_paths["kofam"]}')
        if self.mantis_paths['tcdb'][0:2] != 'NA':
            res.append(f'TCDB reference folder:\n{self.mantis_paths["tcdb"]}')
        res.append('------------------------------------------')
        res='\n'.join(res)
        if custom_res: res+='\n'+custom_res
        ref_weights=', '.join([f'{i}:{self.mantis_ref_weights[i]}' for i in self.mantis_ref_weights if i!='else'])
        if ref_weights:
            res+= f'#  Weights:\n{ref_weights}\n'

        nog_tax=', '.join([i for i in self.mantis_nogt_tax])
        if nog_tax:
            res+= f'#  NOG tax IDs:\n{nog_tax}\n'

        return res

    def requirements_met(self):
        for f in [self.is_conda_available(), self.is_hmmer_available()]:
            if not f:
                kill_switch(RequirementsNotMet)

    def is_conda_available(self):
        process = run_command('conda -V', get_output=True)
        check = re.search('conda', str(process.stdout))
        if not check:
            print_cyan('Conda dependency not met!')
            print_cyan('Install Conda on your system by following the instructions at: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html')
        return check

    def is_hmmer_available(self):
        check_command = ' conda list hmmer '
        process = run_command(check_command, get_output=True)
        check = re.search('hmmer', str(process.stdout))
        if not check:
            print_cyan('HMMER dependency not met!')
            print_cyan('Install HMMER on your conda environment by doing the following:')
            print_cyan('conda activate <conda_environment>')
            print_cyan('conda install -c bioconda hmmer')
        return check


    def check_internet_connection(self):
        try:
            requests.get("http://www.google.com")
            return True
        except requests.ConnectionError:
            print("Could not connect to internet!\nIf you would like to run offline make sure you introduce organism NCBI IDs instead of synonyms!")
            return False

    def get_default_ref_path(self):
        file = open(self.config_file, 'r')
        line = file.readline()
        while line:
            line = line.strip('\n')
            if '#' not in line:
                # data sources configuration
                if 'default_ref_folder=' in line:
                    line_path = add_slash(line.replace('default_ref_folder=', ''))
                    if line_path:
                        default_ref_path = line_path
                        return default_ref_path
            line = file.readline()
        file.close()

    def set_nogt_line(self, line_path):
        if line_path:
            res = set()
            tax_ids = [i for i in line_path.split(',')]
            if tax_ids:
                for t_id in tax_ids:
                    try:
                        ncbi_taxon_id = int(t_id)
                        organism_lineage = self.fetch_ncbi_lineage(ncbi_taxon_id)
                        res.update(organism_lineage)
                    except:
                        ncbi_taxon_id = self.get_taxa_ncbi(t_id)
                        if ncbi_taxon_id:
                            organism_lineage = self.fetch_ncbi_lineage(ncbi_taxon_id)
                            res.update(organism_lineage)
            for i in res:
                self.mantis_nogt_tax.add(str(i))

    def setup_paths_config_file(self):
        self.nog_db = 'dmnd'
        file = open(self.config_file, 'r')
        line = file.readline()
        nogt_line=None
        while line:
            line = line.strip('\n')
            if not line.startswith('#') and line:
                # data sources configuration
                if line.startswith('custom_ref_folder='):
                    line_path = add_slash(line.replace('custom_ref_folder=', ''))
                    if line_path: self.mantis_paths['custom'] = line_path

                elif line.startswith('resources_folder='):
                    line_path = add_slash(line.replace('resources_folder=', ''))
                    if line_path:
                        self.mantis_paths['resources'] = line_path

                elif line.startswith('nog_ref_folder='):
                    line_path = add_slash(line.replace('nog_ref_folder=', ''))
                    if line_path: self.mantis_paths['NOG'] = line_path

                # taxa ids list for only downloading nogt specific to lineage
                elif line.startswith('nog_tax='):
                    nogt_line = line.replace('nog_tax=', '')

                elif line.startswith('pfam_ref_folder='):
                    line_path = add_slash(line.replace('pfam_ref_folder=', ''))
                    if line_path: self.mantis_paths['pfam'] = line_path

                elif line.startswith('kofam_ref_folder='):
                    line_path = add_slash(line.replace('kofam_ref_folder=', ''))
                    if line_path: self.mantis_paths['kofam'] = line_path

                elif line.startswith('ncbi_ref_folder='):
                    line_path = add_slash(line.replace('ncbi_ref_folder=', ''))
                    if line_path: self.mantis_paths['NCBI'] = line_path

                elif line.startswith('tcdb_ref_folder='):
                    line_path = add_slash(line.replace('tcdb_ref_folder=', ''))
                    if line_path: self.mantis_paths['tcdb'] = line_path


                elif line.startswith('_weight='):
                    ref_source, weight = line.split('_weight=')
                    self.mantis_ref_weights[ref_source] = float(weight)

                elif line.startswith('nog_ref='):
                    nog_db = line.replace('nog_ref=', '').split()[0]
                    if nog_db.lower() not in ['dmnd','hmm']:
                        kill_switch(InvalidNOGType)
                    else:
                        self.nog_db=nog_db
            line = file.readline()
        file.close()
        if self.use_taxonomy:
            if nogt_line:
                if self.launch_taxonomy_connector():
                    self.set_nogt_line(nogt_line)


    def read_config_file(self):
        self.mantis_ref_weights = {'else': 0.7}
        self.mantis_nogt_tax = set()

        if self.mantis_config:
            print(f'Using custom MANTIS.config: {self.mantis_config}', flush=True, file=self.redirect_verbose)
            self.config_file = self.mantis_config
        else:
            if not os.path.isdir(MANTIS_FOLDER):
                print('Make sure you are calling the folder to run this package, like so:\n python mantis/ <command>\n ',
                    flush=True, file=self.redirect_verbose)
                raise FileNotFoundError
            self.config_file = MANTIS_FOLDER + 'MANTIS.config'
        try:
            open(self.config_file, 'r')
        except:
            print('MANTIS.config file has been deleted or moved, make sure you keep it in the root of the project!',
                  flush=True, file=self.redirect_verbose)
            raise FileNotFoundError

        default_ref_path = self.get_default_ref_path()

        # if there's no path, we just assume its in the default folder
        if not default_ref_path:  default_ref_path = add_slash(MANTIS_FOLDER + 'References')
        resources_path = add_slash(MANTIS_FOLDER + 'Resources')
        self.mantis_paths = {'default': default_ref_path,
                             'resources': resources_path,
                             'custom': add_slash(default_ref_path + 'Custom_references'),
                             'NOG': add_slash(default_ref_path + 'NOG'),
                             'pfam': add_slash(default_ref_path + 'pfam'),
                             'kofam': add_slash(default_ref_path + 'kofam'),
                             'NCBI': add_slash(default_ref_path + 'NCBI'),
                             'tcdb': add_slash(default_ref_path + 'tcdb'),
                             }
        self.setup_paths_config_file()
        if not self.use_taxonomy:
            self.mantis_paths['NOG']=f'NA{SPLITTER}'
            self.mantis_paths['NCBI']=f'NA{SPLITTER}'
        if not os.path.isdir(self.mantis_paths['custom']):
            Path(self.mantis_paths['custom']).mkdir(parents=True, exist_ok=True)
        if self.verbose: print(self, flush=True, file=self.redirect_verbose)


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

    def check_reference_exists(self, database, taxon_id=None, force_download=False):
        ncbi_resources=add_slash(self.mantis_paths['resources']+'NCBI')

        if database == 'ncbi_res':
            if file_exists(ncbi_resources + 'gc.prt', force_download) and \
                    file_exists(ncbi_resources + 'gc.prt', force_download):
                return True
        elif database == 'taxonomy':
            taxonomy_db=self.mantis_paths['resources'] + 'Taxonomy.db'
            gtdb_resources = add_slash(self.mantis_paths['resources'] + 'GTDB')
            if file_exists(taxonomy_db, force_download):
                return True
        elif database == 'NOGSQL':
            if file_exists(self.mantis_paths['NOG'] + 'eggnog.db', force_download):
                return True
        elif database == 'tcdb':
            if file_exists(self.mantis_paths['tcdb'] + 'tcdb.dmnd', force_download):
                return True
        elif database == 'NOG_DMND':
            if file_exists(self.mantis_paths['NOG'] + 'eggnog_proteins.dmnd', force_download):
                return True
        target_file = self.get_path_default_ref(database, taxon_id)
        if target_file:
            if target_file.endswith('.dmnd'):
                if not file_exists(target_file, force_download=force_download):
                    return False
            else:
                for extension in ['', '.h3f', '.h3i', '.h3m', '.h3p']:
                    if not file_exists(target_file + extension, force_download=force_download):
                        return False
        else:
            return False
        return True

    #####LISTING HMMS DATABASE#####


    def check_installation_extras(self, res, verbose=True):
        ncbi_resources=add_slash(self.mantis_paths['resources']+'NCBI')
        essential_genes = f'{MANTIS_FOLDER}Resources{SPLITTER}essential_genes/essential_genes.txt'
        taxonomy_db=self.mantis_paths['resources']+'Taxonomy.db'

        if verbose: yellow('Checking extra files', flush=True, file=self.redirect_verbose)

        if not file_exists(essential_genes):
            red('Essential genes list is missing, it should be in the github repo!')
            if verbose: red('Failed installation check on [files missing]: ' + essential_genes, flush=True, file=self.redirect_verbose)
            res.append(self.mantis_paths['resources'] + 'essential_genes/')
        else:
            if verbose: green('Passed installation check on: ' + self.mantis_paths['resources'] + 'essential_genes',flush=True, file=self.redirect_verbose)

        if not file_exists(ncbi_resources + 'gc.prt'):
            if verbose: red('Failed installation check on [files missing]: ' + ncbi_resources + 'gc.prt.dmp',flush=True, file=self.redirect_verbose)
            res.append(ncbi_resources)
        else:
            if verbose: green('Passed installation check on: ' + ncbi_resources+'gc.prt.dmp', flush=True,file=self.redirect_verbose)


        if self.use_taxonomy:
            if not file_exists(taxonomy_db):
                if verbose: red(f'Failed installation check on [files missing]: {taxonomy_db}',flush=True, file=self.redirect_verbose)
                res.append(taxonomy_db)
            else:
                if verbose: green(f'Passed installation check on: {taxonomy_db}', flush=True,file=self.redirect_verbose)
        return res

    def check_chunks_dir(self,chunks_dir):
        all_chunks=[]
        for hmm in os.listdir(chunks_dir):
            if hmm.endswith('.hmm'):
                all_chunks.append(hmm)
        for hmm in all_chunks:
            if not self.check_missing_chunk_files(hmm,chunks_dir):
                return False
        return True

    def check_missing_chunk_files(self,hmm,chunks_dir):
        missing_files=['.h3f', '.h3i', '.h3m', '.h3p']
        res=0
        for inner_file in os.listdir(chunks_dir):
            for mf in missing_files:
                if inner_file==f'{hmm}{mf}':
                    res+=1
        if res==len(missing_files): return True
        red(f'Failed installation check on [files missing]: {hmm} in chunks folder: {chunks_dir}',
            flush=True, file=self.redirect_verbose)
        return False


    def check_installation_folder(self, ref_folder_path, res, verbose=True, extra_requirements=[]):
        missing_files = set(extra_requirements)

        try:
            files_dir = os.listdir(ref_folder_path)
        except:
            if verbose: red(f'Failed installation check on [path unavailable]: {ref_folder_path}', flush=True,file=self.redirect_verbose)
            res.append(ref_folder_path)
            self.passed_check = False
            return
        ref_type=None
        for file in files_dir:
            if file.endswith('.dmnd'):
                ref_type='dmnd'
                missing_files.update(['.dmnd'])
            elif file.endswith('.hmm'):
                ref_type='hmm'
                missing_files.update(['.hmm', '.h3f', '.h3i', '.h3m', '.h3p'])

        if not ref_type:
            if verbose: red(f'Failed installation check on [invalid referecence type]: {ref_folder_path}', flush=True,file=self.redirect_verbose)
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
            if ref_type=='hmm':
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
            elif ref_type=='dmnd':
                if file.endswith('.dmnd'):
                    check -= 1
                    missing_files.remove('.dmnd')
            if file in extra_requirements:
                check -= 1
                missing_files.remove(file)
        if check != 0:
            missing_files_str = '; '.join(missing_files)
            red(f'Failed installation check on [files missing]: {ref_folder_path}\n{missing_files_str}',
                flush=True, file=self.redirect_verbose)
            res.append(ref_folder_path)
        else:
            if verbose: green(f'Passed installation check on: {ref_folder_path}', flush=True,
                              file=self.redirect_verbose)

    def compile_sql_metadata(self):
        all_files=set()
        for ref in self.compile_refs_list(folder=True):
            metadata_file=f'{ref}metadata.tsv'
            all_files.add(metadata_file)
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            ncbi_tax=self.get_taxon_refs('NCBI',folder=True)
            for ref in ncbi_tax:
                metadata_file = f'{ref}metadata.tsv'
                all_files.add(metadata_file)
        if self.mantis_paths['NOG'][0:2] != 'NA':
            nog_tax=self.get_taxon_refs('NOG',folder=True)
            for ref in nog_tax:
                metadata_file = f'{ref}metadata.tsv'
                all_files.add(metadata_file)
        for metadata_file in all_files:

            if not file_exists(metadata_file.replace('.tsv','.db')):
                cursor = Metadata_SQLITE_Connector(metadata_file)
                cursor.close_sql_connection()

    def check_sql_databases(self,ref_dbs):
        broken_refs=set()
        broken_ids={}
        for db in ref_dbs:
            yellow(f'Checking {db}metadata.db', flush=True,file=self.redirect_verbose)
            cursor = Metadata_SQLITE_Connector(f'{db}metadata.tsv')
            db_res=cursor.test_database()
            if db_res: broken_refs.add(db)
            if db_res: broken_ids[db]=db_res
            cursor.close_sql_connection()
        for db in broken_ids:
            red(f'Failed SQL check in {db} for the following IDs:\n{broken_ids[db]}', flush=True,file=self.redirect_verbose)
        if not broken_refs:
            green('------------------------------------------', flush=True, file=self.redirect_verbose)
            green('-------------SQL CHECK PASSED-------------', flush=True, file=self.redirect_verbose)
            green('------------------------------------------', flush=True, file=self.redirect_verbose)
        else:
            red('------------------------------------------', flush=True, file=self.redirect_verbose)
            red('-------------SQL CHECK FAILED-------------', flush=True, file=self.redirect_verbose)
            red('------------------------------------------', flush=True, file=self.redirect_verbose)


    def check_installation(self, verbose=True,check_sql=False):
        # we use the verbose mode when running the check_installation directly
        self.compile_sql_metadata()
        self.passed_check = True
        ref_dbs=set()
        if not cython_compiled():
            self.passed_check = False
            if verbose: red('Cython needs to be compiled!', flush=True,
                            file=self.redirect_verbose)
        else:
            if verbose: green('Cython correctly compiled!', flush=True, file=self.redirect_verbose)
        res = []
        res = self.check_installation_extras(res, verbose)

        if verbose: yellow('Checking references installation', flush=True, file=self.redirect_verbose)
        requirements = {
            self.mantis_paths['pfam']: ['metadata.tsv'],
            self.mantis_paths['tcdb']: ['metadata.tsv'],
            self.mantis_paths['kofam']: ['metadata.tsv'],
        }
        # per tax level FOR EGGNOG
        if self.mantis_paths['NOG'][0:2] != 'NA':
            tax_refs = self.get_taxon_refs(db='NOG', folder=True)
            if not tax_refs:
                if verbose: red('Failed installation check on [path unavailable]: ' + self.mantis_paths['NOG'],flush=True, file=self.redirect_verbose)
                res.append(self.mantis_paths['NOG'])
            for tax_ref_folder in tax_refs:
                self.check_installation_folder(tax_ref_folder, res, verbose=False,extra_requirements=['metadata.tsv'])
                ref_dbs.add(tax_ref_folder)
            nogt_check = [i for i in res if self.mantis_paths['NOG'] in i]
            if not nogt_check:
                if verbose: green('Passed installation check on: ' + self.mantis_paths['NOG'], flush=True,file=self.redirect_verbose)

        # per tax level FOR NCBI
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            #checking those already present
            tax_refs = self.get_taxon_refs(db='NCBI', folder=True)
            if not tax_refs:
                if verbose: red('Failed installation check on [path unavailable]: ' + self.mantis_paths['NCBI'],
                                flush=True, file=self.redirect_verbose)
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
            for i in res: fail_res += f'{i}\n'
            if verbose: red(f'Installation check failed on:\n{fail_res}', flush=True, file=self.redirect_verbose)
        if self.passed_check:
            if verbose:
                green('------------------------------------------', flush=True, file=self.redirect_verbose)
                green('--------INSTALLATION CHECK PASSED!--------', flush=True, file=self.redirect_verbose)
                green('------------------------------------------', flush=True, file=self.redirect_verbose)
            else:
                print_cyan('Installation check passed', flush=True, file=self.redirect_verbose)

        else:
            if verbose:
                yellow('------------------------------------------', flush=True, file=self.redirect_verbose)
                red('--------INSTALLATION CHECK FAILED!--------', flush=True, file=self.redirect_verbose)
                yellow('------------------------------------------', flush=True, file=self.redirect_verbose)
            else:
                print_cyan('Installation check failed', flush=True, file=self.redirect_verbose)
        if check_sql: self.check_sql_databases(ref_dbs)


    def get_custom_refs_paths(self, folder=False):
        try:
            custom_refs_folders = os.listdir(self.mantis_paths['custom'])
            for potential_ref_folder in custom_refs_folders:
                try:
                    files = os.listdir(self.mantis_paths['custom'] + potential_ref_folder)
                    for potential_file in files:
                        if potential_file.endswith('.hmm') or potential_file.endswith('.dmnd'):
                            if folder:
                                try:
                                    yield add_slash(self.mantis_paths['custom'] + potential_ref_folder)
                                except GeneratorExit:
                                    return ''
                            else:
                                try:
                                    yield add_slash(self.mantis_paths['custom'] + potential_ref_folder) + potential_file
                                except GeneratorExit:
                                    return ''
                except:
                    pass
        except:
            print('Custom references folder is missing, did you correctly set the path? If path is not set make sure you didn\'t delete the custom_ref folder!',
                flush=True, file=self.redirect_verbose)
            self.passed_check = False
            return
        with open(self.config_file, 'r') as file:
            line = file.readline()
            while line:
                if line[0] != '#':
                    if 'custom_ref=' in line:
                        line = line.strip('\n')
                        ref_path=line.replace('custom_ref=', '')
                        if not (ref_path.endswith('.hmm') or ref_path.endswith('.dmnd')):
                            if os.path.isdir(ref_path):
                                for inner_file in os.listdir(ref_path):
                                    if inner_file.endswith('.hmm') or inner_file.endswith('.dmnd'):
                                        ref_path=add_slash(ref_path)+inner_file
                        if folder:
                            try:
                                yield add_slash(SPLITTER.join(ref_path.split(SPLITTER)[:-1]))
                            except GeneratorExit:
                                return ''
                        else:
                            try:
                                yield ref_path
                            except GeneratorExit:
                                return ''
                line = file.readline()

    def get_taxon_ref_path(self, taxon_id, db):
        tax_refs = self.get_local_ref_taxon_ids(db=db)
        if taxon_id in tax_refs:
            if db=='NOG' and self.nog_db == 'dmnd':
                return add_slash(self.mantis_paths[db] + taxon_id) + f'{taxon_id}.dmnd'
            else:
                return add_slash(self.mantis_paths[db] + taxon_id) + f'{taxon_id}_merged.hmm'
        else:
            return None


    def get_ref_taxon_ids(self, db):
        res = set()
        if not file_exists(self.mantis_paths[db]): return res
        if db=='NOG':
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
                if re.search('\d+',i): res.add(i)
            return res

    def get_local_ref_taxon_ids(self,db):
        res = set()
        if file_exists(self.mantis_paths[db]):
            if db=='NOG':
                if self.mantis_nogt_tax:
                    for i in self.mantis_nogt_tax:
                        res.add(i)
            for i in os.listdir(self.mantis_paths[db]):
                if re.search('\d+', i): res.add(i)
        return res

    def get_taxon_refs(self, db, folder=False):
        #locally available taxon ids
        local_taxon_ids = self.get_local_ref_taxon_ids(db)
        #all taxon ids
        taxon_ids=self.get_ref_taxon_ids(db)
        res = []
        for t in taxon_ids:
            if t in local_taxon_ids:
                if folder:
                    res.append(add_slash(self.mantis_paths[db] + t))
                else:
                    if self.nog_db=='hmm':
                        res.append(add_slash(self.mantis_paths[db] + t) + f'{t}_merged.hmm')
                    else:
                        res.append(add_slash(self.mantis_paths[db] + t) + f'{t}_merged.dmnd')
        global_folder = add_slash(self.mantis_paths[db] + db + 'G')
        if folder:
            if file_exists(global_folder):
                res.append(global_folder)
        else:
            if self.nog_db == 'hmm':
                if file_exists(f'{global_folder}{db}G_merged.hmm'):
                    res.append( f'{global_folder}{db}G_merged.hmm')
            else:
                if file_exists(f'{global_folder}{db}G_merged.dmnd'):

                    res.append( f'{global_folder}{db}G_merged.dmnd')


        return res


    def processes_handler(self, target_worker_function, worker_count, add_sentinels=True):
        '''
        this will first generate one process per worker, then we add sentinels to the end of the list which will basically tell us when the queue is empty
        if we need to add new work (e.g. when doing taxa annotation) we just add the new work to the start of the list
        '''
        # os.getpid to add the master_pid
        processes = [Process(target=target_worker_function, args=(self.queue, os.getpid(),)) for _ in
                     range(worker_count)]
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


if __name__ == '__main__':
    p = Assembler(mantis_config='/media/HDD/data/mantis_references/MANTIS.config')
