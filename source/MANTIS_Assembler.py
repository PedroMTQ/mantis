try:
    from source.MANTIS_DB import *
except:
    from MANTIS_DB import *

'''
This module handles the setup of MANTIS - setting paths and setting up the databases.
The assembler is inherited by the main MANTIS class

check broken installation slurm output
wc -l *.out | sort -rn | awk  '$1 > 18'
grep "Disk quota exceeded" *

'''


def setup_databases(force_download=False, chunk_size=None, mantis_config=None,cores=None):
    print_cyan('Setting up databases')
    if force_download == 'None': force_download = None
    if chunk_size: chunk_size = int(chunk_size)
    if cores: cores=int(cores)
    mantis = MANTIS_Assembler(hmm_chunk_size=chunk_size, mantis_config=mantis_config,user_cores=cores)
    mantis.setup_databases(force_download)


def merge_hmm_folder(target_folder):
    print_cyan('Merging HMM folder')
    s_target_folder = add_slash(target_folder)
    if target_folder == 'None': return
    mantis = MANTIS_Assembler()
    mantis.merge_hmm_folder(target_folder)


def check_installation(mantis_config=None):
    yellow('Checking installation')
    mantis = MANTIS_Assembler(mantis_config=mantis_config)
    mantis.check_installation()


def extract_nog_metadata(metadata_path):
    yellow('Extracting NOG metadata')
    mantis = MANTIS_Assembler()
    mantis.extract_nog_metadata(metadata_path)


class MANTIS_Assembler(MANTIS_DB):
    def __init__(self, verbose=True, redirect_verbose=None, mantis_config=None,
                 hmm_chunk_size=None,keep_files=False,user_cores=None):
        self.redirect_verbose = redirect_verbose
        self.keep_files = keep_files
        self.verbose = verbose
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
            custom_res = 'Custom references:\n' + custom_refs_str
        res= 'Output folder:\n' + self.output_folder if hasattr(self, 'output_folder') else '' + \
                  '#  External data folders:' + '\n' + \
                  '------------------------------------------' + '\n' + \
                  'Default references folder:\n' + \
                  self.mantis_paths['default'] + '\n' + \
                  'Custom references folder:\n' + \
                  self.mantis_paths['custom'] + '\n' + \
                  'TAX NOG HMMs folder:\n' + \
                  self.mantis_paths['NOG'] + '\n' + \
                  'TAX NCBI HMMs folder:\n' + \
                  self.mantis_paths['NCBI'] + '\n' + \
                  'Pfam HMMs folder:\n' + \
                  self.mantis_paths['pfam'] + '\n' + \
                  'KOfam HMMs folder:\n' + \
                  self.mantis_paths['kofam'] + '\n' + \
                  'TIGRFAM HMMs folder:\n' + \
                  self.mantis_paths['tigrfam'] + '\n' + \
                  'TCDB sequences folder:\n' + \
                  self.mantis_paths['tcdb'] + '\n' + \
                '------------------------------------------' + '\n'
        if custom_res: res+='\n'+custom_res
        ref_weights=', '.join([f'{i}:{self.mantis_ref_weights[i]}' for i in self.mantis_ref_weights if i!='else'])
        res+= f'#  Weights: {ref_weights}\n'

        return res

    def print_citation_mantis(self):
        paper_doi=''
        separator='#####################################################################################################################'
        res=f'{separator}\n# Thank you for using Mantis, please make sure you cite the respective paper {paper_doi} #\n{separator}'
        print(res)


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

    def get_taxa_ncbi_url(self,url):
        webpage = None
        c = 0
        while not webpage and c <= 10:
            req = requests.get(url)
            try:
                webpage = req.text
                taxa_id = re.search('<Id>\d+</Id>', webpage)
                return re.search('\d+', taxa_id.group()).group()
            except:
                print('Could not get a response from NCBI, trying again')
                c += 1

    def get_taxa_ncbi(self,organism_name):
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={organism_name}'
        taxa_id = self.get_taxa_ncbi_url(url)
        if not taxa_id:
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=candidatus+{organism_name}'
            taxa_id = self.get_taxa_ncbi_url(url)
        if not taxa_id:    print(f'Could not find taxa ID for {organism_name}')
        return taxa_id

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
                        organism_lineage = self.get_organism_lineage(ncbi_taxon_id)
                        res.update(organism_lineage)
                    except:
                        ncbi_taxon_id = self.get_taxa_ncbi(t_id)
                        if ncbi_taxon_id:
                            organism_lineage = self.get_organism_lineage(ncbi_taxon_id)
                            res.update(organism_lineage)
            for i in res:
                self.mantis_nogt_tax.add(str(i))

    def setup_paths_config_file(self):
        file = open(self.config_file, 'r')
        line = file.readline()
        while line:
            line = line.strip('\n')
            if '#' not in line:
                # data sources configuration
                if 'custom_ref_folder=' in line:
                    line_path = add_slash(line.replace('custom_ref_folder=', ''))
                    if line_path: self.mantis_paths['custom'] = line_path

                elif 'ncbi_resources_folder=' in line:
                    line_path = add_slash(line.replace('ncbi_resources_folder=', ''))
                    if line_path: self.mantis_paths['ncbi_res'] = line_path

                elif 'nog_hmm_folder=' in line:
                    line_path = add_slash(line.replace('nog_hmm_folder=', ''))
                    if line_path: self.mantis_paths['NOG'] = line_path

                # taxa ids list for only downloading nogt specific to lineage
                elif 'nog_tax=' in line:
                    line_path = line.replace('nog_tax=', '')
                    self.set_nogt_line(line_path)

                elif 'pfam_hmm_folder=' in line:
                    line_path = add_slash(line.replace('pfam_hmm_folder=', ''))
                    if line_path: self.mantis_paths['pfam'] = line_path

                elif 'kofam_hmm_folder=' in line:
                    line_path = add_slash(line.replace('kofam_hmm_folder=', ''))
                    if line_path: self.mantis_paths['kofam'] = line_path

                elif 'ncbi_hmm_folder=' in line:
                    line_path = add_slash(line.replace('ncbi_hmm_folder=', ''))
                    if line_path: self.mantis_paths['NCBI'] = line_path

                elif 'tigrfam_hmm_folder=' in line[:len('tigrfam_hmm_folder=')]:
                    line_path = add_slash(line.replace('tigrfam_hmm_folder=', ''))
                    if line_path: self.mantis_paths['tigrfam'] = line_path

                elif 'tcdb_seq_folder=' in line[:len('tcdb_seq_folder=')]:
                    line_path = add_slash(line.replace('tcdb_seq_folder=', ''))
                    if line_path: self.mantis_paths['tcdb'] = line_path


                elif '_weight=' in line:
                    ref_source, weight = line.split('_weight=')
                    self.mantis_ref_weights[ref_source] = float(weight)
            line = file.readline()
        file.close()

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
                             'ncbi_res': add_slash(resources_path + 'NCBI'),
                             'custom': add_slash(default_ref_path + 'Custom_references'),
                             'NOG': add_slash(default_ref_path + 'NOG'),
                             'pfam': add_slash(default_ref_path + 'pfam'),
                             'kofam': add_slash(default_ref_path + 'kofam'),
                             'NCBI': add_slash(default_ref_path + 'NCBI'),
                             'tigrfam': add_slash(default_ref_path + 'tigrfam'),
                             'tcdb': add_slash(default_ref_path + 'tcdb'),
                             }
        self.setup_paths_config_file()
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
            get_ref_in_folder(self.mantis_paths['tigrfam']) if not folder else self.mantis_paths['tigrfam'],
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

    def merge_hmm_folder(self, target_folder):
        self.output_folder = target_folder
        print_cyan(f'Merging hmm folder:\n{target_folder}', flush=True, file=self.redirect_verbose)
        output_file = get_path_level(target_folder)
        print(f'Merging hmm folder: {target_folder}', flush=True, file=self.redirect_verbose)
        run_command(f'[ -f {target_folder}{output_file}_merged.hmm ] && rm {target_folder}{output_file}_merged.hmm*',stdout_file=self.redirect_verbose)
        run_command('for i in ' + target_folder + '*.hmm; do cat ${i} >> ' + target_folder + output_file + '_merged.hmm; done',stdout_file=self.redirect_verbose,shell=True,join_command=True)
        run_command(f'hmmpress {target_folder}{output_file}_merged.hmm', stdout_file=self.redirect_verbose)

    def get_path_default_ref(self, database, taxon_id=None):
        target_file = None
        if 'kofam' in database.lower():
            target_file = get_ref_in_folder(self.mantis_paths['kofam'])
        elif 'pfam' in database.lower():
            target_file = get_ref_in_folder(self.mantis_paths['pfam'])
        elif 'tigrfam' in database.lower():
            target_file = get_ref_in_folder(self.mantis_paths['tigrfam'])
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
        if database == 'ncbi_res':
            if file_exists(self.mantis_paths['ncbi_res'] + 'taxidlineage.dmp', force_download) and \
                    file_exists(self.mantis_paths['ncbi_res'] + 'gc.prt', force_download):
                return True
        elif database == 'NOGSQL':
            if file_exists(self.mantis_paths['default'] + 'eggnog.db', force_download):
                return True
        elif database == 'tcdb':
            if file_exists(self.mantis_paths['tcdb'] + 'tcdb.dmnd', force_download):
                return True
        target_file = self.get_path_default_ref(database, taxon_id)
        if target_file:
            for extension in ['', '.h3f', '.h3i', '.h3m', '.h3p']:
                if not file_exists(target_file + extension, force_download=force_download):
                    return False
        else:
            return False
        return True

    #####LISTING HMMS DATABASE#####


    def check_installation_extras(self, res, verbose=True):
        if verbose: yellow('Checking extra files', flush=True, file=self.redirect_verbose)
        if not file_exists(self.mantis_paths['resources'] + 'essential_genes/essential_genes.txt'):
            red('Essential genes list is missing, it should be in the github repo!')
            if verbose: red('Failed installation check on [files missing]: ' + self.mantis_paths[
                'resources'] + 'essential_genes/essential_genes.txt', flush=True, file=self.redirect_verbose)
            res.append(self.mantis_paths['resources'] + 'essential_genes/')
        else:
            if verbose: green('Passed installation check on: ' + self.mantis_paths['resources'] + 'essential_genes',
                              flush=True, file=self.redirect_verbose)
        if self.mantis_paths['ncbi_res'][0:2] != 'NA':
            if not file_exists(self.mantis_paths['ncbi_res'] + 'taxidlineage.dmp'):
                if verbose: red(
                    'Failed installation check on [files missing]: ' + self.mantis_paths['ncbi_res'] + 'taxidlineage.dmp',flush=True, file=self.redirect_verbose)
                res.append(self.mantis_paths['ncbi_res'])
            else:
                if verbose: green('Passed installation check on: ' + self.mantis_paths['ncbi_res']+ 'taxidlineage.dmp', flush=True,file=self.redirect_verbose)
            if not file_exists(self.mantis_paths['ncbi_res'] + 'gc.prt'):
                if verbose: red(
                    'Failed installation check on [files missing]: ' + self.mantis_paths['ncbi_res'] + 'gc.prt.dmp',flush=True, file=self.redirect_verbose)
                res.append(self.mantis_paths['ncbi_res'])
            else:
                if verbose: green('Passed installation check on: ' + self.mantis_paths['ncbi_res']+'gc.prt.dmp', flush=True,file=self.redirect_verbose)
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

    def check_installation(self, verbose=True):
        # we use the verbose mode when running the check_installation directly
        self.passed_check = True
        if not cython_compiled():
            self.passed_check = False
            if verbose: red('Cython needs to be compiled!', flush=True,
                            file=self.redirect_verbose)
        else:
            if verbose: green('Cython correctly compiled!', flush=True, file=self.redirect_verbose)
        res = []
        res = self.check_installation_extras(res, verbose)

        if verbose: yellow('Checking HMM installation', flush=True, file=self.redirect_verbose)
        requirements = {
            self.mantis_paths['pfam']: ['metadata.tsv'],
            self.mantis_paths['tcdb']: ['metadata.tsv'],
            self.mantis_paths['kofam']: ['metadata.tsv'],
            self.mantis_paths['tigrfam']: ['metadata.tsv'],
        }
        # per tax level FOR EGGNOG
        if self.mantis_paths['NOG'][0:2] != 'NA':
            tax_hmms = self.get_taxon_hmms(db='NOG', folder=True)
            if not tax_hmms:
                if verbose: red('Failed installation check on [path unavailable]: ' + self.mantis_paths['NOG'],
                                flush=True, file=self.redirect_verbose)
                res.append(self.mantis_paths['NOG'])
            for tax_hmm_folder in tax_hmms:
                # we skip the taxon 1 since it has no hmms
                tax_hmm = tax_hmm_folder.split(SPLITTER)[-2]
                if tax_hmm != '1':
                    self.check_installation_folder(tax_hmm_folder, res, verbose=False,
                                                   extra_requirements=[f'{tax_hmm}_sql_annotations.tsv'])
            nogt_check = [i for i in res if self.mantis_paths['NOG'] in i]
            if not nogt_check:
                if verbose: green('Passed installation check on: ' + self.mantis_paths['NOG'], flush=True,
                                  file=self.redirect_verbose)
        # per tax level FOR NCBI
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            #checking those already present
            tax_hmms = self.get_taxon_hmms(db='NCBI', folder=True)
            if not tax_hmms:
                if verbose: red('Failed installation check on [path unavailable]: ' + self.mantis_paths['NCBI'],
                                flush=True, file=self.redirect_verbose)
                res.append(self.mantis_paths['NCBI'])
            for tax_hmm_folder in tax_hmms:
                # we skip the taxon 1 since it has no hmms
                self.check_installation_folder(tax_hmm_folder, res, verbose=False, extra_requirements=['metadata.tsv'])
            ncbi_check = [i for i in res if self.mantis_paths['NCBI'] in i]
            if not ncbi_check:
                if verbose: green('Passed installation check on: ' + self.mantis_paths['NCBI'], flush=True,
                                  file=self.redirect_verbose)
        for ref_folder in self.compile_refs_list(folder=True):
            if ref_folder in requirements:
                self.check_installation_folder(ref_folder, res, verbose, extra_requirements=requirements[ref_folder])
            else:
                self.check_installation_folder(ref_folder, res, verbose)
        if res:
            self.passed_check = False
            fail_res = ''
            for i in res: fail_res += f'{i}\n'
            if verbose: red(f'Installation check failed on:\n{fail_res}', flush=True, file=self.redirect_verbose)
        if self.passed_check:
            if verbose:
                yellow('------------------------------------------', flush=True, file=self.redirect_verbose)
                green('--------INSTALLATION CHECK PASSED!--------', flush=True, file=self.redirect_verbose)
                yellow('------------------------------------------', flush=True, file=self.redirect_verbose)
            else:
                print_cyan('Installation check passed', flush=True, file=self.redirect_verbose)

        else:
            if verbose:
                yellow('------------------------------------------', flush=True, file=self.redirect_verbose)
                red('--------INSTALLATION CHECK FAILED!--------', flush=True, file=self.redirect_verbose)
                yellow('------------------------------------------', flush=True, file=self.redirect_verbose)
            else:
                print_cyan('Installation check failed', flush=True, file=self.redirect_verbose)

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

    def get_hmm_taxon_ids(self, db):
        res=[]
        if not os.path.exists(self.mantis_paths[db]): return res
        if db=='NOG':
            available_taxon_ids = self.get_taxon_ids_NOGT()
            if self.mantis_nogt_tax:
                for tax_id in self.mantis_nogt_tax:
                    res.append(tax_id)
            else:
                res.extend(available_taxon_ids)
            res = set(res)
            if db == 'NOG':
                res = res.intersection(available_taxon_ids)
        else:
            for i in os.listdir(self.mantis_paths[db]):
                if re.search('\d+',i): res.append(i)
        if not os.path.exists(self.mantis_paths[db]): return []
        res = set(res)
        return list(res)

    def get_taxon_hmms(self, db, folder=False):
        taxon_ids = self.get_hmm_taxon_ids(db)
        res = []
        for t in taxon_ids:
            if folder:
                res.append(add_slash(self.mantis_paths[db] + t))
            else:
                res.append(add_slash(self.mantis_paths[db] + t) + f'{t}_merged.hmm')
        if folder:
            res.append(add_slash(self.mantis_paths[db] + db + 'G'))
        else:
            res.append(add_slash(self.mantis_paths[db] + db + 'G') + f'{db}G_merged.hmm')
        return res

    def get_lineage_hmm_path(self, taxon_id, db):
        tax_hmms = self.get_hmm_taxon_ids(db=db)
        if taxon_id in tax_hmms:
            return add_slash(self.mantis_paths[db] + taxon_id) + f'{taxon_id}_merged.hmm'
        else:
            return None

    def get_organism_lineage(self, taxon_id, stdout_file=None):
        lineage_file_path = self.mantis_paths['ncbi_res'] + 'taxidlineage.dmp'
        try:
            lineage_file = open(lineage_file_path, 'r')
        except:
            print_cyan('Lineage dump is not present! If you\'d like to run taxonomic lineage annotation, please run < setup_databases >',
                flush=True, file=stdout_file)
            return []
        line = lineage_file.readline().strip('\n').replace('|', '')
        while line:
            line = line.split()
            if str(taxon_id) == str(line[0]):
                lineage = line[1:]
                lineage.append(taxon_id)
                lineage_file.close()
                return lineage
            line = lineage_file.readline().strip('\n').replace('|', '')
        lineage_file.close()
        return []

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
    p = MANTIS_Assembler()
    print(p.mantis_paths)
    #p.mantis_paths['tcdb']='/home/pedroq/Desktop/test_mantis/tcdb/'
    #test=p.compile_tcdb_metadata()
    p.get_essential_genes_pfam(['P05652','P05649'])
