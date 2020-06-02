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



def setup_databases(force_download=False,chunk_size=None):
    print_cyan('Setting up databases')
    if force_download=='None': force_download=None
    if chunk_size: chunk_size=int(chunk_size)
    mantis = MANTIS_Assembler(hmm_chunk_size=chunk_size)
    mantis.setup_databases(force_download)


def merge_hmm_folder(target_folder):
    print_cyan('Merging HMM folder')
    s_target_folder=add_slash(target_folder)
    if target_folder=='None': return
    mantis = MANTIS_Assembler()
    mantis.merge_hmm_folder(target_folder)


def check_installation(mantis_config=None):
    yellow('Checking installation')
    mantis = MANTIS_Assembler(mantis_config=mantis_config)
    mantis.check_installation()


#hmmer threads:
#http://www.hicomb.org/papers/HICOMB2018-04.pdf

class MANTIS_Assembler(MANTIS_DB):
    def __init__(self,hmmer_threads=None,verbose=True,redirect_verbose=None,job_tracker=None,limit_jobs=None,mantis_config=None,hmm_chunk_size=None):
        self.target_hmm=None
        if hmmer_threads: self.hmmer_threads=hmmer_threads
        else: self.hmmer_threads=5
        self.redirect_verbose=redirect_verbose
        self.verbose=verbose
        self.mantis_config=mantis_config
        self.broken_merged_hmms=set()
        self.clean_merged_hmms=set()
        self.start_time=time()
        #to speed up hmm search we split very big hmms into smaller chunks - better job distribution
        if hmm_chunk_size==0:           self.hmm_chunk_size=None
        elif hmm_chunk_size is None:    self.hmm_chunk_size=5000
        else:                           self.hmm_chunk_size=hmm_chunk_size
        self.read_config_file()
        #self.requirements_met()
        #I use manager instead of queue since I need to be able to add records to the end and start of the 'queue' (actually a list) which is not possible with the multiprocessing.Queue
        self.manager=Manager()
        self.queue = self.manager.list()




    def __str__(self):
        custom_hmms=self.get_custom_hmms_paths(folder=True)
        res=''
        custom_hmms_str=''
        for chmm in custom_hmms:
            custom_hmms_str+=chmm+'\n'
        if custom_hmms_str:
            res = 'Custom hmms folders:\n'+custom_hmms_str
        return 'Output folder:\n' +self.output_folder if hasattr(self,'output_folder') else ''+\
       '#  External data folders:'+'\n'+\
        '------------------------------------------'+'\n'+\
        'Default hmms folder:\n'+self.mantis_paths['default']+'\n'+\
        'Custom hmms folder:\n'+self.mantis_paths['custom']+'\n'+\
        'NCBI dump folder:\n'+self.mantis_paths['ncbi']+'\n'+\
        'Gene ontology folder:\n'+self.mantis_paths['go_obo_nlp']+'\n'+\
        'Uniprot dump folder:\n'+self.mantis_paths['uniprot_nlp']+'\n'+\
        'TAX NOG hmms folder:\n'+self.mantis_paths['NOGT']+'\n'+\
        'Global NOG hmms folder:\n'+self.mantis_paths['NOGG']+'\n'+\
        'Pfam hmms folder:\n'+self.mantis_paths['pfam']+'\n'+\
        'KOfam hmms folder:\n'+self.mantis_paths['kofam']+'\n'+\
        'dbCAN hmms folder:\n'+self.mantis_paths['dbcan']+'\n'+\
        'TIGRFAM hmms folder:\n'+self.mantis_paths['tigrfam']+'\n'+\
        'Resfams hmms folder:\n'+self.mantis_paths['resfams']+'\n'+res+\
        '------------------------------------------'


    def print_citation(self):

        return

    def requirements_met(self):
        for f in [self.is_conda_available(), self.is_hmmer_available()]:
            if not f: raise RequirementsNotMet

    def is_conda_available(self):
        process= self.run_command('conda -V',get_output=True)
        check = re.search('conda', str(process.stdout))
        if not check:
            print_cyan('Conda dependency not met!')
            print_cyan('Install Conda on your system by following the instructions at: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html')
        return check

    def is_hmmer_available(self):
        check_command=' conda list hmmer '
        process=self.run_command(check_command ,get_output=True)
        check = re.search('hmmer', str(process.stdout))
        if not check:
            print_cyan('HMMER dependency not met!')
            print_cyan('Install HMMER on your conda environment by doing the following:')
            print_cyan('conda activate <conda_environment>')
            print_cyan('conda install -c bioconda hmmer')
        return check


    def read_config_file(self):
        #HMMS
        # if there's no path, we just assume its in the default folder
        default_hmm_path=add_slash(mantis_folder+'hmm')
        resources_path=add_slash(mantis_folder + 'Resources')
        self.mantis_paths={'default':default_hmm_path,
                           'resources':resources_path,
                           'go_obo_nlp':add_slash(resources_path + 'Gene_Ontology'),
                           'uniprot_nlp':add_slash(resources_path + 'Uniprot'),
                           'ncbi': add_slash(resources_path + 'NCBI'),
                           'custom':add_slash(default_hmm_path + 'custom_hmms'),
                           'NOGT':add_slash(default_hmm_path + 'NOGT'),
                           'NOGG':add_slash(default_hmm_path+'NOGG'),
                           'pfam':add_slash(default_hmm_path + 'pfam'),
                           'kofam':add_slash(default_hmm_path + 'kofam'),
                           'dbcan':add_slash(default_hmm_path + 'dbcan'),
                           'tigrfam':add_slash(default_hmm_path + 'tigrfam'),
                           'resfams':add_slash(default_hmm_path + 'resfams'),
                           }
        self.mantis_hmm_weights={'else':0.7}
        if self.mantis_config:
            print('Using custom MANTIS.config',flush=True,file=self.redirect_verbose)
            self.config_file=self.mantis_config
        else:
            if not os.path.isdir(mantis_folder):
                print('Make sure you are calling the folder to run this package, like so:\n python mantis/ <command>\n ',flush=True,file=self.redirect_verbose)
                self.cancel_all_jobs()
                raise FileNotFoundError
            self.config_file = mantis_folder+'MANTIS.config'
        try:
            file = open(self.config_file,'r')
        except:
            print('MANTIS.config file has been deleted or moved, make sure you keep it in the root of the project!',flush=True,file=self.redirect_verbose)
            raise FileNotFoundError
        line = file.readline()
        while line:
            line = line.strip('\n')
            if '#' not in line:
                #data sources configuration
                if 'default_hmms_folder=' in line:
                    self.mantis_paths['default'] = add_slash(line.replace('default_hmms_folder=', ''))
                elif 'custom_hmms_folder=' in line:
                    self.mantis_paths['custom'] = add_slash(line.replace('custom_hmms_folder=', ''))
                elif 'ncbi_dmp_path_folder=' in line:
                    self.mantis_paths['ncbi']= add_slash(line.replace('ncbi_dmp_path_folder=', ''))
                elif 'go_obo_nlp_folder=' in line:
                    self.mantis_paths['go_obo_nlp']= add_slash(line.replace('go_obo_nlp_folder=', ''))
                elif 'uniprot_nlp_folder=' in line:
                    self.mantis_paths['uniprot_nlp']= add_slash(line.replace('uniprot_nlp_folder=', ''))
                elif 'NOGT_hmm_folder=' in line:
                    self.mantis_paths['NOGT'] = add_slash(line.replace('NOGT_hmm_folder=', ''))
                elif 'NOGG_hmm_folder=' in line:
                    self.mantis_paths['NOGG'] = add_slash(line.replace('NOGG_hmm_folder=', ''))
                elif 'pfam_hmm_folder=' in line:
                    self.mantis_paths['pfam'] = add_slash(line.replace('pfam_hmm_folder=', ''))
                elif 'kofam_hmm_folder=' in line:
                    self.mantis_paths['kofam'] = add_slash(line.replace('kofam_hmm_folder=', ''))
                elif 'dbcan_hmm_folder=' in line:
                    self.mantis_paths['dbcan'] = add_slash(line.replace('dbcan_hmm_folder=', ''))
                elif 'tigrfam_hmm_folder=' in line[:len('tigrfam_hmm_folder=')]:
                    self.mantis_paths['tigrfam'] = add_slash(line.replace('tigrfam_hmm_folder=', ''))
                elif 'resfams_hmm_folder=' in line:
                    self.mantis_paths['resfams'] = add_slash(line.replace('resfams_hmm_folder=', ''))
                elif '_weight=' in line:
                    hmm_source,weight=line.split('_weight=')
                    self.mantis_hmm_weights[hmm_source]=float(weight)
            line = file.readline()
        file.close()
        if self.verbose: print(self,flush=True,file=self.redirect_verbose)

    def set_path_go_terms_nlp(self):
        self.go_terms_path = self.mantis_paths['go_obo_nlp']+'go.obo'

    def set_path_uniprot_proteins_nlp(self):
        uniprot_folder = self.mantis_paths['uniprot_nlp']
        files_dir=os.listdir(uniprot_folder)
        passed_check=False
        for file in files_dir:
            if 'pickle' not in file and '.gz' not in file:
                self.uniprot_reference=add_slash(uniprot_folder)+file
                passed_check=True
        if not passed_check:
            for file in files_dir:
                if '.gz' in file:
                    print('Uncompressing uniprot reference!')
                    uncompress_archive(add_slash(uniprot_folder)+file)
                    self.uniprot_reference = add_slash(uniprot_folder) + file.replace('.gz','')
                    passed_check=True
        return passed_check


    def order_by_size_descending(self,hmms_list):
        res={}
        for hmm in hmms_list:
            res[hmm]=os.stat(hmm).st_size
        return sorted(res, key=res.get,reverse=True),res


    def compile_hmms_list(self,folder=False):
        #doesnt include NOG
        hmms_list=[]
        default_list = [
            get_hmm_in_folder(self.mantis_paths['pfam']) if not folder else self.mantis_paths['pfam'],
            get_hmm_in_folder(self.mantis_paths['kofam']) if not folder else self.mantis_paths['kofam'],
            get_hmm_in_folder(self.mantis_paths['dbcan']) if not folder else self.mantis_paths['dbcan'],
            get_hmm_in_folder(self.mantis_paths['tigrfam']) if not folder else self.mantis_paths['tigrfam'],
            get_hmm_in_folder(self.mantis_paths['resfams']) if not folder else self.mantis_paths['resfams']]
        for hmm_path in self.get_custom_hmms_paths(folder):
            if hmm_path[0:2] !='NA':
                hmms_list.append(hmm_path)
        for hmm_path in default_list:
            if hmm_path and hmm_path[0:2]!='NA':
                hmms_list.append(hmm_path)
        return hmms_list




    # for executing simple commands in the hpc we run this
    def run_command(self, command, get_output=False,stdout_file=None):
        if get_output:
            process = subprocess.run(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        elif stdout_file:
            process = subprocess.run(command, shell=True, stdout=stdout_file,stderr=stdout_file)
        else:
            process = subprocess.run(command, shell=True)

        return process


    #####SETTING UP DATABASE#####

    def merge_hmm_folder(self, target_folder):
        self.output_folder= target_folder
        print_cyan('Merging hmm folder:\n' + target_folder,flush=True,file=self.redirect_verbose)
        output_file = get_path_level(target_folder)
        print('Merging hmm folder: '+target_folder,flush=True,file=self.redirect_verbose)
        self.run_command('[ -f ' + target_folder + output_file + '_merged.hmm' + ' ] && rm ' + target_folder + output_file + '_merged.hmm*',stdout_file=self.redirect_verbose)
        self.run_command('for i in ' + target_folder + '*.hmm; do cat ${i} >> ' + target_folder + output_file + '_merged.hmm; done',stdout_file=self.redirect_verbose)
        self.run_command('hmmpress '+target_folder + output_file + '_merged.hmm',stdout_file=self.redirect_verbose)

    def file_exists(self,target_file,force_download=False):
        if os.path.exists(target_file) or force_download:
            return True
        return False

    def get_path_default_hmm(self,database,taxon_id):
        target_file=None
        if 'kofam' in database.lower():
            target_file = get_hmm_in_folder(self.mantis_paths['kofam'])
        elif 'dbcan'  in database.lower():
            target_file = get_hmm_in_folder(self.mantis_paths['dbcan'])
        elif 'pfam'  in database.lower():
            target_file = get_hmm_in_folder(self.mantis_paths['pfam'])
        elif 'tigrfam'  in database.lower():
            target_file = get_hmm_in_folder(self.mantis_paths['tigrfam'])
        elif 'resfams'  in database.lower():
            target_file = get_hmm_in_folder(self.mantis_paths['resfams'])
        elif 'NOGG'.lower()  in database.lower():
            target_file = get_hmm_in_folder(self.mantis_paths['NOGG'])
        elif 'NOGT'.lower()  in database.lower():
            target_file = get_hmm_in_folder(add_slash(self.mantis_paths['NOGT'] + taxon_id))
        return target_file



    def check_hmm_exists(self,database,taxon_id=None,force_download=False):
        target_file=self.get_path_default_hmm(database,taxon_id)
        if target_file:
            for extension in ['','.h3f','.h3i','.h3m','.h3p']:
                if not self.file_exists(target_file+extension,force_download=force_download):
                    return False
        return True

    #####LISTING HMMS DATABASE#####

    def check_installation_extras(self,res,verbose=True):
        if verbose: yellow('Checking extra files',flush=True,file=self.redirect_verbose)
        if not self.file_exists(self.mantis_paths['resources'] + 'essential_genes/essential_genes.txt'):
            red('Essential genes list is missing, it should be in the github repo!')
            if verbose: red('Failed installation check on [files missing]: ' + self.mantis_paths['resources'] + 'essential_genes/essential_genes.txt', flush=True,file=self.redirect_verbose)
            res.append(self.mantis_paths['resources'] + 'essential_genes/')
        else:
            if verbose: green('Passed installation check on: ' + self.mantis_paths['resources'] + 'essential_genes', flush=True,file=self.redirect_verbose)


        if not self.file_exists(self.mantis_paths['ncbi'] + 'taxidlineage.dmp'):
            if verbose: red('Failed installation check on [files missing]: ' + self.mantis_paths['ncbi'], flush=True, file=self.redirect_verbose)
            res.append(self.mantis_paths['ncbi'])
        else:
            if verbose: green('Passed installation check on: ' + self.mantis_paths['ncbi'], flush=True, file=self.redirect_verbose)

        if not os.listdir(self.mantis_paths['uniprot_nlp']) or not self.set_path_uniprot_proteins_nlp():
            if verbose: red('Failed installation check on [files missing]: ' + self.mantis_paths['uniprot_nlp'], flush=True, file=self.redirect_verbose)
            res.append(self.mantis_paths['uniprot_nlp'])
        else:
            if verbose: green('Passed installation check on: ' + self.mantis_paths['uniprot_nlp'], flush=True, file=self.redirect_verbose)

        if not self.file_exists(self.mantis_paths['go_obo_nlp'] + 'go.obo'):
            if verbose: red('Failed installation check on [files missing]: ' + self.mantis_paths['go_obo_nlp'], flush=True, file=self.redirect_verbose)
            res.append(self.mantis_paths['go_obo_nlp'])
        else:
            if verbose: green('Passed installation check on: ' + self.mantis_paths['go_obo_nlp'], flush=True, file=self.redirect_verbose)


        return res

    def check_installation_folder(self,hmm_folder_path,res,verbose=True,extra_requirements=[]):
        check = 5+len(extra_requirements)
        try:
            files_dir = os.listdir(hmm_folder_path)
        except:
            if verbose: red('Failed installation check on [path unavailable]: ' + hmm_folder_path,flush=True,file=self.redirect_verbose)
            res.append(hmm_folder_path)
            self.passed_check = False
            return
        for file in files_dir:
            if '.hmm' == file[-4:]:
                check -= 1
            elif '.h3f' == file[-4:]:
                check -= 1
            elif '.h3i' == file[-4:]:
                check -= 1
            elif '.h3m' == file[-4:]:
                check -= 1
            elif '.h3p' == file[-4:]:
                check -= 1
            elif file in extra_requirements:
                check-=1
        if check != 0:
            red('Failed installation check on [files missing]: ' + hmm_folder_path,flush=True,file=self.redirect_verbose)
            res.append(hmm_folder_path)
        else:
            if verbose: green('Passed installation check on: '+hmm_folder_path,flush=True,file=self.redirect_verbose)

    def check_installation(self,verbose=True):
        #we use the verbose mode when running the check_installation directly
        self.passed_check=True
        res=[]
        res=self.check_installation_extras(res,verbose)

        if verbose: yellow('Checking HMM installation',flush=True,file=self.redirect_verbose)
        requirements={
            self.mantis_paths['NOGG']:['NOGG_sql_annotations.tsv'],
            self.mantis_paths['pfam']:['Pfam-A.hmm.dat'],
            self.mantis_paths['kofam']:['ko_list','ko2cog.xl','ko2go.xl','ko2tc.xl','ko2cazy.xl','ko_to_path','map_description'],
            self.mantis_paths['dbcan']:['CAZyDB.07312019.fam.subfam.ec.txt'],
            self.mantis_paths['tigrfam']:['gpl.html','COPYRIGHT','TIGRFAMS_GO_LINK','TIGRFAMS_ROLE_LINK','TIGR_ROLE_NAMES'],
            self.mantis_paths['resfams']:['180102_resfams_metadata_updated_v122.tsv'],

        }
        if self.target_hmm:
            target_hmm_folder=get_folder(self.target_hmm)
            target_hmm_folder=splitter.join(target_hmm_folder)+splitter
            if self.target_hmm in requirements:
                self.check_installation_folder(target_hmm_folder, res, verbose, extra_requirements=requirements[target_hmm_folder])
            else:
                self.check_installation_folder(target_hmm_folder, res, verbose)
        else:
            #per tax level
            if self.mantis_paths['NOGT'][0:2] !='NA':
                nog_tax=self.get_NOG_taxon_hmms(folder=True)
                if not nog_tax:
                    if verbose: red('Failed installation check on [path unavailable]: ' + self.mantis_paths['NOGT'],flush=True,file=self.redirect_verbose)
                    res.append(self.mantis_paths['NOGT'])
                for tax_hmm_folder in nog_tax:
                    #we skip the taxon 1 since it has no hmms
                    tax_hmm=tax_hmm_folder.split(splitter)[-2]
                    if tax_hmm!='1':
                        self.check_installation_folder(tax_hmm_folder,res,verbose=False,extra_requirements=[tax_hmm+'_sql_annotations.tsv'])
                nogt_check=[i for i in res if self.mantis_paths['NOGT'] in i]
                if not nogt_check:
                    if verbose: green('Passed installation check on: ' + self.mantis_paths['NOGT'], flush=True,file=self.redirect_verbose)
            if self.mantis_paths['NOGG'][0:2] !='NA':
                self.check_installation_folder(self.mantis_paths['NOGG'],res,verbose,extra_requirements=requirements[self.mantis_paths['NOGG']])
            for hmm_folder in self.compile_hmms_list(folder=True):
                if hmm_folder in requirements:
                    self.check_installation_folder(hmm_folder, res, verbose,extra_requirements=requirements[hmm_folder])
                else:
                    self.check_installation_folder(hmm_folder, res, verbose)
        if res:
            self.passed_check = False
            if verbose: red('Installation check failed on: '+str(len(res)),flush=True,file=self.redirect_verbose)
        if self.passed_check:
            if verbose:
                yellow('------------------------------------------',flush=True,file=self.redirect_verbose)
                green('--------INSTALLATION CHECK PASSED!--------',flush=True,file=self.redirect_verbose)
                yellow('------------------------------------------',flush=True,file=self.redirect_verbose)
            else: print_cyan('Installation check passed',flush=True,file=self.redirect_verbose)

        else:
            if verbose:
                yellow('------------------------------------------',flush=True,file=self.redirect_verbose)
                red('--------INSTALLATION CHECK FAILED!--------',flush=True,file=self.redirect_verbose)
                yellow('------------------------------------------',flush=True,file=self.redirect_verbose)
            else:
                print_cyan('Installation check failed',flush=True,file=self.redirect_verbose)




    def get_custom_hmms_paths(self,folder=False):
        try:
            custom_hmms_folders = os.listdir(self.mantis_paths['custom'])
            for potential_hmm_folder in custom_hmms_folders:
                try:
                    files = os.listdir(self.mantis_paths['custom'] + potential_hmm_folder)
                    for hmm in files:
                        if '.hmm' in hmm[-4:]:
                            if folder:
                                try:  yield add_slash(self.mantis_paths['custom'] + potential_hmm_folder)
                                except GeneratorExit: return ''
                            else:
                                try: yield add_slash(self.mantis_paths['custom'] + potential_hmm_folder) + splitter + hmm
                                except GeneratorExit: return ''
                except: pass
        except:
            print('Custom hmms folder is missing, did you correctly set the path? If path is not set make sure you didn\'t delete the custom_hmms folder!',flush=True,file=self.redirect_verbose)
            self.passed_check=False
            return
        with open(self.config_file,'r') as file:
            line = file.readline()
            while line:
                if line[0] != '#':
                    if 'custom_hmm=' in line:
                        line = line.strip('\n')
                        if folder:
                            try: yield add_slash(splitter.join(line.replace('custom_hmm=', '').split(splitter)[:-1]))
                            except GeneratorExit:    return ''

                        else:
                            try: yield line.replace('custom_hmm=', '')
                            except GeneratorExit:    return ''
                line = file.readline()


    def get_NOG_taxon_ids(self):
        if not os.path.exists(self.mantis_paths['NOGT']): return []
        NOG_tax_hmms = os.listdir(self.mantis_paths['NOGT'])
        return [i for i in NOG_tax_hmms if re.search('\d+', i)]

    def get_NOG_taxon_hmms(self, folder=False):
        NOG_taxon_ids = self.get_NOG_taxon_ids()
        res = []
        for t in NOG_taxon_ids:
            if folder:
                res.append(self.mantis_paths['NOGT']+t+splitter)
            else:
                res.append(self.mantis_paths['NOGT']+t+splitter+t+'_merged.hmm')
        return res


    def get_lineage_hmm_path(self, taxon_id):
        NOG_tax_hmms = self.get_NOG_taxon_ids()
        if taxon_id in NOG_tax_hmms:
            return self.mantis_paths['NOGT'] + taxon_id + splitter + taxon_id+'_merged.hmm'
        else: return None

    def processes_handler(self,target_worker_function,worker_count,add_sentinels=True):
        processes = [Process(target=target_worker_function, args=(self.queue,)) for _ in range(worker_count)]
        #adding sentinel record since queue can be signaled as empty when its really not
        if add_sentinels:
            for _ in range(worker_count):   self.queue.append(None)
        for process in processes:   process.start()
        for process in processes:   process.join()



if __name__ == '__main__':
    p=MANTIS_Assembler(mantis_config='/home/pedroq/Desktop/test_hmm/t.config')
    hmm_path='/home/pedroq/Desktop/test_hmm/dbcan/dbcan.hmm'
    #a=p.split_hmm_into_chunks(hmm_path)
    #p.setup_databases()