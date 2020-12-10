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
    def setup_databases(self,force_download=False):
        if not cython_compiled():
            compile_cython()
        self.output_folder= mantis_folder+'setup_databases/'
        self.mantis_out=self.output_folder + 'Mantis.out'
        if file_exists(self.mantis_out):
            os.remove(self.mantis_out)
        Path(self.mantis_paths['default']).mkdir(parents=True, exist_ok=True)
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        dbs_list=[
            'go_obo_nlp',
            'uniprot_nlp',
            'ncbi'      if self.mantis_paths['ncbi'][0:2] !='NA' else None,
            'pfam'      if self.mantis_paths['pfam'][0:2] !='NA' else None,
            'kofam'     if self.mantis_paths['kofam'][0:2] !='NA' else None,
            #'dbcan'     if self.mantis_paths['dbcan'][0:2] !='NA' else None,
            'tigrfam'   if self.mantis_paths['tigrfam'][0:2] !='NA' else None,
            #'resfams'   if self.mantis_paths['resfams'][0:2] !='NA' else None,
            #'hamap'     if self.mantis_paths['hamap'][0:2] !='NA' else None,
            'NOGG'      if self.mantis_paths['NOGG'][0:2] !='NA' else None,
        ]
        self.prepare_queue_setup_databases_NOG_sql(force_download)
        self.prepare_queue_setup_databases(dbs_list,force_download)
        #for unzipping tax specific hmms
        self.prepare_queue_setup_databases_tax(force_download)

        #self.queue.append(['NOGT', force_download, '10', self.output_folder + 'Mantis.out'])


        worker_count = estimate_number_workers_setup_database(len(self.queue))
        print('Database will be setup with '+str(worker_count)+' workers!', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_setup_databases, worker_count)
        print_cyan('Finished downloading all data!', flush=True, file=self.redirect_verbose)
        if self.hmm_chunk_size:
            print_cyan('Will now split data into chunks!', flush=True, file=self.redirect_verbose)
            self.prepare_queue_split_hmms()
            worker_count = estimate_number_workers_setup_database(len(self.queue))
            print('Database will be split with ' + str(worker_count) + ' workers!', flush=True,file=self.redirect_verbose)
            self.processes_handler(self.worker_split_hmms, worker_count)
        print_cyan('Will now extract metadata from NOG!', flush=True, file=self.redirect_verbose)
        self.prepare_queue_extract_metadata()
        worker_count = estimate_number_workers_setup_database(len(self.queue),minimum_jobs_per_worker=2)
        print('Metadata will be extracted with ' + str(worker_count) + ' workers!', flush=True, file=self.redirect_verbose)
        self.processes_handler(self.worker_extract_NOG_metadata, worker_count)
        #if os.path.exists(self.mantis_paths['default']+'eggnog.db'):os.remove(self.mantis_paths['default']+"eggnog.db")
        print('Preparing NLP Resources!', flush=True, file=self.redirect_verbose)
        MANTIS_NLP.__init__(self)
        print_cyan('Finished setting up databases!',flush=True,file=self.redirect_verbose)

    ###Filling queue with jobs

    def prepare_queue_setup_databases(self,dbs_list,force_download):
        for database in dbs_list:
            if database:
                self.queue.append([database,force_download,self.mantis_out])

    def prepare_queue_setup_databases_tax(self,force_download):
        stdout_file=open(self.mantis_out,'a+')
        if self.mantis_paths['NOGT'][0:2] !='NA':
            list_taxon_ids = self.get_taxon_for_queue_NOGT(stdout_file=stdout_file)
            for taxon_id in list_taxon_ids:
                self.queue.append(['NOGT',force_download,taxon_id,self.mantis_out])

    def prepare_queue_setup_databases_NOG_sql(self,force_download):
        if self.mantis_paths['NOGG'][0:2] !='NA' or self.mantis_paths['NOGT'][0:2] !='NA':
            self.queue.append(['NOGSQL',force_download,self.mantis_out])

    def prepare_queue_split_hmms(self):
        print('Checking which hmms need to be split, this may take a while...',flush=True,file=self.redirect_verbose)
        taxon_hmms = self.get_taxon_for_queue_NOG_split_hmms()
        general_hmms = self.get_hmm_for_queue_split_hmms()
        print('Will split: ',[get_path_level(i) for i in taxon_hmms+general_hmms],flush=True,file=self.redirect_verbose)
        for hmm_path in taxon_hmms+general_hmms:
            self.queue.append([hmm_path,self.mantis_out])

    def prepare_queue_extract_metadata(self):
        stdout_file = open(self.mantis_out, 'a+')
        print('Checking which NOGs we need to extract metadata from',flush=True,file=stdout_file)
        if self.mantis_paths['NOGG'][0:2] !='NA' or self.mantis_paths['NOGT'][0:2] !='NA':
            self.create_index_eggnog_db()
            self.unpack_NOG_sql(stdout_file=stdout_file)
        if self.mantis_paths['NOGG'][0:2] !='NA':
            target_sql_file = self.mantis_paths['NOGG'] + 'NOGG_sql_annotations.tsv'
            target_annotation_file = self.mantis_paths['NOGG'] + 'NOG.annotations.tsv'
            taxon_id=None
            if file_exists(target_sql_file):
                if len(self.get_hmms_annotation_file(target_sql_file, hmm_col=0)) != len(self.get_hmms_annotation_file(target_annotation_file, hmm_col=1)):
                    self.queue.append([target_sql_file,target_annotation_file,taxon_id, self.mantis_out])
                else:
                    print('Skipping metadata extraction for NOGG',taxon_id, flush=True, file=stdout_file)
            else:
                self.queue.append([target_sql_file, target_annotation_file, taxon_id, self.mantis_out])

        if self.mantis_paths['NOGT'][0:2] !='NA':
            for taxon_id in os.listdir(self.mantis_paths['NOGT']):
                if os.path.isdir(self.mantis_paths['NOGT']+taxon_id):
                    target_annotation_file = self.mantis_paths['NOGT'] + taxon_id + splitter + taxon_id + '_annotations.tsv'
                    target_sql_file = self.mantis_paths['NOGT'] + taxon_id + splitter + taxon_id + '_sql_annotations.tsv'
                    if file_exists(target_sql_file):
                        if len(self.get_hmms_annotation_file(target_sql_file, hmm_col=0)) != len(self.get_hmms_annotation_file(target_annotation_file, hmm_col=1)):
                            self.queue.append([target_sql_file, target_annotation_file, taxon_id,self.mantis_out])
                        else:
                            print('Skipping metadata extraction for NOGT',taxon_id, flush=True, file=stdout_file)
                    else:
                        self.queue.append([target_sql_file, target_annotation_file, taxon_id,self.mantis_out])
        stdout_file.close()

    ###Executing queue jobs

    def worker_setup_databases(self, queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            if len(record)==3:
                database,force_download,stdout_path = record
                taxon_id=None
            else:
                database,force_download,taxon_id,stdout_path = record
            self.run_download_and_unzip(database,force_download,taxon_id,stdout_path)
            if not self.check_reference_exists(database,taxon_id):
                stdout_file=open(stdout_path,'a+')
                print('Setup failed on',database,taxon_id,flush=True,file=stdout_file)
                stdout_file.close()
                queue.insert(0,record)

    def worker_split_hmms(self,queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            hmm_path,stdout_path = record
            self.split_hmm_into_chunks(hmm_path,stdout_path)

    def worker_extract_NOG_metadata(self,queue,master_pid):
        while True:
            record = queue.pop(0)
            if record is None: break
            target_sql_file, target_annotation_file, taxon_id,stdout_path = record
            self.get_metadata_hmms(target_annotation_file=target_annotation_file,target_sql_file=target_sql_file,taxon_id=taxon_id,stdout_path=stdout_path)

    ###Main download function

    def run_download_and_unzip(self,database,force_download=False,taxon_id=None,stdout_path=None):
        stdout_file=open(stdout_path,'a+')
        if database=='ncbi':        self.download_and_unzip_tax_lineage_dmp(force_download=force_download,stdout_file=stdout_file)
        elif database=='go_obo_nlp':    self.download_go_obo(force_download=force_download,stdout_file=stdout_file)
        elif database=='uniprot_nlp':    self.download_uniprot_reference(force_download=force_download,stdout_file=stdout_file)
        elif database=='pfam':      self.download_and_unzip_pfam_hmm(force_download=force_download,stdout_file=stdout_file)
        elif database=='kofam':     self.download_and_unzip_kofam_hmm(force_download=force_download,stdout_file=stdout_file)
        elif database=='dbcan':     self.download_and_unzip_dbcan_hmm(force_download=force_download,stdout_file=stdout_file)
        elif database=='tigrfam':   self.download_and_unzip_tigrfam_hmm(force_download=force_download,stdout_file=stdout_file)
        elif database=='resfams':   self.download_and_unzip_resfams_hmm(force_download=force_download,stdout_file=stdout_file)
        #elif database=='hamap':     self.download_and_unzip_hamap_hmm(force_download=force_download,stdout_file=stdout_file)
        elif database=='NOGT':      self.download_and_unzip_NOGT(taxon_id=taxon_id,stdout_file=stdout_file)
        elif database=='NOGG':      self.download_and_unzip_NOGG(force_download=force_download,stdout_file=stdout_file)
        elif database=='NOGSQL':    self.download_and_unzip_NOGSQL(force_download=force_download,stdout_file=stdout_file)
        if taxon_id:
            print('Finished downloading',database,' with taxon',taxon_id,flush=True,file=stdout_file)
        else:
            print('Finished downloading', database,flush=True,file=stdout_file)
        stdout_file.close()

    def download_and_unzip_tax_lineage_dmp(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['ncbi']).mkdir(parents=True, exist_ok=True)
        if file_exists(self.mantis_paths['ncbi'] + 'taxidlineage.dmp',force_download):
            print('NCBI taxonomic lineage file already exists! Skipping...',flush=True,file=stdout_file)
            return
        try: os.remove(self.mantis_paths['ncbi'] + 'taxidlineage.dmp')
        except: pass
        url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz'
        if not file_exists(self.mantis_paths['ncbi']+'readme.md'):
            with open(self.mantis_paths['ncbi']+'readme.md','w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This file was downloaded on '+datetime_str+' from:')
                file.write(url)
                file.write('It is used to traceback organism lineage for taxonomically resolved annotations')
        print_cyan('Downloading and unzipping ncbi taxon lineage dump ',flush=True,file=stdout_file)
        download_file(url,output_folder=self.mantis_paths['ncbi'],stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['ncbi'] + 'new_taxdump.tar.gz',extract_path=self.mantis_paths['ncbi'],stdout_file=stdout_file,remove_source=True)
        #deleting the extra files that come with the tar dump
        files_dir = os.listdir(self.mantis_paths['ncbi'])
        for file in files_dir:
            if file!='taxidlineage.dmp':
                os.remove(self.mantis_paths['ncbi']+file)

    def download_go_obo(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['go_obo_nlp']).mkdir(parents=True, exist_ok=True)
        if file_exists(self.mantis_paths['go_obo_nlp'] + 'go.obo',force_download):
            print('Gene ontology obo file already exists! Skipping...',flush=True,file=stdout_file)
            return
        try: os.remove(self.mantis_paths['go_obo_nlp'] + 'go.obo')
        except: pass
        url = 'http://purl.obolibrary.org/obo/go.obo'
        if not file_exists(self.mantis_paths['go_obo_nlp']+'readme.md'):
            with open(self.mantis_paths['go_obo_nlp']+'readme.md','w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This file was downloaded on '+datetime_str+' from:')
                file.write(url)
                file.write('It is used to generate the consensus betweent the annotations')
        print_cyan('Downloading Gene ontology obo',flush=True,file=stdout_file)
        download_file(url,output_folder=self.mantis_paths['go_obo_nlp'],stdout_file=stdout_file)

    def download_uniprot_reference(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['uniprot_nlp']).mkdir(parents=True, exist_ok=True)
        files_dir=os.listdir(self.mantis_paths['uniprot_nlp'])
        if not files_dir:
            print('Please download uniprot reference! It should be in the github repo but you can also get it from the uniprot website.'
                  'Go to:  https://www.uniprot.org/uniprot/?query=reviewed \n'
                  '1- Search for all proteins (empty search field)\n'
                  '2- Click in Columns and choose protein_names and annotation\n'
                  '3- Select only the reviewed proteins\n'
                  '4-Download results in tab separated format\n'
                  'Delete old pickled files and set path to new file in the MANTIS.config',flush=True,file=stdout_file)
            return
        else:
            for file in files_dir:
                if '.gz' in file:
                    uncompress_archive(source_filepath=self.mantis_paths['uniprot_nlp'] + file, stdout_file=stdout_file, remove_source=False)

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

    def is_good_description(self,hmm, row_description):
        temp = [i.lower() for i in row_description.split()]
        if hmm.lower() in temp and 'protein' in temp and len(temp) == 2:
            return False
        if re.search(' [uU]nknown [Ff]unction', row_description): return False
        return True

    def get_common_links_pfam(self,string, res):
        ec = find_ecs(string)
        if ec:
            if 'enzyme_ec' not in res: res['enzyme_ec'] = set()
            res['enzyme_ec'].update(ec)
        tc = find_tcdb(string)
        if tc:
            if 'tcdb' not in res: res['tcdb'] = set()
            res['tcdb'].update(tc)
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

    def build_pfam_line(self,hmm, metadata):
        link_line = ''
        pfam_ids=set(metadata['pfam'])
        metadata['pfam']=set()
        for p_id in pfam_ids:
            metadata['pfam'].add(p_id.split('.')[0])
        for link_type in metadata:
            # print('link',link,flush=True)
            for inner_link in metadata[link_type]:
                link_line += '\t' + link_type + ':' + inner_link
        return hmm + '\t|\t' + link_line + '\n'

    def compile_pfam_metadata(self):
        pfam2go = self.read_pfam2go()
        with open(self.mantis_paths['pfam']+'pfam_metadata.tsv','w+') as pfam_metadata_file:
            with open(self.mantis_paths['pfam']+'Pfam-A.hmm.dat') as pfam_dat_file:
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
                        elif row_header == '#=GF DE' and stop:
                            if hmm in pfam2go:
                                current_metadata = pfam2go[hmm]
                            else:
                                current_metadata = {'description': set()}
                            current_metadata['pfam'] = set()
                            current_metadata['pfam'].add(pfam_accession)
                            if self.is_good_description(hmm, row_description):
                                current_metadata['description'].add(row_description)
                            self.get_common_links_pfam(row_description, current_metadata)
                            stop = False
                            metadata_line = self.build_pfam_line(hmm, current_metadata)
                            pfam_metadata_file.write(metadata_line)
                    line = pfam_dat_file.readline()






    def download_and_unzip_pfam_hmm(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['pfam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists(self.mantis_paths['pfam'] + 'Pfam-A.hmm',force_download) and\
            file_exists(self.mantis_paths['pfam'] + 'Pfam-A.hmm.dat', force_download) and \
            file_exists(self.mantis_paths['pfam'] + 'pfam_metadata.tsv', force_download) and \
            file_exists(self.mantis_paths['pfam'] + 'pfam2go', force_download):
            print('Pfam hmm already exists! Skipping...',flush=True,file=stdout_file)
            return
        pfam_hmm = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
        pfam_metadata = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
        pfam2go='http://current.geneontology.org/ontology/external2go/pfam2go'
        if not file_exists(self.mantis_paths['pfam']+'readme.md'):
            with open(self.mantis_paths['pfam']+'readme.md','w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This hmm was downloaded on '+datetime_str+' from:')
                file.write(pfam_hmm)
                file.write('\nMetadata was downloaded from:\n'+pfam_metadata)
                file.write('\nPfam2GO was downloaded from:\n'+pfam2go)
        print_cyan('Downloading and unzipping Pfam hmms ',flush=True,file=stdout_file)
        to_download=[]
        to_unzip=[]
        if not file_exists(self.mantis_paths['pfam'] + 'Pfam-A.hmm.dat', force_download):
            to_download.append(pfam_hmm)
            to_unzip.append('Pfam-A.hmm.gz')
        if not file_exists(self.mantis_paths['pfam'] + 'pfam_metadata.tsv', force_download):
            to_unzip.append('Pfam-A.hmm.dat.gz')
            to_download.append(pfam_metadata)
        if not file_exists(self.mantis_paths['pfam'] + 'pfam2go', force_download): to_download.append(pfam2go)
        for url in to_download:
            download_file(url,output_folder=self.mantis_paths['pfam'],stdout_file=stdout_file)
        for file_to_unzip in to_unzip:
            uncompress_archive(source_filepath=self.mantis_paths['pfam'] + file_to_unzip,stdout_file=stdout_file,remove_source=True)
        if not self.check_reference_exists('pfam'):
            run_command( 'hmmpress ' + self.mantis_paths['pfam'] + 'Pfam-A.hmm',stdout_file=stdout_file)
        self.compile_pfam_metadata()

    def download_and_unzip_kofam_hmm(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['kofam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('kofam',force_download) and\
            file_exists(self.mantis_paths['kofam'] + 'ko_list', force_download) and\
            file_exists(self.mantis_paths['kofam'] + 'ko2cog.xl', force_download) and\
            file_exists(self.mantis_paths['kofam'] + 'ko2go.xl', force_download) and\
            file_exists(self.mantis_paths['kofam'] + 'ko2tc.xl', force_download) and\
            file_exists(self.mantis_paths['kofam'] + 'ko2cazy.xl', force_download) and\
            file_exists(self.mantis_paths['kofam'] + 'ko_to_path', force_download) and\
            file_exists(self.mantis_paths['kofam'] + 'map_description', force_download):
            print('KOfam hmm already exists! Skipping...',flush=True,file=stdout_file)
            return
        kofam_hmm = 'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz'
        ko_list = 'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz'
        ko_to_cog = 'https://www.kegg.jp/kegg/files/ko2cog.xl'
        ko_to_go = 'https://www.kegg.jp/kegg/files/ko2go.xl'
        ko_to_tc = 'https://www.kegg.jp/kegg/files/ko2tc.xl'
        ko_to_cazy = 'https://www.kegg.jp/kegg/files/ko2cazy.xl'
        ko_to_path = 'http://rest.kegg.jp/link/pathway/ko'
        map_description = 'https://www.genome.jp/kegg/pathway.html'
        if not file_exists(self.mantis_paths['kofam']+'readme.md'):
            with open(self.mantis_paths['kofam']+'readme.md','w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This hmm was downloaded on '+datetime_str+' from:')
                file.write(kofam_hmm)
                file.write('\nMetadata was downloaded from:\n'+ko_list+'\n'+ko_to_cog+'\n'+ko_to_go+'\n'+ko_to_tc+'\n'+ko_to_cazy)
        print_cyan('Downloading and unzipping KOfam hmms ',flush=True,file=stdout_file)
        for url in [kofam_hmm,ko_list, ko_to_cog,ko_to_go,ko_to_tc,ko_to_cazy,ko_to_path,map_description]:
            download_file(url, output_folder=self.mantis_paths['kofam'],stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'profiles.tar.gz', extract_path=self.mantis_paths['kofam'],stdout_file=stdout_file,remove_source=True)
        uncompress_archive(source_filepath=self.mantis_paths['kofam']+'ko_list.gz',stdout_file=stdout_file,remove_source=True)
        move_file(self.mantis_paths['kofam']+'ko',self.mantis_paths['kofam']+'ko_to_path')
        move_file(self.mantis_paths['kofam']+'pathway.html',self.mantis_paths['kofam']+'map_description')
        merge_profiles(self.mantis_paths['kofam']+'profiles/',self.mantis_paths['kofam']+'kofam_merged.hmm',stdout_file=stdout_file)
        run_command('hmmpress ' + self.mantis_paths['kofam'] + 'kofam_merged.hmm', stdout_file=stdout_file)

    def download_and_unzip_dbcan_hmm(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['dbcan']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('dbcan',force_download) and\
            file_exists(self.mantis_paths['dbcan'] + 'CAZyDB.07312019.fam.subfam.ec.txt', force_download):
            print('dbCAN hmm already exists! Skipping...',flush=True,file=stdout_file)
            return
        dbcan_hmm = 'http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt'
        dbcan_metadata = 'http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07312019.fam.subfam.ec.txt'
        if not file_exists(self.mantis_paths['dbcan']+'readme.md'):
            with open(self.mantis_paths['dbcan']+'readme.md','w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This hmm was downloaded on '+datetime_str+' from:')
                file.write(dbcan_hmm)
                file.write('\nMetadata was downloaded from:\n'+dbcan_metadata)
        print_cyan('Downloading and unzipping dbCAN hmms ',flush=True,file=stdout_file)
        for url in [dbcan_metadata,dbcan_hmm]:
            download_file(url, output_folder=self.mantis_paths['dbcan'],stdout_file=stdout_file)
        move_file(self.mantis_paths['dbcan']+'dbCAN-HMMdb-V8.txt',self.mantis_paths['dbcan']+'dbCAN-HMMdb-V8.hmm')
        run_command('hmmpress ' + self.mantis_paths['dbcan'] + 'dbCAN-HMMdb-V8.hmm ',stdout_file=stdout_file)

    def download_and_unzip_tigrfam_hmm(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['tigrfam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('tigrfam',force_download) and\
            file_exists(self.mantis_paths['tigrfam'] + 'gpl.html', force_download) and\
            file_exists(self.mantis_paths['tigrfam'] + 'COPYRIGHT', force_download) and\
            file_exists(self.mantis_paths['tigrfam'] + 'TIGRFAMS_GO_LINK', force_download) and\
            file_exists(self.mantis_paths['tigrfam'] + 'TIGRFAMS_ROLE_LINK', force_download) and\
            file_exists(self.mantis_paths['tigrfam'] + 'TIGR_ROLE_NAMES', force_download):
            print('TIGRfam hmm already exists! Skipping...',flush=True,file=stdout_file)
            return
        tigrfam_hmm = 'ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz'
        tigrfam_go_link = 'ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMS_GO_LINK'
        tigrfam_role_link = 'ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMS_ROLE_LINK'
        tigrfam_role_names = 'ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGR_ROLE_NAMES'
        copyright='ftp://ftp.jcvi.org/pub/data/TIGRFAMs/COPYRIGHT'
        license='http://www.gnu.org/licenses/gpl.html'
        if not file_exists(self.mantis_paths['tigrfam']+'readme.md'):
            with open(self.mantis_paths['tigrfam']+'readme.md','w+') as  file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This hmm was downloaded on '+datetime_str+' from:')
                file.write(tigrfam_hmm)
                file.write('\nMetadata was downloaded from:\n'+tigrfam_go_link+'\n'+tigrfam_role_link+'\n'+tigrfam_role_names+'\n'+copyright)
                file.write('License was downloaded from:\n'+license)
        print_cyan('Downloading and unzipping TIGRfam hmms ',flush=True,file=stdout_file)
        for link in [tigrfam_hmm,
                     tigrfam_go_link,
                     tigrfam_role_link,
                     tigrfam_role_names,
                     copyright,
                     license,
                     ]:
            download_file(link,output_folder=self.mantis_paths['tigrfam'],stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['tigrfam'] + 'TIGRFAMs_15.0_HMM.tar.gz',extract_path=self.mantis_paths['tigrfam'] + 'profiles',stdout_file=stdout_file,remove_source=True)
        merge_profiles(self.mantis_paths['tigrfam']+'profiles/',self.mantis_paths['tigrfam']+'tigrfam_merged.hmm',stdout_file=stdout_file)
        run_command('hmmpress ' + self.mantis_paths['tigrfam'] + 'tigrfam_merged.hmm ',stdout_file=stdout_file)

    def download_and_unzip_resfams_hmm(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['resfams']).mkdir(parents=True, exist_ok=True)
        if not file_exists(self.mantis_paths['resfams'] + '180102_resfams_metadata_updated_v122.tsv', force_download):
            print('Resfams metadata tsv missing, it should be in the github repo!',flush=True,file=stdout_file)
            raise RequirementsNotMet('180102_resfams_metadata_updated_v122.tsv')
        if self.check_reference_exists('resfams',force_download):
            print('Resfams hmm already exists! Skipping...',flush=True,file=stdout_file)
            return
        resfams_hmm = 'http://dantaslab.wustl.edu/resfams/Resfams-full.hmm.gz'
        if not file_exists(self.mantis_paths['resfams']+'readme.md'):
            with open(self.mantis_paths['resfams']+'readme.md','w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This hmm was downloaded on '+datetime_str+' from:')
                file.write(resfams_hmm)
                file.write('\nMetadata was first converted from xlsx to tsv and downloaded directly from this tool\'s git repository')
        print_cyan('Downloading and unzipping Resfams hmms ',flush=True,file=stdout_file)
        download_file(resfams_hmm, output_folder=self.mantis_paths['resfams'],stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['resfams'] + 'Resfams-full.hmm.gz',stdout_file=stdout_file,remove_source=True)
        run_command('hmmpress ' + self.mantis_paths['resfams'] + 'Resfams-full.hmm ',stdout_file=stdout_file)

    def download_and_unzip_NOGG(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['NOGG']).mkdir(parents=True, exist_ok=True)
        target_annotation_file = self.mantis_paths['NOGG'] + 'NOG.annotations.tsv'
        target_merged_hmm = self.mantis_paths['NOGG'] + 'NOGG_merged.hmm'
        # this will download all organisms nog hmms (not specific to taxa)
        if self.check_reference_exists('NOGG',force_download) and file_exists(target_annotation_file,force_download):
            profile_count = get_hmm_profile_count(target_merged_hmm)
            annotations_count = len(self.get_hmms_annotation_file(target_annotation_file, 1))
            # should be the same but some HMMs dont have an annotation (probably an error from NOG)
            if profile_count in range(annotations_count - 10, annotations_count + 10):
                print('NOGG already downloaded, skipping download...', flush=True, file=stdout_file)
                return
        eggnog_downloads_page_hmm = 'http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.hmm.tar.gz'
        eggnog_downloads_page_annot = 'http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/NOG.annotations.tsv.gz'
        if not file_exists(self.mantis_paths['NOGG']+'readme.md'):
            with open(self.mantis_paths['NOGG'] + 'readme.md', 'w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This hmms weas downloaded on ' + datetime_str + ' from:')
                file.write(eggnog_downloads_page_hmm)
                file.write('\nMetadata was downloaded from:\n' + eggnog_downloads_page_annot)
        print_cyan('Downloading and unzipping global NOG hmms ', flush=True, file=stdout_file)
        for url in [eggnog_downloads_page_annot, eggnog_downloads_page_hmm]:
            download_file(url, output_folder=self.mantis_paths['NOGG'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['NOGG'] + 'NOG.annotations.tsv.gz',stdout_file=stdout_file,remove_source=True)
        uncompress_archive(source_filepath=self.mantis_paths['NOGG'] + 'NOG.hmm.tar.gz',extract_path=self.mantis_paths['NOGG'],stdout_file=stdout_file,remove_source=True)
        merge_profiles(self.mantis_paths['NOGG'] + 'NOG_hmm/',target_merged_hmm,stdout_file=stdout_file)
        run_command('hmmpress ' + target_merged_hmm,stdout_file=stdout_file)

    def download_and_unzip_NOGT(self,taxon_id,stdout_file=None):
        folder_path =self.mantis_paths['NOGT']+taxon_id+splitter
        Path(folder_path).mkdir(parents=True, exist_ok=True)
        target_annotation_file =self.mantis_paths['NOGT'] + taxon_id +splitter + taxon_id + '_annotations.tsv'
        target_merged_hmm = self.mantis_paths['NOGT'] + taxon_id +splitter + taxon_id + '_merged.hmm'
        eggnog_downloads_page = 'http://eggnog5.embl.de/download/latest/per_tax_level/'+str(taxon_id)+'/'
        for file in ['_annotations.tsv.gz','_hmms.tar']:
            url=eggnog_downloads_page+taxon_id+file
            download_file(url, output_folder=folder_path,stdout_file=stdout_file)
        if os.path.exists(folder_path+'profiles'): shutil.rmtree(folder_path+'profiles')
        uncompress_archive(source_filepath=folder_path+taxon_id+ '_hmms.tar',extract_path=folder_path+'profiles',stdout_file=stdout_file,remove_source=True)
        uncompress_archive(source_filepath=folder_path+taxon_id+ '_annotations.tsv.gz',stdout_file=stdout_file,remove_source=True)
        # removing gz extension
        hmm_files = os.listdir(folder_path+'profiles')
        for hmm_profile in hmm_files:
            if '.hmm' in hmm_profile: move_file(hmm_profile,hmm_profile.strip('.gz'))
        merge_profiles(folder_path+'profiles/'+taxon_id,folder_path+taxon_id+'_merged.hmm',stdout_file=stdout_file)
        if os.path.exists(folder_path+'profiles'):      shutil.rmtree(folder_path+'profiles')
        run_command('hmmpress ' + target_merged_hmm ,stdout_file=stdout_file)

    def download_and_unzip_NOGSQL(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['default']).mkdir(parents=True, exist_ok=True)
        if file_exists(self.mantis_paths['default'] + 'eggnog.db',force_download):
            print('eggnog.db already exists! Skipping...',flush=True,file=stdout_file)
            return
        url = 'http://eggnogdb.embl.de/download/emapperdb-5.0.0/eggnog.db.gz'
        download_file(url, output_folder=self.mantis_paths['default'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['default'] + 'eggnog.db.gz',stdout_file=stdout_file,remove_source=True)

    def unpack_NOG_sql(self,stdout_file=None):
        '''
        splitting into different folders for github pushing (without lfs restrictions)
        x=("NOGT1" "NOGT2" "NOGT3" "NOGT4" "NOGT5")
        c=0
        for f in *
        do
            mv "$f" "${x[c]}"
            c=$(( (c+1)%5 ))
        done

        split into 6 subfolders(5 NOGT + NOGG) which should be compressed with:
        export GZIP=-9
        then in the NOGT folder do:
        for i in *; do tar cvzf $i.tar.gz $i; done
        '''
        resources_path=self.mantis_paths['resources']
        NOG_sql_folder=add_slash(resources_path+'NOG_sql')
        for compressed_sql in os.listdir(NOG_sql_folder):
            if '.tar.gz' in compressed_sql:
                uncompress_archive(source_filepath=NOG_sql_folder+compressed_sql, extract_path=NOG_sql_folder, stdout_file=stdout_file,remove_source=False)
        NOGT_sql_path=add_slash(add_slash(resources_path+'NOG_sql')+'NOGT')
        NOGG_sql_path=add_slash(add_slash(resources_path+'NOG_sql')+'NOGG')
        Path(NOGT_sql_path).mkdir(parents=True, exist_ok=True)
        for sql_folder in os.listdir(NOG_sql_folder):
            if re.search('NOGT\d',sql_folder) and 'tar.gz' not in sql_folder:
                current_nogt_folder=add_slash(add_slash(NOG_sql_folder)+sql_folder)
                for inner_dir in os.listdir(current_nogt_folder):
                    move_file(current_nogt_folder + inner_dir, NOGT_sql_path)
                if os.path.exists(current_nogt_folder):
                    shutil.rmtree(current_nogt_folder)
        #moving the files to their respective folders
        if self.mantis_paths['NOGG']!='NA':
            Path(add_slash(self.mantis_paths['NOGG'])).mkdir(parents=True, exist_ok=True)
            move_file(NOGG_sql_path+'NOGG_sql_annotations.tsv',self.mantis_paths['NOGG']+'NOGG_sql_annotations.tsv')
        if self.mantis_paths['NOGT']!='NA':
            for tax_folder in os.listdir(NOGT_sql_path):
                if os.path.isdir(NOGT_sql_path+tax_folder):
                    for f in os.listdir(add_slash(NOGT_sql_path+tax_folder)):
                        Path(add_slash(self.mantis_paths['NOGT']+tax_folder)).mkdir(parents=True, exist_ok=True)
                        move_file(add_slash(NOGT_sql_path+tax_folder)+f,add_slash(self.mantis_paths['NOGT']+tax_folder)+f)
        #removing empty folders
        if os.path.exists(NOGG_sql_path):
            shutil.rmtree(NOGG_sql_path)
        if os.path.exists(NOGT_sql_path):
            shutil.rmtree(NOGT_sql_path)

    ### Support functions for setting up queue to download NOGT hmms

    def get_taxon_ids_NOGT(self,url):
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

    def get_taxon_for_queue_NOGT(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['NOGT']).mkdir(parents=True, exist_ok=True)
        # we will download tax specific hmms, which can be used for more specific hmmscans
        eggnog_downloads_page = 'http://eggnog5.embl.de/download/latest/per_tax_level/'
        if not file_exists(self.mantis_paths['NOGT']+'readme.md'):
            with open(self.mantis_paths['NOGT']+'readme.md','w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('These hmms were downloaded on '+datetime_str+' from:')
                file.write(eggnog_downloads_page)
        taxon_ids = self.get_taxon_ids_NOGT(eggnog_downloads_page)
        res=[]
        for taxon_id in taxon_ids:
            if taxon_id!='1':
                target_annotation_file = self.mantis_paths['NOGT'] + taxon_id + splitter + taxon_id + '_annotations.tsv'
                if self.check_reference_exists('NOGT',taxon_id,force_download) and file_exists(target_annotation_file,force_download):
                    hmm_path = get_hmm_in_folder(self.mantis_paths['NOGT'] + taxon_id + splitter)
                    profile_count=get_hmm_profile_count(hmm_path)
                    annotations_count=len(self.get_hmms_annotation_file(target_annotation_file, 1))
                    # should be the same but some HMMs dont have an annotation (probably an error from NOG)
                    if profile_count in range(annotations_count-10,annotations_count+10+1):
                        print('Tax ' + taxon_id + ' hmm and annotations.tsv files already exist, and they were merged correctly. Skipping setup...',flush=True, file=stdout_file)
                    else:
                        print('Tax ' + taxon_id + ' hmm already exists but was not merged correctly!',flush=True, file=stdout_file)
                        res.append(taxon_id)
                else:
                    res.append(taxon_id)
        print('Will compile data for the following NOGT:\n' + ','.join(res), flush=True,file=stdout_file)
        return res

    ### Support functions for setting up queue for split hmms

    def get_taxon_for_queue_NOG_split_hmms(self):
        res = []
        stdout_file=open(self.mantis_out,'a+')
        if self.mantis_paths['NOGT'][0:2] != 'NA':
            taxon_ids = os.listdir(self.mantis_paths['NOGT'])
            for taxon_id in taxon_ids:
                if taxon_id != '1' and os.path.isdir(self.mantis_paths['NOGT'] + taxon_id):
                    hmm_path = get_hmm_in_folder(self.mantis_paths['NOGT'] + taxon_id + splitter)
                    print('Checking NOG for splitting:',hmm_path,taxon_id,flush=True,file=stdout_file)
                    if 'chunks' not in os.listdir(self.mantis_paths['NOGT'] + taxon_id):
                        profile_count = get_hmm_profile_count(hmm_path)
                        if profile_count > self.hmm_chunk_size:
                            res.append(hmm_path)
                        else:   print('NOG already split:', hmm_path, taxon_id, flush=True, file=stdout_file)
                    else:   print('NOG already split:', hmm_path, taxon_id, flush=True, file=stdout_file)

        if self.mantis_paths['NOGG'][0:2] != 'NA':
            hmm_path = get_hmm_in_folder(self.mantis_paths['NOGG'])
            print('Checking NOG for splitting:', hmm_path, flush=True, file=stdout_file)
            if 'chunks' not in os.listdir(self.mantis_paths['NOGG']):
                profile_count = get_hmm_profile_count(hmm_path)
                if profile_count > self.hmm_chunk_size:
                    res.append(hmm_path)
                else:    print('NOG already split:', hmm_path, flush=True, file=stdout_file)
            else:   print('NOG already split:', hmm_path, flush=True, file=stdout_file)

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
        profile_count = get_hmm_profile_count(hmm_path)
        hmm_folder = get_folder(hmm_path)
        hmm_chunks_folder = hmm_folder + 'chunks/'
        Path(hmm_chunks_folder).mkdir(parents=True, exist_ok=True)
        load_balanced_chunks = chunk_generator_load_balanced([i for i in range(profile_count)], self.hmm_chunk_size)
        print('Splitting ' + hmm_path + ' into ' + str(len(load_balanced_chunks)) + ' chunks.', flush=True,
              file=stdout_file)
        chunks_name = get_path_level(hmm_path, remove_extension=True) + '_chunk_'
        with open(hmm_path) as hmm_file:
            for chunk_i in range(len(load_balanced_chunks)):
                chunk_file_path = hmm_chunks_folder + chunks_name + str(chunk_i) + '.hmm'
                with open(chunk_file_path, 'w+') as chunk_file:
                    for profile_i in range(len(load_balanced_chunks[chunk_i])):
                        profile = ''.join(read_profile(hmm_file))
                        chunk_file.write(profile)
        chunks_dir = os.listdir(hmm_chunks_folder)
        for chunk in chunks_dir:
            run_command('hmmpress ' + hmm_chunks_folder + chunk, stdout_file=stdout_file)
        stdout_file.close()

    ### Support functions for extracting metadata

    def get_hmms_annotation_file(self,annotations_file,hmm_col):
        res=set()
        if not file_exists(annotations_file): return  res
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

    def get_metadata_hmms_thread(self,hmm,taxon_id,hmm_annotations_path,eggnog_db):

        codes_to_exclude = ['IEA']
        if taxon_id:     taxon_query = '@' + str(taxon_id)
        else:            taxon_query = ''
        connection = sqlite3.connect(eggnog_db)  # open connection to your database
        sql_command = 'SELECT eggnog.groups,gene_ontology.gos,kegg.ec, kegg.ko, kegg.pathway, kegg.module, kegg.reaction, kegg.rclass, kegg.brite, kegg.tc, kegg.cazy, bigg.reaction FROM eggnog LEFT JOIN seq on seq.name = eggnog.name LEFT JOIN gene_ontology on gene_ontology.name = eggnog.name LEFT JOIN kegg on kegg.name = eggnog.name LEFT JOIN bigg on bigg.name = eggnog.name WHERE eggnog.groups LIKE "%' + hmm + taxon_query + '%";'
        cursor = connection.cursor()
        cursor.execute(sql_command)
        rows = cursor.fetchall()
        res = {'go': set(), 'enzyme_ec': set(), 'kegg_ko': set(), 'kegg_pathway': set(), 'kegg_module': set(),
               'kegg_reaction': set(), 'kegg_rclass': set(), 'kegg_brite': set(), 'kegg_cazy': set(),
               'bigg_reaction': set(),'description':set()}
        for row in rows:
            hmms, gene_ontology_gos, kegg_ec, kegg_ko, kegg_pathway, kegg_module, kegg_reaction, kegg_rclass, kegg_brite, kegg_tc, kegg_cazy, bigg_reaction = row
            if kegg_ko: kegg_ko=kegg_ko.replace('ko:','')
            for query_hmm_taxon in hmms.split(','):
                query_hmm,query_taxon = query_hmm_taxon.split('@')
                if not taxon_id: passed_test=True
                elif int(taxon_id) == int(query_taxon): passed_test=True
                else: passed_test=False
                if passed_test:
                    if gene_ontology_gos:
                        gene_ontology_gos_copy = str(gene_ontology_gos)
                        gene_ontology_gos_copy = gene_ontology_gos_copy.split(',')
                        for go_group in gene_ontology_gos_copy:
                            if '|' not in go_group: print(sql_command,'\n',go_group)
                            _, go, evidence_code = go_group.split('|')
                            if evidence_code not in codes_to_exclude:
                                res['go'].add(go.split(':')[1])
                    if kegg_ec: res['enzyme_ec'].update(kegg_ec.split(','))
                    if kegg_ko: res['kegg_ko'].update(kegg_ko.split(','))
                    if kegg_pathway: res['kegg_pathway'].update(kegg_pathway.split(','))
                    if kegg_module: res['kegg_module'].update(kegg_module.split(','))
                    if kegg_reaction: res['kegg_reaction'].update(kegg_reaction.split(','))
                    if kegg_rclass: res['kegg_rclass'].update(kegg_rclass.split(','))
                    if kegg_brite: res['kegg_brite'].update(kegg_brite.split(','))
                    if kegg_cazy: res['kegg_cazy'].update(kegg_cazy.split(','))
                    if bigg_reaction: res['bigg_reaction'].update(bigg_reaction.split(','))
        if taxon_id:
            description = self.get_free_text_NOG(hmm=hmm,hmm_annotations_path= hmm_annotations_path, description_col=3)
        else:
            description = self.get_free_text_NOG(hmm=hmm,hmm_annotations_path=hmm_annotations_path, description_col=5)
        if description:   res['description'].add(description)
        return [hmm,res]

    def launch_pool(self,n_processes, function_to_process, function_args):
        pool = IO_Pool(n_processes)
        res = pool.starmap(function_to_process, function_args)
        pool.close()
        pool.join()
        return res

    def create_index_eggnog_db(self):
        try:
            eggnog_db_path=self.mantis_paths['default'] + "eggnog.db"
            sql_command='create index eggnog_name on eggnog (name);'
            connection = sqlite3.connect(eggnog_db_path)  # open connection to your database
            cursor = connection.cursor()
            cursor.execute(sql_command)
        except:
            pass

    def get_metadata_hmms(self,target_annotation_file,target_sql_file,taxon_id=None,stdout_path=None):
        max_threads=100
        chunk_size=500
        eggnog_db_path=self.mantis_paths['default'] + "eggnog.db"
        hmm_list = self.get_hmms_annotation_file(target_annotation_file, hmm_col=1)
        stdout_file=open(stdout_path,'a+')

        if os.path.exists(target_sql_file):
            already_compiled_hmms_list=self.get_hmms_annotation_file(target_sql_file,hmm_col=0)
            if len(already_compiled_hmms_list)==len(hmm_list):
                print('Target sql annotation file already finished, skipping!',target_sql_file,flush=True, file=stdout_file)
                return
            else:
                for hmm in already_compiled_hmms_list:
                    if hmm in hmm_list:
                        hmm_list.remove(hmm)
                if not hmm_list:
                    print('Target sql annotation file already complete, skipping!', target_sql_file, flush=True,file=stdout_file)
                    return
                else:
                    print('Target sql annotation file already exists but is incomplete, continuing...', target_sql_file,flush=True, file=stdout_file)

        print('Starting metadata extraction into',target_sql_file,flush=True,file=stdout_file)
        # http://geneontology.org/docs/guide-go-evidence-codes/
        while not os.path.exists(self.mantis_paths['default']+"eggnog.db"):sleep(5)
        print('Exporting metadata for NOG',taxon_id,flush=True,file=stdout_file)
        if taxon_id:     hmm_annotations_path = self.mantis_paths['NOGT'] + taxon_id + '/' + taxon_id + '_annotations.tsv'
        else:            hmm_annotations_path = self.mantis_paths['NOGG'] + 'NOG.annotations.tsv'
        pool_args=[]
        for hmm in hmm_list:
            pool_args.append([hmm,taxon_id,hmm_annotations_path,eggnog_db_path])
        if len(pool_args)<max_threads: n_threads= len(pool_args)
        else: n_threads= max_threads
        #to avoid keeping everything in memory, we chunk it
        chunks = chunk_generator(pool_args, chunk_size)
        for chunk in chunks:
            hmm_metadata = self.launch_pool(n_threads,self.get_metadata_hmms_thread,chunk)
            with open(target_sql_file, "a+") as file:
                for hmm_links in hmm_metadata:
                    link_line = ''
                    hmm,metadata=hmm_links
                    #print('hmm_links',hmm_links,flush=True)
                    for link_type in metadata:
                        #print('link',link,flush=True)
                        for inner_link in metadata[link_type]:
                            link_line+='\t'+link_type+':'+inner_link
                    file.write(hmm + '\t|\t' + link_line+'\n')
        print('Finished exporting metadata for NOG', taxon_id,flush=True,file=stdout_file)
        stdout_file.close()


    #LEGACY
    def remove_pattern_hamap(self, string_to_search, pattern, part_to_remove=''):
        patterns_removed = set()
        search = re.search(pattern, string_to_search)
        while search:
            patterns_removed.add(search.group().replace(part_to_remove, ''))
            start = search.span()[0]
            end = search.span()[1]
            if string_to_search[start + 1] == '(': start += 2
            if string_to_search[end - 1] == ')': end -= 1
            string_to_search = list(string_to_search)
            string_to_search[start:end] = ''
            string_to_search = ''.join(string_to_search)
            search = re.search(pattern, string_to_search)
        return patterns_removed

    def parse_metadata_hamap(self, rules_file):
        res = {}
        function_str = ''
        catalytic_str = ''
        with open(rules_file) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.strip(';')
                if 'AC   ' in line:
                    line = line.replace('AC   ', '')
                    hamap_id = line
                    if res:
                        yield res
                    res = {'pfam': set(), 'tigrfam': set(), 'prosite': set(), 'prints': set(), 'pirsf': set(),
                           'smart': set(), 'description': set(), 'go': set(), 'rhea': set(), 'chebi': set(),
                           'enzyme_ec': set()}
                    res['hamap'] = hamap_id
                elif 'GO   ' in line:
                    line = line.replace('GO   ', '')
                    line = line.split(';')[0]
                    line = line.replace('GO:', '')
                    res['go'].add(line)
                elif 'DR   ' in line:
                    line = line.replace('DR   ', '')
                    line = line.split(';')
                    if line[0] == 'Pfam':
                        res['pfam'].add(line[1])
                    elif line[0] == 'TIGRFAMs':
                        res['tigrfam'].add(line[1])
                    elif line[0] == 'PROSITE':
                        res['prosite'].add(line[1])
                    elif line[0] == 'PRINTS':
                        res['prints'].add(line[1])
                    elif line[0] == 'PIRSF':
                        res['pirsf'].add(line[1])
                    elif line[0] == 'SMART':
                        res['smart'].add(line[1])
                elif 'CC   ' in line:
                    line = line.replace('CC   ', '')
                    if '-!- FUNCTION:' in line:
                        line = line.replace('-!- FUNCTION:', '')
                        record_function = True
                        function_str = ''
                    elif '-!-' in line and '-!- FUNCTION:' not in line:
                        if function_str:
                            res['description'].add(function_str)
                        function_str = ''
                        record_function = False
                    if record_function:
                        function_str += ' ' + line.strip()

                    if '-!- CATALYTIC ACTIVITY:' in line:
                        line = line.replace('-!- CATALYTIC ACTIVITY:', '')
                        record_catalyic = True
                        catalytic_str = ''
                    elif '-!-' in line and '-!- CATALYTIC ACTIVITY:' not in line:
                        if catalytic_str:
                            ec_pattern = re.compile('EC=\d(\.(-|\d{1,3}|([a-zA-Z]\d{1,3}))){2,3}')
                            chebi_pattern = re.compile('CHEBI:\d+')
                            rhea_pattern = re.compile('RHEA-COMP:\d+')
                            removed_ids = self.remove_pattern_hamap(catalytic_str, ec_pattern, part_to_remove='EC=')
                            res['enzyme_ec'].update(removed_ids)
                            removed_ids = self.remove_pattern_hamap(catalytic_str, chebi_pattern, part_to_remove='CHEBI:')
                            res['chebi'].update(removed_ids)
                            removed_ids = self.remove_pattern_hamap(catalytic_str, rhea_pattern, part_to_remove='RHEA-COMP:')
                            res['rhea'].update(removed_ids)
                        catalytic_str = ''
                        record_catalyic = False
                    if record_catalyic:
                        catalytic_str += ' ' + line.strip()
                line = file.readline()
        if res: yield res

    def compile_line_hamap(self, metadata):
        headers = ['enzyme_ec', 'pfam', 'tigrfam', 'prosite', 'prints', 'pirsf', 'smart', 'go', 'rhea', 'chebi']
        line = [metadata['hamap']]
        for header in headers:
            line.append(';'.join([i.strip() for i in metadata[header]]))
        line.append(' '.join(metadata['description']))
        return '\t'.join(line)

    def compile_metadata_hamap(self,stdout_file=None):
        hamap_metadata='ftp://ftp.expasy.org/databases/hamap/hamap_rules.dat'
        download_file(hamap_metadata, output_folder=self.mantis_paths['hamap'],stdout_file=stdout_file)
        metadata_file=self.mantis_paths['hamap']+'hamap.tsv'
        with open(metadata_file, 'w+') as file:
            metadata_generator = self.parse_metadata_hamap(self.mantis_paths['hamap']+'hamap_rules.dat')
            headers = ['hamap', 'enzyme_ec', 'pfam', 'tigrfam', 'prosite', 'prints', 'pirsf', 'smart', 'go', 'rhea',
                       'chebi', 'description']
            headers = '\t'.join(headers)
            file.write(headers + '\n')
            for profile in metadata_generator:
                line = self.compile_line_hamap(profile)
                file.write(line + '\n')

    def download_and_unzip_hamap_hmm(self,force_download=False,stdout_file=None):
        Path(self.mantis_paths['hamap']).mkdir(parents=True, exist_ok=True)
        target_merged_hmm=self.mantis_paths['hamap']+'hamap.hmm'
        if not file_exists(self.mantis_paths['hamap'] + 'hamap.tsv', force_download):
            print('HAMAP metadata tsv missing!',flush=True,file=stdout_file)
            self.compile_metadata_hamap(stdout_file=stdout_file)
        if self.check_reference_exists('hamap',force_download):
            print('HAMAP hmm already exists! Skipping...',flush=True,file=stdout_file)
            return
        hamap_alignments='ftp://ftp.expasy.org/databases/hamap/hamap_alignments.tar.gz'
        if not file_exists(self.mantis_paths['hamap']+'readme.md'):
            with open(self.mantis_paths['hamap']+'readme.md','w+') as file:
                datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                file.write('This hmm was downloaded on '+datetime_str+' from:')
                file.write(hamap_alignments)
        print_cyan('Downloading and unzipping HAMAP alignments hmms ',flush=True,file=stdout_file)
        #getting the multiple alignments
        download_file(hamap_alignments, output_folder=self.mantis_paths['hamap'],stdout_file=stdout_file)
        alignments_dir=add_slash(self.mantis_paths['hamap'] + 'hamap_alignments')
        if os.path.exists(alignments_dir):      shutil.rmtree(alignments_dir)
        uncompress_archive(source_filepath=self.mantis_paths['hamap'] + 'hamap_alignments.tar.gz',stdout_file=stdout_file,remove_source=True)
        #building profiles out of multiple alignments
        msa_files = os.listdir(alignments_dir)
        for msa in msa_files:
            profile_name=msa.replace('.msa','')
            run_command('hmmbuild '+alignments_dir+profile_name+'.hmm ' + alignments_dir+msa , stdout_file=stdout_file)
            os.remove(alignments_dir+msa)
        hmm_files = os.listdir(alignments_dir)
        merge_profiles(alignments_dir,target_merged_hmm,stdout_file=stdout_file)
        #pressing merged profiles
        run_command('hmmpress ' +target_merged_hmm,stdout_file=stdout_file)

