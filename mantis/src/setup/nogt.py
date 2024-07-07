import os
import shutil
from pathlib import Path

from mantis.src.metadata.utils import get_common_links_metadata
from mantis.src.utils import print_cyan, remove_file


class Setup():



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
