import os
import shutil
from pathlib import Path

from mantis.src.utils import print_cyan


class Setup():



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
