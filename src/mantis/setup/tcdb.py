import os
import re
from datetime import datetime
from pathlib import Path

from mantis.src.metadata.uniprot_api import UnitprotApi
from mantis.src.setup.base import SetupBase
from mantis.src.utils.downloader import download_file
from mantis.src.utils.logger import logger
from mantis.src.utils.utils import chunk_generator, read_protein_fasta_generator, run_command, write_fasta_generator


class SetupTcdb(SetupBase):
    TCDB_SEQS = 'http://www.tcdb.org/public/tcdb'
    TCDB_GO = 'https://www.tcdb.org/cgi-bin/projectv/public/go.py'
    TCDB_PFAM = 'https://www.tcdb.org/cgi-bin/projectv/public/pfam.py'
    DATABASE = 'tcdb'

    def __init__(self, config_file: str) -> None:
        super().__init__(database=self.DATABASE, config_file=config_file)


    def parse_tsv_tcdb(self, file_path, key_col, val_col, res_key, res, val_clean_function=None):
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.split('\t')
                if line[key_col] not in res:
                    res[line[key_col]] = {}
                if res_key not in res[line[key_col]]:
                    res[line[key_col]][res_key] = set()
                if val_clean_function:
                    line[val_col] = val_clean_function(line[val_col])
                res[line[key_col]][res_key].add(line[val_col])
                line = file.readline()

    def read_tcdb_headers(self) -> dict[dict: [str, str]]:
        file_path = os.path.join(self.output_folder, 'tcdb')
        res = {}
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                if '>' in line:
                    line = line.split('|')
                    uniprot_accession = line[2]
                    tc_id = line[3].split()[0]
                    description = line[3].replace(tc_id, '').strip()
                    if re.search('([Hh]ypothetical|[Uu]ncharacterized|[Uu]ndetermined)', description):
                        description = ''
                    description = description.split('[')[0]
                    if re.search('[A-Z]+=', description):
                        description = description.split(re.search('[A-Z]+=', description).group())[0]
                    if ' - ' in description:
                        description = description.split(' - ')[0]
                    description = description.strip()
                    if description:
                        res[uniprot_accession] = {'tcdb': tc_id, 'description': {description.strip()}}
                    else:
                        res[uniprot_accession] = {'tcdb': tc_id}
                line = file.readline()
        return res

    def remove_bad_entries(self, sequences_list) -> set[str]:
        # some proteins are not in uniprot, but are in some other dbs.... Im not sure if they should be removed, but we will consider uniprot as the central repo, so we will remove them
        chunks_post = chunk_generator(to_chunk=sequences_list, chunk_size=500)
        all_found = set()
        uniprot_api = UnitprotApi()
        for seqs_chunk in chunks_post:
            job_id = uniprot_api.submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=seqs_chunk)
            if uniprot_api.check_id_mapping_results_ready(job_id):
                link = uniprot_api.get_id_mapping_results_link(job_id)
                results = uniprot_api.get_id_mapping_results_search(link)
                if 'failedIds' in results:
                    failed_ids = results['failedIds']
                else:
                    failed_ids = []
                for seq in seqs_chunk:
                    if seq not in failed_ids:
                        all_found.add(seq)
        return all_found

    def yield_tcdb_seqs(self, tcdb_fasta, seqs_found):
        tcdb_seqs = read_protein_fasta_generator(tcdb_fasta)
        for seq_id, seq in tcdb_seqs:
            uniprot_accession = seq_id.split('|')[2]
            if uniprot_accession in seqs_found:
                yield uniprot_accession, seq

    def generate_tcdb_fasta(self, tcdb_fasta, seqs_found):
        write_fasta_generator(seqs_generator=self.yield_tcdb_seqs(tcdb_fasta=tcdb_fasta, seqs_found=seqs_found),
                              fasta_file=os.path.join(self.output_folder, 'tcdb.faa'))

    def compile_tcdb_metadata(self):
        # acession will be the key
        sequences_to_metadata = self.read_tcdb_headers()
        valid_sequences = self.remove_bad_entries(list(sequences_to_metadata.keys()))
        self.generate_tcdb_fasta(tcdb_fasta=os.path.join(self.output_folder, 'tcdb'), sequences=valid_sequences)
        # here tc ids will be keys
        metadata = {}
        self.parse_tsv_tcdb(file_path=os.path.join(self.output_folder, 'go.py'),
                            key_col=1,
                            val_col=0,
                            res_key='go',
                            res=metadata,
                            val_clean_function=lambda a: a.replace('GO:', ''))
        self.parse_tsv_tcdb(file_path=os.path.join(self.output_folder, 'pfam.py'),
                            key_col=1,
                            val_col=0,
                            res_key='pfam',
                            res=metadata)

        # we now add the tc specific metadata to each acession
        metadata_to_write = {}
        for seq in valid_sequences:
            tc_id = sequences_to_metadata[seq]['tcdb']
            if tc_id in metadata:
                metadata_to_write[seq] = metadata[tc_id]
            else:
                metadata_to_write[seq] = {}
            if 'description' in sequences_to_metadata[seq]:
                metadata_to_write[seq]['description'] = sequences_to_metadata[seq]['description']
            metadata_to_write[seq]['tcdb'] = {tc_id}

        self.write_metadata(metadata=metadata_to_write, metadata_file=os.path.join(self.output_folder, 'metadata.tsv'))

    def write_readme(self):
        with open(self.output_folder + 'readme.md', 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'These sequences were downloaded on {datetime_str} from:\n{self.TCDB_SEQS}\n'
                       f'Metadata was downloaded from:\n{self.TCDB_GO}\n{self.TCDB_PFAM}\n')

    def create_diamond_db(self):
        if not os.path.exists(os.path.join(self.output_folder, 'tcdb.dmnd')):
            run_command(commnad=f'diamond makedb --in {self.output_folder} tcdb.faa -d {os.path.join(self.output_folder, 'tcdb')}')


    def run(self):
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('tcdb'):
            logger.info('TCDB sequences already exists! Skipping...', flush=True)
            return
        logger.info('Downloading and unzipping TCDB sequences', flush=True)
        for link in [self.TCDB_SEQS,
                     self.TCDB_GO,
                     self.TCDB_PFAM,
                     ]:
            download_file(link, output_folder=self.output_folder)
        self.write_readme()
        self.compile_tcdb_metadata()
        self.create_diamond_db()
        logger.info(f'Finished downloading {self.DATABASE}')
