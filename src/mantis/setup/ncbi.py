import os
import shutil
from datetime import datetime
from pathlib import Path

from mantis.src.metadata.utils import get_common_links_metadata
from mantis.src.setup.base import SetupBase
from mantis.src.utils.downloader import download_file
from mantis.src.utils.file_processer import copy_file, remove_file, uncompress_archive
from mantis.src.utils.logger import logger
from mantis.src.utils.utils import merge_profiles, run_command


class SetupNcbi(SetupBase):
    HMM_URL = 'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz'
    METADATA_URL = 'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv'
    DATABASE = 'ncbi'
    NCBI_DOMAINS = [
            '1',  # top level in NOG
            '2157',  # Archaea
            '2',  # Bacteria
            '2759',  # Eukaryota
            '10239',  # Viruses
            '28384',  # Others
            '12908',  # Unclassified
        ]

    def __init__(self, config_file: str) -> None:
        super().__init__(database=self.DATABASE, config_file=config_file)


    def write_readme(self):
        with open(os.path.join(self.output_folder, 'readme.md', 'w+')) as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'This hmm was downloaded on {datetime_str} from:\n{self.HMM_URL}\nMetadata was downloaded from:\n{self.METADATA_URL}')

    def run(self,):
        # we cant verify a priori which foulders we should have, so you need to delete the folder to restart
        if self.check_reference_exists('NCBI') and \
                os.path.exists(self.mantis_paths['NCBI'] + 'readme.md'):
            logger.info('NCBI hmm folder already exists! Skipping...')
            return


        logger.info('Downloading and unzipping NCBI hmms')
        for url in [self.HMM_URL, self.METADATA_URL]:
            download_file(url=url, output_folder=self.output_folder)
        compressed_hmm_path = os.path.join(self.output_folder, 'hmm_PGAP.HMM.tgz')
        uncompress_archive(source_filepath=compressed_hmm_path,
                           extract_path=self.output_folder, remove_source=True)
        self.compile_hmms_NCBI()
        remove_file(self.mantis_paths['NCBI'] + 'hmm_PGAP.tsv')
        if os.path.exists(self.mantis_paths['NCBI'] + 'hmm_PGAP'):
            shutil.rmtree(self.mantis_paths['NCBI'] + 'hmm_PGAP')



    # sorting profiles by the taxonomic id from ncbi metadata file
    def compile_hmms_NCBI(self):
        sorted_metadata = self.sort_hmms_NCBI()
        self.write_metadata(sorted_metadata)
        self.assign_hmm_profiles(sorted_metadata)

    def assign_hmm_profiles(self, sorted_metadata):
        for taxa in sorted_metadata:
            folder_path = os.path.join(self.output_folder, taxa, 'to_merge')
            for hmm, _, _, _, _, _ in sorted_metadata[taxa]:
                source_file = os.path.join(self.output_folder, 'hmm_PGAP', f'{hmm}.HMM')
                dest_file = os.path.join(self.output_folder, taxa, 'to_merge', f'{hmm}.HMM')
                try:
                    copy_file(source_file=source_file,
                              dest_file=dest_file)
                except Exception as e:
                    logger.exception(e)
            if os.listdir(folder_path):
                hmm_name = f'{taxa}_merged.hmm'
                hmm_path = os.path.join(self.output_folder, taxa, hmm_name)
                merge_profiles(folder_path=folder_path,
                               output_file=hmm_path)
                run_command(f'hmmpress {hmm_path}')
            else:
                shutil.rmtree(self.mantis_paths['NCBI'] + taxa)

    def write_metadata(self, sorted_metadata):
        for taxa in sorted_metadata:
            Path(self.mantis_paths['NCBI'] + taxa).mkdir(parents=True, exist_ok=True)
            Path(os.path.join(self.output_folder, taxa, 'to_merge')).mkdir(parents=True, exist_ok=True)
            with open(os.path.join(self.output_folder, taxa, + 'metadata.tsv'), 'w+') as metadata_file:
                for _, hmm_label, description, enzyme_ec, go_terms, common_links in sorted_metadata[taxa]:
                    line = [hmm_label, '|', f'description:{description}']

                    for db in common_links:
                        for db_id in common_links[db]:
                            if f'{db}:{db_id}' not in line:
                                line.append(f'{db}:{db_id}')
                    for ec in enzyme_ec:
                        if f'enzyme_ec:{ec}' not in line:
                            line.append(f'enzyme_ec:{ec}')
                    for go_term in go_terms:
                        if f'go:{go_term}' not in line:
                            line.append(f'go:{go_term}')
                    # ncbi also contains tigrfam hmms
                    if hmm_label.startswith('TIGR'):
                        if f'tigrfam:{hmm_label}' not in line:
                            line.append(f'tigrfam:{hmm_label}')

                    line = '\t'.join(line) + '\n'
                    metadata_file.write(line)


    def sort_hmms_NCBI(self):
        res = {'NCBIG': []}
        already_added_NCBIG = set()
        metadata = self.mantis_paths['NCBI'] + 'hmm_PGAP.tsv'
        with open(metadata) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.split('\t')
                hmm, hmm_label, description, enzyme_ec, go_terms, taxa_id = line[0], line[2], line[10], line[12], line[13], line[15]
                common_links = {}
                get_common_links_metadata(input_string=description,
                                          metadata_dict=common_links)
                for db in common_links:
                    for db_id in common_links[db]:
                        description = description.replace(db_id, '').strip()
                description = description.replace('(Provisional)', '')
                description = description.strip()
                enzyme_ec = [i for i in enzyme_ec.split(',') if i]
                go_terms = [i.replace('GO:', '') for i in go_terms.split(',') if i]
                line = file.readline()
                if taxa_id:
                    if taxa_id not in res:
                        res[taxa_id] = []
                    res[taxa_id].append([hmm, hmm_label, description, enzyme_ec, go_terms, common_links])
                    if taxa_id in self.NCBI_DOMAINS:
                        if hmm not in already_added_NCBIG:
                            res['NCBIG'].append([hmm, hmm_label, description, enzyme_ec, go_terms, common_links])
                            already_added_NCBIG.add(hmm)
                else:
                    if hmm not in already_added_NCBIG:
                        res['NCBIG'].append([hmm, hmm_label, description, enzyme_ec, go_terms, common_links])
                        already_added_NCBIG.add(hmm)
        return res
