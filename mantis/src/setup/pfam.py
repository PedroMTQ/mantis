
import os
import re
from datetime import datetime

from mantis.src.metadata.utils import get_common_links_metadata
from mantis.src.setup.base import SetupBase
from mantis.src.utils.downloader import download_file
from mantis.src.utils.file_processer import remove_file, uncompress_archive
from mantis.src.utils.logger import logger
from mantis.src.utils.utils import run_command


class c(SetupBase):
    HMM_URL = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
    METADATA_URL = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
    PFAM2GO_URL= 'http://current.geneontology.org/ontology/external2go/pfam2go'
    DATABASE = 'pfam'

    def __init__(self, config_file: str) -> None:
        super().__init__(database=self.DATABASE, config_file=config_file)
        self.pfam_metadata_path = os.path.join(self.output_folder, 'metadata.tsv')
        self.pfam_hmm_path = os.path.join(self.output_folder, 'Pfam-A.hmm.dat')
        self.pfam2go_path = os.path.join(self.output_folder, 'pfam2go')


    def read_pfam2go(self):
        res = {}
        with open(self.pfam2go_path) as pfam2go_file:
            line = pfam2go_file.readline()
            while line:
                line = line.strip('\n')
                if '!' not in line[0]:
                    line = line.split('>')
                    line = [i.strip() for i in line]
                    pfam_id = line[0].split()[0].replace('Pfam:', '')
                    go_annots = line[1].split(';')
                    go_description = [i.replace('GO:', '').strip() for i in go_annots if not re.search(r'GO:\d{3,}', i)]
                    go_ids = [i.replace('GO:', '').strip() for i in go_annots if re.search(r'GO:\d{3,}', i)]
                    if pfam_id not in res:
                        res[pfam_id] = {'go': set(go_ids), 'description': set(go_description)}
                    else:
                        res[pfam_id]['go'].update(go_ids)
                        res[pfam_id]['description'].update(go_description)
                line = pfam2go_file.readline()
        return res

    def is_good_description(self, hmm, row_description):
        temp = [i.lower() for i in row_description.split()]
        if hmm.lower() in temp and 'protein' in temp and len(temp) == 2:
            return False
        if re.search(' [uU]nknown [Ff]unction', row_description): return False
        return True

    def build_pfam_line(self, hmm, metadata):
        link_line = []
        pfam_ids = set(metadata['pfam'])
        metadata['pfam'] = set()
        for p_id in pfam_ids:
            metadata['pfam'].add(p_id.split('.')[0])
        for link_type in metadata:
            for inner_link in metadata[link_type]:
                link_line.append(f'{link_type}:{inner_link}')
        link_line = '\t'.join(link_line)
        return f'{hmm}\t|{link_line}\n'

    def get_hmm_info_pfam(self, line, file, pfam2go):
        if line.startswith('#=GF ID'):
            hmm = line.replace('#=GF ID', '').strip('\n').strip()
        if hmm:
            line = file.readline()
            if line.startswith('#=GF AC'):
                pfam_accession = line.replace('#=GF AC', '').strip('\n').strip().split('.')[0]
            line = file.readline()
            if line.startswith('#=GF DE'):
                hmm_description = line.replace('#=GF DE', '').strip('\n').strip()

            if pfam_accession in pfam2go:
                current_metadata = pfam2go[pfam_accession]
            else:
                current_metadata = {'description': set()}

            current_metadata['pfam'] = set()
            current_metadata['pfam'].add(pfam_accession)
            current_metadata['pfam'].add(hmm)
            if self.is_good_description(hmm, hmm_description):
                current_metadata['description'].add(hmm_description)
            get_common_links_metadata(input_string=hmm_description,
                                      metadata_dict=current_metadata)
            metadata_line = self.build_pfam_line(hmm, current_metadata)
            return metadata_line

    def compile_pfam_metadata(self):
        pfam2go = self.read_pfam2go()

        with open(self.pfam_metadata_path, 'w+') as pfam_metadata_file:
            with open(self.pfam_hmm_path) as pfam_dat_file:
                for line in pfam_dat_file:
                    line = line.strip('\n')
                    if line.startswith('#=GF ID'):
                        metadata_line = self.get_hmm_info_pfam(line, pfam_dat_file, pfam2go)
                        pfam_metadata_file.write(metadata_line)
        remove_file(self.pfam_hmm_path)
        pfam2go_path = os.path.join(self.output_folder, 'pfam2go')
        remove_file(pfam2go_path)

    def write_readme(self):
        with open(os.path.join(self.output_folder, 'readme.md'), 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'This hmm was downloaded on {datetime_str} from:\n{self.HMM_URL}\nMetadata was downloaded from:\n{self.METADATA_URL}\nPfam2GO was downloaded from:\n{self.PFAM2GO_URL}')

    def run(self):
        if self.check_reference_exists('pfam') and \
                os.path.exists(self.metadata_path):
            logger.info('Pfam hmm already exists! Skipping...')
            return

        logger.info('Downloading and unzipping Pfam hmms ')
        to_download = []
        to_unzip = []
        if not os.path.exists(self.hmm_path):
            to_download.append(self.HMM_URL)
            to_unzip.append('Pfam-A.hmm.gz')
        if not os.path.exists(self.metadata_path):
            to_unzip.append('Pfam-A.hmm.dat.gz')
            to_download.append(self.METADATA_URL)
            to_download.append(self.PFAM2GO_URL)
        for url in to_download:
            download_file(url=url, output_folder=self.output_folder)
        self.write_readme()
        for file_to_unzip in to_unzip:
            uncompress_archive(source_filepath=os.path.join(self.output_folder, file_to_unzip),
                               remove_source=True)
        if not self.check_reference_exists('pfam'):
            run_command(f'hmmpress {self.hmm_path}')
        self.compile_pfam_metadata()
