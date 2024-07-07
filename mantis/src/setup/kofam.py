import os
from pathlib import Path

from mantis.src.metadata.utils import get_common_links_metadata
from mantis.src.utils import print_cyan, remove_file


class SetupKofam():

    #####################   KOFAM

    def compile_kofam_metadata(self):
        metadata_to_write = {}
        self.get_link_kofam_ko_list(metadata_to_write)
        self.get_link_kofam_ko_to_binary(metadata_to_write, target_file='ko2cog.xl')
        self.get_link_kofam_ko_to_binary(metadata_to_write, target_file='ko2go.xl')
        self.get_link_kofam_ko_to_binary(metadata_to_write, target_file='ko2tc.xl')
        self.get_link_kofam_ko_to_binary(metadata_to_write, target_file='ko2cazy.xl')
        self.write_metadata(metadata_to_write, self.mantis_paths['kofam'] + 'metadata.tsv')
        remove_file(self.mantis_paths['kofam'] + 'ko_list')
        remove_file(self.mantis_paths['kofam'] + 'ko2cazy.xl')
        remove_file(self.mantis_paths['kofam'] + 'ko2cog.xl')
        remove_file(self.mantis_paths['kofam'] + 'ko2go.xl')
        remove_file(self.mantis_paths['kofam'] + 'ko2tc.xl')

    def get_link_kofam_ko_list(self, res):
        file_path = self.mantis_paths['kofam'] + 'ko_list'
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                ko, description = line[0], line[-1]
                if ko not in res: res[ko] = {}
                if '[EC:' in description:
                    description, temp_links = description.split('[EC:')
                else:
                    temp_links = description
                get_common_links_metadata(input_string=temp_links,
                                          metadata_dict=res[ko])
                if 'kegg_ko' not in res[ko]: res[ko]['kegg_ko'] = set()
                res[ko]['kegg_ko'].add(ko)
                if 'description' not in res[ko]:
                    res[ko]['description'] = set()
                res[ko]['description'].add(description)
                line = file.readline()

    def get_link_kofam_ko_to_binary(self, res, target_file):
        file_path = self.mantis_paths['kofam'] + target_file
        if 'ko2tc' in target_file:
            target_link = 'tcdb'
        else:
            target_link = target_file.replace('ko2', '').replace('.xl', '')
        with open(file_path) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                ko, link = line
                link = link.strip('[]').split(':')[1].split()
                if ko not in res: res[ko] = {}
                if target_link not in res[ko]: res[ko][target_link] = set()
                res[ko][target_link].update(link)
                line = file.readline()

    def download_kofam(self, stdout_file=None):
        Path(self.mantis_paths['kofam']).mkdir(parents=True, exist_ok=True)
        if self.check_reference_exists('kofam') and \
                os.path.exists(self.mantis_paths['kofam'] + 'metadata.tsv'):
            print('KOfam HMM already exists! Skipping...', flush=True, file=stdout_file)
            return
        kofam_hmm = 'https://www.genome.jp/ftp/db/kofam/profiles.tar.gz'
        ko_list = 'https://www.genome.jp/ftp/db/kofam/ko_list.gz'
        ko_to_cog = 'https://www.kegg.jp/kegg/files/ko2cog.xl'
        ko_to_go = 'https://www.kegg.jp/kegg/files/ko2go.xl'
        ko_to_tc = 'https://www.kegg.jp/kegg/files/ko2tc.xl'
        ko_to_cazy = 'https://www.kegg.jp/kegg/files/ko2cazy.xl'
        with open(self.mantis_paths['kofam'] + 'readme.md', 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(
                f'This hmm was downloaded on {datetime_str} from:\n{kofam_hmm}\nMetadata was downloaded from:\n{ko_list}\n{ko_to_cog}\n{ko_to_go}\n{ko_to_tc}\n{ko_to_cazy}')
        print_cyan('Downloading and unzipping KOfam hmms ', flush=True, file=stdout_file)
        for url in [kofam_hmm, ko_list, ko_to_cog, ko_to_go, ko_to_tc, ko_to_cazy]:
            download_file(url, output_folder=self.mantis_paths['kofam'], stdout_file=stdout_file)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'profiles.tar.gz',
                           extract_path=self.mantis_paths['kofam'], stdout_file=stdout_file, remove_source=True)
        uncompress_archive(source_filepath=self.mantis_paths['kofam'] + 'ko_list.gz', stdout_file=stdout_file,
                           remove_source=True)
        merge_profiles(self.mantis_paths['kofam'] + 'profiles/', self.mantis_paths['kofam'] + 'kofam_merged.hmm',
                       stdout_file=stdout_file)
        run_command('hmmpress ' + self.mantis_paths['kofam'] + 'kofam_merged.hmm', stdout_file=stdout_file)
        self.compile_kofam_metadata()
