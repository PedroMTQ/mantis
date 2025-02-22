import os
from mantis.src.utils.config_reader import ConfigReader
from mantis.src.utils.logger import logger
from pathlib import Path

class SetupBase():
    def __init__(self, datatabase: str, config_file: str) -> None:
        # TODO finish config reading
        self.config = ConfigReader(config_file=config_file)
        self.output_folder = self.config.get_folder(datatabase))
        Path(self.output_folder).mkdir(parents=True, exist_ok=True)


    def write_metadata(self, metadata, metadata_file):
        with open(metadata_file, 'w+') as file:
            for seq in metadata:
                link_line = f'{seq}\t|'
                for link_type in metadata[seq]:
                    for inner_link in metadata[seq][link_type]:
                        if inner_link:
                            link_line += f'\t{link_type}:{inner_link}'
                link_line += '\n'
                file.write(link_line)


    def check_reference_exists(self, database, taxon_id=None):
        ncbi_resources = os.path.join(self.mantis_paths['resources'], 'NCBI')

        if database == 'ncbi_res':
            if os.path.exists(ncbi_resources + 'gc.prt') and \
                    os.path.exists(ncbi_resources + 'gc.prt'):
                return True
        elif database == 'taxonomy':
            taxonomy_db = self.mantis_paths['resources'] + 'taxonomy.db'
            if os.path.exists(taxonomy_db):
                return True
        elif database == 'NOGSQL':
            if os.path.exists(self.mantis_paths['NOG'] + 'eggnog.db'):
                return True
        elif database == 'tcdb':
            if os.path.exists(self.mantis_paths['tcdb'] + 'tcdb.dmnd') and  os.path.exists(self.mantis_paths['tcdb'] + 'metadata.tsv'):
                return True
        elif database == 'NOG_DMND':
            if os.path.exists(self.mantis_paths['NOG'] + 'eggnog_proteins.dmnd'):
                return True
        target_file = self.get_path_default_ref(database, taxon_id)
        if target_file:
            if target_file.endswith('.dmnd'):
                if not os.path.exists(target_file):
                    return False
            else:
                for extension in ['', '.h3f', '.h3i', '.h3m', '.h3p']:
                    if not os.path.exists(target_file + extension):
                        return False
        else:
            return False
        return True
