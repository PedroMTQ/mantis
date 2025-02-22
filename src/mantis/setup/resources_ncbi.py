import os
from datetime import datetime

from mantis.src.setup.base import SetupBase
from mantis.src.utils.downloader import download_file
from mantis.src.utils.logger import logger


class SetupResourcesNcbi(SetupBase):
    TRANSLATION_TABLES_URL = 'https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt'
    DATABASE = 'resources'

    def __init__(self, config_file: str) -> None:
        super().__init__(database=self.DATABASE, config_file=config_file)
        self.ncbi_resources = os.path.join(self.output_folder, 'ncbi')


    def write_readme(self):
        with open(self.ncbi_resources + 'readme.md', 'w+') as file:
            datetime_str = str(datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            file.write(f'These files were downloaded on {datetime_str} from:\n{self.TRANSLATION_TABLES_URL}\nThey are used to translate CDS')

    def run(self):
        translation_table_path = os.path.join(self.ncbi_resources, 'gc.prt')
        if os.path.exists(translation_table_path):
            logger.info('Translation tables already exist! Skipping...')
            return
        try:
            os.remove(translation_table_path)
        except Exception:
            pass
        logger.info('Downloading and unzipping NCBI resources')
        download_file(url=self.TRANSLATION_TABLES_URL, output_folder=self.ncbi_resources)
        logger.info(f'Finished downloading {self.DATABASE}')
