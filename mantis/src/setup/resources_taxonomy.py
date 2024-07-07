import os

from mantis.src.setup.base import SetupBase
from mantis.src.taxonomy.taxonomy_sqlite_connector import TaxonomySqliteConnector
from mantis.src.utils.logger import logger


class SetupResourcesTaxonomy(SetupBase):
    def __init__(self) -> None:
        super().__init__(database='resources')

    def run(self):
        taxonomy_connector = TaxonomySqliteConnector()
        db_path = os.path.join(self.output_folder, 'taxonomy.db')
        if not os.path.exists(db_path):
            taxonomy_connector.launch_taxonomy_connector()
            taxonomy_connector.create_taxonomy_db()
            taxonomy_connector.close_taxonomy_connection()
        logger.info(f'Finished downloading {self.DATABASE}')
