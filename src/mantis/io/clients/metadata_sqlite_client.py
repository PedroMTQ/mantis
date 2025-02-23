import os
import sqlite3
from pathlib import Path
from typing import Generator


from mantis.settings import DATA, BATCH_SIZE
from mantis.io.logger import logger
from mantis.core.metadata.metadata_document import get_extended_metadata_document, BaseMetadataDocument



class MetadataSqliteClient():
    def __init__(self,
                 metadata_file: str,
                 data_folder: str=DATA):
        self.metadata_file = metadata_file
        if not os.path.exists(self.metadata_file):
            logger.info(f'Metadata file missing: {metadata_file}')
            return
        self.db_file = os.path.join(data_folder, f'{Path(self.metadata_file).stem}.db')
        if os.path.exists(self.db_file):
            self.sqlite_connection = sqlite3.connect(self.db_file)
            self.metadata_types = self.get_metadata_types_database()
            self.set_extended_metadata_document()
        else:
            self.metadata_types = self.get_metadata_types_file()
            self.set_extended_metadata_document()
            self.sqlite_connection = sqlite3.connect(self.db_file)
            logger.info(f'Creating SQL database {self.db_file}')
            self.create_sql_table()
            self.insert_metadata()

    def set_extended_metadata_document(self):
        self.document_factory: BaseMetadataDocument = get_extended_metadata_document(metadata_types=self.metadata_types)

    def generate_insert_command(self):
        headers_str = ', '.join(self.metadata_types)
        headers_str = f'(REFERENCE, {headers_str})'.upper()
        n_values = ['?' for i in range(len(self.metadata_types) + 1)]
        n_values_str = ', '.join(n_values)
        n_values_str = f'({n_values_str})'
        insert_command = f'INSERT INTO METADATA {headers_str} values {n_values_str}'
        return insert_command


    def yield_metadata_file(self):
        for line in open(self.metadata_file, 'r'):
            line = line.strip('\n')
            if not line:
                continue
            line = line.split('\t')
            reference_id = line[0]
            if '|' in line:
                line.remove('|')
            annotations = line[1:]
            metadata_document: BaseMetadataDocument = self.document_factory()
            for annotation in annotations:
                if not annotation:
                    continue
                annotation = annotation.split(':')
                metadata_type = annotation[0]
                metadata = ':'.join(annotation[1:])
                metadata = metadata.strip()
                metadata_document.add(metadata_type=metadata_type,
                                      metadata=metadata)
                if metadata_type == 'description' and (metadata != 'NA' and metadata is not None):
                    metadata_document.add_by_search(metadata)
            yield [reference_id] + metadata_document.to_sql(metadata_types=self.metadata_types)

    def get_metadata_types_database(self):
        res = set()
        cursor = self.sqlite_connection.cursor()
        schema_command = 'PRAGMA table_info(METADATA);'
        res_fetch = cursor.execute(schema_command).fetchall()
        res_fetch.pop(0)
        for line in res_fetch:
            metadata_type = line[1]
            res.add(metadata_type.lower())
        cursor.close()
        return res

    def get_metadata_types_file(self):
        res = set()
        for line in open(self.metadata_file, 'r'):
            line = line.strip('\n')
            line = line.split('\t')
            annotations = line[2:]
            for link in annotations:
                if not link:
                    continue
                temp_link = link.split(':')
                metadata_type = temp_link[0]
                res.add(metadata_type.lower())
        return res

    def create_sql_table(self):
        cursor = self.sqlite_connection.cursor()
        create_table_command = 'CREATE TABLE METADATA (REFERENCE TEXT, '
        for header in self.metadata_types:
            create_table_command += f'{header.upper()} JSON, '
        create_table_command = create_table_command.rstrip(', ')
        create_table_command += ')'
        cursor.execute(create_table_command)
        self.sqlite_connection.commit()
        create_index_command = 'CREATE INDEX REFERENCE_IDX ON METADATA (REFERENCE)'
        cursor.execute(create_index_command)
        self.sqlite_connection.commit()
        cursor.close()

    @staticmethod
    def generate_inserts(metadata_yielder):
        temp = []
        for i in metadata_yielder:
            temp.append(i)
            if len(temp) >= BATCH_SIZE:
                yield temp
                temp = []
        yield temp

    def insert_metadata(self):
        insert_command = self.generate_insert_command()
        metadata_yielder = self.yield_metadata_file()
        cursor = self.sqlite_connection.cursor()
        generator_insert = self.generate_inserts(metadata_yielder)
        for table_chunk in generator_insert:
            if table_chunk:
                cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()
        cursor.close()

    def get_metadata(self, reference_id: str) -> BaseMetadataDocument:
        if not os.path.exists(self.db_file):
            raise FileNotFoundError(f'Metadata database {self.db_file} does not exist!')
        headers_str = ', '.join(self.metadata_types)
        headers_str = f'{headers_str}'.upper()
        fetch_command = f'''
        SELECT {headers_str} FROM METADATA WHERE REFERENCE="{reference_id}"
        '''
        cursor = self.sqlite_connection.cursor()
        try:
            metadata_list = cursor.execute(fetch_command).fetchone()
            metadata_document = self.document_factory.from_sql(metadata_types=self.metadata_types,
                                                               metadata_list=metadata_list)
        except Exception as e:
            logger.exception(e)
            logger.error(f'Failed retrieving {reference_id} in metadata database {self.db_file}')
            metadata_document = BaseMetadataDocument()
        cursor.close()
        return metadata_document

    def get_metadata_batch(self, reference_ids: list[str]) -> dict[str, BaseMetadataDocument]:
        if not os.path.exists(self.db_file):
            raise FileNotFoundError(f'Metadata database {self.db_file} does not exist!')
        placeholders = ','.join('?' for _ in reference_ids)
        headers_str = ', '.join(self.metadata_types)
        headers_str = f'{headers_str}'.upper()
        fetch_command = f'''
        SELECT REFERENCE, {headers_str} FROM METADATA WHERE REFERENCE IN ({placeholders})
        '''
        res = {reference_id: self.document_factory() for reference_id in reference_ids}
        cursor = self.sqlite_connection.cursor()
        for metadata_list in cursor.execute(fetch_command, reference_ids):
            reference_id = metadata_list[0]
            metadata_list = metadata_list[1:]
            metadata_document = self.document_factory.from_sql(metadata_types=self.metadata_types,
                                                               metadata_list=metadata_list)
            res[reference_id] = metadata_document
        cursor.close()
        return res

    def is_valid(self):
        if not os.path.exists(self.db_file):
            raise FileNotFoundError(f'Metadata database {self.db_file} does not exist!')
        res = set()
        with open(self.metadata_file) as file:
            for line in file:
                reference_id = line.split('\t')[0]
                try:
                    self.get_metadata(reference_id=reference_id)
                except Exception:
                    logger.error(f'Failed retrieving {reference_id} in {self.db_file}')
                    res.add(reference_id)
        return res

if __name__ == '__main__':
    client = MetadataSqliteClient(metadata_file='tests/test_hmm/test1/metadata.tsv')
    doc = client.get_metadata('K00036')
    print(doc)
    docs: dict = client.get_metadata_batch(reference_ids=['K00036', 'K00003'])
    print(docs)