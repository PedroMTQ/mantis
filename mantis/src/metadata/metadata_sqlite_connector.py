import os
import sqlite3

from mantis.src.metadata.utils import get_common_links_metadata
from mantis.src.settings import SQLITE_INFO_SPLITTER
from mantis.src.utils.logger import logger


class MetadataSqliteConnector():
    def __init__(self, metadata_file):
        self.metadata_file_tsv = metadata_file
        self.db_file = metadata_file.replace('.tsv', '') + '.db'
        if not os.path.exists(self.metadata_file_tsv):
            logger.info(f'Metadata file missing: {metadata_file}')
            return
        self.insert_step = 5000
        if not os.path.exists(self.db_file):
            logger.info(f'Creating SQL database {self.db_file}')
            self.create_sql_table()
        self.start_sqlite_cursor()

    def start_sqlite_cursor(self):
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()
        self.get_db_headers()

    def commit_and_close_sqlite_cursor(self):
        self.sqlite_connection.commit()
        self.sqlite_connection.close()

    def close_sql_connection(self):
        try:
            self.sqlite_connection.close()
        except Exception:
            return

    def check_all_tables(self):
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        all_tables = self.cursor.fetchall()
        print(all_tables)

    def convert_row_to_sql(self, ref, row_info):
        res = [ref]
        for db in self.db_headers:
            if db in row_info:
                res.append(SQLITE_INFO_SPLITTER.join(row_info[db]))
            else:
                res.append(None)
        return res

    def generate_insert_command(self):
        headers_str = ', '.join(self.db_headers)
        headers_str = f'(REF, {headers_str})'.upper()
        n_values = ['?' for i in range(len(self.db_headers) + 1)]
        n_values_str = ', '.join(n_values)
        n_values_str = f'({n_values_str})'
        insert_command = f'INSERT INTO METADATA {headers_str} values {n_values_str}'
        return insert_command

    def generate_fetch_command(self, ref_id):
        headers_str = ', '.join(self.db_headers)
        headers_str = f'REF, {headers_str}'.upper()
        fetch_command = f'SELECT {headers_str} FROM METADATA WHERE REF="{ref_id}"'
        return fetch_command

    def yield_metadata(self):
        with open(self.metadata_file_tsv, 'r') as file:
            for line in file:
                row_info = {}
                line = line.strip('\n')
                if not line:
                    continue
                line = line.split('\t')
                current_ref = line[0]
                if '|' in line:  line.remove('|')
                annotations = line[1:]
                for link in annotations:
                    if link:
                        temp_link = link.split(':')
                        link_type = temp_link[0]
                        link_text = ':'.join(temp_link[1:])
                        link_text = link_text.strip()
                        if link_type not in row_info: row_info[link_type] = set()
                        row_info[link_type].add(link_text)
                        if link_type == 'description' and link_text == 'NA':
                            link_text = None
                        if link_text and link_type == 'description':
                            get_common_links_metadata(input_string=link_text,
                                                      metadata_dict=row_info)
                yield self.convert_row_to_sql(current_ref, row_info)

    def get_db_headers(self):
        res = set()
        try:
            schema_command = 'PRAGMA table_info(METADATA);'
            res_fetch = self.cursor.execute(schema_command).fetchall()
            res_fetch.pop(0)
            for line in res_fetch:
                link_type = line[1]
                res.add(link_type)
        except Exception:
            with open(self.metadata_file_tsv, 'r') as file:
                for line in file:
                    line = line.strip('\n')
                    line = line.split('\t')
                    annotations = line[2:]
                    for link in annotations:
                        if not link:
                            continue
                        temp_link = link.split(':')
                        link_type = temp_link[0]
                        res.add(link_type)
        self.db_headers = sorted(list(res))

    def create_sql_table(self):
        self.get_db_headers()
        if os.path.exists(self.db_file):
            os.remove(self.db_file)
        self.start_sqlite_cursor()

        create_table_command = 'CREATE TABLE METADATA (REF TEXT, '
        for header in self.db_headers:
            create_table_command += f'{header.upper()} TEXT, '
        create_table_command = create_table_command.rstrip(', ')
        create_table_command += ')'
        self.cursor.execute(create_table_command)
        self.sqlite_connection.commit()
        create_index_command = 'CREATE INDEX REF_IDX ON METADATA (REF)'
        self.cursor.execute(create_index_command)

        self.store_metadata()

        self.commit_and_close_sqlite_cursor()

    def generate_inserts(self, metadata_yielder):
        step = self.insert_step
        temp = []
        for i in metadata_yielder:
            temp.append(i)
            if len(temp) >= step:
                yield temp
                temp = []
        yield temp

    def store_metadata(self):
        insert_command = self.generate_insert_command()
        metadata_yielder = self.yield_metadata()
        if metadata_yielder:
            generator_insert = self.generate_inserts(metadata_yielder)
            for table_chunk in generator_insert:
                if table_chunk:
                    self.cursor.executemany(insert_command, table_chunk)
            self.sqlite_connection.commit()

    def convert_sql_to_dict(self, sql_result):
        sql_result = sql_result[1:]
        res = {}
        for i in range(len(self.db_headers)):
            db = self.db_headers[i].lower()
            db_res = sql_result[i]
            if db_res:
                db_res = db_res.split(SQLITE_INFO_SPLITTER)
                if db not in res: res[db] = set()
                res[db].update(db_res)
        return res

    def fetch_metadata(self, ref_id):
        if not os.path.exists(self.db_file):
            return {}
        fetch_command = self.generate_fetch_command(ref_id)
        try:
            res_fetch = self.cursor.execute(fetch_command).fetchone()
            res = self.convert_sql_to_dict(res_fetch)
            return res
        except Exception:
            logger.error(f'Failed retrieving {ref_id} in {self.db_file}')
            return {}

    def test_database(self):
        res = set()
        if not os.path.exists(self.metadata_file_tsv): return res
        with open(self.metadata_file_tsv) as file:
            for line in file:
                ref = line.split('\t')[0]
                try:
                    self.fetch_metadata(ref)
                except Exception:
                    logger.error(f'Failed retrieving {ref} in {self.db_file}')
                    res.add(ref)
        return res

