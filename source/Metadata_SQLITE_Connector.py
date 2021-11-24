try:
    from source.utils import *
except:
    from utils import *


class Metadata_SQLITE_Connector():
    def __init__(self,metadata_file):
        self.metadata_file_tsv=metadata_file
        self.db_file=metadata_file.replace('.tsv','')+'.db'
        self.db_headers=self.get_db_headers()
        self.insert_step=5000
        if not file_exists(self.db_file):
            print('Creating SQL database',self.db_file)
            self.create_sql_table()
        self.start_sqlite_cursor()




    def start_sqlite_cursor(self):
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()

    def commit_and_close_sqlite_cursor(self):
        self.sqlite_connection.commit()
        self.sqlite_connection.close()

    def close_sql_connection(self):
        self.sqlite_connection.close()

    def check_all_tables(self):
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        all_tables = self.cursor.fetchall()
        print(all_tables)

    def convert_row_to_sql(self,ref,row_info):
        res=[ref]
        for db in self.db_headers:
            if db in row_info:
                res.append(','.join(row_info[db]))
            else:
                res.append(None)
        return res


    def generate_insert_command(self):
        headers_str=', '.join(self.db_headers)
        headers_str=f'(REF, {headers_str})'.upper()
        n_values=['?' for i in range(len(self.db_headers)+1)]
        n_values_str=', '.join(n_values)
        n_values_str=f'({n_values_str})'
        insert_command = f'INSERT INTO METADATA {headers_str} values {n_values_str}'
        return insert_command

    def generate_fetch_command(self,ref_id):
        headers_str=', '.join(self.db_headers)
        headers_str=f'REF, {headers_str}'.upper()
        fetch_command = f'SELECT {headers_str} FROM METADATA WHERE REF="{ref_id}"'
        return fetch_command

    def yield_metadata(self):
        res=[]
        with open(self.metadata_file_tsv, 'r') as file:
            for line in file:
                row_info={}
                line = line.strip('\n')
                line = line.split('\t')
                current_ref = line[0]
                annotations = line[2:]
                for link in annotations:
                    if link:
                        temp_link = link.split(':')
                        link_type = temp_link[0]
                        link_text = ':'.join(temp_link[1:])
                        link_text=link_text.strip()
                        if link_type not in row_info: row_info[link_type]=set()
                        row_info[link_type].add(link_text)
                        if link_type == 'description' and link_text == 'NA':
                            link_text = None
                        if link_text and link_type == 'description':
                            get_common_links_metadata(link_text, res=row_info)
                res.append(self.convert_row_to_sql(current_ref,row_info))
        return res

    def get_db_headers(self):
        res=set()
        with open(self.metadata_file_tsv, 'r') as file:
            for line in file:
                line = line.strip('\n')
                line = line.split('\t')
                annotations = line[2:]
                for link in annotations:
                    if link:
                        temp_link = link.split(':')
                        link_type = temp_link[0]
                        res.add(link_type)
        return sorted(list(res))


    def create_sql_table(self):
        if os.path.exists(self.db_file):
            os.remove(self.db_file)
        self.start_sqlite_cursor()

        create_table_command = f'CREATE TABLE METADATA (REF TEXT, '
        for header in self.db_headers:
            create_table_command+=f'{header.upper()} TEXT, '
        create_table_command=create_table_command.rstrip(', ')
        create_table_command+=')'
        self.cursor.execute(create_table_command)
        self.sqlite_connection.commit()
        self.store_metadata()

        self.commit_and_close_sqlite_cursor()

    def generate_inserts(self, metadata):
        step=self.insert_step
        for i in range(0, len(metadata), step):
            yield metadata[i:i + step]


    def store_metadata(self):
        insert_command=self.generate_insert_command()
        metadata_yielder=self.yield_metadata()
        generator_insert = self.generate_inserts(metadata_yielder)
        for table_chunk in generator_insert:
            self.cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()

    def convert_sql_to_dict(self,sql_result):
        sql_result=sql_result[1:]
        res={}
        for i in range(len(self.db_headers)):
            db=self.db_headers[i]
            db_res=sql_result[i]
            if db_res:
                db_res=db_res.split(',')
                if db not in res: res[db]=set()
                res[db].update(db_res)
        return res


    def fetch_metadata(self,ref_id):
        fetch_command=self.generate_fetch_command(ref_id)
        res_fetch=self.cursor.execute(fetch_command).fetchone()
        res=self.convert_sql_to_dict(res_fetch)
        return res



if __name__ == '__main__':
    metadata_connector=Metadata_SQLITE_Connector('/home/pedroq/Desktop/test_cr/metadata.tsv')
    res=metadata_connector.fetch_metadata('D0I9C4')
    print(res)