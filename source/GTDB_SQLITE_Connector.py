try:
    from source.utils import *
except:
    from utils import *


class GTDB_SQLITE_Connector():
    '''
    this just creates an sql database from two gtdb files to convert gtdb to ncbi. first we download them and create the db
    then anytime we need to fetch info we just open the db, fetch the info, and close the connection
    '''


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

    def process_gtdb_taxonomy(self,gtdb_lineage):
        res = None
        temp=gtdb_lineage.split(';')
        temp = ['_'.join(i.split('_')[2:]) for i in temp]
        while temp and not res:
            res = temp.pop(-1).strip()
        return res

    def read_gtdb_tsv(self,gtdb_tsv):
        res=[]
        with open(gtdb_tsv) as file:
            file.readline()
            for line in file:
                line = line.strip('\n')
                line = line.split('\t')
                gtdb_taxonomy = line[16]
                most_resolved = self.process_gtdb_taxonomy(gtdb_taxonomy)
                ncbi_id = line[77]
                res.append([most_resolved, ncbi_id])
        return res

    def download_data(self):
        url = 'https://data.gtdb.ecogenomic.org/releases/latest/'
        ar_url = f'{url}ar122_metadata.tar.gz'
        bac_url = f'{url}bac120_metadata.tar.gz'

        ar_file = f'{self.gtdb_folder}ar122_metadata.tar.gz'
        bac_file = f'{self.gtdb_folder}bac120_metadata.tar.gz'

        try:
            ar_file_unc = [i for i in os.listdir(self.gtdb_folder) if i.startswith('ar122') and i.endswith('.tsv')][0]
            ar_file_unc = f'{self.gtdb_folder}{ar_file_unc}'
        except:
            ar_file_unc = None
        try:
            bac_file_unc = [i for i in os.listdir(self.gtdb_folder) if i.startswith('bac120') and i.endswith('.tsv')][0]
            bac_file_unc = f'{self.gtdb_folder}{bac_file_unc}'
        except:
            bac_file_unc = None

        if file_exists(ar_file) or file_exists(ar_file_unc):
            pass
        else:
            download_file(ar_url, output_folder=self.gtdb_folder)

        if file_exists(bac_file) or file_exists(bac_file_unc):
            pass
        else:
            download_file(bac_url, output_folder=self.gtdb_folder)

        if file_exists(ar_file):         uncompress_archive(ar_file, remove_source=True)
        if file_exists(bac_file):        uncompress_archive(bac_file, remove_source=True)
        self.bac_file = [i for i in os.listdir(self.gtdb_folder) if i.startswith('bac120') and i.endswith('.tsv')][0]
        self.ar_file = [i for i in os.listdir(self.gtdb_folder) if i.startswith('ar122') and i.endswith('.tsv')][0]
        self.bac_file=f'{self.gtdb_folder}{self.bac_file}'
        self.ar_file=f'{self.gtdb_folder}{self.ar_file}'

    def create_sql_table(self):
        self.download_data()
        #this will probably need to be changed to an output_folder provided by the user
        if os.path.exists(self.db_file):
            os.remove(self.db_file)
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()

        create_table_command = f'CREATE TABLE GTDB2NCBI (' \
                            f'GTDB TEXT,' \
                            f'NCBI  INTEGER )'
        self.cursor.execute(create_table_command)
        self.sqlite_connection.commit()
        self.store_gtdb2ncbi()
        if file_exists(self.bac_file):
            os.remove(self.bac_file)
        if file_exists(self.ar_file):
            os.remove(self.ar_file)


    def generate_inserts(self, chebi2others):
        step=self.insert_step
        for i in range(0, len(chebi2others), step):
            yield chebi2others[i:i + step]


    def store_gtdb2ncbi(self):

        gtdb2ncbi=self.read_gtdb_tsv(self.bac_file)+self.read_gtdb_tsv(self.ar_file)
        generator_insert = self.generate_inserts(gtdb2ncbi)
        for table_chunk in generator_insert:
            insert_command = f'INSERT INTO GTDB2NCBI (GTDB, NCBI) values (?,?)'
            self.cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()

    def fetch_ncbi_id(self,gtdb_lineage):
        res=set()
        gtdb_id=self.process_gtdb_taxonomy(gtdb_lineage)
        fetch_command = f'SELECT GTDB,NCBI FROM GTDB2NCBI WHERE GTDB = "{gtdb_id}"'
        res_fetch=self.cursor.execute(fetch_command).fetchall()
        for i in res_fetch:
            res.add(i[1])
        return res

    def fetch_gtdb_id(self,ncbi_id):
        res=set()
        fetch_command = f"SELECT GTDB,NCBI FROM GTDB2NCBI WHERE NCBI = {ncbi_id}"
        res_fetch=self.cursor.execute(fetch_command).fetchall()
        for i in res_fetch:
            res.add(i[0])
        return res



class GTDB_SQLITE_Connector(GTDB_SQLITE_Connector):
    def launch_gtdb_connector(self,resources_folder):

        self.gtdb_folder=resources_folder
        Path(self.gtdb_folder).mkdir(parents=True, exist_ok=True)
        self.db_file = f'{self.gtdb_folder}gtdb_to_ncbi.db'
        self.insert_step=50000
        if os.path.exists(self.db_file):
            self.start_sqlite_cursor()
        else:
            self.create_sql_table()


    def get_organism_lineage(self, taxon_id, stdout_file=None):
        ncbi_resources=add_slash(self.mantis_paths['resources']+'NCBI')

        lineage_file_path = f'{ncbi_resources}taxidlineage.dmp'
        try:
            lineage_file = open(lineage_file_path, 'r')
        except:
            print_cyan('Lineage dump is not present! If you\'d like to run taxonomic lineage annotation, please run < setup_databases >',flush=True, file=stdout_file)
            return []
        line = lineage_file.readline().strip('\n').replace('|', '')
        while line:
            line = line.split()
            if str(taxon_id) == str(line[0]):
                lineage = line[1:]
                lineage.append(taxon_id)
                lineage_file.close()
                return lineage
            line = lineage_file.readline().strip('\n').replace('|', '')
        lineage_file.close()
        return []

    def get_taxa_ncbi_url(self,url):
        webpage = None
        c = 0
        while not webpage and c <= 10:
            req = requests.get(url)
            try:
                webpage = req.text
                taxa_id = re.search('<Id>\d+</Id>', webpage)
                return re.search('\d+', taxa_id.group()).group()
            except:
                print('Could not get a response from NCBI, trying again')
                c += 1

    def get_taxa_ncbi(self,organism_name):
        taxa_id=self.fetch_ncbi_id(organism_name)
        if taxa_id: return str(list(taxa_id)[0])
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={organism_name}'
        taxa_id = self.get_taxa_ncbi_url(url)
        if not taxa_id:
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=candidatus+{organism_name}'
            taxa_id = self.get_taxa_ncbi_url(url)
        if not taxa_id:    print(f'Could not find taxa ID for {organism_name}')
        return str(taxa_id)


if __name__ == '__main__':
    gtdb_connector=GTDB_SQLITE_Connector()
    gtdb_connector.launch_gtdb_connector(resources_folder='/home/pedroq/Desktop/test_cr/')
    a=gtdb_connector.fetch_ncbi_id('d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanolobus;s__Methanolobus psychrophilus')
    print(a)
    a=gtdb_connector.fetch_ncbi_id('s__Methanolobus psychrophilus')
    print(a)
    a=gtdb_connector.fetch_gtdb_id('1094980')
    print(a)
