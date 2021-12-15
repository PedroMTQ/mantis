try:
    from source.utils import *
except:
    from utils import *




class Taxonomy_SQLITE_Connector():

    def __init__(self,resources_folder):
        self.info_splitter='##'
        self.insert_step=50000
        self.resources_folder=resources_folder
        self.taxonomy_db_file = f'{self.resources_folder}Taxonomy.db'
        Path(self.resources_folder).mkdir(parents=True, exist_ok=True)



    '''
    this just creates an sql database from two gtdb files to convert gtdb to ncbi. first we download them and create the db
    then anytime we need to fetch info we just open the db, fetch the info, and close the connection
    '''
    def start_taxonomy_connection(self):
        self.sqlite_connection = sqlite3.connect(self.taxonomy_db_file)
        self.cursor = self.sqlite_connection.cursor()

    def commit_and_close_sqlite_cursor(self):
        self.sqlite_connection.commit()
        self.sqlite_connection.close()

    def close_taxonomy_connection(self):
        self.sqlite_connection.close()

    def check_table(self):
        self.cursor.execute("SELECT * FROM GTDB2NCBI limit 10")
        res_fetch = self.cursor.fetchall()
        print(res_fetch)

    def process_gtdb_taxonomy(self,gtdb_lineage):
        res = None
        temp=gtdb_lineage.split(';')
        if len(temp)>1:
            temp = ['_'.join(i.split('_')[2:]) for i in temp]
            for i in range(len(temp)):
                current_str=temp[i].split('_')
                current_str=[i for i in current_str if len(i)>1]
                temp[i]='_'.join(current_str)

        lineage_str=self.info_splitter.join(temp)
        while temp and not res:
            res = temp.pop(-1).strip()
        return res,lineage_str

    def read_gtdb_tsv(self,gtdb_tsv):
        #to avoid duplicate entries in sql
        already_added=set()
        with open(gtdb_tsv) as file:
            file.readline()
            for line in file:
                line = line.strip('\n')
                line = line.split('\t')
                gtdb_taxonomy = line[16]
                most_resolved,lineage_str = self.process_gtdb_taxonomy(gtdb_taxonomy)
                ncbi_id = line[77]
                yield_str=self.info_splitter.join([most_resolved,ncbi_id])
                if yield_str not in already_added:
                    yield most_resolved,lineage_str,ncbi_id
                already_added.add(yield_str)

    def download_data(self):
        url = 'https://data.gtdb.ecogenomic.org/releases/latest/'
        ar_url = f'{url}ar122_metadata.tar.gz'
        bac_url = f'{url}bac120_metadata.tar.gz'
        ar_file = f'{self.temp_folder}ar122_metadata.tar.gz'
        bac_file = f'{self.temp_folder}bac120_metadata.tar.gz'
        taxonomy_file = f'{self.temp_folder}new_taxdump.tar.gz'

        taxonomy_url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz'
        download_file(taxonomy_url, output_folder=self.temp_folder)

        try:
            ar_file_unc = [i for i in os.listdir(self.temp_folder) if i.startswith('ar122') and i.endswith('.tsv')][0]
            ar_file_unc = f'{self.temp_folder}{ar_file_unc}'
        except:
            ar_file_unc = None
        try:
            bac_file_unc = [i for i in os.listdir(self.temp_folder) if i.startswith('bac120') and i.endswith('.tsv')][0]
            bac_file_unc = f'{self.temp_folder}{bac_file_unc}'
        except:
            bac_file_unc = None

        if file_exists(ar_file) or file_exists(ar_file_unc):
            pass
        else:
            download_file(ar_url, output_folder=self.temp_folder)

        if file_exists(bac_file) or file_exists(bac_file_unc):
            pass
        else:
            download_file(bac_url, output_folder=self.temp_folder)

        if file_exists(ar_file):         uncompress_archive(ar_file, remove_source=True)
        if file_exists(bac_file):        uncompress_archive(bac_file, remove_source=True)
        if file_exists(taxonomy_file):        uncompress_archive(taxonomy_file, remove_source=True)

        self.bac_file = [i for i in os.listdir(self.temp_folder) if i.startswith('bac120') and i.endswith('.tsv')][0]
        self.ar_file = [i for i in os.listdir(self.temp_folder) if i.startswith('ar122') and i.endswith('.tsv')][0]
        self.bac_file=f'{self.temp_folder}{self.bac_file}'
        self.ar_file=f'{self.temp_folder}{self.ar_file}'



    def create_taxonomy_db(self):
        if not os.path.exists(self.taxonomy_db_file):
            self.temp_folder = f'{self.resources_folder}Taxonomy_temp{SPLITTER}'

            Path(self.temp_folder).mkdir(parents=True, exist_ok=True)
            self.download_data()
            #this will probably need to be changed to an output_folder provided by the user
            if os.path.exists(self.taxonomy_db_file):
                os.remove(self.taxonomy_db_file)
            self.sqlite_connection = sqlite3.connect(self.taxonomy_db_file)
            self.cursor = self.sqlite_connection.cursor()

            create_table_command = f'CREATE TABLE GTDB2NCBI (' \
                                f'GTDB TEXT,' \
                                f'GTDBLINEAGE TEXT,' \
                                f'NCBI  INTEGER )'
            self.cursor.execute(create_table_command)
            create_index_command = f'CREATE INDEX GTDB2NCBI_IDX ON GTDB2NCBI (GTDB,GTDBLINEAGE,NCBI)'
            self.cursor.execute(create_index_command)


            create_table_command = f'CREATE TABLE NCBILINEAGE (' \
                                f'NCBI INTEGER,' \
                                f'LINEAGE  TEXT )'
            self.cursor.execute(create_table_command)
            create_index_command = f'CREATE INDEX NCBILINEAGEI_IDX ON NCBILINEAGE (NCBI)'
            self.cursor.execute(create_index_command)
            self.sqlite_connection.commit()
            self.store_gtdb2ncbi()
            self.store_ncbi_lineage()
            shutil.rmtree(self.temp_folder)





    def generate_inserts(self, input_generator):
        step=self.insert_step
        temp=[]
        for i in input_generator:
            if len(temp)<step:
                temp.append(i)
            elif len(temp)==step:
                yield temp
                temp=[]
        yield temp

    def chain_generators(self):
        for i in self.read_gtdb_tsv(self.bac_file):
            yield i
        for i in self.read_gtdb_tsv(self.ar_file):
            yield i

    def store_gtdb2ncbi(self):
        gtdb2ncbi=self.chain_generators()
        generator_insert = self.generate_inserts(gtdb2ncbi)
        for table_chunk in generator_insert:
            insert_command = f'INSERT INTO GTDB2NCBI (GTDB,GTDBLINEAGE, NCBI) values (?,?,?)'
            self.cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()

    def read_ncbi_lineage(self):
        with open(f'{self.temp_folder}taxidlineage.dmp') as file:
            for line in file:
                line = line.strip('\n').replace('|', '')
                line = line.split()
                taxon_id=line[0]
                lineage = line[1:]
                yield taxon_id,','.join(lineage)



    def store_ncbi_lineage(self):
        ncbi_lineage=self.read_ncbi_lineage()
        generator_insert = self.generate_inserts(ncbi_lineage)
        for table_chunk in generator_insert:
            insert_command = f'INSERT INTO NCBILINEAGE (NCBI, LINEAGE) values (?,?)'
            self.cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()


    def fetch_ncbi_id(self,gtdb_lineage):
        res=set()
        gtdb_id,_=self.process_gtdb_taxonomy(gtdb_lineage)
        fetch_command = f'SELECT GTDB,NCBI FROM GTDB2NCBI WHERE GTDB = "{gtdb_id}"'
        res_fetch=self.cursor.execute(fetch_command).fetchall()
        for i in res_fetch:
            res.add(i[1])
        if len(res)==1: res=list(res)[0]
        elif len(res)==0: res=None
        elif len(res)>1: res= self.get_lowest_common_ancestor_ncbi(res)
        return res

    def fetch_gtdb_id(self,ncbi_id):
        res=set()
        fetch_command = f'SELECT GTDB,GTDBLINEAGE,NCBI FROM GTDB2NCBI WHERE NCBI = "{ncbi_id}"'
        res_fetch=self.cursor.execute(fetch_command).fetchall()
        for i in res_fetch:
            res.add(i)
        if len(res)==1: res=list(res)[0][1]
        elif len(res)==0: res=None
        #if the ncbi matches with more than one gtdb id
        elif len(res)>1:
            all_lineages=[]
            for gtdb_id,lineage,ncbi_id in res:
                lineage=lineage.split(self.info_splitter)
                all_lineages.append(lineage)
            lca=self.get_lowest_common_ancestor_gtdb(all_lineages)
            res=lca
        return res

    def fetch_ncbi_lineage(self,ncbi_id):
        fetch_command = f'SELECT NCBI,LINEAGE FROM NCBILINEAGE WHERE NCBI = "{ncbi_id}"'
        res_fetch=self.cursor.execute(fetch_command).fetchone()
        if res_fetch:
            taxon_id,lineage=res_fetch
            lineage=lineage.split(',')
            lineage.append(str(taxon_id))
            return lineage
        return []


    def launch_taxonomy_connector(self):
        if os.path.exists(self.taxonomy_db_file):
            self.start_taxonomy_connection()
            return True
        else:
            print(f'Database file {self.taxonomy_db_file} does not exist')
            return False

    def get_lowest_common_ancestor_ncbi(self,list_taxons):
        print('Retrieved multiple NCBI entries, retrieving lowest common ancestry NCBI ID')
        all_lineages=[]
        min_size=None
        for taxon_id in list_taxons:
            taxon_lineage=self.fetch_ncbi_lineage(taxon_id)
            if taxon_lineage:
                all_lineages.append(taxon_lineage)
                if not min_size: min_size=len(taxon_lineage)
                if len(taxon_lineage)<min_size: min_size=len(taxon_lineage)
        i=0
        for i in range(min_size):
            temp=set()
            for taxon_lineage in all_lineages:
                temp.add(taxon_lineage[i])
            if len(temp)>1:
                break
        lca=all_lineages[0][0:i+1]
        return lca[-1]

    def get_lowest_common_ancestor_gtdb(self,list_taxon_lineages):
        print('Retrieved multiple GTDB entries, retrieving lowest common ancestry GTDB ID')
        min_size=None
        for taxon_lineage in list_taxon_lineages:
            if not min_size: min_size=len(taxon_lineage)
            if len(taxon_lineage)<min_size: min_size=len(taxon_lineage)
        i=0
        for i in range(min_size):
            temp=set()
            for taxon_lineage in list_taxon_lineages:
                temp.add(taxon_lineage[i])
            if len(temp)>1:
                break
        lca=list_taxon_lineages[0][0:i+1]
        return lca[-1]

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
        if taxa_id: return str(taxa_id)
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={organism_name}'
        taxa_id = self.get_taxa_ncbi_url(url)
        if not taxa_id:
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=candidatus+{organism_name}'
            taxa_id = self.get_taxa_ncbi_url(url)
        if not taxa_id:    print(f'Could not find taxa ID for {organism_name}')
        return str(taxa_id)


if __name__ == '__main__':
    gtdb_connector=Taxonomy_SQLITE_Connector(resources_folder='/home/pedroq/Desktop/test_cr/')
    gtdb_connector.launch_taxonomy_connector()
    gtdb_connector.create_taxonomy_db()
    #gtdb_connector.process_gtdb_taxonomy('d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae_C;g__JACEMX01;s__JACEMX01 sp011773785')
    a=gtdb_connector.fetch_ncbi_id('d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanolobus;s__Methanolobus psychrophilus')
    print(a)
    a=gtdb_connector.fetch_ncbi_id('Clostridium_P perfringens')
    print(a)
    a=gtdb_connector.fetch_gtdb_id('1423')
    print(a)
