try:
    from source.MANTIS_Assembler import *
except:
    from MANTIS_Assembler import *

class MANTIS_Interpreter():

    def get_target_custom_hmms_paths(self,target,folder):
        for custom_hmm in self.get_custom_hmms_paths(folder=folder):
            if target in custom_hmm: return custom_hmm

    def add_to_dict(self,target_hmm_dict,dict_key,to_add):
        if not to_add: return
        if dict_key not in target_hmm_dict['link']:
            target_hmm_dict['link'][dict_key] = []
        if isinstance(to_add,str): list_to_add=[to_add]
        else: list_to_add=to_add
        for i in list_to_add:
            if 'unknown' not in i.lower():
                if i not in target_hmm_dict['link'][dict_key]:
                    target_hmm_dict['link'][dict_key].append(i)

    #This is the default interpreter since we should always have NOG annotations, the others interpreters are built according to the available hmms
    def get_link_NOG(self,dict_hmms,hmm_file_path):
        #we can extract ECs,KOs and pfam ids from NOG hmmss
        with open(hmm_file_path,'r') as file:
            line=file.readline()
            while line:
                line=line.strip('\n')
                line=line.split('\t')
                current_hmm=line[0]
                if current_hmm in dict_hmms:
                    annotations=line[2:]
                    for link in annotations:
                        if link:
                            temp_link=link.split(':')
                            link_type=temp_link[0]
                            if link_type=='kegg_cazy': link_type='cazy'
                            if link_type=='kegg_ec': link_type='enzyme_ec'
                            link_text=':'.join(temp_link[1:])
                            if link_type=='description' and link_text=='NA':
                                link_text=''
                            if link_text and link_type=='description':
                                self.get_common_links(link_text, res=dict_hmms[current_hmm])
                            if link_text:
                                self.add_to_dict(dict_hmms[current_hmm],link_type,link_text)
                line=file.readline()

    def is_redundant_description_pfam(self,hmm,row_description):
        temp=[i.lower() for i in row_description.split()]
        if hmm.lower() in temp and 'protein' in temp and len(temp)==2:
            return True
        return False

    def get_link_pfam(self,dict_hmms):
        file_path=self.mantis_paths['pfam']+'Pfam-A.hmm.dat'
        with open(file_path) as file:
            line=file.readline()
            stop=False
            while line:
                line=line.strip('\n').split('   ')
                if len(line)==2:
                    row_header, row_description=line
                    if row_header=='#=GF ID':
                        if row_description in dict_hmms:
                            stop = True
                            hmm = str(row_description)
                    if row_header=='#=GF DE' and stop:
                        self.add_to_dict(dict_hmms[hmm], 'pfam', hmm)
                        pfam_accession=dict_hmms[hmm]['link']['accession'].split('.')[0]
                        self.add_to_dict(dict_hmms[hmm], 'pfam', pfam_accession)
                        if not self.is_redundant_description_pfam(hmm,row_description):
                            self.add_to_dict(dict_hmms[hmm], 'description', row_description)
                        self.get_common_links(row_description,dict_hmms[hmm])
                        stop = False
                line=file.readline()

    def get_link_hamap(self,dict_hmms):
        '''
        article:
        https://www.ncbi.nlm.nih.gov/pubmed/24642063
        Data:
        https://github.com/tseemann/prokka/tree/master/db/hmm
        '''
        self.get_direct_link(dict_hmms,'hamap')

    def get_link_cas(self,dict_hmms):
        '''
        article:
        https://www.nature.com/articles/nature21059
        Data:
        http://www.nature.com.proxy.bnl.lu/nature/journal/v542/n7640/full/nature21059.html
        http://www.nature.com/nature/journal/v542/n7640/extref/nature21059-s3.zip
        '''
        self.get_direct_link(dict_hmms,'cas')

    def get_link_metacyc(self,dict_hmms):
        '''
        Custom built hmm
        '''
        #gives metacyc RXN ids
        self.get_direct_link(dict_hmms,'biocyc_rxn')

    def get_direct_link(self,dict_hmms,link_type):
        for hmm in dict_hmms:
            self.add_to_dict(dict_hmms[hmm],link_type,hmm)


    def get_link_custom_hmm(self,dict_hmms,custom_hmm_path):
        if not custom_hmm_path: return
        file_path = custom_hmm_path.replace('.hmm','.tsv')
        headers={}
        if not os.path.exists(file_path): return
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                #first line should be the metadata headers, it should always start with the hmm id and then come the corresponding data
                if not headers:
                    headers={1:'hmm_id'}
                    for header_i in range(len(line[1:])):
                        headers[header_i]=line[header_i].lower()
                else:
                    hmm_id=line[0]
                    if hmm_id in dict_hmms:
                        for annot_i in range(len(line[1:])):
                            annot_i+=1
                            if annot_i in headers:
                                self.add_to_dict(dict_hmms[hmm_id],headers[annot_i],line[annot_i])
                            else:
                                self.get_common_links(line[annot_i],dict_hmms[hmm_id])
                                self.add_to_dict(dict_hmms[hmm_id], 'description', line[annot_i])
                line = file.readline()

    def get_common_links(self,string,res={}):
        ec=find_ecs(string)
        if ec:
            self.add_to_dict(res,'enzyme_ec',ec)
        tc=find_tcdb(string)
        if tc:
            self.add_to_dict(res,'tcdb',tc)
        ko=find_ko(string)
        if ko:
            self.add_to_dict(res,'kegg_ko',ko)
        pfam=find_pfam(string)
        if pfam:
            self.add_to_dict(res,'pfam',pfam)
        cog=find_cog(string)
        if cog:
            self.add_to_dict(res,'cog',cog)
        go=find_go(string)
        if go:
            self.add_to_dict(res,'go',cog)
        return res


    def get_link_kofam(self,dict_hmms):
        self.get_link_kofam_ko_list(dict_hmms)
        self.get_link_kofam_ko_to_binary(dict_hmms,target_file='ko2cog.xl')
        self.get_link_kofam_ko_to_binary(dict_hmms,target_file='ko2go.xl')
        self.get_link_kofam_ko_to_binary(dict_hmms,target_file='ko2tc.xl')
        self.get_link_kofam_ko_to_binary(dict_hmms,target_file='ko2cazy.xl')
        self.get_link_kofam_ko_to_pathway(dict_hmms)

    def get_link_kofam_ko_list(self,dict_hmms):
        '''
        Data:
        ftp://ftp.genome.jp/pub/db/kofam/
        '''
        file_path = self.mantis_paths['kofam'] + 'ko_list'
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                ko, description = line[0], line[-1]
                if ko in dict_hmms:
                    if '[EC:' in description:
                        description, temp_links = description.split('[EC:')
                    else:
                        temp_links = description
                    self.get_common_links(temp_links, dict_hmms[ko])
                    self.add_to_dict(dict_hmms[ko],'kegg_ko',ko)
                    self.add_to_dict(dict_hmms[ko],'description',description)
                line = file.readline()

    def get_link_kofam_ko_to_binary(self, dict_hmms,target_file):
        '''
        Data:
        ftp://ftp.genome.jp/pub/db/kofam/
        '''
        file_path = self.mantis_paths['kofam'] + target_file
        if 'ko2tc' in target_file: target_link='tcdb'
        else: target_link=target_file.replace('ko2','').replace('.xl','')
        with open(file_path) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                ko, link = line
                link=link.strip('[]').split(':')[1].split()
                if ko in dict_hmms:
                    self.add_to_dict(dict_hmms[ko],target_link,link)
                line = file.readline()

    def get_link_kofam_ko_to_pathway(self,dict_hmms):
        file_path = self.mantis_paths['kofam'] + 'ko_to_path'
        map_description=self.get_kofam_pathway_description()
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                ko, link = line
                ko=ko.split(':')[1]
                link=link.split(':')[1]

                if ko in dict_hmms:
                    if 'map' in link:
                        self.add_to_dict(dict_hmms[ko],'kegg_map',link)
                        if link in map_description: link_id=str(link)
                        else: link_id='kegg_ko'+str(link)[3:]
                        if link_id in map_description:
                            self.add_to_dict(dict_hmms[ko],'kegg_map_lineage',
                                             map_description[link_id]['grand_parent'] + ' -> '+
                                             map_description[link_id]['parent'] +'-> '+
                                             map_description[link_id]['description'] +
                                             ' ('+link+')')
                line = file.readline()



    def get_kofam_pathway_description(self):
        file_path = self.mantis_paths['kofam'] + 'map_description'
        res = {}
        with open(file_path) as file:
            line = file.readline()
            while line:
                main_search = re.search('<h4>', line)
                if main_search:
                    main_tile = line.replace('<h4>', '').replace('</h4>', '').split()
                    main_tile = ' '.join(main_tile[1:])
                sub_search = re.search('<b>', line)
                if sub_search:
                    sub_title = line.replace('<b>', '').replace('</b>', '').split()
                    sub_title = ' '.join(sub_title[1:])
                map_search = re.search('show_pathway\?map', line)
                if map_search:
                    kegg_map = re.search('show_pathway\?map.*\d+', line).group()
                    kegg_map = kegg_map.split('&amp')[0]
                    if 'map=' in kegg_map: kegg_map = kegg_map.split('map=')[-1]
                    else: kegg_map = kegg_map.split('show_pathway?')[-1]
                    description = re.search('<a href=.*>.*<\/a>', line).group()
                    description = description.split('>')[1].split('<')[0]
                    res[kegg_map] = {'description': description, 'parent': sub_title, 'grand_parent': main_tile}
                line = file.readline()
        return res

    def get_link_resfams(self,dict_hmms):
        '''
        Data:
        http://www.dantaslab.org/resfams
        Interesting for antibiotic resistance research
        '''
        file_path=self.mantis_paths['resfams']+'180102_resfams_metadata_updated_v122.tsv'
        hmm_hit_accession={dict_hmms[i]['link']['accession']:i for i in dict_hmms}
        #just adding ids but metadata file has some more info, not sure if it's relevant
        with open(file_path) as file:
            line=file.readline()
            while line:
                line=line.strip('\n').split('\t')
                resfam_hmm_id,family_name,description,aro,hmm_source=line[0],line[1],line[2],line[4],line[5]
                if resfam_hmm_id in hmm_hit_accession:
                    if hmm_source=='Pfam':
                        self.add_to_dict(dict_hmms[hmm_hit_accession[resfam_hmm_id]],'pfam',family_name)
                    else:
                        self.add_to_dict(dict_hmms[hmm_hit_accession[resfam_hmm_id]],'resfams',family_name)
                        if aro !='Not Available':
                            aro=aro.replace('ARO:','')
                            aro=aro.split('; ')
                            self.add_to_dict(dict_hmms[hmm_hit_accession[resfam_hmm_id]],'aro',aro)
                        self.add_to_dict(dict_hmms[hmm_hit_accession[resfam_hmm_id]],'description',description)
                line=file.readline()

    def get_link_dbcan(self,dict_hmms):
        '''
        Data:
        http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt
        http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07312019.fam.subfam.ec.txt
        '''
        #dbcan has ec numbers specific to ena ids (european nucleotide archive- organism specific) but then merges these into a singular hmm, so we will be using all the ecs split between one HMM different accessions
        file_path=self.mantis_paths['dbcan']+'CAZyDB.07312019.fam.subfam.ec.txt'
        with open(file_path) as file:
            line=file.readline()
            while line:
                cazy_id,ena,ecs=line.strip('\n').split('\t')
                hmm_target=cazy_id+'.hmm'
                if hmm_target in dict_hmms:
                    ecs=ecs.split('|')
                    res=set()
                    for e in ecs:
                        if e not in res:
                            res.add(e)
                    self.add_to_dict(dict_hmms[hmm_target], 'enzyme_ec', res)
                    self.add_to_dict(dict_hmms[hmm_target], 'cazy', cazy_id)
                line=file.readline()

    def add_tigrfam_go_link(self,dict_hmms):
        go_link_path=self.mantis_paths['tigrfam']+'TIGRFAMS_GO_LINK'
        with open(go_link_path) as file:
            line=file.readline()
            while line:
                line=line.strip('\n').split('\t')
                if line[0] in dict_hmms:
                    self.add_to_dict(dict_hmms[line[0]], 'go', line[1].split(':')[-1])
                    self.add_to_dict(dict_hmms[line[0]], 'tigrfam', line[0])
                line=file.readline()

    def get_tigrfam_role_link_id(self,dict_hmms):
        role_link_path=self.mantis_paths['tigrfam']+'TIGRFAMS_ROLE_LINK'
        role_links={}
        with open(role_link_path) as file:
            line=file.readline()
            while line:
                line=line.strip('\n').split('\t')
                #line[0] is the hmm, line[1] is the role link
                if line[0] in dict_hmms:
                    if line[1] not in role_links: role_links[line[1]]=[]
                    role_links[line[1]].append(line[0])
                line=file.readline()
        return role_links


    def add_tigrfam_role_names(self,dict_hmms):
        role_links = self.get_tigrfam_role_link_id(dict_hmms)
        role_names_path=self.mantis_paths['tigrfam']+'TIGR_ROLE_NAMES'
        with open(role_names_path) as file:
            line=file.readline()
            while line:
                line=line.strip('\n').split('\t')
                #line[1] is the role link
                if line[1] in role_links:
                    for hmm in dict_hmms:
                        if hmm in role_links[line[1]]:
                            if line[3] not in ['Unknown','Other','General']:
                                self.add_to_dict(dict_hmms[hmm], 'description', line[3]+' ('+line[2][:-1]+')')
                line=file.readline()

    def get_link_tigrfam(self,dict_hmms):
        '''
        Data:
        ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz
        ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMS_GO_LINK
        ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMS_ROLE_LINK
        ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGR_ROLE_NAMES

        '''
        self.add_tigrfam_go_link(dict_hmms)
        self.add_tigrfam_role_names(dict_hmms)


    def get_link_transportdb(self,dict_hmms):
        hmm_path=self.get_target_custom_hmms_paths('transportdb',folder=True)
        folder_path=get_folder(hmm_path)
        folder_path=splitter.join(folder_path)
        if not os.path.exists(add_slash(folder_path)+'transportdb_metadata'):
            return
        with open(add_slash(folder_path)+'transportdb_metadata') as file:
            line=file.readline()
            while line:
                family,tcdb_id=line.split()
                if family in dict_hmms:
                    self.add_to_dict(dict_hmms[family], 'tcdb', tcdb_id)
                line=file.readline()

    def get_essential_genes_list(self):
        essential_genes=self.mantis_paths['default']+'essential_genes/essential_genes.txt'
        if os.path.exists(essential_genes):
            with open(essential_genes) as file: lines=file.readlines()
            lines=[l.strip('\n') for l in lines]
            return lines

    def is_essential(self,dict_hmms):
        essential_genes_list=self.get_essential_genes_list()
        if essential_genes_list:
            for essential_gene in essential_genes_list:
                if essential_gene in dict_hmms:
                    self.add_to_dict(dict_hmms[essential_gene], 'is_essential_gene', 'True')



    def get_hmm_links(self,dict_hmms,hmm_file):
        #here we need to customize the data we want to extract from each HMM. EGGNOG will be present by default
        #default linking:
        if 'NOGG' in hmm_file:
            target_sql_file = self.mantis_paths['NOGG'] + 'NOGG_sql_annotations.tsv'
            self.get_link_NOG(dict_hmms=dict_hmms,hmm_file_path=target_sql_file)
        elif 'NOGT' in hmm_file:
            taxon_id=re.search('NOGT\d+',hmm_file).group().replace('NOGT','')
            target_sql_file = self.mantis_paths['NOGT'] + taxon_id + splitter + taxon_id + '_sql_annotations.tsv'
            self.get_link_NOG(dict_hmms=dict_hmms,hmm_file_path=target_sql_file)
        elif 'Pfam' in hmm_file:
            self.get_link_pfam(dict_hmms)
            self.is_essential(dict_hmms)
        #direct linking:
        elif 'HAMAP' in hmm_file:
            self.get_link_hamap(dict_hmms)
        elif 'metacyc' in hmm_file:
            self.get_link_metacyc(dict_hmms)
        elif 'Burstein2016' in hmm_file:
            self.get_link_cas(dict_hmms)
        #indirect linking
        elif 'Resfams' in hmm_file:
            self.get_link_resfams(dict_hmms)
        elif 'dbCAN' in hmm_file:
            self.get_link_dbcan(dict_hmms)
        elif 'kofam' in hmm_file:
            self.get_link_kofam(dict_hmms)
        elif 'tigrfam' in hmm_file:
            self.get_link_tigrfam(dict_hmms)
            self.is_essential(dict_hmms)
        else:
            self.get_link_custom_hmm(dict_hmms,self.get_target_custom_hmms_paths(hmm_file,folder=False))
        for hmm in dict_hmms:
            self.get_common_links(hmm,dict_hmms[hmm])
        return dict_hmms



    def generate_interpreted_line(self,query,hmm_file,link,evalue,query_len,query_start,query_end,hmm_start,hmm_end):
        temp_link=dict(link)
        hmm=temp_link.pop('hmm')
        hmm_accession='-'
        if 'accession' in temp_link: hmm_accession=temp_link.pop('accession')
        row_start = [query, hmm_file, hmm, hmm_accession, evalue,query_len,query_start,query_end,hmm_start,hmm_end, '|']
        res = list(row_start)
        sorted_keys = sorted(temp_link.keys())
        # so that description always comes in the end
        if 'kegg_map_lineage' in sorted_keys:
            sorted_keys.remove('kegg_map_lineage')
            sorted_keys.append('kegg_map_lineage')
        if 'description' in sorted_keys:
            sorted_keys.remove('description')
            sorted_keys.append('description')
        for link_key in sorted_keys:
            if isinstance(temp_link[link_key], str): temp_link[link_key] = [temp_link[link_key]]
            for inner_l in temp_link[link_key]:
                res.append(link_key + ':' + inner_l)
        return res


    def read_and_interpret_output_annotation(self,output_annotation_tsv):
        c=0
        links_to_get={}
        lines_info={}
        with open(output_annotation_tsv) as file:
            line=file.readline()
            line=file.readline()
            while line:
                line=line.strip('\n').split('\t')
                if len(line)==10:
                    query,hmm_file,hmm_hit,hmm_hit_accession,evalue,query_len,query_start,query_end,hmm_start,hmm_end=line
                    hmm_file=hmm_file.replace('domtblout_annotation_','')
                    if 'NOGG' in hmm_file: hmm_hit=hmm_hit.split('.')[1]
                    elif 'NOGT' in hmm_file: hmm_hit=hmm_hit.split('.')[0]
                    if hmm_file not in links_to_get: links_to_get[hmm_file]={}
                    if hmm_hit not in links_to_get[hmm_file]: links_to_get[hmm_file][hmm_hit]={'link':{'hmm':hmm_hit},'lines':[]}
                    if hmm_hit_accession!='-':links_to_get[hmm_file][hmm_hit]['link']['accession']=hmm_hit_accession
                    links_to_get[hmm_file][hmm_hit]['lines'].append(c)
                    lines_info[c]={'query':query,'evalue':evalue,'query_len':query_len,'query_start':query_start,'query_end':query_end,'hmm_start':hmm_start,'hmm_end':hmm_end}
                    c+=1
                line=file.readline()
        res={}
        for hmm_file in links_to_get:
            hmm_file_links=self.get_hmm_links(links_to_get[hmm_file],hmm_file)
            for hmm in hmm_file_links:
                for line in hmm_file_links[hmm]['lines']:
                    res[line]=self.generate_interpreted_line(query=lines_info[line]['query'],
                                                             hmm_file=hmm_file,
                                                             link=hmm_file_links[hmm]['link'],
                                                             evalue=lines_info[line]['evalue'],
                                                             query_len=lines_info[line]['query_len'],
                                                             query_start=lines_info[line]['query_start'],
                                                             query_end=lines_info[line]['query_end'],
                                                             hmm_start=lines_info[line]['hmm_start'],
                                                             hmm_end=lines_info[line]['hmm_end'],
                                                             )
        return res



    def generate_interpreted_output(self,output_annotation_tsv,interpreted_annotation_tsv):
        first_line = ['Query',
                      'HMM_file',
                      'HMM_hit',
                      'HMM_hit_accession',
                      'evalue',
                      'Query_length',
                      'Query_hit_start',
                      'Query_hit_end',
                      'HMM_hit_start',
                      'HMM_hit_end',
                      '|',
                      'Links']
        with open(interpreted_annotation_tsv,'w+') as file:
            first_line='\t'.join(first_line)
            file.write(first_line+'\n')
            output_annotation = self.read_and_interpret_output_annotation(output_annotation_tsv)
            for line in range(len(output_annotation)):
                if output_annotation[line]:
                    out_line='\t'.join(output_annotation[line])
                    file.write(out_line+'\n')

if __name__ == '__main__':

    f='/home/pedroq/Desktop/test_inter/output_annotation.tsv'
    f2='/home/pedroq/Desktop/test_inter/interpreted_annotation.tsv'
    custom_hmm='/home/pedroq/Desktop/test_inter/custom.hmm'

    m=MANTIS_Interpreter()
    m.annotation_output_file=f
    m.interpreted_output_file=f2
    m.custom=custom_hmm
    a=m.read_and_interpret_output_annotation(f)
    for i in a:
        print(a[i])