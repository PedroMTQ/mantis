try:
    from source.MANTIS_Assembler import *
except:
    from MANTIS_Assembler import *


class MANTIS_Metadata():

    def get_target_custom_hmms_paths(self, target, folder):
        for custom_hmm in self.get_custom_hmms_paths(folder=folder):
            if target in custom_hmm: return custom_hmm

    def is_good_annotation(self, to_add):
        if 'unknown' in to_add:
            return False
        elif 'hypothetical protein' in to_add:
            return False
        elif 'hypothetical enzyme' in to_add:
            return False
        elif 'putative protein' in to_add:
            return False
        elif 'putative enzyme' in to_add:
            return False

        return True

    def add_to_dict(self, target_hmm_dict, dict_key, to_add):
        if not to_add: return
        if dict_key not in target_hmm_dict['link']:
            target_hmm_dict['link'][dict_key] = []
        if isinstance(to_add, str):
            list_to_add = [to_add]
        else:
            list_to_add = to_add
        for i in list_to_add:
            if self.is_good_annotation(i.lower()):
                if i not in target_hmm_dict['link'][dict_key]:
                    i = i.strip()
                    if i:
                        target_hmm_dict['link'][dict_key].append(i)

    # This is the default interpreter since we should always have NOG annotations, the others interpreters are built according to the available hmms
    def get_link_compiled_metadata(self, dict_hmms, hmm_file_path):
        # we can extract ECs,KOs and pfam ids from NOG hmmss
        with open(hmm_file_path, 'r') as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.split('\t')
                current_hmm = line[0]
                if current_hmm in dict_hmms:
                    annotations = line[2:]
                    for link in annotations:
                        if link:
                            temp_link = link.split(':')
                            link_type = temp_link[0]
                            if link_type == 'kegg_cazy': link_type = 'cazy'
                            if link_type == 'kegg_ec': link_type = 'enzyme_ec'
                            link_text = ':'.join(temp_link[1:])
                            link_text=link_text.strip()
                            if link_type == 'description' and link_text == 'NA':
                                link_text = ''
                            if link_text and link_type == 'description':
                                self.get_common_links(link_text, res=dict_hmms[current_hmm])
                            if link_text:
                                self.add_to_dict(dict_hmms[current_hmm], link_type, link_text)
                    if '_sql_annotations.tsv' in hmm_file_path:
                        self.add_to_dict(dict_hmms[current_hmm], 'eggnog', current_hmm)
                line = file.readline()

    def get_direct_link(self, dict_hmms, link_type):
        for hmm in dict_hmms:
            self.add_to_dict(dict_hmms[hmm], link_type, hmm)

    def get_link_custom_hmm(self, dict_hmms, custom_hmm_path):
        if not custom_hmm_path: return
        file_path = custom_hmm_path.replace('.hmm', '.tsv')
        headers = {}
        if not os.path.exists(file_path): return
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                # first line should be the metadata headers, it should always start with the hmm id and then come the corresponding data
                if not headers:
                    headers = {0: 'hmm_id'}
                    for header_i in range(len(line[1:])):
                        headers[header_i+1] = line[header_i+1].lower()
                else:
                    hmm_id = line[0]
                    if hmm_id in dict_hmms:
                        for annot_i in range(len(line[1:])):
                            annot_i += 1
                            if annot_i in headers:
                                self.add_to_dict(dict_hmms[hmm_id], headers[annot_i], line[annot_i].split(';'))
                            else:
                                self.get_common_links(line[annot_i], dict_hmms[hmm_id])
                                self.add_to_dict(dict_hmms[hmm_id], 'description', line[annot_i])
                line = file.readline()

    def get_common_links(self, string, res={}):
        ec = find_ecs(string)
        if ec:
            self.add_to_dict(res, 'enzyme_ec', ec)
        tc = find_tcdb(string)
        if tc:
            self.add_to_dict(res, 'tcdb', tc)
        tigr = find_tigrfam(string)
        if tigr:
            self.add_to_dict(res, 'tigrfam', tigr)
        ko = find_ko(string)
        if ko:
            self.add_to_dict(res, 'kegg_ko', ko)
        pfam = find_pfam(string)
        if pfam:
            self.add_to_dict(res, 'pfam', pfam)
        cog = find_cog(string)
        if cog:
            self.add_to_dict(res, 'cog', cog)
        go = find_go(string)
        if go:
            self.add_to_dict(res, 'go', cog)
        return res

    def get_link_kofam(self, dict_hmms):
        self.get_link_kofam_ko_list(dict_hmms)
        self.get_link_kofam_ko_to_binary(dict_hmms, target_file='ko2cog.xl')
        self.get_link_kofam_ko_to_binary(dict_hmms, target_file='ko2go.xl')
        self.get_link_kofam_ko_to_binary(dict_hmms, target_file='ko2tc.xl')
        self.get_link_kofam_ko_to_binary(dict_hmms, target_file='ko2cazy.xl')
        self.get_link_kofam_ko_to_pathway(dict_hmms)

    def get_link_kofam_ko_list(self, dict_hmms):
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
                    self.add_to_dict(dict_hmms[ko], 'kegg_ko', ko)
                    self.add_to_dict(dict_hmms[ko], 'description', description)
                line = file.readline()

    def get_link_kofam_ko_to_binary(self, dict_hmms, target_file):
        '''
        Data:
        ftp://ftp.genome.jp/pub/db/kofam/
        '''
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
                if ko in dict_hmms:
                    self.add_to_dict(dict_hmms[ko], target_link, link)
                line = file.readline()

    def get_link_kofam_ko_to_pathway(self, dict_hmms):
        file_path = self.mantis_paths['kofam'] + 'ko_to_path'
        map_description = self.get_kofam_pathway_description()
        with open(file_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                ko, link = line
                ko = ko.split(':')[1]
                link = link.split(':')[1]

                if ko in dict_hmms:
                    if 'map' in link:
                        self.add_to_dict(dict_hmms[ko], 'kegg_map', link)
                        if link in map_description:
                            link_id = str(link)
                        else:
                            link_id = 'kegg_ko' + str(link)[3:]
                        if link_id in map_description:
                            self.add_to_dict(dict_hmms[ko], 'kegg_map_lineage',
                                             map_description[link_id]['grand_parent'] + ' -> ' +
                                             map_description[link_id]['parent'] + '-> ' +
                                             map_description[link_id]['description'] +
                                             ' (' + link + ')')
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
                    if 'map=' in kegg_map:
                        kegg_map = kegg_map.split('map=')[-1]
                    else:
                        kegg_map = kegg_map.split('show_pathway?')[-1]
                    description = re.search('<a href=.*>.*<\/a>', line).group()
                    description = description.split('>')[1].split('<')[0]
                    res[kegg_map] = {'description': description, 'parent': sub_title, 'grand_parent': main_tile}
                line = file.readline()
        return res

    def add_tigrfam_go_link(self, dict_hmms):
        go_link_path = self.mantis_paths['tigrfam'] + 'TIGRFAMS_GO_LINK'
        with open(go_link_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                if line[0] in dict_hmms:
                    self.add_to_dict(dict_hmms[line[0]], 'go', line[1].split(':')[-1])
                    self.add_to_dict(dict_hmms[line[0]], 'tigrfam', line[0])
                line = file.readline()

    def get_tigrfam_role_link_id(self, dict_hmms):
        role_link_path = self.mantis_paths['tigrfam'] + 'TIGRFAMS_ROLE_LINK'
        role_links = {}
        with open(role_link_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                # line[0] is the hmm, line[1] is the role link
                if line[0] in dict_hmms:
                    if line[1] not in role_links: role_links[line[1]] = []
                    role_links[line[1]].append(line[0])
                line = file.readline()
        return role_links

    def add_tigrfam_role_names(self, dict_hmms):
        role_links = self.get_tigrfam_role_link_id(dict_hmms)
        role_names_path = self.mantis_paths['tigrfam'] + 'TIGR_ROLE_NAMES'
        with open(role_names_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                # line[1] is the role link
                if line[1] in role_links:
                    for hmm in dict_hmms:
                        if hmm in role_links[line[1]]:
                            if line[3] not in ['Unknown', 'Other', 'General']:
                                self.add_to_dict(dict_hmms[hmm], 'description', line[3] + ' (' + line[2][:-1] + ')')
                line = file.readline()

    def get_link_tigrfam(self, dict_hmms):
        '''
        Data:
        ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.tar.gz
        ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMS_GO_LINK
        ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGRFAMS_ROLE_LINK
        ftp://ftp.tigr.org/pub/data/TIGRFAMs/TIGR_ROLE_NAMES

        '''
        self.add_tigrfam_go_link(dict_hmms)
        self.add_tigrfam_role_names(dict_hmms)

    def get_essential_genes_list(self) -> object:
        essential_genes = self.mantis_paths['resources'] + 'essential_genes/essential_genes.txt'
        if os.path.exists(essential_genes):
            with open(essential_genes) as file: lines = file.readlines()
            lines = [l.strip('\n') for l in lines]
            return lines

    def is_essential(self, dict_hmms):
        essential_genes_list = self.get_essential_genes_list()
        if essential_genes_list:
            for essential_gene in essential_genes_list:
                if essential_gene in dict_hmms:
                    self.add_to_dict(dict_hmms[essential_gene], 'is_essential_gene', 'True')

    def get_hmm_links(self, dict_hmms, hmm_file):
        if 'NOG' in hmm_file:
            if 'NOGG' in hmm_file:
                taxon_id = 'NOGG'
            else:
                taxon_id = re.search('NOGT\d+', hmm_file).group().replace('NOGT', '')
            target_sql_file = add_slash(self.mantis_paths['NOG'] + taxon_id) + f'{taxon_id}_sql_annotations.tsv'
            self.get_link_compiled_metadata(dict_hmms=dict_hmms, hmm_file_path=target_sql_file)
        elif 'Pfam' in hmm_file:
            self.get_link_compiled_metadata(dict_hmms=dict_hmms,
                                            hmm_file_path=self.mantis_paths['pfam'] + 'pfam_metadata.tsv')
            self.is_essential(dict_hmms)
        # direct linking:
        elif 'metacyc' in hmm_file:
            self.get_direct_link(dict_hmms, 'biocyc_rxn')
        elif 'Burstein2016' in hmm_file:
            # article: https://www.nature.com/articles/nature21059
            # Data:
            # http://www.nature.com.proxy.bnl.lu/nature/journal/v542/n7640/full/nature21059.html
            # http://www.nature.com/nature/journal/v542/n7640/extref/nature21059-s3.zip
            self.get_direct_link(dict_hmms, 'cas')

        elif 'kofam' in hmm_file:
            self.get_link_kofam(dict_hmms)
        elif 'NCBI' in hmm_file:
            if 'NCBIG' in hmm_file:
                taxon_id = 'NCBIG'
            else:
                taxon_id = re.search('NCBIT\d+', hmm_file).group().replace('NCBIT', '')
            metadata_file = add_slash(self.mantis_paths['NCBI'] + taxon_id) + 'metadata.tsv'
            self.get_link_compiled_metadata(dict_hmms=dict_hmms, hmm_file_path=metadata_file)
        elif 'tigrfam' in hmm_file:
            self.get_link_tigrfam(dict_hmms)
            self.is_essential(dict_hmms)
        else:
            self.get_link_custom_hmm(dict_hmms, self.get_target_custom_hmms_paths(hmm_file, folder=False))

        for hmm in dict_hmms:
            self.get_common_links(hmm, dict_hmms[hmm])
            if 'accession' in dict_hmms[hmm]['link']:   self.get_common_links(dict_hmms[hmm]['link']['accession'], dict_hmms[hmm])

        return dict_hmms

    def remove_ids_text(self, sorted_keys, temp_link, target_removal):
        for link_key in sorted_keys:
            if isinstance(temp_link[link_key], str): temp_link[link_key] = [temp_link[link_key]]
            for inner_l in temp_link[link_key]:
                if target_removal == 'description':
                    for i in range(len(temp_link['description'])):
                        temp_link['description'][i] = temp_link['description'][i].replace(inner_l, '').replace('()','').strip()
                elif target_removal == 'kegg_map_lineage':
                    if 'kegg_map' == link_key and 'kegg_map_lineage' in temp_link:
                        for i in range(len(temp_link['kegg_map_lineage'])):
                            temp_link['kegg_map_lineage'][i] = temp_link['kegg_map_lineage'][i].replace(inner_l,'').replace('()', '').strip()

    def generate_interpreted_line(self, query, hmm_file, link, evalue, bitscore, direction, query_len, query_start, query_end,
                                  hmm_start, hmm_end, hmm_len):
        temp_link = dict(link)
        hmm = temp_link.pop('hmm')
        hmm_accession = '-'
        if 'accession' in temp_link: hmm_accession = temp_link.pop('accession')
        row_start = [query, hmm_file, hmm, hmm_accession, evalue, bitscore,direction, query_len, query_start, query_end,
                     hmm_start, hmm_end, hmm_len, '|']
        res = list(row_start)
        sorted_keys = sorted(temp_link.keys())
        if 'enzyme_ec' in sorted_keys:
            sorted_keys.remove('enzyme_ec')
            sorted_keys.insert(0, 'enzyme_ec')
        # so that description always comes in the end
        if 'kegg_map_lineage' in sorted_keys:
            sorted_keys.remove('kegg_map_lineage')
            self.remove_ids_text(sorted_keys, temp_link, target_removal='kegg_map_lineage')
            sorted_keys.append('kegg_map_lineage')
        if 'description' in sorted_keys:
            sorted_keys.remove('description')
            #self.remove_ids_text(sorted_keys, temp_link, target_removal='description')
            sorted_keys.append('description')
        for link_key in sorted_keys:
            if isinstance(temp_link[link_key], str): temp_link[link_key] = [temp_link[link_key]]
            for inner_l in temp_link[link_key]:
                res.append(link_key + ':' + inner_l)
        return res

    def read_and_interpret_output_annotation(self, output_annotation_tsv):
        c = 0
        links_to_get = {}
        lines_info = {}
        with open(output_annotation_tsv) as file:
            line = file.readline()
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                query, hmm_file, hmm_hit, hmm_hit_accession, evalue, bitscore,direction, query_len, query_start, query_end, hmm_start, hmm_end, hmm_len = line
                hmm_file = hmm_file.replace('domtblout_annotation_', '')
                if 'NOG' in hmm_file: hmm_hit = hmm_hit.split('.')[0]
                if hmm_file not in links_to_get: links_to_get[hmm_file] = {}
                if hmm_hit not in links_to_get[hmm_file]: links_to_get[hmm_file][hmm_hit] = {'link': {'hmm': hmm_hit},
                                                                                             'lines': []}
                if hmm_hit_accession != '-': links_to_get[hmm_file][hmm_hit]['link']['accession'] = hmm_hit_accession
                links_to_get[hmm_file][hmm_hit]['lines'].append(c)
                lines_info[c] = {'query': query, 'evalue': evalue, 'bitscore': bitscore,'direction': direction, 'query_len': query_len,
                                 'query_start': query_start, 'query_end': query_end, 'hmm_start': hmm_start,
                                 'hmm_end': hmm_end, 'hmm_len': hmm_len}
                c += 1
                line = file.readline()
        res = {}
        for hmm_file in links_to_get:
            hmm_file_links = self.get_hmm_links(links_to_get[hmm_file], hmm_file)
            for hmm in hmm_file_links:
                for line in hmm_file_links[hmm]['lines']:
                    res[line] = self.generate_interpreted_line(query=lines_info[line]['query'],
                                                               hmm_file=hmm_file,
                                                               link=hmm_file_links[hmm]['link'],
                                                               evalue=lines_info[line]['evalue'],
                                                               bitscore=lines_info[line]['bitscore'],
                                                               direction=lines_info[line]['direction'],
                                                               query_len=lines_info[line]['query_len'],
                                                               query_start=lines_info[line]['query_start'],
                                                               query_end=lines_info[line]['query_end'],
                                                               hmm_start=lines_info[line]['hmm_start'],
                                                               hmm_end=lines_info[line]['hmm_end'],
                                                               hmm_len=lines_info[line]['hmm_len'],
                                                               )
        return res

    def generate_interpreted_output(self, output_annotation_tsv, interpreted_annotation_tsv):
        first_line = ['Query',
                      'HMM_file',
                      'HMM_hit',
                      'HMM_hit_accession',
                      'evalue',
                      'bitscore',
                      'Direction',
                      'Query_length',
                      'Query_hit_start',
                      'Query_hit_end',
                      'HMM_hit_start',
                      'HMM_hit_end',
                      'HMM_length',
                      '|',
                      'Links']
        with open(interpreted_annotation_tsv, 'w+') as file:
            first_line = '\t'.join(first_line)
            file.write(first_line + '\n')
            output_annotation = self.read_and_interpret_output_annotation(output_annotation_tsv)
            for line in range(len(output_annotation)):
                if output_annotation[line]:
                    out_line = '\t'.join(output_annotation[line])
                    file.write(out_line + '\n')





    ######for KEGG module matrix#####
    def generate_module_col(self, tree_modules):
        sorted_keys = [
            'Carbohydrate metabolism',
            'Energy metabolism',
            'Lipid metabolism',
            'Nucleotide metabolism',
            'Amino acid metabolism',
            'Glycan metabolism',
            'Metabolism of cofactors and vitamins',
            'Biosynthesis of terpenoids and polyketides',
            'Biosynthesis of other secondary metabolites',
            'Xenobiotics biodegradation',
        ]
        res = []
        for sk in sorted_keys:
            sorted_sub_paths = sorted(tree_modules[sk].keys())
            for ssp in sorted_sub_paths:
                sorted_modules = sorted(tree_modules[sk][ssp])
                for sm in sorted_modules:
                    res.append(sm)
        return res

    def generate_sample_col(self, sample_kos, tree_modules):
        sorted_keys = [
            'Carbohydrate metabolism',
            'Energy metabolism',
            'Lipid metabolism',
            'Nucleotide metabolism',
            'Amino acid metabolism',
            'Glycan metabolism',
            'Metabolism of cofactors and vitamins',
            'Biosynthesis of terpenoids and polyketides',
            'Biosynthesis of other secondary metabolites',
            'Xenobiotics biodegradation',
        ]
        res = []
        for sk in sorted_keys:
            sorted_sub_paths = sorted(tree_modules[sk].keys())
            for ssp in sorted_sub_paths:
                sorted_modules = sorted(tree_modules[sk][ssp])
                for sm in sorted_modules:
                    sm_kos = tree_modules[sk][ssp][sm]
                    module_perc = len(sm_kos.intersection(sample_kos)) / len(sm_kos)
                    res.append(str(round(module_perc, 3)))
        return res

    def get_sample_kos(self, annotation_file):
        res = set()
        sample_name = annotation_file.split(SPLITTER)[-2]
        with open(annotation_file) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.split('\t')
                kegg_kos = [i.replace('kegg_ko:', '') for i in line if 'kegg_ko:' in i]
                res.update(kegg_kos)
                line = file.readline()
        return {'sample': sample_name, 'kegg_ko': res}

    def generate_matrix(self):
        sample_paths = []
        out_file = self.output_folder + 'kegg_modules.tsv'
        for i in self.fastas_to_annotate:
            file_path, sample_output_path, organism_details, count_seqs_original_file = i
            sample_paths.append(sample_output_path + 'consensus_annotation.tsv')

        samples_info = [self.get_sample_kos(i) for i in sample_paths]
        tree = load_metrics(add_slash(self.mantis_paths['resources'] + 'KEGG') + 'modules.pickle')
        if tree:
            module_col = self.generate_module_col(tree)
            with open(out_file, 'w+') as file:
                top_line = [s['sample'] for s in samples_info]
                top_line.insert(0, 'Module')
                top_line = '\t'.join(top_line) + '\n'
                file.write(top_line)
                samples_perc = [self.generate_sample_col(s['kegg_ko'], tree) for s in samples_info]
                for i in range(len(module_col)):
                    module_line = [module_col[i]]
                    for s in samples_perc:
                        module_line.append(s[i])
                    if i == len(module_col) - 1:
                        module_line = '\t'.join(module_line)
                    else:
                        module_line = '\t'.join(module_line) + '\n'
                    file.write(module_line)
        else:
            print('KEGG modules pickle is not present, so Mantis cannot create the KEGG matrix', flush=True,file=self.redirect_verbose)


if __name__ == '__main__':
    m = MANTIS_Metadata()
