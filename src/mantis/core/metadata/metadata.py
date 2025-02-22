import os
import re

from io.clients.metadata_sqlite_client import MetadataSqliteClient
from mantis.src.metadata.utils import find_arcog, find_cog, find_ecs, find_go, find_ko, find_pfam, find_tcdb, find_tigrfam
from mantis.src.settings import ROOT
from mantis.src.utils import load_metrics


class Metadata():
    def __init__(self):
        self.__metadata_clientt = MetadataSqliteClient


    def get_target_custom_ref_paths(self, target, folder):
        for custom_ref in self.get_custom_refs_paths(folder=folder):
            if target in custom_ref:
                return custom_ref

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

    def add_to_dict(self, dict_hits, dict_key, to_add):
        if not to_add: return
        if dict_key not in dict_hits['link']:
            dict_hits['link'][dict_key] = set()
        if isinstance(to_add, str):
            list_to_add = [to_add]
        else:
            list_to_add = to_add
        for i in list_to_add:
            if self.is_good_annotation(i.lower()):
                if i not in dict_hits['link'][dict_key]:
                    i = i.strip()
                    if i:
                        dict_hits['link'][dict_key].add(i)

    def get_link_compiled_metadata(self, dict_hits, ref_file_path):
        cursor = MetadataSqliteConnector(ref_file_path)
        for hit in dict_hits:
            hit_dict = dict_hits[hit]
            hit_info = cursor.get_metadata(hit)
            for db in hit_info:
                if db not in hit_dict['link']:  hit_dict['link'][db] = set()
                hit_dict['link'][db].update(hit_info[db])
        cursor.close_sql_connection()

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
        arcog = find_arcog(string)
        if arcog:
            self.add_to_dict(res, 'arcog', arcog)
        go = find_go(string)
        if go:
            self.add_to_dict(res, 'go', cog)
        return res

    def get_essential_genes(self):
        essential_genes = os.path.join(ROOT, 'Resources', 'essential_genes', 'essential_genes.txt')
        if os.path.exists(essential_genes):
            with open(essential_genes) as file:
                lines = file.readlines()
            lines = [line.strip('\n') for line in lines]
            return set(lines)

    def is_essential(self, dict_hits):
        essential_genes = self.get_essential_genes()
        if essential_genes:
            for hit in dict_hits:
                valid_ids = set()
                if 'pfam' in dict_hits[hit]['link']:
                    valid_ids.update(dict_hits[hit]['link']['pfam'])
                if 'tigrfam' in dict_hits[hit]['link']:
                    valid_ids.update(dict_hits[hit]['link']['tigrfam'])
                valid_ids.update(find_tigrfam(dict_hits[hit]['link']['hit']))
                valid_ids.update(find_tigrfam(dict_hits[hit]['link']['accession']))
                valid_ids.update(find_pfam(dict_hits[hit]['link']['hit']))
                valid_ids.update(find_pfam(dict_hits[hit]['link']['accession']))
                if valid_ids.intersection(essential_genes):
                    self.add_to_dict(dict_hits[hit], 'is_essential_gene', 'True')

    def get_hit_links(self, dict_hits, ref_file):
        if re.search('NOG[GT]', ref_file):
            if 'NOGG' in ref_file:
                taxon_id = 'NOGG'
            else:
                taxon_id = re.search(r'NOGT\d+', ref_file).group().replace('NOGT', '')
            metadata_file = os.path.join(self.mantis_paths['NOG'] + taxon_id, 'metadata.tsv')
            self.get_link_compiled_metadata(dict_hits=dict_hits, ref_file_path=metadata_file)
        elif re.search('NCBI[GT]', ref_file):
            if 'NCBIG' in ref_file:
                taxon_id = 'NCBIG'
            else:
                taxon_id = re.search(r'NCBIT\d+', ref_file).group().replace('NCBIT', '')
            metadata_file = os.path.join(self.mantis_paths['NCBI'] + taxon_id, 'metadata.tsv')
            self.get_link_compiled_metadata(dict_hits=dict_hits, ref_file_path=metadata_file)
            self.is_essential(dict_hits)

        elif ref_file == 'Pfam-A':
            self.get_link_compiled_metadata(dict_hits=dict_hits,
                                            ref_file_path=self.mantis_paths['pfam'] + 'metadata.tsv')
            self.is_essential(dict_hits)
        elif ref_file == 'kofam_merged':
            self.get_link_compiled_metadata(dict_hits=dict_hits,
                                            ref_file_path=self.mantis_paths['kofam'] + 'metadata.tsv')
        elif ref_file == 'tcdb':
            self.get_link_compiled_metadata(dict_hits=dict_hits,
                                            ref_file_path=self.mantis_paths['tcdb'] + 'metadata.tsv')
        else:
            custom_ref_path = self.get_target_custom_ref_paths(ref_file, folder=True)
            if custom_ref_path:
                file_path = os.path.join(custom_ref_path, 'metadata.tsv')
                if os.path.exists(file_path):
                    self.get_link_compiled_metadata(dict_hits=dict_hits, ref_file_path=file_path)
        for hit in dict_hits:
            self.get_common_links(hit, dict_hits[hit])
            if 'accession' in dict_hits[hit]['link']:
                self.get_common_links(dict_hits[hit]['link']['accession'],
                                                                              dict_hits[hit])
        return dict_hits

    def remove_ids_text(self, sorted_keys, temp_link, target_removal):
        for link_key in sorted_keys:
            if isinstance(temp_link[link_key], str): temp_link[link_key] = [temp_link[link_key]]
            for inner_l in temp_link[link_key]:
                if target_removal == 'description':
                    for i in range(len(temp_link['description'])):
                        temp_link['description'][i] = temp_link['description'][i].replace(inner_l, '').replace('()', '').strip()
                elif target_removal == 'kegg_map_lineage':
                    if 'kegg_map' == link_key and 'kegg_map_lineage' in temp_link:
                        for i in range(len(temp_link['kegg_map_lineage'])):
                            temp_link['kegg_map_lineage'][i] = temp_link['kegg_map_lineage'][i].replace(inner_l, '').replace('()', '').strip()

    def generate_interpreted_line(self, query, ref_file, link, evalue, bitscore, direction, query_len, query_start,
                                  query_end,
                                  ref_start, ref_end, ref_len):
        temp_link = dict(link)
        hit = temp_link.pop('hit')
        hit_accession = '-'
        if 'accession' in temp_link:
            hit_accession = temp_link.pop('accession')
        row_start = [query,
                     ref_file,
                     hit,
                     hit_accession,
                     evalue,
                     bitscore,
                     direction,
                     query_len,
                     query_start,
                     query_end,
                     ref_start,
                     ref_end,
                     ref_len,
                     '|']
        res = []
        sorted_keys = sorted(temp_link.keys())
        if 'enzyme_ec' in sorted_keys:
            sorted_keys.remove('enzyme_ec')
            sorted_keys.insert(0, 'enzyme_ec')
        # so that description always comes in the end
        if 'kegg_map_lineage' in sorted_keys:
            sorted_keys.remove('kegg_map_lineage')
            self.remove_ids_text(sorted_keys=sorted_keys, temp_link=temp_link, target_removal='kegg_map_lineage')
            sorted_keys.append('kegg_map_lineage')
        if 'description' in sorted_keys:
            sorted_keys.remove('description')
            # self.remove_ids_text(sorted_keys, temp_link, target_removal='description')
            sorted_keys.append('description')
        for link_key in sorted_keys:
            if isinstance(temp_link[link_key], str):
                temp_link[link_key] = [temp_link[link_key]]
            for inner_l in temp_link[link_key]:
                res.append(link_key + ':' + inner_l)
        res = sorted(res)
        res = row_start + res
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
                query, ref_file, ref_hit, ref_hit_accession, evalue, bitscore, direction, query_len, query_start, query_end, ref_start, ref_end, ref_len = line
                if self.nog_db == 'hmm' and 'NOG' in ref_file:
                    ref_hit = ref_hit.split('.')[0]
                if ref_file not in links_to_get:
                    links_to_get[ref_file] = {}
                if ref_hit not in links_to_get[ref_file]:
                    links_to_get[ref_file][ref_hit] = {'link': {'hit': ref_hit}, 'lines': []}
                if ref_hit_accession != '-':
                    links_to_get[ref_file][ref_hit]['link']['accession'] = ref_hit_accession
                links_to_get[ref_file][ref_hit]['lines'].append(c)
                lines_info[c] = {'query': query,
                                 'evalue': evalue,
                                 'bitscore': bitscore,
                                 'direction': direction,
                                 'query_len': query_len,
                                 'query_start': query_start,
                                 'query_end': query_end,
                                 'ref_start': ref_start,
                                 'ref_end': ref_end,
                                 'ref_len': ref_len}
                c += 1
                line = file.readline()
        res = {}
        for ref_file in links_to_get:
            ref_file_links = self.get_hit_links(links_to_get[ref_file], ref_file)
            for ref_hit in ref_file_links:
                for line in ref_file_links[ref_hit]['lines']:
                    res[line] = self.generate_interpreted_line(query=lines_info[line]['query'],
                                                               ref_file=ref_file,
                                                               link=ref_file_links[ref_hit]['link'],
                                                               evalue=lines_info[line]['evalue'],
                                                               bitscore=lines_info[line]['bitscore'],
                                                               direction=lines_info[line]['direction'],
                                                               query_len=lines_info[line]['query_len'],
                                                               query_start=lines_info[line]['query_start'],
                                                               query_end=lines_info[line]['query_end'],
                                                               ref_start=lines_info[line]['ref_start'],
                                                               ref_end=lines_info[line]['ref_end'],
                                                               ref_len=lines_info[line]['ref_len'],
                                                               )
        return res

    def generate_gff_line_integrated(self, integrated_line):
        # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        # verified with http://genometools.org/cgi-bin/gff3validator.cgi
        line_split = integrated_line.index('|')
        line_data, annotation_data = integrated_line[:line_split], integrated_line[line_split + 1:]
        query, ref_file, hit, hit_accession, evalue, bitscore, direction, query_len, query_start, query_end, ref_start, ref_end, ref_len = line_data
        # attributes of gff line:

        if hit_accession != '-':
            attributes = f'Name={hit};Target={hit} {ref_start} {ref_end};Alias={hit_accession}'
        else:
            attributes = f'Name={hit};Target={hit} {ref_start} {ref_end}'
        notes = f'Note=ref_file:{ref_file},ref_len:{ref_len}'
        dbxref = []
        ontology_terms = []
        descriptions = []

        for i in annotation_data:
            if not i.startswith('go:'):
                dbxref.append(i)
            elif i.startswith('description:'):
                descriptions.append(i)
            else:
                ontology_terms.append(i)
        if descriptions:
            notes += ',' + ','.join(descriptions)
        all_annotations = None
        if dbxref and ontology_terms:
            all_annotations = 'Dbxref=' + ','.join(dbxref) + ';' + 'Ontology_term=' + ','.join(ontology_terms)
        elif dbxref and not ontology_terms:
            all_annotations = 'Dbxref=' + ','.join(dbxref)
        elif not dbxref and ontology_terms:
            all_annotations = 'Ontology_term=' + ','.join(ontology_terms)
        gff_line = '\t'.join([query, 'Mantis', 'CDS', query_start, query_end, evalue, direction, '0']) + '\t'
        if all_annotations:
            gff_line += ';'.join([attributes, notes, all_annotations])
        else:
            gff_line += ';'.join([attributes, notes])
        sequence_region_line = f'##sequence-region {query} 1 {query_len}'
        return query, sequence_region_line, gff_line

    def generate_integrated_output(self, output_annotation_tsv, interpreted_annotation_tsv):
        first_line = ['Query',
                      'Ref_file',
                      'Ref_hit',
                      'Ref_hit_accession',
                      'evalue',
                      'bitscore',
                      'Direction',
                      'Query_length',
                      'Query_hit_start',
                      'Query_hit_end',
                      'Ref_hit_start',
                      'Ref_hit_end',
                      'Ref_length',
                      '|',
                      'Links']
        first_line = '\t'.join(first_line)
        output_file = open(interpreted_annotation_tsv, 'w+')
        output_file.write(first_line + '\n')

        gff_file = None
        gff_file_path = interpreted_annotation_tsv.replace('.tsv', '.gff')
        gff_dict = {}
        if self.output_gff:
            gff_file = open(gff_file_path, 'w+')
            gff_file.write('##gff-version 3' + '\n')

        # generating output
        output_annotation = self.read_and_interpret_output_annotation(output_annotation_tsv)
        for line in range(len(output_annotation)):
            current_output_line = output_annotation[line]
            if current_output_line:
                out_line = '\t'.join(current_output_line)
                output_file.write(out_line + '\n')

                if gff_file:
                    seq_id, sequence_region_line, gff_line = self.generate_gff_line_integrated(current_output_line)
                    if seq_id not in gff_dict:
                        gff_dict[seq_id] = {'seq_region': sequence_region_line, 'lines': []}
                    gff_dict[seq_id]['lines'].append(gff_line)
        # writing gff
        for seq_id in gff_dict:
            gff_file.write(gff_dict[seq_id]['seq_region'] + '\n')
        for seq_id in gff_dict:
            for annot_line in gff_dict[seq_id]['lines']:
                gff_file.write(annot_line + '\n')

        if gff_file:
            gff_file.close()
        output_file.close()

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
                    module_name, module_paths = tree_modules[sk][ssp][sm]
                    res.append([sm, module_name])
        return res

    def get_best_sample_module_path(self, sample_kos, all_paths):
        best_score = 0
        best_path = set()

        for current_path in all_paths:
            available_kos = current_path.intersection(sample_kos)
            current_score = len(available_kos) / len(current_path)

            if current_score > best_score:
                best_score = current_score
                best_path = current_path

        if best_score:
            available_kos = best_path.intersection(sample_kos)
            missing_kos = best_path.difference(available_kos)
            available = ','.join(available_kos)
            missing = ','.join(missing_kos)
        else:
            available = 'NA'
            missing = 'NA'
        return best_score, available, missing

    def generate_sample_col_verbose(self, sample_kos, tree_modules):
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
        res = {}
        for sk in sorted_keys:
            sorted_sub_paths = sorted(tree_modules[sk].keys())
            for ssp in sorted_sub_paths:
                sorted_modules = sorted(tree_modules[sk][ssp])
                for sm in sorted_modules:
                    module_name, module_paths = tree_modules[sk][ssp][sm]
                    module_perc, available, missing = self.get_best_sample_module_path(sample_kos, module_paths)
                    module_perc = str(round(module_perc, 3))
                    res[sm] = [module_perc, available, missing]
        return res

    def generate_sample_col_non_verbose(self, sample_kos, tree_modules):
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
        res = {}
        for sk in sorted_keys:
            sorted_sub_paths = sorted(tree_modules[sk].keys())
            for ssp in sorted_sub_paths:
                sorted_modules = sorted(tree_modules[sk][ssp])
                for sm in sorted_modules:
                    _, module_paths = tree_modules[sk][ssp][sm]
                    module_perc, _, _ = self.get_best_sample_module_path(sample_kos, module_paths)
                    module_perc = str(round(module_perc, 3))
                    res[sm] = [module_perc]

        return res

    def generate_sample_col(self, sample_kos, tree_modules):
        if self.verbose_kegg_matrix:
            return self.generate_sample_col_verbose(sample_kos, tree_modules)
        else:
            return self.generate_sample_col_non_verbose(sample_kos, tree_modules)

    def get_sample_kos(self, annotation_file):
        res = set()
        sample_name = os.path.split(annotation_file)[-2]
        with open(annotation_file) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.split('\t')
                kegg_kos = [i.replace('kegg_ko:', '') for i in line if 'kegg_ko:' in i]
                res.update(kegg_kos)
                line = file.readline()
        return {'sample': sample_name, 'kegg_ko': res}

    def get_ids_consensus(self, input_file, wanted_db):
        res = {}
        with open(input_file) as file:
            file.readline()
            for line in file:
                line = line.strip('\n')
                line = line.split('\t')
                seq_id = line[0]
                if seq_id not in res:
                    res[seq_id] = set()
                separator = line.index('|')
                annotations = line[separator + 1:]
                for db_annot in annotations:
                    db = db_annot.split(':')[0]
                    if wanted_db == db:
                        # to avoid bad splitting when dealing with descriptions
                        annot = db_annot[len(db) + 1:]
                        res[seq_id].add(annot)
        return res

    def export_sample_kos(self, samples_info):
        out_file = self.output_folder + 'sample_kos.tsv'
        with open(out_file, 'w+') as file:
            firstline = 'Sample\tKOs\n'
            file.write(firstline)
            for sample_path in samples_info:
                sample_name = os.path.split(sample_path)[-2]
                seq_kos = self.get_ids_consensus(sample_path, wanted_db='kegg_ko')
                for seq in seq_kos:
                    seq_name = f'{sample_name}@{seq}'
                    for ko in seq_kos[seq]:
                        line = f'{seq_name}\t{ko}\n'
                        file.write(line)

    def generate_matrix(self):
        sample_paths = []
        figure_info = {}
        out_file = self.output_folder + 'kegg_modules.tsv'
        for i in self.fastas_to_annotate:
            _, sample_output_path, _, _, _, _ = i
            sample_paths.append(os.path.join(sample_output_path, 'consensus_annotation.tsv'))

        samples_info = [self.get_sample_kos(i) for i in sample_paths]
        tree = load_metrics(os.path.join(self.mantis_paths['resources'] + 'KEGG'), 'modules.pickle')
        self.export_sample_kos(sample_paths)
        if tree:
            module_col = self.generate_module_col(tree)
            with open(out_file, 'w+') as file:
                if self.verbose_kegg_matrix:
                    top_line = [f'Completeness_score_{s["sample"]}\tKOs_{s["sample"]}\tMissing_KOs_{s["sample"]}' for s in samples_info]
                    top_line.insert(0, 'Module_ID\tModule_name')

                else:
                    top_line = [s['sample'] for s in samples_info]
                    top_line.insert(0, 'Module_ID')
                top_line = '\t'.join(top_line) + '\n'
                file.write(top_line)
                samples_perc = [[s['sample'], self.generate_sample_col(s['kegg_ko'], tree)] for s in samples_info]
                for i in range(len(module_col)):
                    module_id, module_name = module_col[i]
                    if self.verbose_kegg_matrix:
                        module_key = f'{module_id}\t{module_name}'
                    else:
                        module_key = module_id
                    module_line = [module_key]
                    for sample_names, s in samples_perc:
                        module_line.extend(s[module_id])
                        if sample_names not in figure_info: figure_info[sample_names] = {}
                        figure_info[sample_names][module_id] = s[module_id][0]

                    if i == len(module_col) - 1:
                        module_line = '\t'.join(module_line)
                    else:
                        module_line = '\t'.join(module_line) + '\n'
                    file.write(module_line)
        else:
            print('KEGG modules pickle is not present, so Mantis cannot create the KEGG matrix', flush=True,
                  file=self.redirect_verbose)


if __name__ == '__main__':
    m = Metadata()
