import os
import json
import re
from typing import Iterator
import pickle

from mantis.io.clients.metadata_sqlite_client import MetadataSqliteClient
from mantis.core.metadata.metadata_document import BaseMetadataDocument
from mantis.core.metadata.metadata_line_document import MetadataLineDocument
from mantis.core.metadata.metadata_file_writers import MetadataFileWriter, MetadataGffFileWriter
from mantis.io.logger import logger
from mantis.settings import STATIC_DATA
from mantis.core.utils.utils import BATCH_SIZE
from mantis.io.config_reader import ConfigReader

class Metadata():
    def __init__(self):
        self.__metadata_clients = {}
        self.__invalid_annotations = json.loads(open(os.path.join(STATIC_DATA, 'invalid_annotations.json')).read())
        self.__metabolic_modules = json.loads(open(os.path.join(STATIC_DATA, 'metabolic_modules1.json')).read())
        self.__essential_genes = json.loads(open(os.path.join(STATIC_DATA, 'essential_genes', 'essential_genes.json')).read())
        self.__config = ConfigReader()

    def get_metadata_client(self, metadata_file: str):
        if metadata_file not in self.__metadata_clients:
            self.__metadata_clients[metadata_file] = MetadataSqliteClient(metadata_file=metadata_file)
        return self.__metadata_clients[metadata_file]


    def is_invalid_metadata(self, metadata: str):
        normalized_metadata = metadata.lower()
        for invalid_str in self.__invalid_annotations:
            if invalid_str in normalized_metadata:
                return True
        return False

    def get_custom_reference_path(self, target, folder):
        for custom_reference in self.__config.get_custom_reference_paths(folder=folder):
            if target in custom_reference:
                return custom_reference

    def set_is_essential(self, line_documents: list[MetadataLineDocument]):
        for line_document in line_documents:
            metadata: BaseMetadataDocument = line_document.metadata
            valid_ids = set()
            valid_ids.update(metadata.tigrfam)
            valid_ids.update(metadata.pfam)
            if valid_ids.intersection(self.__essential_genes):
                line_document.is_essential_gene = True

    def get_metadata_file(self, reference_file: str):
        metadata_file = None
        if re.search('NOG[GT]', reference_file):
            if 'NOGG' in reference_file:
                taxon_id = 'NOGG'
            else:
                taxon_id = re.search(r'NOGT\d+', reference_file).group().replace('NOGT', '')
            metadata_file = f"{self.mantis_paths['NOG']}{taxon_id}"
        elif re.search('NCBI[GT]', reference_file):
            if 'NCBIG' in reference_file:
                taxon_id = 'NCBIG'
            else:
                taxon_id = re.search(r'NCBIT\d+', reference_file).group().replace('NCBIT', '')
            metadata_file = f"{self.mantis_paths['NCBI']}{taxon_id}"
        elif reference_file == 'Pfam-A':
            metadata_file = self.mantis_paths['pfam']
        elif reference_file == 'kofam_merged':
            metadata_file = self.mantis_paths['kofam']
        elif reference_file == 'tcdb':
            reference_path = self.mantis_paths['tcdb']
        else:
            reference_path = self.get_custom_reference_path(reference_file=reference_file, folder=True)
        metadata_file = os.path.join(reference_path, 'metadata.tsv')
        if not os.path.exists(reference_path):
            logger.debug(f'Metadata file does not exist, skipping... {metadata_file}')
            return None
        return metadata_file

    def get_lines_metadata(self, reference_file: str, line_documents: list[MetadataLineDocument]):
        metadata_file = self.get_metadata_file(reference_file=reference_file)
        if not metadata_file:
            return
        if not os.path.exists(metadata_file):
            logger.info(f'Metadata file does not exist for {reference_file}. Skipping metadata imputating...')
            return
        metadata_client: MetadataSqliteClient = self.get_metadata_client(metadata_file)
        # since multiple lines can share the same reference ID, we first map them by reference_id
        reference_ids_to_lines = []
        for line_document in line_documents:
            if line_document.reference_id not in reference_ids_to_lines:
                reference_ids_to_lines[line_document.reference_id] = []
            reference_ids_to_lines[line_document.reference_id].append(line_document)

        # now we get metadata by reference_id
        reference_ids_to_metadata_documents: dict[str, BaseMetadataDocument] = metadata_client.get_metadata_batch(reference_ids=reference_ids_to_lines.keys())
        # and then we simply assign the same metadata document to each line that shares the same reference id
        for reference_id, metdata_document in reference_ids_to_metadata_documents.items():
            line_document: MetadataLineDocument
            for line_document in reference_ids_to_lines[reference_id]:
                line_document.metadata_documents.append(metdata_document)
        # we also extract metadata from the reference ids or accession numbers
        for line_document in line_documents:
            line_document.metadata.append(BaseMetadataDocument.from_string(line_document.reference_id))
            if line_document.reference_id_accession:
                line_document.metadata.append(BaseMetadataDocument.from_string(line_document.reference_id_accession))
        # when all the metadata is extracted, we merge all the metadata documents
        for line_document in line_documents:
            line_document.merge_metadata()
        # and we finally check if the line has an essential gene or not
        self.set_is_essential(line_documents=line_documents)

    def yield_output_annotation(self, input_file: str) -> Iterator[MetadataLineDocument]:
        with open(input_file) as file:
            file.readline()
            for line_idx, line in enumerate(file):
                line = line.strip('\n').split('\t')
                line.insert(0, line_idx)
                line_document = MetadataLineDocument(*line)
                if self.nog_db == 'hmm' and 'NOG' in line_document.reference_file:
                    line_document.reference_id = line.reference_id.split('.')[0]
                yield line_document

    def yield_output_annotation_by_reference_file(self, input_file: str) -> Iterator[str, list[MetadataLineDocument]]:
        res = {}
        line_document: MetadataLineDocument
        for line_document in self.yield_output_annotation(input_file=input_file):
            if line_document.reference_file not in res:
                res[line_document.reference_file] = []
            if res[line_document.reference_file] > BATCH_SIZE:
                yield line_document.reference_file, res[line_document.reference_file]
                res[line_document.reference_file] = []
            res[line_document.reference_file].append(line_document)
        for reference_file, line_documents in res.items():
            yield reference_file, line_documents

    def yield_annotated_lines(self, input_file) -> Iterator[MetadataLineDocument]:
        for reference_file, line_documents in self.yield_output_annotation_by_reference_file(input_file=input_file):
            self.get_lines_metadata(reference_file=reference_file,
                                    line_documents=line_documents)
            line_document: MetadataLineDocument
            for line_document in line_documents:
                yield line_document.to_file()

    def generate_integrated_output(self, input_file, output_file):
        file_writer = MetadataFileWriter(file_path=output_file)
        if self.output_gff:
            gff_file_writer = MetadataGffFileWriter(file_path=output_file)
        line_document: MetadataLineDocument
        for line_document in self.yield_annotated_lines(input_file=input_file):
            file_writer.write(line_document=line_document)
            if self.output_gff:
                gff_file_writer.write(line_document=line_document)


    ##### TODO continue here

    # TODO check if this is necessary
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




    ######for KEGG module matrix#####
    def generate_module_col(self, tree_modules):
        res = []
        for sk in self.__metabolic_modules:
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

        res = {}
        for sk in self.__metabolic_modules:
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
        res = {}
        for sk in self.__metabolic_modules:
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
        modules_tree = pickle.load(open(os.path.join(self.mantis_paths['resources'], 'kegg', 'modules.json'), 'rb'))
        self.export_sample_kos(sample_paths)
        if modules_tree:
            module_col = self.generate_module_col(modules_tree)
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
            logger.error('KEGG modules pickle is not present, so Mantis cannot create the KEGG matrix')


if __name__ == '__main__':
    m = Metadata()
