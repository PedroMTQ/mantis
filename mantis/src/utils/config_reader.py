import os
from pathlib import Path

from mantis.src.settings import ROOT, DEFAULT_CONFIG
from mantis.src.taxonomy.taxonomy_sqlite_connector import Taxonomy_SQLITE_Connector
from mantis.src.utils.exceptions import InvalidNOGType
from mantis.src.utils.logger import logger


class ConfigReader():
    def __init__(self, config_file: str='', use_taxonomy: bool=False) -> None:
        self.use_taxonomy = use_taxonomy
        self.config_file = config_file
        if self.config_file:
            logger.info(f'Using custom MANTIS.cfg: {self.mantis_config}')
        else:
            self.config_file = DEFAULT_CONFIG
            if not os.path.isdir(self.config_file):
                raise FileNotFoundError(self.config_file)
            logger.info(f'Using default MANTIS.cfg: {self.config_file}')
        self.mantis_ref_weights = {'else': 0.7}
        self.mantis_nogt_tax = set()
        self.mantis_paths = {'default': None,
                             'resources': None,
                             'custom': None,
                             'NOG': None,
                             'pfam': None,
                             'kofam': None,
                             'NCBI': None,
                             'tcdb': None,
                             }

        self.nog_db = 'dmnd'
        self.nogt_line = None




    def read(self):



        config_to_dict_mapping = {
            'default_ref_folder=': 'default',
            'resources_folder=': 'resources',
            'custom_ref_folder=': 'custom',
            'nog_ref_folder=': 'NOG',
            'pfam_ref_folder=': 'pfam',
            'kofam_ref_folder=': 'kofam',
            'ncbi_ref_folder=': 'NCBI',
            'tcdb_ref_folder=': 'tcdb',
            }
        with open(self.config_file) as file:
            for line in file:
                line = line.strip('\n')
                if not line.startswith('#') and line:
                    for config_str, dict_str in config_to_dict_mapping.items():
                        if line.startswith(config_str):
                            line_path = line.replace(config_str, '')
                            self.mantis_paths[dict_str] = line_path
                            break
                    # taxa ids list for only downloading nogt specific to lineage
                    if line.startswith('nog_tax='):
                        self.nogt_line = line.replace('nog_tax=', '')

                    elif '_weight=' in line:
                        ref_source, weight = line.split('_weight=')
                        self.mantis_ref_weights[ref_source] = float(weight)

                    elif line.startswith('nog_ref='):
                        nog_db = line.replace('nog_ref=', '').split()[0]
                        if nog_db.lower() not in ['dmnd', 'hmm']:
                            raise InvalidNOGType()
                        else:
                            self.nog_db = nog_db

        if not self.mantis_paths['default']:
            self.mantis_paths['default'] = os.path.join(ROOT, 'references')

        if not self.mantis_paths['resources']:
            self.mantis_paths['resources'] = os.path.join(ROOT, 'resources')

        for ref in ['custom',
                     'NOG',
                     'pfam',
                     'kofam',
                     'NCBI',
                     'tcdb',
                     ]:
            self.mantis_paths[ref] = self.mantis_paths.get(ref, os.paht.join(self.mantis_paths['default'], ref))



    def __str__(self):
        custom_refs = self.get_custom_refs_paths(folder=True)
        custom_refs_str = ''
        custom_res = ''
        for cref in custom_refs:
            custom_refs_str += cref + '\n'
        if custom_refs_str:
            custom_res = f'# Custom references:\n{custom_refs_str}'

        res = []
        if hasattr(self, 'output_folder'):
            res.append('Output folder:\nself.output_folder')
        res.append(f'Default references folder:\n{self.mantis_paths["default"]}')
        res.append(f'Resources folder:\n{self.mantis_paths["resources"]}')
        res.append(f'Custom references folder:\n{self.mantis_paths["custom"]}')
        if self.mantis_paths['NOG'][0:2] != 'NA':
            res.append(f'TAX NOG references folder:\n{self.mantis_paths["NOG"]}')
        if self.mantis_paths['NCBI'][0:2] != 'NA':
            res.append(f'TAX NCBI references folder:\n{self.mantis_paths["NCBI"]}')
        if self.mantis_paths['pfam'][0:2] != 'NA':
            res.append(f'Pfam reference folder:\n{self.mantis_paths["pfam"]}')
        if self.mantis_paths['kofam'][0:2] != 'NA':
            res.append(f'KOfam reference folder:\n{self.mantis_paths["kofam"]}')
        if self.mantis_paths['tcdb'][0:2] != 'NA':
            res.append(f'TCDB reference folder:\n{self.mantis_paths["tcdb"]}')
        res.append('------------------------------------------')
        res = '\n'.join(res)
        if custom_res: res += '\n' + custom_res
        ref_weights = ', '.join([f'{i}:{self.mantis_ref_weights[i]}' for i in self.mantis_ref_weights if i != 'else'])
        if ref_weights:
            res += f'#  Weights:\n{ref_weights}\n'
        nog_tax = ', '.join([i for i in self.mantis_nogt_tax])
        if nog_tax:
            res += f'\n#  NOG tax IDs:\n{nog_tax}\n'

        return res

    def read_config_file(self):

        self.setup_paths_config_file()
        if not self.use_taxonomy:
            self.mantis_paths.pop('NOG')
            self.mantis_paths.pop('NCBI')
        if not os.path.isdir(self.mantis_paths['custom']):
            Path(self.mantis_paths['custom']).mkdir(parents=True, exist_ok=True)
        logger.debug(self)


    def get_custom_refs_paths(self, folder=False):
        try:
            custom_refs_folders = os.listdir(self.mantis_paths['custom'])
            for potential_ref_folder in custom_refs_folders:
                try:
                    files = os.listdir(self.mantis_paths['custom'] + potential_ref_folder)
                    for potential_file in files:
                        if potential_file.endswith('.hmm') or potential_file.endswith('.dmnd'):
                            if folder:
                                try:
                                    yield os.path.join(self.mantis_paths['custom'], potential_ref_folder)
                                except GeneratorExit:
                                    return ''
                            else:
                                try:
                                    yield os.path.join(self.mantis_paths['custom'], potential_ref_folder, potential_file)
                                except GeneratorExit:
                                    return ''
                except Exception:
                    pass
        except Exception:
            logger.error('Custom references folder is missing, did you correctly set the path? If path is not set make sure you didn\'t delete the custom_ref folder!')
            self.passed_check = False
            return
        with open(self.config_file, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                if 'custom_ref=' in line:
                    line = line.strip('\n')
                    ref_path = line.replace('custom_ref=', '')
                    if not (ref_path.endswith('.hmm') or ref_path.endswith('.dmnd')):
                        if os.path.isdir(ref_path):
                            for inner_file in os.listdir(ref_path):
                                if inner_file.endswith('.hmm') or inner_file.endswith('.dmnd'):
                                    ref_path = os.path.join(ref_path, inner_file)
                    if folder:
                        try:
                            yield os.path.dirname(ref_path)
                        except GeneratorExit:
                            return ''
                    else:
                        try:
                            yield ref_path
                        except GeneratorExit:
                            return ''


    def set_nogt_line(self, line_path):
        if line_path:
            res = set()
            tax_ids = [i for i in line_path.split(',')]
            if tax_ids:
                for t_id in tax_ids:
                    try:
                        ncbi_taxon_id = int(t_id)
                        organism_lineage = self.fetch_ncbi_lineage(ncbi_taxon_id)
                        res.update(organism_lineage)
                    except:
                        ncbi_taxon_id = self.get_taxa_ncbi(t_id)
                        if ncbi_taxon_id:
                            organism_lineage = self.fetch_ncbi_lineage(ncbi_taxon_id)
                            res.update(organism_lineage)
            for i in res:
                self.mantis_nogt_tax.add(str(i))

    def setup_paths_config_file(self):
        self.read()

        # setting up which taxa we need to have references for
        Taxonomy_SQLITE_Connector.__init__(self, resources_folder=self.mantis_paths['resources'])
        if self.use_taxonomy:
            if self.nogt_line:
                if self.launch_taxonomy_connector():
                    self.set_nogt_line(self.nogt_line)
