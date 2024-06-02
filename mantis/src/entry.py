from mantis.src.utils.logger import logger
from mantis.src.setup.installer import MantisInstaller


def setup_databases(chunk_size=None, no_taxonomy=False, mantis_config=None, cores=None):
    logger.info('Setting up databases')
    if chunk_size:
        chunk_size = int(chunk_size)
    if cores:
        cores = int(cores)
    mantis = MantisInstaller(hmm_chunk_size=chunk_size,
                       mantis_config=mantis_config,
                       user_cores=cores,
                       no_taxonomy=no_taxonomy)
    try:
        mantis.setup_databases()
    except Exception as e:
        logger.exception(e)
        raise e

def check_installation(mantis_config=None, no_taxonomy=False, check_sql=False):
    logger.info('Checking installation')
    mantis = MantisInstaller(mantis_config=mantis_config, no_taxonomy=no_taxonomy)
    mantis.check_installation(check_sql=check_sql)



def print_citation_mantis():
    paper_doi = 'https://doi.org/10.1093/gigascience/giab042'
    separator = '##########################################################################################################################'
    res = f'{separator}\n# Thank you for using Mantis, please make sure you cite the respective paper {paper_doi} #\n{separator}'
    print(res)


def print_version(user, project):
    import requests
    response = requests.get(f"https://api.github.com/repos/{user}/{project}/releases/latest")
    json = response.json()
    if 'name' in json:
        print(f'{project}\'s latest release is:', json['name'])
    else:
        print('No release available')


def run_mantis(input_path,
               output_folder,
               mantis_config=None,
               evalue_threshold=None,
               overlap_value=None,
               minimum_consensus_overlap=None,
               organism_details=None,
               genetic_code=None,
               domain_algorithm=None,
               best_combo_formula=None,
               sorting_type=None,
               keep_files=False,
               skip_consensus=False,
               skip_managed_memory=False,
               force_evalue=False,
               no_consensus_expansion=False,
               no_taxonomy=False,
               no_unifunc=False,
               kegg_matrix=False,
               verbose_kegg_matrix=False,
               output_gff=False,
               verbose=True,
               default_workers=None,
               chunk_size=None,
               time_limit=None,
               hmmer_threads=None,
               cores=None,
               memory=None,
               ):
    if evalue_threshold:
        if evalue_threshold != 'dynamic':   evalue_threshold = float(evalue_threshold)
    if overlap_value:                       overlap_value = float(overlap_value)
    if minimum_consensus_overlap:           minimum_consensus_overlap = float(minimum_consensus_overlap)
    if best_combo_formula:                  best_combo_formula = int(best_combo_formula)
    if default_workers:                     default_workers = int(default_workers)
    if chunk_size:                          chunk_size = int(chunk_size)
    if time_limit:                          time_limit = int(time_limit)
    if hmmer_threads:                       hmmer_threads = int(hmmer_threads)
    if cores:                               cores = int(cores)
    if memory:                              memory = int(memory)
    if genetic_code:                        genetic_code = int(genetic_code)
    mantis = MANTIS(
        input_path=input_path,
        output_folder=output_folder,
        mantis_config=mantis_config,
        evalue_threshold=evalue_threshold,
        overlap_value=overlap_value,
        minimum_consensus_overlap=minimum_consensus_overlap,
        organism_details=organism_details,
        genetic_code=genetic_code,
        domain_algorithm=domain_algorithm,
        best_combo_formula=best_combo_formula,
        sorting_type=sorting_type,
        keep_files=keep_files,
        skip_consensus=skip_consensus,
        skip_managed_memory=skip_managed_memory,
        force_evalue=force_evalue,
        no_consensus_expansion=no_consensus_expansion,
        no_taxonomy=no_taxonomy,
        no_unifunc=no_unifunc,
        kegg_matrix=kegg_matrix,
        verbose_kegg_matrix=verbose_kegg_matrix,
        output_gff=output_gff,
        verbose=verbose,
        default_workers=default_workers,
        chunk_size=chunk_size,
        time_limit=time_limit,
        hmmer_threads=hmmer_threads,
        user_cores=cores,
        user_memory=memory,
    )
    mantis.run_mantis()


def run_mantis_test(input_path,
                    output_folder,
                    mantis_config,
                    ):
    mantis = MANTIS(
        input_path=input_path,
        output_folder=output_folder,
        mantis_config=mantis_config,
        keep_files=True)
    mantis.run_mantis_test()