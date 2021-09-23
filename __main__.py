try:
    import argparse
    import os
    from datetime import datetime
    import sys
    import uuid

    from source.MANTIS import run_mantis, run_mantis_test,print_citation_mantis
    from source.MANTIS_NLP import test_nlp
    from source.MANTIS_Assembler import add_slash, \
        get_path_level, \
        check_installation, \
        extract_nog_metadata, \
        setup_databases, \
        merge_hmm_folder
    from source.utils import MANTIS_FOLDER
except ImportError as e:
    import signal
    master_pid = os.getpid()
    print('Import Error!')
    os.kill(master_pid, signal.SIGKILL)


if __name__ == '__main__':
    print('Executing command:\n', ' '.join(sys.argv))
    parser = argparse.ArgumentParser(description='___  ___               _    _      \n'
                                                 '|  \\/  |              | |  (_)     \n'
                                                 '| .  . |  __ _  _ __  | |_  _  ___ \n'
                                                 '| |\\/| | / _` || \'_ \\ | __|| |/ __|\n'
                                                 '| |  | || (_| || | | || |_ | |\\__ \\\n'
                                                 '\\_|  |_/ \\__,_||_| |_| \\__||_||___/, a consensus driven protein function annotation tool\n'
                                     , formatter_class=argparse.RawTextHelpFormatter)
    #run mantis
    parser.add_argument('execution_type',
                        help='Please choose from :\n\trun_mantis\n\tsetup_databases\n\tmerge_hmm_folder\n\textract_nog_metadata\n\tcheck_installation\n\trun_test\n\n' +
                             'If this is your first time running this software, please run <setup_databases> to download and unzip the necessary files.\n'
                             'If you have custom hmms, please include them in the <custom_hmms> folder.\n' +
                             'If your custom hmms are split 1 file/1 hmm please use <merge_hmm_folder> followed by the hmm folder path. These will be automatically pressed\n' +
                             'Custom hmms need to be pressed, to do so just run HMMER\'s hmmpress.' +
                             'To check recognized hmms please run <check_installation>\n\n' +
                             'If you have a taxonomic classification of this sample, include <-od> followed by the organism name or NCBI taxon ID\n' +
                             'For multiple protein fastas annotations, use <run_mantis>, with a tsv file path.\n' +
                             'This file should have the following structure:\n' +
                             '\tQuery name\tQuery path\tOrganism details\n' +
                             '\tquery_name_1\ttarget_path_1\t561\n' +
                             '\tquery_name_2\ttarget_path_2\tProteobacteria\n' +
                             '\tquery_name_3\ttarget_path_3\t\n' +
                             '\tquery_name_4\ttarget_path_4\tEscherichia coli\n',
                        choices=['run_mantis', 'setup_databases', 'merge_hmm_folder', 'check_installation', 'run_test',
                                 'extract_nog_metadata','test_nlp','citation'])
    parser.add_argument('-t', '--target', help='[required]\tAnnotation target file path. Required when using <run_mantis>.')
    parser.add_argument('-o', '--output_folder', help='[optional]\tOutput folder path')
    parser.add_argument('-mc', '--mantis_config',
                        help='Custom MANTIS.config file. Default is in Mantis\' folder')
    parser.add_argument('-et', '--evalue_threshold',
                        help='[optional]\tCustom e-value threshold. Default is 1e-3. You can use <dynamic> to take into account sequence length.')
    parser.add_argument('-ov', '--overlap_value',
                        help='[optional]\tcustom value for the allowed overlap between hits! Default is 0.1, maximum is 0.3')
    parser.add_argument('-mco', '--minimum_consensus_overlap',
                        help='[optional]\tcustom value for the minimum overlap between hits when generating the consensus annotation. Default is 0.7, 0 to accept any consistent hit, regardless of coordinate overlap.')
    parser.add_argument('-da', '--domain_algorithm', choices=['dfs', 'heuristic', 'bpo'],
                        help='[optional]\tChoose how multiple domains should be processed. Default is dfs, more information on the algorithms in the wiki page.')
    parser.add_argument('-tl', '--time_limit',
                        help='[optional]\ttime limit in seconds when running Mantis\' DFS algorithm. Default is 60 seconds')
    parser.add_argument('-od', '--organism_details',
                        help='[optional]\tIf your target fasta has been taxonimically classified please introduce details.\n'
                             '\t\tTwo formats are allowed:\n'
                             '\t\t\ttaxon name, e.g. "Proteobacteria" or "Escherichia coli"\n'
                             '\t\t\tNCBI taxon ID, e.g.: 561 for Escherichia coli\n'
                             'Providing NCBI IDs is faster and safer.')
    parser.add_argument('-gc', '--genetic_code',
                        help='[optional]\tIf you want Mantis to translate your target fasta, please provide a genetic code. Default is 11. \n'
                             '\t\tFor further details please see https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes\n')
    parser.add_argument('-k', '--keep_files', action='store_true',
                        help='[optional]\tKeep intermediary output files')
    parser.add_argument('-sc', '--skip_consensus', action='store_true',
                        help='[optional]\tSkip consensus generation.')
    parser.add_argument('-nuf', '--no_unifunc', action='store_true',
                        help='[optional]\tdo not use UniFunc for similarity analysis during consensus generation.')
    parser.add_argument('-nce', '--no_consensus_expansion', action='store_true',
                        help='[optional]\tdo not expand hits during consensus generation.')
    parser.add_argument('-km', '--kegg_matrix', action='store_true',
                        help='[optional]\tgenerate KEGG modules completeness matrix.')
    parser.add_argument('-vkm', '--verbose_kegg_matrix', action='store_true',
                        help='[optional]\tgenerate KEGG modules completeness matrix in verbose mode. Verbose mode gives, in addition to the default matrix, complete module name and missing KOs; it also exports a summary figure.')
    parser.add_argument('-fo', '--force_output', action='store_true',
                        help='[optional]\tIf you would like to force the output to the folder you specified. This may result in errrors!')
    #setup databases


    parser.add_argument('-f', '--force_download', action='store_true',
                        help='[optional]\tIf you would like to force the download of the databases when running <setup_databases>')


    #general args
    parser.add_argument('-dw', '--default_workers',
                        help='[optional]\tnumber of virtual workers used by Mantis. This is different from the physical <cores>. Default number of workers corresponds to the number of physical cores.')
    parser.add_argument('-cs', '--chunk_size', help='[optional]\tchunk size when running Mantis')
    parser.add_argument('-ht', '--hmmer_threads', help='[optional]\tnumber of threads used by HMMER. Default is 1.')
    parser.add_argument('-c', '--cores',
                        help='[optional]\tset the number of physical cores used by Mantis. Mantis uses all available physical cores by default.')
    parser.add_argument('-m', '--memory',
                        help='[optional]\tset the amount of RAM used by Mantis (in GB). Mantis uses all available RAM by default.')

    #developers only / testing tools
    parser.add_argument('-bcf', '--best_combo_formula', choices=[str(i) for i in range(1,13)],
                        help='[developers_only]\tChoose which scoring formula to use to determine the best combination.')
    parser.add_argument('-st', '--sorting_type', choices=['evalue', 'bitscore'],
                        help='[developers_only]\tPlease choose the score to sort hits.')
    parser.add_argument('-smm', '--skip_managed_memory', action='store_true',
                        help='[developers_only]\tskip memory management. No HMMER memory management (less stable), but may allow for runs to finish in low memory environments')
    parser.add_argument('-fe', '--force_evalue', action='store_true',
                        help='[developers_only]\tE-value is scaled to chunk size to avoid introduction of FPs, with this you can force it to be independent of chunk size. Use at your own risk.')




    args = parser.parse_args()
    # if no input is given , arg is of class None. If it's a store_true or store_false , arg is bool
    # otherwise it's a str
    if args.execution_type == 'run_mantis':
        target_path = args.target
        output_folder = args.output_folder
        mantis_config = args.mantis_config
        evalue_threshold = args.evalue_threshold
        overlap_value = args.overlap_value
        minimum_consensus_overlap = args.minimum_consensus_overlap
        organism_details = args.organism_details
        genetic_code = args.genetic_code
        domain_algorithm = args.domain_algorithm
        best_combo_formula = args.best_combo_formula
        sorting_type = args.sorting_type
        keep_files = args.keep_files
        skip_consensus = args.skip_consensus
        skip_managed_memory = args.skip_managed_memory
        force_evalue = args.force_evalue
        no_consensus_expansion = args.no_consensus_expansion
        no_unifunc = args.no_unifunc
        kegg_matrix = args.kegg_matrix
        verbose_kegg_matrix = args.verbose_kegg_matrix
        force_output = args.force_output
        default_workers = args.default_workers
        chunk_size = args.chunk_size
        time_limit = args.time_limit
        hmmer_threads = args.hmmer_threads
        cores = args.cores
        memory = args.memory
        if target_path:
            if os.path.exists(target_path):
                if not output_folder:
                    datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                    output_folder = add_slash(os.getcwd()) + get_path_level(target_path,remove_extension=True) + '_' + datetime_str
                    print(f'No output folder provided! Saving data to: {output_folder}')
                if os.path.exists(output_folder):
                    if not force_output and os.listdir(output_folder):
                        datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                        hex_random= '_hex_'+uuid.uuid4().hex[:10]
                        output_folder += '_' + datetime_str+hex_random
                        print(f'The output folder already contains something! New output folder will be: {output_folder}')
                output_folder = add_slash(output_folder)

                run_mantis(target_path=target_path,
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
                           no_unifunc=no_unifunc,
                           kegg_matrix=kegg_matrix,
                           verbose_kegg_matrix=verbose_kegg_matrix,
                           default_workers=default_workers,
                           chunk_size=chunk_size,
                           time_limit=time_limit,
                           hmmer_threads=hmmer_threads,
                           cores=cores,
                           memory=memory,
                           )
                print_citation_mantis()

            else:
                print('Target path not found, quitting now!')
        else:
            print("Missing target, quitting now!")
    elif args.execution_type == 'setup_databases':
        mantis_config = args.mantis_config
        force_download = args.force_download
        chunk_size = args.chunk_size
        cores = args.cores
        setup_databases(force_download=force_download, chunk_size=chunk_size, mantis_config=mantis_config,cores=cores)
        print_citation_mantis()
    elif args.execution_type == 'citation':
        print_citation_mantis()
    elif args.execution_type == 'merge_hmm_folder':
        target = args.target
        merge_hmm_folder(target_folder=target)
        print_citation_mantis()
    elif args.execution_type == 'check_installation':
        mantis_config = args.mantis_config
        check_installation(mantis_config=mantis_config)
    elif args.execution_type == 'extract_nog_metadata':
        output_folder = args.output_folder
        if not output_folder:
            output_folder = add_slash(os.getcwd()) + 'metadata_extraction'
            print(f'No output folder provided! Saving data to: {output_folder}')
        output_folder = add_slash(output_folder)
        extract_nog_metadata(metadata_path=output_folder)
        print_citation_mantis()

    elif args.execution_type == 'test_nlp':
        test_nlp()
    elif args.execution_type == 'run_test':
        output_folder = args.output_folder
        if not output_folder:
            output_folder = add_slash(os.getcwd()) + 'test_run'
            print(f'No output folder provided! Saving data to: {output_folder}')
        if os.path.exists(output_folder):
            if os.listdir(output_folder):
                datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                hex_random = '_hex_' + uuid.uuid4().hex[:10]
                output_folder += '_' + datetime_str+hex_random
                print(f'The output folder already contains something! New output folder will be: {output_folder}')
        output_folder = add_slash(output_folder)
        run_mantis_test(target_path=add_slash(MANTIS_FOLDER + 'tests')+ 'test_sample.faa',
                        output_folder=output_folder,
                        mantis_config=add_slash(MANTIS_FOLDER + 'tests')+ 'test_MANTIS.config',
                        )
        print_citation_mantis()

