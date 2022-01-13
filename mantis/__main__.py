try:
    import argparse
    import os
    from datetime import datetime
    import sys
    import uuid
    from mantis.MANTIS import run_mantis, run_mantis_test,print_citation_mantis,print_version
    from mantis.UniFunc_wrapper import test_nlp
    from mantis.Assembler import add_slash, \
        get_path_level, \
        check_installation, \
        setup_databases
    from mantis.utils import MANTIS_FOLDER,SPLITTER

except ImportError as e:
    import signal
    master_pid = os.getpid()
    print('Import Error:\n',e)
    os.kill(master_pid, signal.SIGKILL)

def main():
    default_config_path=f'{MANTIS_FOLDER}config{SPLITTER}MANTIS.cfg'
    print('Executing command:\n', ' '.join(sys.argv))
    parser = argparse.ArgumentParser(description='___  ___               _    _      \n'
                                                 '|  \\/  |              | |  (_)     \n'
                                                 '| .  . |  __ _  _ __  | |_  _  ___ \n'
                                                 '| |\\/| | / _` || \'_ \\ | __|| |/ __|\n'
                                                 '| |  | || (_| || | | || |_ | |\\__ \\\n'
                                                 '\\_|  |_/ \\__,_||_| |_| \\__||_||___/, a consensus driven protein function annotation tool\n'
                                                 'Documentation is available at https://github.com/PedroMTQ/mantis/wiki'
                                     , formatter_class=argparse.RawTextHelpFormatter)
    # run mantis
    parser.add_argument('execution_type',
                        help='[required]\tExecution mode',
                        choices=['run', 'setup', 'check', 'run_test','citation','version','test_nlp', 'check_sql'])
    parser.add_argument('-i', '--input',
                        help='[required]\tInput file path. Required when using <run>.')
    parser.add_argument('-o', '--output_folder',
                        help='[optional]\tOutput folder path')
    parser.add_argument('-mc', '--mantis_config',
                        help=f'Custom MANTIS.cfg file. Default is in:{default_config_path}')
    parser.add_argument('-et', '--evalue_threshold',
                        help='[optional]\tCustom e-value threshold. Default is 1e-3.')
    parser.add_argument('-ov', '--overlap_value',
                        help='[optional]\tcustom value for the allowed overlap between hits! Default is 0.1, maximum is 0.3')
    parser.add_argument('-mco', '--minimum_consensus_overlap',
                        help='[optional]\tcustom value for the minimum overlap between hits when generating the consensus annotation. Default is 0.7, 0 to accept any consistent hit, regardless of coordinate overlap.')
    parser.add_argument('-da', '--domain_algorithm', choices=['dfs', 'heuristic', 'bpo'],
                        help='[optional]\tChoose how multiple domains should be processed. Default is dfs.')
    parser.add_argument('-tl', '--time_limit',
                        help='[optional]\ttime limit in seconds when running Mantis\' DFS algorithm. Default is 60 seconds')
    parser.add_argument('-od', '--organism_details',
                        help='[optional]\tIf your input fasta has been taxonimically classified please introduce details.\n')
    parser.add_argument('-gc', '--genetic_code',
                        help='[optional]\tIf you want Mantis to translate your input fasta, please provide a genetic code. Default is 11.')
    parser.add_argument('-k', '--keep_files', action='store_true',
                        help='[optional]\tKeep intermediary output files')
    parser.add_argument('-sc', '--skip_consensus', action='store_true',
                        help='[optional]\tSkip consensus generation.')
    parser.add_argument('-nuf', '--no_unifunc', action='store_true',
                        help='[optional]\tdo not use UniFunc for descriptions similarity analysis during consensus generation.')
    parser.add_argument('-nce', '--no_consensus_expansion', action='store_true',
                        help='[optional]\tdo not expand hits during consensus generation.')
    parser.add_argument('-notax', '--no_taxonomy', action='store_true',
                        help='[optional]\tDo not download and use taxonomy databases.')
    parser.add_argument('-km', '--kegg_matrix', action='store_true',
                        help='[optional]\tgenerate KEGG modules completeness matrix.')
    parser.add_argument('-vkm', '--verbose_kegg_matrix', action='store_true',
                        help='[optional]\tgenerate KEGG modules completeness matrix in verbose mode. Verbose mode additionally outputs the complete module name and missing KOs.')
    parser.add_argument('-gff', '--output_gff', action='store_true',
                        help='[optional]\tgenerate GFF-formatted output files.')

    parser.add_argument('-fo', '--force_output', action='store_true',
                        help='[optional]\tIf you would like to force the output to the folder you specified.')

    # general args
    parser.add_argument('-dw', '--default_workers',
                        help='[optional]\tnumber of virtual workers used by Mantis. This is different from the physical <cores>. By default, the number of workers is the same as <cores>.')
    parser.add_argument('-cs', '--chunk_size',
                        help='[optional]\tchunk size when running Mantis')
    parser.add_argument('-ht', '--hmmer_threads',
                        help='[optional]\tnumber of threads used by HMMER. Default is 1.')
    parser.add_argument('-c', '--cores',
                        help='[optional]\tset the number of physical cores used by Mantis (all are used by default).')
    parser.add_argument('-m', '--memory',
                        help='[optional]\tset the amount of RAM used by Mantis (in GB) (all is used by default).')

    # developers only / testing tools
    parser.add_argument('-bcf', '--best_combo_formula', choices=[str(i) for i in range(1, 13)],
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
    if args.execution_type == 'run':
        input_path = args.input
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
        no_taxonomy = args.no_taxonomy
        no_unifunc = args.no_unifunc
        kegg_matrix = args.kegg_matrix
        verbose_kegg_matrix = args.verbose_kegg_matrix
        output_gff = args.output_gff
        force_output = args.force_output
        default_workers = args.default_workers
        chunk_size = args.chunk_size
        time_limit = args.time_limit
        hmmer_threads = args.hmmer_threads
        cores = args.cores
        memory = args.memory
        if input_path:
            if os.path.exists(input_path):
                if not output_folder:
                    datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                    output_folder = add_slash(os.getcwd()) + get_path_level(input_path,
                                                                            remove_extension=True) + '_' + datetime_str
                    print(f'No output folder provided! Saving data to: {output_folder}')
                if os.path.exists(output_folder):
                    if not force_output and os.listdir(output_folder):
                        datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                        hex_random = '_hex_' + uuid.uuid4().hex[:10]
                        output_folder += '_' + datetime_str + hex_random
                        print(
                            f'The output folder already contains something! New output folder will be: {output_folder}')
                output_folder = add_slash(output_folder)

                run_mantis(input_path=input_path,
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
                           default_workers=default_workers,
                           chunk_size=chunk_size,
                           time_limit=time_limit,
                           hmmer_threads=hmmer_threads,
                           cores=cores,
                           memory=memory,
                           )
                print_citation_mantis()

            else:
                print('Input path not found, quitting now!')
        else:
            print("Missing input file, quitting now!")

    elif args.execution_type == 'setup':
        mantis_config = args.mantis_config
        chunk_size = args.chunk_size
        no_taxonomy = args.no_taxonomy
        cores = args.cores
        setup_databases(chunk_size=chunk_size, no_taxonomy=no_taxonomy,
                        mantis_config=mantis_config, cores=cores)
        print_citation_mantis()
    elif args.execution_type == 'check':
        mantis_config = args.mantis_config
        no_taxonomy = args.no_taxonomy
        check_installation(mantis_config=mantis_config, no_taxonomy=no_taxonomy)

    elif args.execution_type == 'run_test':
        output_folder = args.output_folder
        if not output_folder:
            output_folder = add_slash(os.getcwd()) + 'test_run'
            print(f'No output folder provided! Saving data to: {output_folder}')
        if os.path.exists(output_folder):
            if os.listdir(output_folder):
                datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                hex_random = '_hex_' + uuid.uuid4().hex[:10]
                output_folder += '_' + datetime_str + hex_random
                print(f'The output folder already contains something! New output folder will be: {output_folder}')
        output_folder = add_slash(output_folder)
        run_mantis_test(input_path=add_slash(MANTIS_FOLDER + 'tests') + 'test_sample.faa',
                        output_folder=output_folder,
                        mantis_config=add_slash(MANTIS_FOLDER + 'tests') + 'test_MANTIS.cfg',
                        )
        print_citation_mantis()

    elif args.execution_type == 'citation':
        print_citation_mantis()

    elif args.execution_type == 'version':
        print_version('pedromtq', 'mantis')

    elif args.execution_type == 'test_nlp':
        test_nlp()
    elif args.execution_type == 'check_sql':
        mantis_config = args.mantis_config
        no_taxonomy = args.no_taxonomy
        check_installation(mantis_config=mantis_config, check_sql=True, no_taxonomy=no_taxonomy)



if __name__ == '__main__':
    main()