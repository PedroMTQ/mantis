import argparse
import os
from datetime import datetime
import sys

from source.MANTIS import run_mantis
from source.MANTIS_Assembler import add_slash,\
                                    splitter,\
                                    get_path_level,\
                                    check_installation,\
                                    setup_databases,\
                                    merge_hmm_folder


if __name__ == '__main__':
    print('Executing command:\n',' '.join(sys.argv))
    parser = argparse.ArgumentParser(description='MANTIS is a simple standalone workflow to annotate with HMMER.\n'+
                                                 'The main goal of this tool is to both be easy to setup (minimal requirements and automatic download and compilation of several HMMs sources) as well as easy to customize (possibility to include additional HMMs).\n'+
                                                 'After running, the hits (HMMER\'s domtblout) will be analyzed and processed. The processing of hits ensures linking of hits with an HMM\'s metadata; ensuring the results are more easily integrated and interpreted by the end user (human or machine).\n'+
                                                 'This tool takes as input an aminoacids sequences fasta; if you are starting out with a genome, you can use gene prediction tools (e.g. prodigal) to "convert" it into a protein fasta.'
                                     ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('execution_type',help='Please choose from :\n\trun_mantis\n\tsetup_databases\n\tmerge_hmm_folder\n\tcheck_installation\n\n'+
                             'If this is your first time running this software, please run <setup_databases> to download and unzip the necessary files.\n'
                             'If you have custom hmms, please include them in the <custom_hmms> folder.\n'+
                             'If your custom hmms are split 1 file/1 hmm please use <merge_hmm_folder> followed by the hmm folder path. These will be automatically pressed\n'+
                             'Custom hmms need to be pressed, to do so just run HMMER\'s hmmpress.'+
                             'To check recognized hmms please run <check_installation>\n\n'+
                             'For a single protein fasta annotation, use <run_mantis> with a fasta as target.\n'+
                             'If you have a taxonomic classification of this sample, include <-od> followed by the organism name or NCBI taxon ID\n'+
                             'For serialized protein fastas annotations, use <run_mantis>, with a tsv file path.\n'+
                             'This file should have the following structure:\n'+
                             '\tQuery name\tQuery path\tOrganism details\n'+
                             '\tquery_name_1\ttarget_path_1\t561\n'+
                             '\tquery_name_2\ttarget_path_2\tProteobacteria\n'+
                             '\tquery_name_3\ttarget_path_3\t\n'+
                             '\tquery_name_4\ttarget_path_4\tEscherichia coli\n',
                        choices=['run_mantis','setup_databases','merge_hmm_folder','check_installation'])
    parser.add_argument('-t','--target',help='Please provide the target file path')
    parser.add_argument('-o','--output_folder',help='Please provide the output folder path')
    parser.add_argument('-mc','--mantis_config',help='If you would like to provide your own custom MANTIS.config file.')
    parser.add_argument('-et','--evalue_threshold',help='Please define the e-value threshold! Default is 1e-9')
    parser.add_argument('-ar','--acceptable_range',help='Please define the acceptable range when getting the best hits! Default is 0.05')
    parser.add_argument('-ov','--overlap_value',help='Please define the allowed overlap between hits! Default is 0.1. Max is 0.3')
    parser.add_argument('-da','--domain_algorithm',choices=['exact','approximation','lowest'],help='Please choose one how multiple domains should be processed. Default is exact')
    parser.add_argument('-od','--organism_details',help='If your target fasta has been taxonimically classified please introduce details.\nTwo formats are allowed:\n'
                                                        '\ttaxon name, e.g. "Proteobacteria" or "Escherichia coli"\n'
                                                        '\tNCBI taxon ID, e.g.: 561 for Escherichia coli\n')
    parser.add_argument('-f','--force_download',action='store_true',help='If you would like to force the download of the databases when running <setup_databases>')
    parser.add_argument('-k','--keep_files',action='store_true',help='If you would like to keep intermediary output files')
    parser.add_argument('-thmm','--target_hmm',help='If you would like to restrict your sequence search to a single hmm. For the most complete annotation, don\'t change this!')
    parser.add_argument('-dw','--default_workers',help='If you would like to set the number of virtual workers used by Mantis. This is different from the physical <cores>')
    parser.add_argument('-cs','--chunk_size',help='If you would like to set the chunk size when running Mantis')
    parser.add_argument('-ht','--hmmer_threads',help='If you would like to set the number of threads used by HMMER.')
    parser.add_argument('-c','--cores',help='If you would like to set the number of physical cores used by Mantis. This is the different from the virtual <default_workers>')


    args=parser.parse_args()
    #if no input is given , arg is of class None. If it's a store_true or store_false , arg is bool
    # otherwise it's a str
    if args.execution_type=='run_mantis':
        target_path=args.target
        output_folder = args.output_folder
        mantis_config = args.mantis_config
        evalue_threshold = args.evalue_threshold
        acceptable_range = args.acceptable_range
        overlap_value = args.overlap_value
        organism_details = args.organism_details
        domain_algorithm = args.domain_algorithm
        keep_files = args.keep_files
        target_hmm = args.target_hmm
        default_workers = args.default_workers
        chunk_size = args.chunk_size
        hmmer_threads = args.hmmer_threads
        cores = args.cores
        if target_path:
            if os.path.exists(target_path):
                if not output_folder:
                    datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                    output_folder=add_slash(os.getcwd())+ get_path_level(target_path,remove_extension=True)+'_'+datetime_str
                    print('No output folder provided! Saving data to: '+output_folder)
                output_folder=add_slash(output_folder)

                run_mantis(target_path=target_path,
                                        output_folder=output_folder,
                                        mantis_config=mantis_config,
                                        evalue_threshold=evalue_threshold,
                                        acceptable_range=acceptable_range,
                                        overlap_value=overlap_value,
                                        organism_details=organism_details,
                                        domain_algorithm=domain_algorithm,
                                        keep_files=keep_files,
                                        target_hmm=target_hmm,
                                        default_workers=default_workers,
                                        chunk_size=chunk_size,
                                        hmmer_threads=hmmer_threads,
                                        cores=cores,
                                        )

            else:
                print('Target path not found, quitting now!')
        else:
            print("Missing target, quitting now!")
    elif args.execution_type=='setup_databases':
        force_download=args.force_download
        chunk_size = args.chunk_size
        setup_databases(force_download=force_download,chunk_size=chunk_size)
    elif args.execution_type=='merge_hmm_folder':
        target=args.target
        merge_hmm_folder(target_folder=target)
    elif args.execution_type=='check_installation':
        mantis_config = args.mantis_config
        check_installation(mantis_config=mantis_config)