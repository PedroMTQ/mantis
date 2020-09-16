
# Functionalities

To run this pipeline don't call the py files (**calling the code from inside the folder won't work**) instead call on the project folder, like so:
````
python  mantis/  "function"
````
Documentation is also available via console:
````
python  mantis/  -h
````
## Setup databases:
````
python  mantis/  setup_databases

Optional arguments: --force_download / -f
                    --chunk_size / -cs
````
This method will download and unzip all the default data into their respective path folders. Please don't move any of this data around during execution!  
This execution mode takes into account the existing data it has downloaded and generated, if for some reason you'd like to download all the data again, use  `-f` or `--force_download` .  
You can use `-cs` or `--chunk_size` to set the amount of HMM profiles per chunk. By default each HMM chunk will have 5000 profiles, this allows better throughput as it ensures resources saturation.

**Note:** It's usually good practice to check your installation with `python  mantis/  check_installation`, if something has failed you should run `python  mantis/  setup_databases` again.



## Check installation:
````
python  mantis/  check_installation
````
This method will check the paths from the **MANTIS.config** and check the installed custom HMMs.

## Merge HMM folder:
````
python  mantis/  merge_hmm_folder -t target

Mandatory arguments: --target / -t
````

This method will concatenate all the **.hmms** files in a specified target folder **target**. This is useful when you have downloaded a lot of different HMM profiles from the same data source and they are all split into one HMM/file . This can not only cause issues in storage (due to limits in number of files) but it also reduces efficiency with HMMER due to the lack of proper HMM indexation.

## Annotate one sample:
````
python mantis/ run_mantis -t target.faa -o output_folder -od organism_details

Mandatory arguments: --target / -t
Optional arguments:  --output_folder / -o
                     --organism_details / -od
                     --evalue_threshold / -et
                     --acceptable_range / -ar
                     --overlap_value / -ov
                     --target_hmm / -thmm
                     --delete_files / -d
                     --mantis_config / -mc
                     --default_workers / -dw
                     --chunk_size / -cs
                     --hmmer_threads / -ht
````


- **target** Mantis will run on the `target` fasta path.  
- **output_folder** If no output folder is provided, data will be saved to the `current_path/target_date_time`.  
- **organism_details** this variable is used for determining which hmms to use during taxonomic specific annotation. It accepts an NCBI ID or an organism name; if a string contains a blank space please include the string within quotes (e.g. "genus species"). Any taxon level can be provided. If none are provided, only the general HMMs are used.   
- **evalue_threshold** see [what is the default e-value?](https://github.com/PedroMTQ/mantis/wiki/Additional-information#what-is-the-e-value-threshold)  
- **acceptable_range**  see [wow are hits generated?](https://github.com/PedroMTQ/mantis/wiki/Additional-information#how-are-hits-generated)  
- **overlap_value** If you would like to allow partial overlap between hits, please use `--overlap_value <overlap_value>` or `-ov <overlap_value>`. Default is 0.1, maximum is 0.3.  
- **target_hmm** If you would like to annotate your sample with only one hmm, please use `--target_hmm <absolute path to hmm>` or `-thmm <absolute path to hmm>`. This is mostly used for testing purposes.  
- **delete_files** use this option to delete extra output files  
- **mantis_config** use this option to use a custom MANTIS.config file  
- **default_workers** use this to set the amount of processes that will run HMMER (calculated automatically).This is mostly used for testing purposes. If you set it to 1, Mantis will run HMMER without parallelization.  
- **chunk_size** use this to set the size of the chunks that your sample files will be divided in. 1000 sequences per chunk by default.  
- **hmmer_threads** use this to set the number of threads used by HMMER. 5 by default.






***Example***
````
python mantis run_mantis -t mantis/tests/test_sample.faa -od "Escherichia coli"
````
## Annotate multiple samples:
````
python mantis/ run_mantis -t target.tsv -o output_folder

Mandatory arguments: --target / -t
Optional arguments:  --output_folder / -o
                     --evalue_threshold / -et
                     --acceptable_range / -ar
                     --overlap_value / -ov
                     --target_hmm / -thmm
                     --delete_files / -d
                     --mantis_config / -mc
                     --default_workers / -dw
                     --chunk_size / -cs
                     --hmmer_threads / -ht
````

- **target** Mantis will parse `target` and annotate the respective fasta file paths.  			 
- **output_folder** If no output folder is provided, data will be saved to the `current_path/target_date_time/query_name`.  
- **evalue_threshold** see [what is the default e-value?](https://github.com/PedroMTQ/mantis/wiki/Additional-information#what-is-the-e-value-threshold)  
- **acceptable_range**  see [wow are hits generated?](https://github.com/PedroMTQ/mantis/wiki/Additional-information#how-are-hits-generated)  
- **overlap_value** If you would like to allow partial overlap between hits, please use `--overlap_value <overlap_value>` or `-ov <overlap_value>`. Default is 0.1, maximum is 0.3.  
- **target_hmm** If you would like to annotate your sample with only one hmm, please use `--target_hmm <absolute path to hmm>` or `-thmm <absolute path to hmm>`. This is mostly used for testing purposes.  
- **delete_files** use this option to delete extra output files  
- **mantis_config** use this option to use a custom MANTIS.config file  
- **default_workers** use this to set the amount of processes that will run HMMER (calculated automatically).This is mostly used for testing purposes. If you set it to 1, Mantis will run HMMER without parallelization.  
- **chunk_size** use this to set the size of the chunks that your sample files will be divided in. 1000 sequences per chunk by default.  
- **hmmer_threads** use this to set the number of threads used by HMMER. 5 by default.


The `target` tsv file should have the following format:

| Query name  | Absolute sample path | Organism details |
| :-------------: | :-------------: | :-------------: |
| query_name_1  | target_path_1  | 561  |
| query_name_2  | target_path_2  | Proteobacteria  |
| query_name_3  | target_path_3  |   |
| query_name_4  | target_path_4  | Escherichia coli  |

An example file is provided `example_file.tsv`. The query name and the sample path are mandatory, the organism details column is optional.

## Other input formats

Mantis also accepts directory paths or compressed files (`.gz`,`.zip`,`.tar.gz`). Samples will be uncompressed and deleted after execution. Keep in mind that with this method it's not possible to input the taxonomical classification of each sample; for taxa-resolved annotations use the previous input methods.

## Annotating metagenomes

Mantis scales well with Metagenomes, since it automatically splits fasta files into evenly sized chunks, ensuring parallelization without the potential idle time you'd get due to iterating over sample sequences or the differently sized HMM profiles references.
