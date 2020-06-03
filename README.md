# MANTIS

![mantis_icon](Images/mantis_icon_small.png)

This tool can be used for protein function annotation, it is a standalone tool that uses HMMER to match sequences against multiple reference datasets. It accepts as input an aminoacids sequence fasta.  
The main goals of this tool are to:
- consider multiple protein domains
- annotate with taxonomy resolution
- use different reference datasets and provide a consensus annotation
- be easy to setup and/or customize
- scale well with multiple samples and/or metagenomes


If you have only loose reads, you need to assemble them first; when you have assembled reads/genomes you need to predict the protein coding regions (gene prediction - e.g. prodigal) to convert your data into a protein fasta that Mantis can then use.

**Mantis is compatible with genomes and metagenomes.**

- [Requirements](#requirements)
- [Quick configuration](#quick-configuration)
- [Functions](#functions)
- [Further details](#further-details)




# Straight to the point
### Requirements

- Python>=v3   
- Conda  
- HMMER>=v3

### Quick configuration
1. `git clone ssh://git@git-r3lab-server.uni.lu:8022/pedro.queiros/mantis.git`  
2. Run `conda env create -f mantis_env.yml`
3. Run `conda activate mantis_env`
4. Run `python mantis setup_databases`
5. Run `python mantis run_mantis -t target_faa`

**Custom hmms**  

        custom_hmms_folder=/path/to/mantis/hmm/custom_hmms/   
        custom_hmm=/path/to/HMM_folder/file.hmm

Custom hmms can be added in **MANTIS.config** by adding their absolute path, alternatively you may add them to the **custom_hmms** folder.
Mantis will read the folders within the custom hmms folder and use the .hmm stored in each of those folders.  

### Functions

**1. Help**  
````
python  mantis/  -h
````
**2. Setup databases**  
````
python  mantis/  setup_databases
````

**3. Check installation**  
````
python  mantis/  check_installation
````
**4. Merge HMM folder**  
````
python  mantis/  merge_hmm_folder -t target
````
**5. Annotate one sample**  
````
python mantis/ run_mantis -t target.faa -o output_folder-od organism_details -et evalue_threshold -ar acceptable_range -ov overlap_value -mc custom_MANTIS.config    
````
*example*: `python mantis run_mantis -t mantis/tests/test_sample.faa -od "Escherichia coli"`

**6. Annotate multiple samples**  
````
python mantis/ run_mantis -t target.tsv -o output_folder -et evalue_threshold -ar acceptable_range -ov overlap_value -mc custom_MANTIS.config
````
*example*: `python mantis run_mantis -t mantis/tests/test_file.tsv`

### Output files  

There are 3 output files:
 - `output_annotation.tsv`, which has all hits and their coordinates and e-values;
 - `interpreted_annotation.tsv` which has all hits, their coordinates and e-value, as well as the respective hit metadata;
- `consensus_annotation.tsv` which has all hits and their respective metadata from the best hmm sources consensus.  

The first two files can have the same query sequence in several lines (query sequence/hmm source) while the `consensus_annotation.tsv` will only have one line per query sequence (consensus/query).

# Further details

1. [Configuration](https://github.com/PedroMTQ/mantis/wiki/Configuration)  
  1.1 [Requirements](https://github.com/PedroMTQ/mantis/wiki/Configuration#requirements)  
  1.2 [Installation](https://github.com/PedroMTQ/mantis/wiki/Configuration#installation)  
  1.2 [Conda environment](https://github.com/PedroMTQ/mantis/wiki/Configuration#conda-environment)  
  1.3 [Setting your own paths](https://github.com/PedroMTQ/mantis/wiki/Configuration#setting-your-own-paths)  
  1.4 [Custom HMMs](https://github.com/PedroMTQ/mantis/wiki/Configuration#custom-hmms)  
  1.5 [Setting HMMs weight](https://github.com/PedroMTQ/mantis/wiki/Configuration#setting-hmms-weight)  
2. [Functionalities](https://github.com/PedroMTQ/mantis/wiki/Functionalities)  
  2.1 [Setup databases](https://github.com/PedroMTQ/mantis/wiki/Functionalities#setup-databases)  
  2.2 [Check installation](https://github.com/PedroMTQ/mantis/wiki/Functionalities#check-installation)  
  2.3 [Merge HMM folder](https://github.com/PedroMTQ/mantis/wiki/Functionalities#merge-hmm-folder)  
  2.4 [Annotating one sample](https://github.com/PedroMTQ/mantis/wiki/Functionalities#annotate-one-sample)  
  2.5 [Annotating multiple samples](https://github.com/PedroMTQ/mantis/wiki/Functionalities#annotate-multiple-samples)  
  2.6 [Annotating Metagenomes](https://github.com/PedroMTQ/mantis/wiki/Functionalities#annotating-metagenomes)   
3. [Output](https://github.com/PedroMTQ/mantis/wiki/Output)  
4. [Additional information](https://github.com/PedroMTQ/mantis/wiki/Additional-information)  
  4.1 [Reference data](https://github.com/PedroMTQ/mantis/wiki/Additional-information#reference-data)  
  4.2 [What is the default e-value?](https://github.com/PedroMTQ/mantis/wiki/Additional-information#what-is-the-e-value-threshold)  
  4.3 [How are hits generated?](https://github.com/PedroMTQ/mantis/wiki/Additional-information#how-are-hits-generated)  
  4.4 [How are hit consensus generated?](https://github.com/PedroMTQ/mantis/wiki/Additional-information#how-are-hit-consensus-generated)  
  4.5 [A note on efficiency](https://github.com/PedroMTQ/mantis/wiki/Additional-information#notes-on-efficiency)  
5. [Project structure and architecture](https://github.com/PedroMTQ/mantis/wiki/Project-structure-and-architecture)  
6. [Copyright](https://github.com/PedroMTQ/mantis/wiki/Copyright)  

# Copyright

This project is available under the [MIT license](https://github.com/PedroMTQ/mantis/wiki/Copyright).

# References and acknowledgements

>1. S. R. Eddy. HMMER: biosequence analysis using profile hidden Markov models. HMMER v.3.2.1 www.hmmer.org
>2. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Jaime Huerta-Cepas, Damian Szklarczyk, Davide Heller, Ana Hernández-Plaza, Sofia K Forslund, Helen Cook, Daniel R Mende, Ivica Letunic, Thomas Rattei, Lars J Jensen, Christian von Mering, Peer Bork Nucleic Acids Res. 2019 Jan 8; 47(Database issue): D309–D314. https://doi.org/10.1093/nar/gky1085
>3. The Pfam protein families database in 2019: S. El-Gebali, J. Mistry, A. Bateman, S.R. Eddy, A. Luciani, S.C. Potter, M. Qureshi, L.J. Richardson, G.A. Salazar, A. Smart, E.L.L. Sonnhammer, L. Hirsh, L. Paladin, D. Piovesan, S.C.E. Tosatto, R.D. Finn Nucleic Acids Research (2019)  https://doi.org/10.1093/nar/gky995
>4. Haft DH, Loftus BJ, Richardson DL, et al. TIGRFAMs: a protein family resource for the functional identification of proteins. Nucleic Acids Res. 2001;29(1):41–43. https://doi.org/10.1093/nar/29.1.41
>5. Aramaki T., Blanc-Mathieu R., Endo H., Ohkubo K., Kanehisa M., Goto S., Ogata H. KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold. Bioinformatics. 2019 Nov 19. pii: btz859. https://doi.org/10.1093/bioinformatics/btz859.
>6. Han Zhang, Tanner Yohe, Le Huang, Sarah Entwistle, Peizhi Wu, Zhenglu Yang, Peter K Busk, Ying Xu, Yanbin Yin, dbCAN2: a meta server for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 46, Issue W1, 2 July 2018, Pages W95–W101, https://doi.org/10.1093/nar/gky418
>7. Yanbin Yin, Xizeng Mao, Jincai Yang, Xin Chen, Fenglou Mao, Ying Xu, dbCAN: a web resource for automated carbohydrate-active enzyme annotation, Nucleic Acids Research, Volume 40, Issue W1, 1 July 2012, Pages W445–W451, https://doi.org/10.1093/nar/gks479
>8. Gibson MK, Forsberg KJ, Dantas G. Improved annotation of antibiotic resistance functions reveals microbial resistomes cluster by ecology. The ISME Journal. 2014, https://doi.org/ISMEJ.2014.106
>9. Albertsen, M., Hugenholtz, P., Skarshewski, A. et al. Genome sequences of rare, uncultured bacteria obtained by differential coverage binning of multiple metagenomes. Nat Biotechnol 31, 533–538 (2013). https://doi.org/10.1038/nbt.2579
>10. W. Arndt, "Modifying HMMER3 to Run Efficiently on the Cori Supercomputer Using OpenMP Tasking," 2018 IEEE International Parallel and Distributed Processing Symposium Workshops (IPDPSW), Vancouver, BC, 2018, pp. 239-246. https://doi.org/10.1109/IPDPSW.2018.00048
