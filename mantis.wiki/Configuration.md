# Requirements

## Software requirements
This tool requires a **Conda** environment with the following packages:
- **Python**, tested with v3.7.3 but anything above v3 should be fine
- **requests**, tested with v2.22.0
- **numpy**, tested with v1.18.1
- **nltk**, tested with v3.4.4
- **sqlite**, tested with v3.30.1
- **psutil**, tested with 5.6.7
- [HMMER](http://hmmer.org/), tested with v3.2.1

A conda environment is available - `mantis_env.yml`.  

## How much space do we need?
The lineage annotation requires quite a lot of space since NOG's HMM database is quite extensive. For the taxonomy you will need around 2.8 terabytes. The rest of the HMMs only take up around 200 gibabytes. To check default datasets see [Reference data](https://github.com/PedroMTQ/mantis/wiki/Additional-information#reference-data)  
You don't need to use all of this data though!


# Installation

Mantis is easy to setup, simply:  

1. Clone the repository with git  
2. Edit MANTIS.config with desired paths
3. Create a Conda environment for Mantis with `conda env create -f mantis/mantis_env.yml`
4. Activate the previously created Conda environment
5. Setup all default databases for Mantis with`python mantis setup_databases`

To check your installation run:
````
python mantis check_installation
````
Keep in mind the installation will take a while as a lot of data is downloaded. If NOG's hmms are not used it can finish within a couple of hours, otherwise it may take a few days.  
To customize your installation (setting installation paths or removing certain HMMs) please refer to [configuration](#setting-your-own-paths).


# Configuration
The **MANTIS.config** allows the user to edit and add custom HMMs. An example config file is included, please use the same syntax, otherwise configuration won't be taken into account.

Mantis comes with a MANTIS.config file which serves as the default to all the users in the system. You can configure your own MANTIS.config file by copying this file and editing it as you wish. Afterwards you can just add `-mc <path/to/edited_MANTIS.config>`.  

## Conda environment  

It's preferable to use a self contained environment, avoiding compatibility issues, but you can run Mantis in whichever Conda environment you'd like, simply active it and run Mantis.     

This is not necessary, but if you'd like to share your Mantis environment across multiple users do the following:
1. Create the Mantis environment in a group folder location, by running `conda env create -f mantis_env.yml -p <path/to/group/folder/>`
Future Mantis users now need to tot hef following:
2. Run `conda config` to generate the `.condarc` file
3. Edit `.condarc` file (usually located in your root folder) and add:
````
envs_dirs:  
    - path/to/group/folder/  
````

## Setting your own paths

After running `setup_databases` you may wish to move data around, if so, make sure you change all these paths:
````
uniprot_folder=/path/to/mantis/Resources/Uniprot/  
go_obo_folder=/path/to/mantis/Resources/Gene_Ontology/  
ncbi_dmp_path_folder=/path/to/mantis/Resources/NCBI/  
default_hmms_folder=/path/to/mantis/hmm/  
NOGT_hmm_folder=/path/to/mantis/hmm/eggnogdb.embl.de/NOGT/    
NOGG_hmm_folder=/path/to/mantis/hmm/eggnogdb.embl.de/NOGG/  
pfam_hmm_folder=/path/to/mantis/hmm/pfam/  
kofam_hmm_folder=/path/to/mantis/hmm/kofam/  
dbcan_hmm_folder=/path/to/mantis/hmm/dbcan/  
tigrfam_hmm_folder=/path/to/mantis/hmm/tigrfam/  
resfams_hmm_folder=/path/to/mantis/hmm/resfams/  
````
If you don't move any of these folders, don't worry about configuring this.  
If you don't want all the hmm files to be used, you can change the path to 'NA', for example: `NOGT_hmm_folder=NA`

**Important**: All of the default hmms belong to their respective authors, I haven't compiled any of this data, I'm merely distributing it in a more automated manner! Make sure you cite them when using this tool/their data.

`NOGT` is the collection of taxon specific HMMs, `NOGG` the collection of all HMMs.

## Custom hmms
````
custom_hmms_folder=/path/to/mantis/hmm/custom_hmms/  
custom_hmm=/path/to/HMM_folder/file.hmm
````
Custom hmms can be added in **MANTIS.config** by adding their absolute path, alternatively you may add them to the **custom_hmms** folder.
This tool will read the folders within the custom hmms folder and use the .hmm stored in each of those folders.  

**Important:**  
Remember to use HMMER's **hmmpress** on the custom hmms!  
If custom hmms are divided 1 hmm/hmm file make sure you merge them together using the `merge_hmm_folder` function.    
If hmms from the same source are not merged, hits processing won't take into account potential hmm hits overlaps.  

Most metadata is formatted differently, therefore, for custom hmms this tool requires the metadata to be formatted in a specific manner, otherwise only the hmm name will be extracted as "metadata".
To see an example please go to `hmm/custom_hmms/` where you will find two files `custom.hmm` and `custom.tsv`.  
In the `custom.tsv` you can see how the metadata should be formatted.
In the first column there should be the HMM name, in the columns that come after any kind of metadata can be added. To specify the type of metadata simply add the type to the headers of the `.tsv` file. Columns without any headers will be assumed to be a free-text description. Some identifiers will still be searched for in this free text (EC, KO, TCDB, DUF, GO, and COG).
For the custom metadata to be recognized please place the custom metadata in the same folder as the custom hmm file and use the same name but with a `.tsv` extension, for example: `path/to/custom_hmm/custom.hmm` and `path/to/custom_hmm/custom.tsv`.  
The metadata tsv files should have the following format:

| HMM_profile  | Metadata_type_1 | Metadata_type_2 | Metadata_type_3 |
| :-------------: | :-------------: | :-------------: | :-------------: |
| HMM_1  | 2.1.15.64  |  | this is a description  |
| HMM_1  | 3.2.9.13  | KO0002  | this is a description  |

An experienced user can add their linking method to the `MANTIS_Interpreter.py` . You can also post an issue in this repository and I can try to write a linking method for your custom hmm.

## Setting HMMs weight

When generating the consensus, some HMMs can be given more weight, this is important because some HMMs are more specific than others.
By default NOG has the most weight since its HMMs are specific to taxon.  
To configure the weight of an HMM simply change the MANTIS.config file:  
- example: **NOGT**_hmm_folder should be `NOGT_weight=X` where X is the weight of the HMM (0-1)  
- example: custom_hmm=path/to/**customHMM**.hmm should be `customHMM_weight=X` where X is the weight of the HMM (0-1)

In essence make sure the  names of the weights correspond to the path of the HMMs.  
If no weight is given to the HMM, it will default to **0.7**.
