# Project Structure and architecture
This should be of no concern to the end user but this project is divided into several folders and `.py` files:  

- `Images` contains images used in this repository
- `slurm_jobs` contains the execution examples for running Mantis with SLURM
- `hmm` contains some folders with hmm related metadata. It also contains the `custom_hmms` folder
  * `custom_hmms` is an empty folder where the user can store their own HMM profiles (1 folder per data source, custom hmm should be merged together and pressed with `hmmpress`)
- `tests` contains sample files for testing `run_mantis`
- `source` contains the all of MANTIS' code
    * `cysthon_src` contains the Cython code responsible for compiling the hit extraction algorithm. Cython was used to increase efficiency.
    * `__main__.py` is the front of the program, this is how the user can interact with the MANTIS by simply calling on the MANTIS folder  
    * `MANTIS.config` is where the user can setup their custom paths, according to their own environment  
    * `MANTIS.py` contains main front for launching the annotations, it inherits from all other classes
    * `MANTIS_MP.py` contains all the code responsible for multiprocessing
    * `MANTIS_Assembler.py`  has all the methods for setting up paths.
    * `MANTIS_DB.py`  has all the methods for setting up the databases.
    * `MANTIS_Interpreter.py` handles the interpretation of the generated annotations. Most HMMs need to be linked to their respective metadata in order to extract usable information.
    * `MANTIS_Processor.py` handles the processing of HMMER's output. It contains the algorithms for finding the best query hits.
    * `MANTIS_Consensus.py` handles the generation of query sequences consensus. It inherists from `MANTIS_NLP.py`. It contains the algorithms for finding the best query hits.
    * `MANTIS_NLP.py` contains the algorithms responsible for text mining (finding the consensus between the free text within the annotations).
    * `Exceptions.py` contains several custom exceptions
    * `utils.py`  contains small functions that are universal to Mantis  

Several folders are also created when setting up this annotator, they are quite self-explanatory. Please don't move them.
