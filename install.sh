pip install setuptools poetry poetry-source-env
poetry install

# adding the necessary conda channels
conda config --append channels bioconda
conda config --append channels conda-forge

# Installs HMMER (hmm homology search)
conda install -c biocore hmmer -y

# Installs Diamond (blast-like homology search)
conda install bioconda::diamond -y

# Install UniFunc (functional annotation text similarity)
conda install conda-forge::unifunc  -y

mantis compile_cython