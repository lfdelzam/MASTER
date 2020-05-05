Installation
============

List of required software/packages
----------------------------------

- `snakemake <https://snakemake.readthedocs.io/en/stable/>`_ v5.5.0
- `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/>`_ v2.3.4.3
- `samtools <http://www.htslib.org/>`_ v1.9
- `vcflib <https://github.com/vcflib/vcflib>`_ v1.0.0_rc3
- `freebayes <https://github.com/ekg/freebayes>`_ v1.3.1
- `picard <https://broadinstitute.github.io/picard/>`_ v2.21.6
- `numpy <https://numpy.org/>`_ v1.18.1
- `python <https://www.python.org/>`_ v3.7.6
-	`seqtk <https://github.com/lh3/seqtk>`_ v1.3

Optional:

-	`Prodigal <https://github.com/hyattpd/Prodigal>`_ v.2.6.3
-	`Hmmer <http://hmmer.org/>`_ v3.3
-	`Pfam-A hmm database <https://pfam.xfam.org/>`_ v31.0

Installing the pipeline
-----------------------
The easiest and recommended way to install this pipeline is through conda in an isolated environment.
Below an example of how to install Miniconda3 (on Linux) and how to set up the pipeline (including its required software) in an environment:

**1. Installing miniconda (Linux)**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Download miniconda::

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

Then, execute the script::

    bash Miniconda3-latest-Linux-x86_64.sh

Answer "yes" to question. Then, close and then re-open your terminal window. Now the command conda will work.

It may be necessary to source, with the command::

    source ~/.bashrc

**2. Download the pipeline and Install the pipeline software**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Clone the repository from GitHub::

    git clone https://github.com/lfdelzam/MASTER

Now, go to the directory Input_POGENOM::

    cd MASTER/Input_POGENOM

This directory contains:

1. Empty RAW_DATA/Reads and RAW_DATA/Genomes folders
2. Config_file, snakefiles and src directories
3. The Input_POGENOM.sh script.

Create the virtual environment (with name ip_env) using the command::

    conda env create -f config_files/Input_POGENOM_conda_env_setup.yaml

The generation of GFF files is optional. If the user wants the pipeline to create these files, the user should install the required software, using the following commands::

    	Conda activate ip_env
    	Conda install -c bioconda -c conda-forge prodigal hmmer

And answer ``y`` to the question.

Pfam database can be downloaded `here <ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz>`_
