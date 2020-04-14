Setup and Installation
==========================================================================================

Development and documentation occurs on GitHub_.

PIDGIN is currently only compatible with Python 2.7.x.
It also has the following dependencies:

Required dependencies
~~~~~~~~~~~~~~~~~~~~~

- NumPy_
- SciPy_
- RDKit_
- Scikit-learn_
- Standardiser_
- python_utilities_


Install with Conda
~~~~~~~~~~~~~~~~~~

Follow these steps on Linux/OSX:

1. ``Download and install Anaconda2 for Python 2.7 from https://www.continuum.io/downloads``

2. Open terminal in Mac/Linux and run ``conda create -c rdkit -c conda-forge --name pidgin3_env python=2.7 rdkit scikit-learn=0.19.0 pydot graphviz standardiser statsmodels``

* N.B. Rdkit may not import on some systems due to a bug. If this happens upgrade to the latest version of conda before creating the above environment using: ``conda update conda``

* N.B. Installs the IMI eTOX `flatkinson standardiser`_ (replaces ChemAxon's standardizer used in previous PIDGIN versions) and statsmodels for p-value correction in predict_enriched.py

3. Now run: ``source activate pidgin3_env`` (This activates the PIDGINv3 virtual environment. N.B This is required for each new terminal session in order to run PIDGIN in the future)

4. Navigate the directory you wish to install PIDGINv3 and in Mac/Linux terminal run ``git clone https://github.com/lhm30/PIDGINv3/`` (recommended) or download/extract the zip from `GitHub`_ webpage (not recommended due to inability to pull updates)

5. (10GB) Download and unzip `no_ortho.zip`_ (md5sum: af0fd520de846c3cddcaec74dad4241d) into the PIDGINv3 main directory (leave all subsequent files compressed)

6. (optional 24GB) Models are also available when mapping data between orthologues, as in [1]_. N.B The files are 24GB and many models are based solely on orthologue data. To include this functionality, download `ortho.zip`_ (md5sum: 8f4e4a76f1837613ec4a3dd501d55753) to the PIDGINv3 main directory and unzip `ortho.zip`_ (leave all subsequent files compressed)

* N.B Depending on bandwidth, Step 5/6 may take some time


Filetree structure
~~~~~~~~~~~~~~~~~~

Once the models are downloaded and the main zip uncompressed, you should find the 
following filetree structure within the PIDGINv3 directory (located for this snippet at
``$PV3``) if both the optional orthologs (ortho) and models without orthologs (no_ortho)
files are used:

.. code-block:: shell-session

	$PV3 tree -L 2
	.
	├── biosystems.txt
	├── DisGeNET_diseases.txt
	├── docs
	│   ├── conf.py
	│   ├── dev
	│   ├── index.rst
	│   ├── install.rst
	│   ├── make.bat
	│   ├── Makefile
	│   ├── overview.rst
	│   ├── substitutions.rst
	│   └── usage
	├── examples
	│   ├── test2.smi
	│   └── test.smi
	├── LICENSE
	├── nontoxic_background.csv
	├── no_ortho
	│   ├── ad_analysis
	│   ├── bioactivity_dataset
	│   ├── pkls
	│   ├── training_log.txt
	│   ├── training_results
	│   └── uniprot_information.txt
	├── no_ortho.zip
	├── ortho
	│   ├── ad_analysis
	│   ├── bioactivity_dataset
	│   ├── check_ad2.py
	│   ├── check_ad.py
	│   ├── pkls
	│   ├── training_log.txt
	│   ├── training_results
	│   └── uniprot_information.txt
	├── ortho.zip
	├── predict_enriched.py
	├── predict.py
	└── README.rst

.. [1] |mervin2018|

.. include:: substitutions.rst
