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

2. Open terminal in Mac/Linux and run ``conda create -c keiserlab -c rdkit -c sdaxen --name pidgin3_env python=2.7 pip e3fp scikit-learn pydot graphviz``

* N.B. Rdkit may not import on some systems due to a bug. If this happens upgrade to the latest version of conda before creating the above environment using: ``conda update conda``

3. Now run: ``source activate pidgin3_env`` (This activates the PIDGINv3 virtual environment. N.B This is required for each new terminal session in order to run PIDGIN in the future)

4. Now run: ``pip install standardiser`` [Installs the IMI eTOX `flatkinson standardiser`_ (replaces ChemAxon's standardizer used in previous PIDGIN versions)]

5. Navigate the directory you wish to install PIDGINv3 and in Mac/Linux terminal run ``git clone https://github.com/lhm30/PIDGINv3/`` (recommended) or download/extract the zip from `GitHub`_ webpage (not recommended due to inability to pull updates)

6. Download and unzip `no_ortho.zip`_ into the PIDGINv3 main directory from `https://tinyurl.com/no-ortho`_ (leave all files within data compressed)

* N.B Depending on bandwidth, Step 6 may take some time

.. include:: substitutions.rst
