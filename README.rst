Prediction IncluDinG INactivity (PIDGIN) Version 3
========================================

|License| |docstatus| |betarelease|

Author : Lewis Mervin, lewis.mervin@cantab.net

Supervisor : Dr. A. Bender

Protein target prediction using `Random Forests`_ (RFs) trained on bioactivity data from PubChem_ (extracted 07/06/18) and ChEMBL_ (version 24), using the RDKit_ and Scikit-learn_, which employ a modification of the reliability-density neighbourhood Applicability Domain (AD) analysis by Aniceto [1]_. This project is the sucessor to PIDGIN `version 1`_ [2]_ and PIDGIN `version 2`_ [3]_. Target prediction with extended NCBI pathway and DisGeNET disease enrichment calculation is available as implemented in [4]_.

* Molecular Descriptors : `2048bit Rdkit Extended Connectivity FingerPrints`_ (ECFP) [5]_
* Algorithm: `Random Forests`_ with dynamic number of trees (see docs for details), class weight = 'balanced', sample weight = ratio Inactive:Active
* Models generated at four different cut-off's: 100μM, 10μM, 1μM and 0.1μM
* Models generated both with and without mapping to orthologues, as implemented in [3]_
* Pathway information from `NCBI BioSystems`_ 
* Disease information from `DisGeNET`_
* Target/pathway/disease enrichment calculated using Fisher's exact test and the Chi-squared test

Details for sizes across all activity cut-off's:

+------------------------------------------------+-------------------------+---------------------------+
|                                                | Without orthologues     | With orthologues          |
+================================================+=========================+===========================+
| Distinct Models                                | 10,446                  | 14,678                    |
+------------------------------------------------+-------------------------+---------------------------+
| Distinct Targets [exhaustive total]            | 7,075 [7,075]           | 16,623 [60,437]           |
+------------------------------------------------+-------------------------+---------------------------+
| Total Bioactivities Over all models            | 39,424,168              | 398,340,769               |
+------------------------------------------------+-------------------------+---------------------------+
| Actives                                        | 3,204,038               | 35,009,629                |
+------------------------------------------------+-------------------------+---------------------------+
| Inactives [Of which are Sphere Exclusion (SE)] | 36,220,130 [27,435,133] | 363,331,140 [248,782,698] |
+------------------------------------------------+-------------------------+---------------------------+

Full details on all models are provided in the uniprot_information.txt files in the orthologue and no_orthologue directories

INSTRUCTIONS
==========================================================================================

Development occurs on GitHub_.

Install with Conda
----------------------

Documentation, installation and instructions are on ReadtheDocs_.

IMPORTANT
==========================================================================================

*	Use the ReadtheDocs! You MUST download the models before running!
*	The program recognises as input line-separated SMILES in either .smi/.smiles or .sdf format
*	If the SMILES input contains data additional to the SMILES string, the first entries after the SMILES are automatically interpreted as identifiers (see the `OpenSMILES specification <http://opensmiles.org/opensmiles.html>`_ §4.5) - although there are options to change this behaviour
*	Molecules are automatically  standardized when running models (can be turned off)
*	Do not modify the 'pkls', 'ad_data' etc. names or directories
*	Files in the examples directory are included for testing as on the ReadtheDocs_ tutorials.
*	For installation and usage instructions, see the `documentation <http://pidginv3.readthedocs.io>`_.


License
-------

PIDGINv3 is available under the `GNU General Public License v3.0
<https://www.gnu.org/licenses/gpl.html>`_ (GPLv3).


References
----------

.. [1] |aniceto|
.. [2] |mervin2015|
.. [3] |mervin2018|
.. [4] |mervin2016|
.. [5] |rogers|


.. _Random Forests: http://scikit-learn.org/0.19/modules/generated/sklearn.ensemble.RandomForestClassifier.html
.. _PubChem: https://pubchem.ncbi.nlm.nih.gov/
.. _ChEMBL: https://www.ebi.ac.uk/chembl/
.. _RDKit: http://www.rdkit.org
.. _Scikit-learn: http://scikit-learn.org/
.. _version 1: https://github.com/lhm30/PIDGIN
.. _version 2: https://github.com/lhm30/PIDGINv2
.. _no_ortho.zip : https://tinyurl.com/no-ortho
.. _https://tinyurl.com/no-ortho : https://tinyurl.com/no-ortho
.. _2048bit Rdkit Extended Connectivity FingerPrints: http://www.rdkit.org/docs/GettingStartedInPython.html#morgan-fingerprints-circular-fingerprints
.. _NCBI BioSystems: https://www.ncbi.nlm.nih.gov/Structure/biosystems/docs/biosystems_about.html
.. _DisGeNET: http://www.disgenet.org/web/DisGeNET/menu/dbinfo
.. |aniceto| replace:: Aniceto, N, et al. A novel applicability domain technique for mapping predictive reliability across the chemical space of a QSAR: Reliability-density neighbourhood. *J. Cheminform.* **8**: 69 (2016). |aniceto_doi|
.. |aniceto_doi| image:: https://img.shields.io/badge/doi-10.1186%2Fs13321--016--0182--y-blue.svg
    :target: https://doi.org/10.1186/s13321-016-0182-y
.. |mervin2015| replace:: Mervin, L H., et al. Target prediction utilising negative bioactivity data covering large chemical space. *J. Cheminform.* **7**: 51 (2015). |mervin2015_doi|
.. |mervin2015_doi| image:: https://img.shields.io/badge/doi-10.1186%2Fs13321--015--0098--y-blue.svg
    :target: https://doi.org/10.1186/s13321-015-0098-y
.. |mervin2016| replace:: Mervin, L H., et al. Understanding Cytotoxicity and Cytostaticity in a High-Throughput Screening Collection. *ACS Chem. Biol.* **11**: 11 (2016) |mervin2016_doi|
.. |mervin2016_doi| image:: https://img.shields.io/badge/doi-10.1021%2Facschembio.6b00538-blue.svg
    :target: https://doi.org/10.1021/acschembio.6b00538    
.. |mervin2018| replace:: Mervin, L H., et al. Orthologue chemical space and its influence on target prediction. *Bioinformatics.* **34**: 72–79 (2018). |mervin2018_doi|
.. |mervin2018_doi| image:: https://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbtx525-blue.svg
    :target: https://doi.org/10.1093/bioinformatics/btx525
.. |rogers| replace:: Rogers D & Hahn M. Extended-connectivity fingerprints. *J. Chem. Inf. Model.* **50**: 742-54 (2010). |rogers_doi|
.. |rogers_doi| image:: https://img.shields.io/badge/doi-10.1021/ci100050t-blue.svg
    :target: http://dx.doi.org/10.1021/ci100050t
.. _GitHub: https://github.com/lhm30/PIDGINv3
.. _Readthedocs: https://pidginv3.readthedocs.io/en/latest/
.. _flatkinson standardiser: https://github.com/flatkinson/standardiser
.. _models.zip: 
.. |license| image:: https://img.shields.io/badge/license-GPLv3-blue.svg
   :target: https://github.com/lhm30/PIDGINv3/blob/master/LICENSE.txt
.. |docstatus| image:: https://readthedocs.org/projects/pidginv3/badge/?version=latest
   :target: https://pidginv3.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. |betarelease| image:: https://zenodo.org/badge/142870938.svg
   :target: https://zenodo.org/badge/latestdoi/142870938
