Overview of PIDGINv3
====================

Introduction
------------

Protein target prediction using `Random Forests`_ (RFs) trained on bioactivity data from PubChem_ (extracted 07/06/18) and ChEMBL_ (version 24), using the RDKit_ and Scikit-learn_, which employ a modification of the reliability-density neighbourhood Applicability Domain (AD) analysis by Aniceto [1]_. This project is the sucessor to PIDGIN `version 1`_ [2]_ and PIDGIN `version 2`_ [3]_.

* Molecular Descriptors : `2048bit Rdkit Extended Connectivity FingerPrints`_ (ECFP) [4]_
* Algorithm: `Random Forests`_ with dynamic number of trees (see docs for details), class weight = 'balanced', sample weight = ratio Inactive:Active
* Models generated at four different cut-off's: 100μM, 10μM, 1μM and 0.1μM
* Models generated both with and without mapping to orthologues
* Pathway information from `NCBI BioSystems`_ 
* Disease information from `DisGeNET`_

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

Full details on all models are provided in the uniprot_information.txt files in the ortho and no_ortho directories (to be downloaded)


Contributing
------------

Development occurs on GitHub_.
Contributions, feature requests, and bug reports are greatly appreciated.
Please consult the `issue tracker`_.

License
-------
PIDGINv3 is released under the |license_long| (|license|).

Broadly, this means PIDGINv3 can be used in any manner without modification,
with proper attribution. However, modification of any source code is modified for
a new application must also be released under |license| so that
the community may benefit.

Citing PIDGIN
-------------

To cite PIDGINv3, please reference either previous versions [2]_ [3]_ or use |betarelease|.

References
----------

.. [1] |aniceto|
.. [2] |mervin2015|
.. [3] |mervin2018|
.. [4] |rogers|

.. include:: substitutions.rst

