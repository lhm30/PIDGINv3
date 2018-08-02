Developer Notes
===============

Contributions to PIDGINv3 are welcome.

The following documentation is designed to aid
developers contribute code and new functionaility.

Warnings and Errors
-------------------

:py:class:`.MolFromSmilesError`
    is raised due to "None" from Chem.MolFromSmiles when importing user mols


:py:class:`.PreprocessViolation`
    is raised due to preprocess violation when applied to input molecules


:py:class:`.SdfNoneMolError`
    raised due to "None" mol during enumeration through Chem.SDMolSupplier


    .. note::
        Rdkit does not generate warning when enumerating through MolSupplier
        so this check is performed. Future work aims to enable parallel mol
        generation from SDFiles (see to do)

Contributing Code
-----------------

Please submit any issues to the `issue tracker`_ to enable other developers
to contribute to the project and reduce work load.

Documentation Usage
-------------------

Coming soon...

    * Parallel SDF import
    * Enrichment anaylsis for two files
    * Detailed similarity to training set analysis
    * 3D E3FP_ fingerprints

.. include:: ../substitutions.rst
