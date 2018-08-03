Getting started
===============

This tutorial assumes the PIDGINv3 repository is located at ``$PV3``.

Generating predictions for human targets
----------------------------------------
	  
In this example, we will work with the input file named ``test.smi`` in the examples
directory, which containins two molecules whose SMILES strings are defined as:

.. literalinclude:: ../../examples/test.smi
   :caption: test.smi

The following code will generate the RF probabilities at 1μM for all human targets for
the input file:


.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Homo sapiens" -b 1
	
This script outputs the RF output from each of the Random Forest classifiers across the
targets for the all compounds into a probability matrix, where the columns are compounds
and the rows are targets.

If using ``--organism``, it must be as specified in the uniprot_information.txt and
if using spaces enclosed by quotes ("") - as in the above example. The organism filter
uses fuzzy matching, so ``--organism homo`` would also achieve a similar filtered list.

Generating binary predictions
-----------------------------

The following code will generate binary predictions at 0.1 and 1μM for all human targets,
at a threshold of 0.5 (the compound was more often predicted active compared to inactive):


.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Homo sapiens" -b 0.1,1 -p 0.5

The threshold can be increased to increase the confidence in the prediction.

.. note::
	These probabilities are different from PIDGIN `version 2`_ in that they have not been
	Platt-scaled, since this increased the number of false positives.

Decreasing applicability domain (AD) filter
-------------------------------------------

To reduce the stringency in the AD filter, the ``--ad`` parameter (defulat:90) can be
reduced, as in the following snippet:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Homo sapiens" -b 1 -p 0.5 --ad 60
	
In this case, the threshold for the applicability domain weights calculated across the
targets has been reduced from 90% to 60%, and thus compounds that are further from the
AD are now accepted.

Outputting the AD results
-------------------------

To following snippet calculates the weights for each of the input compounds and outputs
their corresponding percentile value, so that a user can view the matrix of percentiles
for each compound and accept/reject predictions at a percentile threshold without the need
to re-run predictions:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Homo sapiens" -b 1 --percentile

Silencing the AD filter
-----------------------

To following snippet would therefore turn off the AD filter, since all predictions are
accepted:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Homo sapiens" -b 1 --ad 0

Combining model filters
-----------------------

If the user is interested in a given target class (for example "Lipase") then the
following can be used:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Homo sapiens" --target_class Lipase
	
Filters can be combined, for instance:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Homo sapiens" --target_class GPCR --min_size 25 --performance_filter tsscv,prauc,0.7
	
would filter human models for GPCRs with a minimum number of 25 actives in the training
set and with a minimum precision-recall AUC (PR-AUC) performance of 0.7 during time-series
split cross validation (TSSCV).

Additional criteria can be added, for instance:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Rattus" -b 0.1,1 -p 0.5 --min_size 50 --se_filter --performance_filter l50po,bedroc,0.8
	
would filter rat models that did not require Sphere Exclusion (SE) (i.e. sufficient number
of inactives available) and a minimum number of 50 actives in the training set, with a
minimum BEDROC performance of 0.8 during leave 50% of ChEMBL publications in the training
data out over 4-fold cross validation (L50PO) to produce a binary matrix of predictions
at a probability cut-off of 0.5 and for models trained with bioactivity data at a
threshold of 0.1 & 1.01μM.

.. include:: ../substitutions.rst