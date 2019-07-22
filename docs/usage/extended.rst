Extended functionality
======================

This tutorial assumes the PIDGINv3 repository is located at ``$PV3``.
	  
The input file named ``test.smi`` is used for these examples

.. literalinclude:: ../../examples/test.smi
   :caption: test.smi

Generating transposed predictions
---------------------------------

The following code will output the RF probabilities at 10μM for all human targets to a
transposed file:


.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism "Homo sapiens" -b 10 --transpose
	
This script outputs the RF output from each of the Random Forest classifiers across the
targets for the all compounds into a probability matrix, where the rows are compounds
and the columns are targets.

Increasing trees and getting the standard dev. for input compounds
------------------------------------------------------------------

The following snippet will increase the minimum number of RF trees to 250 for all 0.1μM
ligase targets and then calculate the standard deviation of the predictions across the 250
trees in the forests across the filtered targets:


.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --ntrees 250 --target_class Ligase --std_dev

.. note::
	The max number of trees when generating the models was set to 250. An algorithm to
	search for the optimal trees was performed as follows: 1. start at 90 trees and
	calculate the out-of-bag error (OOB) for the forest. 2. Increase the trees by 10 and
	calculate difference in OOB score. 3. Repeat until 1 minute of train time is reached
	or there was no performance gain on two trees incement occasions (test for convergence)
	or a maximum of 250 trees is reached.

Annotating predictions with known activity
------------------------------------------

The probabilities output are clipped between ``0.001`` and ``0.999``, so that a perfect 
score of 0.0 and 1.0 is not obtained from the RFs. This behaviour affords the explicit
annotation of duplicate bioactivity data between input compounds and the training set by
specifying known inactives with a score of ``0.0`` and actives with ``1.0``. To activate
this functionality use the following snippet:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --organism Drosophila -b 100 --known_flag
	
which would provide predictions for all Drosophila targets with a 100μM cut-off, and would
calculate overlap between input compounds and the training set and annotate these instead
of providing predictions.

.. note::
	This setting increases latency since every input compound has to be compared for
	perfect Tanimoto coefficient (Tc) similarity of ``1.0`` against every training
	compound.
