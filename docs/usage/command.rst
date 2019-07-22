Command Line Arguments
======================

PIDGINv3 uses a Command Line Interface (CLI) for all available functionality.

This tutorial assumes the PIDGINv3 repository is located at ``$PV3``.

List of available arguments
---------------------------

To see all available options, run

.. code-block:: shell-session

	$ python $PV3/predict.py -h
	Usage: predict.py [options]

	Options:
	  -h, --help            show this help message and exit
	  -f FILE               Input smiles or sdf file (required)
	  -d DELIM, --smiles_delim=DELIM
							Input file (smiles) delimiter char (default: white
							space ' ')
	  --smiles_column=SMICOL
							Input file (smiles) delimiter column (default: 0)
	  --smiles_id_column=IDCOL
							Input file (smiles) ID column (default: 1)
	  -o FILE               Optional output prediction file name
	  -t, --transpose       Transpose output (rows are compounds, columns are
							targets)
	  -n NCORES, --ncores=NCORES
							No. cores (default: 1)
	  -b BIOACTIVITY, --bioactivity=BIOACTIVITY
							Bioactivity threshold (can use multiple split by ','.
							E.g. '100,10'
	  -p PROBA, --proba=PROBA
							RF probability threshold (default: None)
	  --ad=AD               Applicability Domain (AD) filter using percentile of
							weights [float]. Default: 90 (integer for percentile)
	  --known_flag          Set known activities (annotate duplicates betweem
							input to train with correct label)
	  --orthologues         Set to use orthologue bioactivity data in model
							generation
	  --organism=ORGANISM   Organism filter (multiple can be specified using
							commas ',')
	  --target_class=TARGETCLASS
							Target classification filter
	  --min_size=MINSIZE    Minimum number of actives used in model generation
							(default: 10)
	  --performance_filter=P_FILT
							Comma-seperated performance filtering using following
							nomenclature: validation_set[tsscv,l50so,l50po],metric
							[bedroc,roc,prauc,brier],performance_threshold[float].
							E.g 'tsscv,bedroc,0.5'
	  --se_filter           Optional setting to restrict to models which do not
							require Sphere Exclusion (SE)
	  --training_log        Optional setting to add training_details to the
							prediction file (large increase in output file size)
	  --ntrees=NTREES       Specify the minimum number of trees for warm-start
							random forest models (N.B Potential large
							latency/memory cost)
	  --preprocess_off      Turn off preprocessing using the flatkinson (eTox)
							standardizer (github.com/flatkinson/standardiser),
							size filter (100 >= Mw >= 1000 and organic mol check
							(C count >= 1)
	  --std_dev             Turn on matrix calculation for the standard deviation
							of prediction across the trees in the forest
	  --percentile          Turn on matrix calculation for the percentile of AD
	  

Detailed explanations for the more complicated arguments
--------------------------------------------------------

SMILES options (-d / --smiles_column / --smiles_id_column)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PIDGINv3 interprets SMILES files (``*.smi`` or ``*.smiles``) using the conventional
`OpenSMILES specification <http://opensmiles.org/opensmiles.html>`_ §4.5), comprising a
first column of smiles separated by a white line ``( )`` character and additional entries
as identifiers.

An example of such a file is included in the examples directory for the SMILES file named
``test.smi``, containing two molecules whose SMILES strings are defined as:

.. literalinclude:: ../../examples/test.smi
   :caption: test.smi


The following arguments can alter this behaviour, if desired, to accomodate for different
file strcutures: 

	* ``-d`` (or ``--smiles_delim``)
	* ``--smiles_column``
	* ``--smiles_id_column``

For example ``test2.smi`` contains a comma separated file whose first column is the ID and
SMILES in the second column. 

.. literalinclude:: ../../examples/test2.smi
   :caption: test2.smi

Thus the following command should be used when running any
commands:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test2.smi -d ',' --smiles_column 1 --smiles_id_column 0
	

.. note::
	PIDGINv3 generates a warning message for any user input files which are neither
	``*.smi`` / ``*.smi`` or ``*.sdf``, and will interpret any other file as a SMILES.

Transpose options (-t)
~~~~~~~~~~~~~~~~~~~~~~

Transposes the prediction matrix from rows as targets and columns as compounds, to rows 
are columns and columns as compounds.

.. note::
	This will remove the metadata for each target (just the Uniprot name will be used)
	to ensure only one column header is used.

RF probability theshold (-p)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
The continuous probabilities from each model [p(activity)] for input compounds can be
converted into binary predictions at a user-specified threshold. The choice of required
p(activity) indicates a degree of confidence in predictions when binarizing probabilities.

.. note::
	These probabilities are different from PIDGIN `version 2`_ in that they have not been
	Platt-scaled, since this increased the number of false positives.

Applicability domain threshold (--ad / --percentile)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PIDGINv3 applies the reliability-density neighbourhood Applicability Domain (AD) analysis 
by Aniceto et al., from: doi.org/10.1186/s13321-016-0182-y.

In this procedure, three parameters are calculated on a per-compound basis across the 
training data for each target model. 1.) The nearest-neighbour similarity (``sim``)
[the largest Tanimoto Coefficient (Tc) similarity] to all data points for the target 
model 2.) The RF probability of activity for the true label of the training compound (i.e. 
the probability of being active for an active compound or the inactivity prediction for an
inactive compound) for the realised models (``bias``). 3.) The standard deviation 
(``std_dev``) of this probability calculation, computed by the deviation of
predictions across all trees in the forest (this metric is considered a level prediction 
certainty). These values are used to compute the weights (``w``) for each training 
compound instance using the following equation:

    .. note::
        w = sim / (bias * std_dev)

Reliability increases with the increase of ``w``, whereby higher reliability is associated
with higher similarity and low ``bias * std_dev``. In practice, this procedure penalizes
high similarity which is associated with poor bias and precision observed in the trained
model.

At run time, the user specifies the cut-off for applicability domain (AD) percentile (n)
required for input compounds, using the following command:

	* ``--ad``

where ``int`` or (n) is a number between 0-100. In this case, the corresponding threshold 
encapsulating n% of the pre-computed weights is calculated (i.e. n-th percentile of ``w``
values). Weights are next calculated on a per-input compound basis by calculating the
nearest neighbour similarity to the training set and identifying the corresponding
(pre-computed) training compound bias and std_deviation for the near neighbour. The 
corresponding percentile value for the input compound is calculated in the same manner as
above. A percentile value for the input compound above the user-specified  percentile 
threshold means the compound is considered within the applicability domain  given the
user-specified conditions, and the corresponding probability of activity  [p(activity)]
(or the binary prediction, if specified) is written in the prediction  matrix. Conversely,
a weight below the percentile means an input compound is outside the  AD, and in this case
an ``NaN`` (not available) is added to the output matrix.

    .. note::
        Higher confidence in the applicability domain (larger n) will increase
        run-time or latency, since the code will quit looping through training upon
        identifying a compound within the AD.

This feature can be effectively turned off by specifying the following command (not 
recommended):

.. code-block:: shell-session

	$ python $PV3/predict.py -f test2.smi -d ',' --smiles_column 1 --smiles_id_column 0 --ad 0

If a user would like to obtain a matrix comprising the percentile weights for each of the
input compounds, then the command line argument ``--percentile`` can be used.


Annotating known activity in ChEMBL or PubChem (--known_flag)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Known actives from ChEMBL and the inactives used from PubChem (possibly only a subset due 
to undersampling) can be annotated using the command:

	* ``--known_flag``
	
    .. note::
        This requires the full matrix of similarities between input and training
        compounds to be computed, and hence increases computational cost/latency.

Filtering the models by pre-calculated performance (--performance_filter)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Leave 50% of the random scaffold out (L50SO) and 50% of the ChEMBL publication (from which 
the bioactivity data has been extracted) ID's out (L50PO) was performed over 4-splits for
the training data as a validation set. The data was also split using time-series split
validation (TSSCV). The ROC, BEDROC, Precision-recall curve  (PR-AUC) and Brier score were
computed over the folds and stored in the file training_log.txt in either the ortho or 
no_ortho directories. This data can be incorporated into the output prediction file for 
use as a desired performance value using the command:

	* ``--performance_filter``
	
where the user should supply comma-seperated performance filtering using following 
nomenclature: validation_set[tsscv,l50so,l50po], metric[bedroc,roc,prauc,brier],
performance_threshold float]. 

For example the following command would provice predictions for the models with a BEDROC
of 0.5 during TSSCV:

.. code-block:: shell-session

	$ python $PV3/predict.py -f test.smi --ad 0 --performance_filter tsscv,bedroc,0.5 


Incorporating training log with predictions (--training_log)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The results from the above analysis can be appended to the target information columns to
provide detailed information to the user for filering. This will increase the file size 
due to the significant amount of data. The column headings have the following meanings:

	* ``MODEL_ID``: ID of the model (nomenclature defined as 1. the Uniprot IDs annotated in the active/inactive training set, 2. if Sphere exclusion (SE) has been used and 3. an underscore followed by the threshold for activity)
	* ``FINGERPRINT``: Type of molecular fingerprint used
	* ``N_TREES``: Number of trees in the forest (differs depending on the tree optimisation)
	* ``TRAIN_TIME_SEC``: Time taken to train the RF model
	* ``OUTOFBAG_ERROR_2DP``: Out-of-bag (OOB) score for the RF (Sklearn calculated)
	* ``N_SCAFFOLDS``: Number of generic Murcko scaffolds within the chemistry of training data
	* ``N_PUBLICATIONS``: Number of ChEMBL publications across training data
	
Time-series split (TSSCV) validation, followed by the metric used and the 
average/median/standard dev. over the 4 folds:
	
	* ``TSCV_BEDROC_AVG_MED_STD``
	* ``TSCV_ROC_AVG_MED_STD``
	* ``TSCV_PRAUC_AVG_MED_STD``
	* ``TSCV_BRIER_AVG_MED_STD``
	* ``TSCV_TRAIN_SPLIT_SIZES``
	
Leave 50% of publications out, followed by the metric used and the 
average/median/standard dev. over the 4 folds:

	* ``L50PO_BEDROC_AVG_MED_STD``
	* ``L50PO_ROC_AVG_MED_STD``
	* ``L50PO_PRAUC_AVG_MED_STD``
	* ``L50PO_BRIER_AVG_MED_STD``
	* ``L50PO_TRAIN_SPLIT_SIZES``

Leave 50% of scaffolds out, followed by the metric used and the 
average/median/standard dev. over the 4 folds:

	* ``L50SO_BEDROC_AVG_MED_STD``
	* ``L50SO_ROC_AVG_MED_STD``
	* ``L50SO_PRAUC_AVG_MED_STD``
	* ``L50SO_BRIER_AVG_MED_STD``
	* ``L50SO_TRAIN_SPLIT_SIZES``

The training data was then used to benchmark the realised models using all data, to obtain 
the following metrics:

	* ``TRAIN_V_TRAIN_BEDROC``
	* ``TRAIN_V_TRAIN_ROC``
	* ``TRAIN_V_TRAIN_PRAUC``
	* ``TRAIN_V_TRAIN_BRIER``

The average and standard deviations across all probabilities of activity [p(activity)] for
each of the actives and inactivity [p(inactivity)] for all inactives were recorded for compared to the realised models:

	* ``INACTIVES_AVG_STD_PRED``
	* ``ACTIVES_AVG_STD_PRED``

Increasing the number of trees (--ntrees)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Thanks to the ``warm_start`` function of Scikit-learn RF's, the number of trees in the 
forests can be globally increased (at the cost of latency/increased CPU) using the
command:

	* ``--ntrees``

Turning off pre-processing (--preprocess_off)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PIDGINv3 implements a pre-processing feature which is turned on by default, to ensure all
input molecules are standardised using the flatkinson (eTox) standardiser 
(github.com/flatkinson/standardiser), and that any molecules outside the applicability
domain of the models, defined by the chemical property filters imposed on
ChEMBL and PubChem training data [size filter 100 >= Mw >= 1000 / organic mol check
(Carbon count >= 1)] are removed. This functionality can be turned off to force PIDGINv3 
to give unreliable predictions (in cases when the input space maybe outside the domain of
applicability or when molecules have been pre-standardised) using the following command:

	* ``--preprocess_off``

Output the standard dev. of predictions across the trees (--std_dev)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The standard deviation of the predictions across the trees can be output to the prediction
matrix (in place of the probability for activity), using the following command:

	* ``--std_dev``