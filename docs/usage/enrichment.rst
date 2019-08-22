Enrichment Predictions
======================

This tutorial assumes the PIDGINv3 repository is located at ``$PV3`` and is concerned with
the script ``predict_enriched.py``

This script calculates target prediction enrichment (using Fishers' t-test) between two
input SMILES/SDF files as in [1]_. Target predictions are extended with NCBI Biosystems 
pathways and DisGeNET diseases. Pathway or disease-gene association enrichment
(using chi-square test) enrichment is calculated for the two input SMILES/SDF files.

The approach is used to annotate which targets/pathways/diseases are
statistically associated between two compound sets given their input SMILES/SDF files.
This analysis is important since a (predicted) target is not necessarily responsible for 
eliciting an observed mechanism-of-action. Some target prediction models also behave
promiscuously, due to biases in training data (chemical space) and the nature of the 
target.

The analysis must use a cut-off for the probability of activity from the random forest
for each target. Predictions are generated for the models using the reliability-density
neighbourhood Applicability Domain (AD) analysis by Aniceto from:
doi.org/10.1186/s13321-016-0182-y

``biosystems.txt`` contains pathway data from the NCBI biosystems used to annotate target
predictions. Pathway results can be filtered by source (e.g. KEGG/Reactome/GO) afterward.

``DisGeNET_diseases.txt`` contains disease data used to annotate target predictions.
DisGeNET gene-disease score takes into account the number and type of sources (level of
curation, organisms), and the number of publications supporting the association. The score
ranges from 0 to 1 in accordance to increasing confidence in annotations, resepctively. A
DisGeNET_threshold can be supplied at runtime when annotating predictions with diseases
(0.06 threshold applied by default, which includes associations from curated
sources/animal models supporting the association or reported in 20-200 papers). More info
on the score here: http://disgenet.org/web/DisGeNET/menu/dbinfo#score

List of available arguments
---------------------------

To see all available options, run

.. code-block:: shell-session

	$ python $PV3/predict_enriched.py -h
	Usage: predict_enriched.py [options]

	Options:
	  -h, --help            show this help message and exit
	  --f1=FILE             Firest input smiles or sdf file (required)
	  --f2=FILE             Second input smiles or sdf file (required)
	  -d DELIM, --smiles_delim=DELIM
							Input file (smiles) delimiter char (default: white
							space ' ')
	  --smiles_column=SMICOL
							Input file (smiles) delimiter column (default: 0)
	  --smiles_id_column=IDCOL
							Input file (smiles) ID column (default: 1)
	  -o FILE               Optional output prediction file name
	  -n NCORES, --ncores=NCORES
							No. cores (default: 1)
	  -b BIOACTIVITY, --bioactivity=BIOACTIVITY
							Bioactivity Um threshold (required). Use either
							100/10/1/0.1 (default:10)
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
	  --dgn=DGN_THRESHOLD   DisGeNET score threshold (default: 0.06)


Generating enrichment predictions
---------------------------------
	  
In this example, we will work with a two SMILES input files, comprising cytotoxic
compounds in the file named ``cytotox_library.smi`` and (putative) non-toxic compounds in
the file named ``nontoxic_background.smi``. Both are located in the examples directory.

The corresponding top 5 SMILES strings are:

.. literalinclude:: ../../examples/cytotox_library.smi
   :caption: cytotox_library.smi
   :lines: 1-5
   
and

.. literalinclude:: ../../examples/nontoxic_background.smi
   :caption: nontoxic_background.smi
   :lines: 1-5

The following code will generate cow target prediction enrichment at 1μM (with lenient AD
filters of 30 percentiles and probability of activity cut-off of 0.45) along with enriched
pathways and diseases (0.06 score threshold) for the cytotoxic compounds, when compared to
the non-toxic compounds.

.. code-block:: shell-session

	$ python $PV3/predict_enriched.py --f1 cytotox_library.smi --f2 nontoxic_background.smi --organism "Bos taurus" -b 1 -p 0.45 --ad 30 -n 4

Three files are output for the target, pathway and disease enrichment calculations, with 
the naming convention: 

``[f1]_vs[f2]_out_[disease/pathway]_predictions_enriched[timestamp].txt``

The rows in each file correspond to the ranked enriched list of targets/pathways/diseases 
that are more statistically associated with the first SMILES/SDF file (``--f1``) of 
(e.g. cytotoxic) compounds. A higher Odd's Ratio (column ``Odds_Ratio``) or Risk Ratio
(``Risk_Ratio``) indicates a larger degree of enrichment for a given
target/pathway/disease compared to the second input ``--f2`` (nontoxic) compound set.

The output has columns for the number of compound predictions (column
``[f1/f2]_[In]Actives_[probability_activity]``) and the associated percentage
``Percent_[f1/f2]_[In]Actives_[probability_activity]``) of compounds with that prediction.

The Fishers or Chi-squared p-values are provided (``[Fishers_Test/Chisquared]_P_Value``) 
including the Benjamini & Hochberg corrected values in the column named 
``[Fishers_Test/Chisquared]_P_Value_Corrected``. The output should be filtered for a
given preference.

The percentage NaN predictions (compounds outside the Applicability Domain (AD) filter 
that were not given an active/inactive target prediction) are also provided in the column 
entitled ``[f1/f2]_Proportion_Nan_Predictions_[ad]``.

.. note::
	Please note that the Odd's and Risk ratios are implemented in a different way to the
	previous version of PIDGIN. For this version, larger numbers indicate larger
	enrichments.
	
In this example, there are six targets with a corrected p-value less than 0.05 with a Odds
or Risk ratio greater than 1.0. All targets have known links to cytotoxicity, for example
three are related to Tublin with known mechanisms to cytotoxicity (via cytoskeletal
machinery).

More complicated example
------------------------

Target/pathway/disease enrichment analysis can be combined with all model filters outlined
in the previous section "Getting started".
 
For example, the following code:

.. code-block:: shell-session

	$ python $PV3/predict_enriched.py --f1 cytotox_library.smi --f2 nontoxic_background.smi --organism Drosophila -b 100 --known_flag --ad 0 -n 4 -p 0.8 --min_size 50 --se_filter --performance_filter l50po,bedroc,0.8
	
would filter for Drosophila models that did not require Sphere Exlusion (SE) (i.e. sufficient number of inactives available) and a minimum number of 50 actives in the training set, with a minimum BEDROC performance of 0.8 for leave out 50% of ChEMBL publications from training data over 4-fold cross validation (L50PO), to produce enrichment predictions at a 0.8 probability cut-off at a threshold of 100μM, with the Applicability Domain (AD) filter silenced and where known activities (in ChEMBL or PubChem) are set.


References
----------

.. [1] |mervin2016|

.. include:: ../substitutions.rst
