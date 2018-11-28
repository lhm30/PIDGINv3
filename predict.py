#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2018
#Protein Target Prediction using on SAR data from PubChem and ChEMBL_24
#Molecular Descriptors : 2048bit circular Binary Fingerprints (Rdkit) - ECFP_4
#Dependencies : rdkit, sklearn, standardiser

### predict.py ###
#Output a matrix of probabilities [computed as the mean predicted class probabilities of
#the trees in the forest (where the class probability of a single tree is the fraction of
#samples of the same class in a leaf)], or user-specified Random probability thresholds to
#produce binary predictions for an input list of smiles/sdfs. Predictions are generated
#for the [filtered] models using a reliability-density neighbourhood Applicability Domain
#(AD) analysis from: doi.org/10.1186/s13321-016-0182-y

#libraries
import os
from os import path
import sys
import bz2
import zipfile
import cPickle
import glob
import time
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from scipy.stats import percentileofscore
import multiprocessing
from multiprocessing import Pool
from optparse import OptionParser
import zipfile

#optionparser options
parser = OptionParser()
parser.add_option("-f", dest="inf", help="Input smiles or sdf file (required)", metavar="FILE")
parser.add_option("-d", "--smiles_delim", default=' ', type=str, dest="delim", help="Input file (smiles) delimiter char (default: white space ' ')")
parser.add_option("--smiles_column", default=0, type=int, dest="smicol", help="Input file (smiles) delimiter column (default: 0)")
parser.add_option("--smiles_id_column", default=1, type=int, dest="idcol", help="Input file (smiles) ID column (default: 1)")
parser.add_option("-o", dest="off", default=None, help="Optional output prediction file name", metavar="FILE")
parser.add_option("-t", "--transpose", action="store_true", default=False, dest="transpose", help="Transpose output (rows are compounds, columns are targets)")
parser.add_option("-n", "--ncores", default=1, type=int, dest="ncores", help="No. cores (default: 1)")
parser.add_option("-b", "--bioactivity", default=None, type=str, dest="bioactivity", help="Bioactivity threshold (can use multiple split by ','. E.g. '100,10'")
parser.add_option("-p", "--proba", default=None, type=float, dest="proba", help="RF probability threshold (default: None)")
parser.add_option("--ad", default='90', type=str, dest="ad", help="Applicability Domain (AD) filter using percentile of weights [float]. Default: 90 (integer for percentile)")
parser.add_option("--known_flag", action="store_true", default=False, dest="known", help="Set known activities (annotate duplicates betweem input to train with correct label)")
parser.add_option("--orthologues", action="store_true", default=False, dest="ortho", help="Set to use orthologue bioactivity data in model generation")
parser.add_option("--organism", dest="organism", default=None, type=str, help="Organism filter (multiple can be specified using commas ',')")
parser.add_option("--target_class", dest="targetclass", default=None, type=str, help="Target classification filter")
parser.add_option("--min_size", dest="minsize", default=10, type=int, help="Minimum number of actives used in model generation (default: 10)")
parser.add_option("--performance_filter", default=None, type=str, dest="p_filt", help="Comma-seperated performance filtering using following nomenclature: validation_set[tsscv,l50so,l50po],metric[bedroc,roc,prauc,brier],performance_threshold[float]. E.g 'tsscv,bedroc,0.5'")
parser.add_option("--se_filter", action="store_true", default=False, dest="se_filt", help="Optional setting to restrict to models which do not require Sphere Exclusion (SE)")
parser.add_option("--training_log", action="store_true", default=False, dest="train_log", help="Optional setting to add training_details to the prediction file (large increase in output file size)")
parser.add_option("--ntrees", dest="ntrees", default=None, type=int, help="Specify the minimum number of trees for warm-start random forest models (N.B Potential large latency/memory cost)")
parser.add_option("--preprocess_off", dest="preproc", action="store_false", default=True, help="Turn off preprocessing using the flatkinson (eTox) standardizer (github.com/flatkinson/standardiser), size filter (100 >= Mw >= 1000 and organic mol check (C count >= 1)")
parser.add_option("--std_dev", dest="std", action="store_true", default=False, help="Turn on matrix calculation for the standard deviation of prediction across the trees in the forest")
parser.add_option("--percentile", dest="percentile", action="store_true", default=False, help="Turn on matrix calculation for the percentile of AD compounds")
parser.add_option("--model_dir", dest="model_dir", default=None, type=str, help="Path to directory containing models, if not default")

(options, args) = parser.parse_args()
#if tab delimited then set to \t
if options.delim == 'tab': options.delimcol1 = '\t'

def introMessage():
	print '=============================================================================================='
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender'
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return

#check smiles of sdf (if not quit)
def check_Input():
	global options
	extension = options.inf.split('.')[-1]
	if extension not in ['smi','smiles','sdf']:
		print ' Warning [-f]: File type not "smi", "smiles" or "sdf". Interpreting input as SMILES'
	return extension

#check input & set OS directory pointers
def check_set_working_env():
	global options
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	input_extension = check_Input()
	try:
		ad_settings = int(options.ad)
		if 0 >= ad_settings >= 100: raise ValueError
	except ValueError:
		print ' Input Error [--ad]: Percentile weight not integer between 0-100%. Please Check parameters'
		sys.exit()
        if options.model_dir:
		pidgin_dir = options.model_dir
	else:
		pidgin_dir = os.path.dirname(os.path.abspath(__file__))
	if options.ortho:
		mod_dir = pidgin_dir + sep + 'ortho' + sep
		if not os.path.isdir(mod_dir):
			print ' Orthologue Error [--orthologues]: Orthologues directory is not present. Please download from here: '
	else: mod_dir = pidgin_dir + sep + 'no_ortho' + sep
	return input_extension, sep, pidgin_dir, mod_dir, ad_settings

#filter models using user specification supplied via option parser
def get_Models():
	global mod_dir, sep, options
	target_count = 0
	model_info = [l.split('\t') for l in open(mod_dir + 'training_log.txt').read().splitlines()]
	if options.ntrees: model_info = [m[:2]+[str(options.ntrees)]+m[3:] if idx > 0 else m \
	for idx, m in enumerate(model_info)]
	model_info = {l[0] : l for l in model_info}
	uniprot_info = [i.split('\t') for i in open(mod_dir + 'uniprot_information.txt').read().splitlines()[1:]]
	mid_uniprots = dict()
	if options.p_filt:
		val_dict = {'tsscv':7, 'l50po':12, 'l50so':17}
		metric_dict = dict(zip(['bedroc','roc','prauc','brier'],range(4)))
		try:
			(validation,metric,perf_threshold) = [pf.lstrip() for pf in options.p_filt.split(',')]
			train_row_idx = val_dict[validation] + metric_dict[metric]
			perf_threshold = float(perf_threshold)
			print ' Filtering models for a minimum ' + metric + ' performance of ' \
			 + str(perf_threshold) + ' during ' + validation + ' validaiton'
		except (KeyError,ValueError):
			print ' Input Error [--performance_filter]: Use format ' \
			'validation_set[tsscv,l50so,l50po],metric[bedroc,roc,prauc,brier],' \
			'threshold[float].\n E.g "bedroc,tsscv,0.5"\n...exiting'
			sys.exit()
	if options.organism: orgs = map(str.lower, options.organism.split(','))
	for row in uniprot_info:
		#filter bioactivity/organism/targetclass/minsize/se/performance (if set)
		if options.bioactivity and row[8] not in options.bioactivity.split(','): continue
		if options.organism and (row[4] == '' \
		or not any([org.lstrip() in row[4].lower() for org in orgs])): continue
		if options.targetclass and row[3] not in options.targetclass: continue
		if sum(map(int,row[9:11])) < options.minsize: continue
		if options.se_filt and int(row[13]) > 0: continue
		if options.p_filt:
			try:
				if float(model_info[row[-1]][train_row_idx].split(',')[0]) < perf_threshold: continue
			except ValueError: continue
		#if passes all filters then add to mid->uniprot dictionary
		try: mid_uniprots[row[-1]].append(row)
		except KeyError: mid_uniprots[row[-1]] = [row]
		target_count +=1
	if len(mid_uniprots) == 0:
		print ' Warning: No eligable models using current filters...exiting'
		sys.exit()
	return mid_uniprots, model_info, target_count

#preprocess rdkit molecule
def preprocessMolecule(inp):
	def checkC(mm):
		mwt = Descriptors.MolWt(mm)
		for atom in mm.GetAtoms():
			if atom.GetAtomicNum() == 6 and 100 <= mwt <= 1000: return True
		return False
	def checkHm(mm):
		for atom in mm.GetAtoms():
			if atom.GetAtomicNum() in [2,10,13,18]: return False
			if 21 <= atom.GetAtomicNum() <= 32: return False
			if 36 <= atom.GetAtomicNum() <= 52: return False
			if atom.GetAtomicNum() >= 54: return False
		return True
	try: std_mol = standardise.run(inp)
	except standardise.StandardiseException: return None
	if not std_mol or checkHm(std_mol) == False or checkC(std_mol) == False: return None
	else: return std_mol

#preprocess exception to catch
class MolFromSmilesError(Exception):
    'raise due to "None" from Chem.MolFromSmiles'

#preprocess exception to catch
class PreprocessViolation(Exception):
    'raise due to preprocess violation'

#calculate 2048bit morgan fingerprints, radius 2, for smiles or sdf input
def calcFingerprints(input,qtype='smiles'):
	if qtype == 'smiles': m = Chem.MolFromSmiles(input)
	else: m = input
	if not m: raise MolFromSmilesError(' None mol in function')
	if options.preproc:
		m = preprocessMolecule(m)
		if not m: raise PreprocessViolation(' Molecule preprocessing violation')
	fp = AllChem.GetMorganFingerprintAsBitVect(m,2, nBits=2048)
	binary = fp.ToBitString()
	if qtype == 'sdf': return Chem.MolToSmiles(m), map(int,list(binary)), fp
	else: return map(int,list(binary)), fp

#calculate fingerprints for chunked array of smiles
def arrayFP(inp):
	outfp = []
	outmol = []
	outsmi_id = []
	for idx, i in inp:
		i = i.split(options.delim)
		try:
			fp, mol = calcFingerprints(i[options.smicol])
			outfp.append(fp)
			outmol.append(mol)
			try: outsmi_id.append(i[options.idcol])
			except IndexError: outsmi_id.append(i[options.smicol])
		except PreprocessViolation: print ' SMILES preprocessing violation (line ' + str(idx+1) + ' will be removed): ' + i[options.smicol]
		except MolFromSmilesError: print ' SMILES parse error (line ' + str(idx+1) + '): ' + i[options.smicol]
	return outfp, outmol, outsmi_id

#import user smiles query
def importQuerySmiles(in_file):
	query = open(in_file).read().splitlines()
	query = zip(range(len(query)),query)
	matrix = np.empty((len(query), 2048), dtype=np.uint8)
	num_chunks = min(len(query), options.ncores)
	chunked_smiles = np.array_split(query, num_chunks)
	pool = Pool(processes=num_chunks)  # set up resources
	jobs = pool.imap(arrayFP, chunked_smiles)
	current_end = 0
	processed_mol = []
	processed_id = []
	for i, result in enumerate(jobs):
		percent = (float(i)/float(num_chunks))*100 + 1
		sys.stdout.write(' Processing Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		matrix[current_end:current_end+len(result[0]), :] = result[0]
		current_end += len(result[0])
		processed_mol += result[1]
		processed_id += result[2]
	pool.close()
	pool.join()
	print
	return matrix[:current_end], processed_mol, processed_id

#preprocess exception to catch
class SdfNoneMolError(Exception):
    'raised due to "None" mol during enumeration through Chem.SDMolSupplier'

#import user query for sdf
def importQuerySDF(in_file):
	outfp = []
	outid= []
	outmol = []
	query = Chem.SDMolSupplier(in_file)
	for idx, m in enumerate(suppl):
		sys.stdout.write(' Importing SDF file. Compound number: %s\r' % idx)
		sys.stdout.flush()
		try:
			if not m: raise SdfNoneMolError('None mol')
			smi, fp, mol = calcFingerprints(m,qtype='sdf')
			try: outid.append(m.GetProp("_Name"))
			except KeyError: outid.append(smi)
			outfp.append(fp)
			outmol.append(mol)
		except SdfNoneMolError: print ' SDF parse error (compound index: ' + str(idx) + ')'
	print
	return np.array(outfp,dtype=np.uint8),outmol,outid

#unzip a pkl model
def open_Model(mod):
	global mod_dir, options
	with bz2.BZ2File(mod_dir + 'pkls' + sep + mod + '.pkl.bz2', 'rb') as bzfile:
		clf = cPickle.load(bzfile)
		#if set, increase number of trees in forest
		if options.ntrees and clf.n_estimators < options.ntrees:
			clf.set_params(n_estimators=options.ntrees)
	return clf

#import the training data similarity, bias and standard deviation file for given model
def getAdData(model_name):
	global mod_dir, sep
	actual_mid = model_name.split('/')[-1].split('.pkl.bz2')[0]
	ad_file = mod_dir + 'ad_analysis' + sep + actual_mid + '.pkl.bz2'
	with bz2.BZ2File(ad_file, 'rb') as bzfile:
		ad_data = cPickle.load(bzfile)
	return ad_data

#perform AD analysis using full similarity weighted threshold [similarity/(bias * std_dev]
#adapted from DOI:10.1186/s13321-016-0182-y
def doSimilarityWeightedAdAnalysis(model_name):
	global rdkit_mols, ad_settings
	ad_idx = []
	known = []
	ad_data = getAdData(model_name)
	required_threshold = np.percentile(ad_data[:,5],ad_settings)
	for mol_idx, m in enumerate(rdkit_mols):
		ad_flag = False
		#only check for known compounds if set in options (True means dont check)
		if options.known: k_flag = False
		else: k_flag = True
		for training_instance in ad_data:
			sim = DataStructs.TanimotoSimilarity(m,training_instance[0])
			if sim == 1.0 and k_flag == False:
				known.append([mol_idx,training_instance[1]])
				k_flag = True
			weight = sim/(training_instance[2]*training_instance[3])
			if weight >= required_threshold and ad_flag != True:
				ad_idx.append(mol_idx)
				ad_flag = True
			#if compound is in AD and no need to check accross all comps for known then break
			if k_flag == True and ad_flag == True: break
	return ad_idx, np.array(known)

#return percentile AD analysis for [similarity/(bias * std_dev] vs training data
# def doPercentileSimilarityWeightedAdAnalysis(model_name):
# 	global rdkit_mols,ad_settings
# 	ret = []
# 	ad_data = getAdData(model_name)
# 	for mol_idx, m in enumerate(rdkit_mols):
# 		percentiles = []
# 		for training_instance in ad_data:
# 			sim = DataStructs.TanimotoSimilarity(m,training_instance[0])
# 			weight = sim/(training_instance[2]*training_instance[3])
# 			percentiles.append(percentileofscore(ad_data[:,5],weight))
# 		ret.append(max(percentiles))
# 	return ret

#return percentile AD analysis for [similarity/(bias * std_dev] vs training data
def doPercentileCalculation(model_name):
	global rdkit_mols
	ad_data = getAdData(model_name)
	def calcPercentile(rdkit_mol):
		sims = DataStructs.BulkTanimotoSimilarity(rdkit_mol,ad_data[:,0])
		bias = ad_data[:,2].astype(float)
		std_dev = ad_data[:,3].astype(float)
		scores = ad_data[:,5].astype(float)
		weights = sims / (bias * std_dev)
		critical_weight = weights.max()
		percentile = percentileofscore(scores,critical_weight)
		return percentile
	ret = [calcPercentile(x) for x in rdkit_mols]
	return model_name, ret

#prediction runner for percentile calculation
def performPercentileCalculation(models):
	print " Starting percentile calculation"
	input_len = len(models)
	percentile_results = np.empty(input_len, dtype=object)
	pool = Pool(processes=options.ncores)
	chunksize = max(1, int(input_len / (10 * options.ncores)))
	jobs = pool.imap_unordered(doPercentileCalculation, models, chunksize)
	for i, result in enumerate(jobs):
		progress = str(i+1) + '/' + str(input_len)
		percent = '%3d%%\r' % (float(i)/float(input_len)*100 + 1)
		sys.stdout.write(' Performing percentile calculation: ' + progress + ', ' + percent)
		sys.stdout.flush()
		if result is not None: percentile_results[i] = result
	pool.close()
	pool.join()
	sys.stdout.write(' Performing percentile calculation: ' + progress + ', 100%')
	print "\n Percentile calculation completed"
	return percentile_results

#calculate standard deviation for an input compound
def getStdDev(clf, querymatrix):
	std_dev = []
	for tree in range(len(clf.estimators_)):
		std_dev.append(clf.estimators_[tree].predict_proba(querymatrix)[:,1])
	std_dev = np.clip(np.std(std_dev,axis=0),0.001,None)
	return std_dev

#raw or binary prediction worker
def doTargetPrediction(model_name):
	global ad_settings
	try:
		#percentile ad analysis calculated from [near_neighbor_similarity/(bias*deviation)]
		clf = open_Model(model_name)
		ret = np.zeros(len(querymatrix))
		ret.fill(np.nan)
		try:
			ad_idx, known = doSimilarityWeightedAdAnalysis(model_name)
		except: return model_name, ret
		#if no mols in AD then return
		if len(ad_idx) == 0: return model_name, ret
		probs = clf.predict_proba(querymatrix[ad_idx])[:,1].clip(0.001,0.999)
		#return the standard deviation if turned on
		if options.std:
			std_dev = getStdDev(clf,querymatrix)
			ret[ad_idx] = std_dev[ad_idx]
			return model_name, ret
		#will only have known if was set on
		if len(known) > 0: probs[known[:,0]] = known[:,1]
		if threshold: ret[ad_idx] = map(int,probs > threshold)
		else: ret[ad_idx] = probs
	except IOError: return None
	return model_name, ret

#prediction runner for prediction or standard deviation calculation
def performTargetPrediction(models):
	print " Starting classification"
	input_len = len(models)
	prediction_results = np.empty(input_len, dtype=object)
	pool = Pool(processes=options.ncores, initializer=initPool, initargs=(querymatrix,options.proba,))
	chunksize = max(1, int(input_len / (10 * options.ncores)))
	jobs = pool.imap_unordered(doTargetPrediction, models, chunksize)
	for i, result in enumerate(jobs):
		progress = str(i+1) + '/' + str(input_len)
		percent = '%3d%%\r' % ((float(i)/float(input_len))*100 + 1)
		sys.stdout.write(' Performing classification on query molecules: ' + progress + ', ' + percent)
		sys.stdout.flush()
		if result is not None: prediction_results[i] = result
	pool.close()
	pool.join()
	sys.stdout.write(' Performing classification on query molecules: ' + progress + ', 100%')
	print "\n Classification completed"
	return prediction_results

#write out normal results (rows are targets, columns are compounds)
def writeOutResults(results):
	global query_id, mid_uniprots, model_info, options
	out_data = []
	for mid, preds in results:
		for uniprot_rows in mid_uniprots[mid]:
			if options.train_log: out_data += [uniprot_rows + model_info[mid][1:] + list(preds)]
			else: out_data += [uniprot_rows + list(preds)]
	timestr = time.strftime("%Y%m%d-%H%M%S")
	if not options.off: output_name = options.inf + '_out_predictions_' + timestr + '.txt'
	else: output_name = options.off + '_out_predictions_' + timestr + '.txt'
	print ' Writing predictions to file: ' + output_name
	out_file = open(output_name, 'w')
	header = 'Uniprot\tName\tGene_ID\tProtein_Classification\tOrganism' \
	'\tPDB_IDs\tDisGeNET_0.06\tChEMBL_First_Publication\tActivity_Threshold' \
	'\tActives\tOrthologue_Actives\tInactives\tOrthologue_Inactives' \
	'\tSphere_Exclusion_Inactives\tRatio\tModel_ID\t'
	if options.train_log: out_file.write(header + \
	'\t'.join(map(str,model_info['MODEL_ID'][1:])) + '\t' + '\t'.join(map(str,query_id)) + '\n')
	else:
		out_file.write(header + '\t'.join(map(str,query_id)) + '\n')
	for row in out_data: out_file.write('\t'.join(map(str,row)) + '\n')
	out_file.close()
	return

#write out transposed results (columns are targets, rows are compounds)
def writeOutTransposed(results):
	global query_id, mid_uniprots, model_info, options
	out_data = [['Uniprot'] + query_id]
	for mid, preds in results:
		for uniprot_rows in mid_uniprots[mid]:
			out_data += [[uniprot_rows[0]] + list(preds)]
	timestr = time.strftime("%Y%m%d-%H%M%S")
	if not options.off: output_name = options.inf + '_out_predictions_' + timestr + '.txt'
	else: output_name = options.off + '_out_predictions_' + timestr + '.txt'
	print ' Writing out transposed predictions to file: ' + output_name
	out_file = open(output_name, 'w')
	out_data = np.transpose(out_data)
	for row in out_data: out_file.write('\t'.join(map(str,row)) + '\n')
	out_file.close()
	return

#nt (Windows) compatibility initializer for the pool
def initPool(querymatrix_, threshold_=None):
	global querymatrix, threshold
	querymatrix = querymatrix_
	threshold = threshold_

#main
if __name__ == '__main__':
	introMessage()
	print ' Predicting Targets for input: ' + options.inf
	print ' Using ' + str(options.ncores) + ' core(s)'
	if options.ntrees: print 'Latency warning: Number of trees will be increased to minimum: ' + str(options.ntrees)
	print ' Using bioactivity thresholds (IC50/EC50/Ki/Kd) of: ' + str(options.bioactivity)
	print ' Using orthologues: ' + str(options.ortho)
	if options.organism: print ' Organism filter: ' + options.organism
	if options.targetclass: print ' Target class filter: ' + options.targetclass
	if options.minsize: print ' Minimum number of actives in training: ' + str(options.minsize)
	if options.se_filt: print ' Filtering out models with Sphere Exclusion (SE)'
	if options.train_log: print ' Training log will be added to output'
	#set up environment
	input_extension, sep, pidgin_dir, mod_dir, ad_settings = check_set_working_env()
	#if preprocessing import standardizer
	if options.preproc: from standardiser import standardise
	#gather the models required and their information
	mid_uniprots, model_info, target_count = get_Models()
	print ' Total number of protein targets: ' + str(target_count)
	print ' Total number of distinct models: ' + str(len(mid_uniprots))
	print ' Using p(activity) threshold of: ' + str(options.proba)
	print ' Importing query (calculating ECFP_4 fingerprints)'
	#import user query files
	if input_extension == 'sdf': querymatrix,rdkit_mols,query_id = importQuerySDF(options.inf)
	else:  querymatrix,rdkit_mols,query_id = importQuerySmiles(options.inf)
	print ' Total number of query molecules: ' + str(len(querymatrix))
	#perform target prediction on (filtered) models (using model ids)
	if options.percentile:
		prediction_results = performPercentileCalculation(mid_uniprots.keys())
	else:
		prediction_results = performTargetPrediction(mid_uniprots.keys())
	#write out
	if options.transpose: writeOutTransposed(prediction_results)
	else: writeOutResults(prediction_results)
