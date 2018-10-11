#Author : Lewis Mervin lhm30@cam.ac.uk
#Supervisor : Dr. A. Bender
#All rights reserved 2018
#Protein Target Prediction using on SAR data from PubChem and ChEMBL_24
#Molecular Descriptors : 2048bit circular Binary Fingerprints (Rdkit) - ECFP_4
#Dependencies : rdkit, sklearn, standardiser

### predict_enriched.py ###
#Calculate target, pathway and disease-gene association enrichment (Fishers' exact t-test)
#for one (vs. a background) or two smiles/sdf files. This analysis must use a cut-off
#for probability [computed as the mean predicted class probabilities of
#the trees in the forest (where the class probability of a single tree is the fraction of 
#samples of the same class in a leaf)]. Predictions are generated 
#for the [filtered] models using a reliability-density neighbourhood Applicability Domain
#(AD) analysis from: doi.org/10.1186/s13321-016-0182-y

#libraries
import os
import sys
import bz2
import zipfile
import cPickle
import glob
import time
import math
import operator
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Descriptors
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
from scipy.stats import percentileofscore
from scipy.stats import chi2_contingency
import multiprocessing
from multiprocessing import Pool
from optparse import OptionParser

#optionparser options
parser = OptionParser()
parser.add_option("--f1", dest="inf1", help="Firest input smiles or sdf file (required)", metavar="FILE")
parser.add_option("--f2", dest="inf2", default=None, help="Second input smiles or sdf file (optional)", metavar="FILE")
parser.add_option("-d", "--smiles_delim", default=' ', type=str, dest="delim", help="Input file (smiles) delimiter char (default: white space ' ')")
parser.add_option("--smiles_column", default=0, type=int, dest="smicol", help="Input file (smiles) delimiter column (default: 0)")
parser.add_option("--smiles_id_column", default=1, type=int, dest="idcol", help="Input file (smiles) ID column (default: 1)")
parser.add_option("-o", dest="off", default=None, help="Optional output prediction file name", metavar="FILE")
parser.add_option("-n", "--ncores", default=1, type=int, dest="ncores", help="No. cores (default: 1)")
parser.add_option("-b", "--bioactivity", default='10', type=str, dest="bioactivity", help="Bioactivity Um threshold (required). Use either 100/10/1/0.1 (default:10)")
parser.add_option("-p", "--proba", default=0.5, type=float, dest="proba", help="RF probability threshold (default: None)")
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
parser.add_option("--dgn", default=0.06, type=float, dest="dgn_threshold", help="DisGeNET score threshold (default: 0.06)")

(options, args) = parser.parse_args()
#if tab delimited then set to \t
if options.delim == 'tab': options.delimcol1 = '\t'

def introMessage():
	print '=============================================================================================='
	print ' Author: Lewis Mervin\n Email:  lhm30@cam.ac.uk\n Supervisor: Dr. A. Bender'
	print ' Address: Centre For Molecular Informatics, Dept. Chemistry, Lensfield Road, Cambridge CB2 1EW'
	print '==============================================================================================\n'
	return

#check smiles of sdf (if not warn)	
def check_Input(inp):
	global options
	extension = inp.split('.')[-1]
	if extension not in ['smi','smiles','sdf']:
		print ' Warning [--f1/f2]: File' + inp + ' not "smi", "smiles" or "sdf". Interpreting input as SMILES'
	return extension

#check input & set OS directory pointers
def check_set_working_env():
	global options
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	input_extension1 = check_Input(options.inf1)
	try: input_extension2 = check_Input(options.inf2)
	except AttributeError: input_extension2 = None
	try:
		ad_settings = int(options.ad)
		if 0 >= ad_settings >= 100: raise ValueError
	except ValueError:
		print ' Input Error [--ad]: Percentile weight not integer between 0-100%. Please Check parameters'
		sys.exit()
	pidgin_dir = os.path.dirname(os.path.abspath(__file__))
	if options.ortho:
		mod_dir = pidgin_dir + sep + 'ortho' + sep
		if not os.path.isdir(mod_dir):
			print ' Orthologue Error [--orthologues]: Orthologues directory is not present. Please download from here: '
	else: mod_dir = pidgin_dir + sep + 'no_ortho' + sep
	return input_extension1, input_extension2, sep, pidgin_dir, mod_dir, ad_settings

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

#get info for diseases
def getDisgenetInfo():
	global pidgin_dir
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	return_dict1 = dict()
	return_dict2 = dict()
	disease_file = [l.split('\t') for l in open(pidgin_dir + sep + 'DisGeNET_diseases.txt').read().splitlines()]
	for l in disease_file:
		try:
			return_dict1[l[0]].append(l[1])
		except KeyError:
			return_dict1[l[0]] = [l[1]]
		try:
			return_dict2[(l[1],l[0])] = float(l[2])
		except ValueError: pass
	return return_dict1, return_dict2 
	
#get info for biosystems pathways
def getPathwayInfo():
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	return_dict1 = dict()
	return_dict2 = dict()
	pathway_info = [l.split('\t') for l in open(pidgin_dir + sep + 'biosystems.txt').read().splitlines()]
	for l in pathway_info:
		try:
			return_dict1[l[0]].append(l[1])
		except KeyError:
			return_dict1[l[0]] = [l[1]]
		return_dict2[l[1]] = l[2:]
	return return_dict1, return_dict2

#get pre-calculated bg hits from PubChem
def getBGhits(threshold):
	if os.name == 'nt': sep = '\\'
	else: sep = '/'
	bg_column = int((threshold*100)+1)
	bg_file = [l.split('\t') for l in open(os.path.dirname(pidgin_dir + sep + 'bg_predictions.txt').read().splitlines())]
	bg_file.pop(0)
	bg_predictions = {l[0] : int(l[bg_column]) for l in bg_file}
	return bg_predictions

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
	for idx, i in inp:
		i = i.split(options.delim)
		try:
			fp, mol = calcFingerprints(i[options.smicol])
			outfp.append(fp)
			outmol.append(mol)
		except PreprocessViolation: print ' SMILES preprocessing violation (line ' + str(idx+1) + ' will be removed): ' + i[options.smicol]
		except MolFromSmilesError: print ' SMILES parse error (line ' + str(idx+1) + '): ' + i[options.smicol]
	return outfp, outmol

#import user smiles query
def importQuerySmiles(in_file):
	query = open(in_file).read().splitlines()
	query = zip(range(len(query)),query)
	matrix = np.empty((len(query), 2048), dtype=np.uint8)
	smiles_per_core = int(math.ceil(len(query) / options.ncores)+1)
	chunked_smiles = [query[x:x+smiles_per_core] for x in xrange(0, len(query), smiles_per_core)]
	pool = Pool(processes=options.ncores)  # set up resources
	jobs = pool.imap(arrayFP, chunked_smiles)
	current_end = 0
	processed_mol = []
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(chunked_smiles)))*100 + 1
		sys.stdout.write(' Processing Molecules for ' + in_file + ': %3d%%\r' % percent)
		sys.stdout.flush()
		matrix[current_end:current_end+len(result[0]), :] = result[0]
		current_end += len(result[0])
		processed_mol += result[1]
	sys.stdout.write(' Processing Molecules for ' + in_file + ': %3d%%\r' % float(100))
	sys.stdout.flush()
	pool.close()
	pool.join()
	print
	return matrix[:current_end], processed_mol
	
#preprocess exception to catch
class SdfNoneMolError(Exception):
    'raised due to "None" mol during enumeration through Chem.SDMolSupplier'

#import user query for sdf
def importQuerySDF(in_file):
	outfp = []
	outmol = []
	query = Chem.SDMolSupplier(in_file)
	for idx, m in enumerate(suppl):
		sys.stdout.write(' Importing SDF file. Compound number: %s\r' % idx)
		sys.stdout.flush()
		try:
			if not m: raise SdfNoneMolError('None mol')
			smi, fp, mol = calcFingerprints(m,qtype='sdf')
			outfp.append(fp)
			outmol.append(mol)
		except SdfNoneMolError: print ' SDF parse error (compound index: ' + str(idx) + ')'
	print
	return np.array(outfp,dtype=np.uint8),outmol

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

#perform AD analysis using full similarity weighted threshold [similarity/(bias * std_dev)]
#adapted from DOI:10.1186/s13321-016-0182-y
def doSimilarityWeightedAdAnalysis(model_name, rdkit_mols):
	global ad_settings
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

#return target prediction information for first and second files (p1,p2)
def calcProportionPredsPercentageFisher(model_name,p1,p2):
	#calculate proportion of nan (outside ad) compounds in each set
	propnanp1 = round((float(sum(np.isnan(p1))) / float(len(p1))),3)
	propnanp2 = round((float(sum(np.isnan(p2))) / float(len(p2))),3)
	#filter out nans
	p1 = p1[np.isfinite(p1)]
	p2 = p2[np.isfinite(p2)]
	#calculate active and inactive predictions for first file
	p1_0, p1_1 = len(p1)-sum(p1), sum(p1)
	#calculate active and inactive predictions for second file
	p2_0, p2_1 = len(p2)-sum(p2), sum(p2)
	#return if no compounds in either set due to ad filter
	if p1_0+p1_1 ==0 or p2_0+p2_1 ==0: return None
	#if no actives in either set return none
	if p1_1 == 0 and p2_1 == 0: return None
	oddsratio, pvalue = fisher_exact([[p1_1,p1_0],[p2_1,p2_0]])
	#if no inactives in either set then set risk & odds to 1.0 and pval to NA
	if p1_0 == 0 and p2_0 == 0: rr, oddsratio, pvalue = 1.0, 'NA'
	else:
		#calculate risk ratio
		try: rr = (float(p1_1)/(float(p1_1)+float(p1_0)))/(float(p2_1)/(float(p2_1)+float(p2_0)))
		except ZeroDivisionError: rr = np.inf
	#calculate percentages
	pcp1_0, pcp1_1 = round(float(p1_0)/float(len(p1)),3), round(float(p1_1)/float(len(p1)),3)
	pcp2_0, pcp2_1 = round(float(p2_0)/float(len(p2)),3), round(float(p2_1)/float(len(p2)),3)
	return oddsratio, model_name, p1_0, pcp1_0, p1_1, pcp1_1, p2_0, pcp2_0, \
	p2_1, pcp2_1, propnanp1, propnanp2, rr, pvalue

def get_bg_preds(inp):
	return 50,50,50,50

#raw or binary prediction worker
def doTargetPrediction(model_name):
	try:
		clf = open_Model(model_name)
		ret = np.zeros(len(querymatrix))
		ret.fill(np.nan)
		#percentile ad analysis calculated from [near_neighbor_similarity/(bias*deviation)]
		ad_idx, known = doSimilarityWeightedAdAnalysis(model_name,rdkit_mols)
		#if no mols in AD then return
		if len(ad_idx) == 0: return model_name, ret
		probs = clf.predict_proba(querymatrix[ad_idx])[:,1]
		#will only have known if was set on
		if len(known) > 0: probs[known[:,0]] = known[:,1]
		ret[ad_idx] = map(int,probs >= threshold)
	except IOError: return None
	#if have second file
	if n_inf2 > 0:
		preds1 = ret[:len(ret)-n_inf2]
		preds2 = ret[-n_inf2:]
		pred_info = calcProportionPredsPercentageFisher(model_name,preds1,preds2)
	#else get bg preds
	else:
		pred_info = calcProportionPredsPercentageFisher(model_name,ret,get_bg_preds(model_name))
	if pred_info is None: return None
	return pred_info

#prediction runner for prediction or standard deviation calculation
def performTargetPrediction(models):
	global disease_links, disease_hits, pathway_links, pathway_hits
	prediction_results = []
	pool = Pool(processes=options.ncores, initializer=initPool, initargs=(querymatrix,rdkit_mols,options.proba,ad_settings,n_inf2))
	jobs = pool.imap_unordered(doTargetPrediction, models)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(models)))*100 + 1
		sys.stdout.write(' Performing Classification on Query Molecules: %3d%%\r' % percent)
		sys.stdout.flush()
		if result is not None:
			prediction_results.append(result)
			#update hits for each of the uniprots linked to mid
			for uniprot_row in mid_uniprots[result[1]]:
				updateHits(disease_links,disease_hits,uniprot_row[0],result[4],result[8])
				updateHits(pathway_links,pathway_hits,uniprot_row[0],result[4],result[8])
	pool.close()
	pool.join()
	if len(prediction_results) == 0:
		print '\n Warning: No active target predictions...quitting'
		sys.exit()
	return sorted(prediction_results,reverse=True)

def initPool(querymatrix_, rdkit_mols_, threshold_, ad_settings_, n_inf2_):
	global querymatrix, rdkit_mols, threshold, ad_settings, n_inf2
	querymatrix = querymatrix_
	threshold = threshold_
	ad_settings = ad_settings_
	rdkit_mols = rdkit_mols_
	n_inf2 = n_inf2_

#update counts for each pathway/disease that is hit by predictions	
def updateHits(links,hits,uniprot,hit1,hit2):
	try:
		for idx in links[uniprot]:
			#try checks if pw or dnet
			try:
				if disease_score[(idx,uniprot)] < options.dgn_threshold: continue
			except KeyError: pass
			try:
				hits[idx] = hits[idx] + np.array([hit1,hit2])
			except KeyError:
				hits[idx] = np.array([hit1,hit2])
	except KeyError: return
	return

#worker for the processHits to calculate the chi2 in parallel
def doHitProcess(inp):
	idx, hits, n_f1_hits, n_f2_hits = inp
	p1_0, p1_1 = n_f1_hits-hits[0], hits[0]
	p2_0, p2_1 = n_f2_hits-hits[1], hits[1]
	#if no actives in either set return
	if p1_1 == 0 and p2_1 == 0: return
	#calculate percentage of hits for file1 and file2
	pcp1_1 = float(p1_1)/float(p1_0)
	pcp2_1 = float(p2_1)/float(p2_0)
	#if no inactives in either set, set chi2 to 1.0 and pvalue to 0
	if p1_0 == 0 and p2_0 == 0: return 1.0, idx, p1_1, pcp1_1, p2_1, pcp2_1, 1.0, 'NA'
	chi, pvalue = chi2_contingency([[p1_1,p1_0],[p2_1,p2_0]])[:2]
	#calculate odds ratio
	try: odr = (float(p1_1)/float(p1_0))/(float(p2_1)/float(p2_0))
	except ZeroDivisionError: odr = np.inf
	#calculate risk ratio
	try: rr = (float(p1_1)/(float(p1_1)+float(p1_0)))/(float(p2_1)/(float(p2_1)+float(p2_0)))
	except ZeroDivisionError: rr = np.inf
	return odr, idx, p1_1, pcp1_1, p2_1, pcp2_1, rr, pvalue
	
#calculate the chi2 and odds ratio between pathway and disease predictions
def processHits(inp_dict):
	out_data = []
	total_hits = np.array(inp_dict.values()).sum(axis=0)
	if total_hits.shape is (): return out_data, 0, 0
	n_f1_hits = total_hits[0]
	n_f2_hits = total_hits[1]
	tasks = [[idx,hits,n_f1_hits,n_f2_hits] for idx, hits in inp_dict.iteritems()]
	pool = Pool(processes=options.ncores)  # set up resources
	jobs = pool.imap_unordered(doHitProcess, tasks)
	for i, result in enumerate(jobs):
		percent = (float(i)/float(len(tasks)))*100 + 1
		sys.stdout.write(' Calculating Chi-square test: %3d%%\r' % percent)
		sys.stdout.flush()
		if result is None: continue
		out_data.append(result)
	out_data=np.array(sorted(out_data,reverse=True))
	pvals = np.copy(out_data[:,-1])
	if sum(pvals!='NA') > 2:
		pvals[np.where(pvals!='NA')] = \
		multipletests(pvals[np.where(pvals!='NA')].astype(float),alpha=0.05,method="fdr_bh")[1]
	out_data = np.append(out_data,pvals.reshape((-1,1)),axis=1)
	return out_data, n_f1_hits, n_f2_hits

#write out TP results (rows are targets, columns are compounds)
def writeOutTPResults(results):
	global mid_uniprots, model_info, options, timestr
	out_data = []
	for res in results:
		mid = res[1]
		for uniprot_rows in mid_uniprots[mid]:
			if options.train_log: out_data += [uniprot_rows + model_info[mid][1:] + list(res[2:])+[res[0]]]
			else: out_data += [uniprot_rows + list(res[2:])+[res[0]]]
	out_data = np.array(out_data)
	if not options.off:
		if options.inf2:
			output_name = options.inf1 + '_vs_' + options.inf2.split(sep)[-1] + '_out_predictions_enriched' + timestr + '.txt'
		else: options.inf1 + '_out_predictions_enriched' + timestr + '.txt'
	else: output_name = options.off + '_out_predictions_enriched' + timestr + '.txt'
	print ' Writing predictions to file: ' + output_name
	f1 = options.inf1
	if options.inf2: f2 = options.inf2
	else: f2 = 'Background'
	out_file = open(output_name, 'w')
	p=str(options.proba)
	a=str(options.ad)
	header = 'Uniprot\tName\tGene_ID\tProtein_Classification\tOrganism' \
	'\tPDB_IDs\tDisGeNET_0.06\tChEMBL_First_Publication\tActivity_Threshold' \
	'\tActives\tOrthologue_Actives\tInactives\tOrthologue_Inactives' \
	'\tSphere_Exclusion_Inactives\tRatio\tModel_ID\t'+f1+'_Inactives_'+p+ \
	'\tPercent_'+f1+'_Inactives_'+p+'\t'+f1+'_Actives_'+p+'\tPercent_'+f1+'_Actives_'+p+ \
	'\t'+f2+'_Inactives_'+p+'\tPercent_'+f2+'_Inactives_'+p+'\t'+f2+'_Actives_'+p+ \
	'\tPercent_'+f2+'_Actives_'+p+'\t'+f1+'_Proportion_Nan_Predictions_'+a+ \
	'\t'+f2+'_Proportion_Nan_Predictions_'+a+'\tRisk_Ratio\tFishers_Test_P_Value' \
	'\tFishers_Test_P_Value_Corrected\tOdds_Ratio\n'
	if options.train_log: out_file.write(header + \
	'\t'.join(map(str,model_info['MODEL_ID'][1:])) + '\n')
	else:
		out_file.write(header)
	#perform FDR Benjamini-Hochberg p-value correction on non-nan pvalues
	pvals = np.copy(out_data[:,-2])
	if sum(pvals!='NA') > 2:
		pvals[np.where(pvals!='NA')] = \
		multipletests(pvals[np.where(pvals!='NA')].astype(float),alpha=0.05,method="fdr_bh")[1]
	for idx, row in enumerate(out_data): \
	out_file.write('\t'.join(map(str,list(row[:-1]) + [pvals[idx]] + [row[-1]])) + '\n')
	out_file.close()
	return

#write out disgenet or pathway predictions
def writeOutResults(inp,name):
	global options, timestr, pathway_info
	if not options.off:
		if options.inf2:
			output_name = options.inf1 + '_vs_' + options.inf2.split(sep)[-1] + '_out_'+name+'_predictions_enriched' + timestr + '.txt'
		else: options.inf1 + '_out_'+name+'_predictions_enriched' + timestr + '.txt'
	else: output_name = options.off + '_out_'+name+'_predictions_enriched' + timestr + '.txt'
	print ' Writing '+name+' predictions to file: ' + output_name
	f1 = options.inf1
	if options.inf2: f2 = options.inf2
	else: f2 = 'Background'
	out_file = open(output_name, 'w')
	processed, inp1_total, inp2_total = processHits(inp)
	if name == 'disease':
		out_file.write('Disease_Name\t'+f1+'_Hits_p'+str(options.proba)+'\t' \
		+f1+'_Precent_Hits_p' + str(options.proba)+'\t'+f2+'_Hits_p'+str(options.proba) \
		+'\t'+f2+'_Precent_Hits_p'+str(options.proba)+'\tRisk_Ratio\tChi2_pValue' \
		'\tChi2_pValue_Corrected\tOdds_ratio\n')
		for row in processed:
			out_file.write('\t'.join(map(str,row[1:])) + '\t' + str(row[0]) + '\n')
	if name == 'pathway':
		out_file.write('Pathway_ID\tPathway_Name\tSource\tClass' \
		'\t'+f1+'_Hits_p'+str(options.proba)+'\t'+f1+'_Precent_Hits_p'+str(options.proba) \
		+'\t'+f2+'_Hits\t'+f2+'_Precent_Hits'\
		'\tRisk_Ratio\tChi2_pValue\tCorrected_Chi2_pValue\tChi2_Test_Statistic\n')
		for row in processed:
			out_file.write(row[1] + '\t' + '\t'.join(map(str,pathway_info[row[1]])) \
			+ '\t' + '\t'.join(map(str,row[2:])) + '\t' + str(row[0]) + '\n')
	print '\n Wrote '+name+' results to: ' + output_name
	out_file.close()

#main
if __name__ == '__main__':
	introMessage()
	if options.inf2: print ' Predicting targets for input: ' + options.inf1 + ' and ' + options.inf2
	else:
		print ' Please provide second file [--f2]: Background predictions not yet supported...quiting'
		sys.exit()
		print ' Predicting targets for input : ' + options.inf1 + ' and a background'
	print ' Using ' + str(options.ncores) + ' core(s)'
	if options.ntrees: print ' Warning: Number of trees will be increased to minimum: ' + str(options.ntrees)
	print ' Using bioactivity thresholds (IC50/EC50/Ki/Kd) of : ' + str(options.bioactivity)
	print ' Using orthologues: ' + str(options.ortho)
	if options.organism: print ' Organism filter : ' + options.organism
	if options.targetclass: print ' Target class filter : ' + options.targetclass
	if options.minsize: print ' Minimum number of actives in training : ' + str(options.minsize)
	if options.se_filt: print ' Filtering out models with Sphere Exclusion (SE)'
	if options.train_log: print ' Training log will be added to output'
	#set up environment
	input_extension1, input_extension2, sep, pidgin_dir, mod_dir, ad_settings = check_set_working_env()
	disease_links, disease_score = getDisgenetInfo()
	pathway_links, pathway_info = getPathwayInfo()
	disease_hits, pathway_hits = dict(), dict()
	#if preprocessing import standardizer
	if options.preproc: from standardiser import standardise
	#gather the models required and their information
	mid_uniprots, model_info, target_count = get_Models()
	print ' Total number of protein targets: ' + str(target_count)
	print ' Total number of distinct models: ' + str(len(mid_uniprots))
	print ' Using p(activity) threshold of: ' + str(options.proba)
	print ' Importing query (calculating ECFP_4 fingerprints)'
	#import user query files
	if input_extension1 == 'sdf': querymatrix,rdkit_mols = importQuerySDF(options.inf1)
	else:  querymatrix,rdkit_mols = importQuerySmiles(options.inf1)
	print ' Total number of ' + options.inf1 + ' query molecules : ' + str(len(querymatrix))
	if options.inf2:
		if input_extension2 == 'sdf': querymatrix2,rdkit_mols2 = importQuerySDF(options.inf2)
		else:  querymatrix2,rdkit_mols2 = importQuerySmiles(options.inf2)
		n_inf2 = len(querymatrix2)
		print ' Total number of ' + options.inf1 + ' query molecules : ' + str(len(querymatrix2))
		querymatrix = np.vstack((querymatrix,querymatrix2))
		rdkit_mols = rdkit_mols + rdkit_mols2
	else: n_inf2 = 0
	timestr = time.strftime("%Y%m%d-%H%M%S")
	#perform target prediction on (filtered) models (using model ids)
	prediction_results = performTargetPrediction(mid_uniprots.keys())
	#write out target predictions
	writeOutTPResults(prediction_results)
	#write out disease and pathway results
	writeOutResults(disease_hits,'disease')
	writeOutResults(pathway_hits,'pathway')
