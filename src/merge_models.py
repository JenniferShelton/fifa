## Need to change data + model paths to relative paths / 
## Add them to the github repo 

import os
import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score, precision_score, f1_score, recall_score
import multiprocessing as mp
from multiprocessing import Pool
from json import dump
from scipy.stats import trim_mean
from sklearn.preprocessing import FunctionTransformer, StandardScaler
import sys
import glob
import logging
from typing import Any, List, Optional

from interpret.glassbox import merge_ebms as interpret_merge_ebms
from classify_with_scaling import get_zero_score_names, remove_zero_score_terms_by_name # in future can move these functions to a utils file
from helper_funcs import align_model_features, convert_hot_encodings, scale_features

np.random.seed(42)
label_encoder = LabelEncoder() 
logger = logging.getLogger(__name__)

def merge_ebms(
    input_models: List[str], output_path: Optional[str] = None
) -> Any:
    if not input_models:
        raise ValueError("At least one model path must be provided.")

    fifa_install_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bundled_models = {
        "NYGC1": "ebm_hyperparams_NYGC1.pkl",
        "NYGC2": "ebm_hyperparams_NYGC2.pkl",
        "BLGSP": "ebm_hyperparams_CGCI-BLGSP.pkl",
        "CGCI-BLGSP": "ebm_hyperparams_CGCI-BLGSP.pkl",
        "HTMCP": "ebm_hyperparams_CGCI-HTMCP.pkl",
        "CGCI-HTMCP": "ebm_hyperparams_CGCI-HTMCP.pkl",
    }
    models: List[Any] = []
    for model_path in input_models:
        filename = bundled_models.get(model_path)
        path = (
            os.path.join(fifa_install_dir, "models", filename)
            if filename
            else model_path
        )
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Model file {path} does not exist.")
        with open(path, "rb") as model_file:
            models.append(pickle.load(model_file))

    reference_feature_names = list(models[0].feature_names_in_)
    for model_path, model in zip(input_models, models):
        align_model_features(model, reference_feature_names, model_label=model_path)

    if len(models) > 1:
        zero_score_names = get_zero_score_names(models)
        models = [
            remove_zero_score_terms_by_name(model, zero_score_names)
            for model in models
        ]
        ebm = interpret_merge_ebms(models)
    else:
        ebm = models[0]

    if output_path: 
        logger.info(f'Output path provided; Merged model saved to: {output_path}')
        with open(output_path, 'wb') as file:
            pickle.dump(ebm, file)


    ## Test in HCC1395
    directory = os.path.join(os.path.dirname(os.getcwd()), 'data/')
    pattern = "FFG*_extracted_features.csv"

    matching_files = glob.glob(directory + pattern, recursive=True)
    try: 
        features = pd.concat((pd.read_csv(f).fillna(0) for f in matching_files), ignore_index=True)
        features = features.drop(['Cohort_x', 'Cohort_y', '0'], axis=1)
        features['Cohort'] = 'HCC1395'
        features['tumor_depth'] = 0
    except Exception as e:
        logger.error(f"An error occurred: {e}")
    labels_path=os.path.join(directory, "real_labels.csv")
    labels = pd.read_csv(labels_path, sep=',').fillna(0)
    labels['Label'] = labels['Label'].apply(lambda x: 1 if x == 'Real' else 0)
    variants = pd.merge(features, labels, on=['Sample', 'Variant'], how='inner')

    variants = convert_hot_encodings(variants)
    variants = scale_features(variants)

    logger.info('Uploaded cell line test samples')

    ## Generate PDF plots of the performance of the merged model on HCC1395
    ## Do these need to be plots ? Or can they be just summary metrics ? 
    ## Could probably also generate a confusion matrix here too

    ## Testing on chr1 variants from HCC1395 with Mutect2
    logger.info(f'Metrics in HCC1395 Cell Line Replicates (chr1 variants)')
    for sample in variants[variants['Cohort'] == 'HCC1395']['Sample'].unique():
        testing_set = variants[(variants['Sample'] == sample) & (variants['Variant'].str.startswith('chr1:'))].drop(['Sample', 'Variant', 'Cohort'], axis=1)

        X_test = testing_set.drop(['Label'], axis=1)
        y_test = label_encoder.fit_transform(testing_set['Label'])
        
        y_pred = ebm.predict(X_test)
        
        logger.info(f'Mutect2 calls from {sample}')
        f1 = f1_score(y_test, y_pred)
        logger.info("F1: {:.3f}".format(f1))
        
        accuracy = accuracy_score(y_test, y_pred)
        logger.info("Accuracy: {:.3f}".format(accuracy))
        
        precision = precision_score(y_test, y_pred)
        logger.info("Precision: {:.3f}".format(precision))
        
        recall = recall_score(y_test, y_pred)
        logger.info("Recall: {:.3f}".format(recall))

    ## Testing on DeepSomatic's chr 1 variants
    ## Again this is code that I can clean up
    file_path = os.path.join(directory, "1395_tumor_ffpe_wgs_extracted_features.csv")
    features=pd.read_csv(file_path).fillna(0)
    labels_path=os.path.join(directory, "labels.txt")

    column_names=['CHROM', 'POS', 'REF', 'ALT', 'Label']

    labels = pd.read_csv(labels_path, sep='\t', names=column_names).fillna(0)
    labels['Variant'] = labels.apply(lambda row: f"{row['CHROM']}:{row['POS']}_{row['REF']}>{row['ALT']}", axis=1)
    labels = labels.drop(['CHROM', 'POS', 'REF', 'ALT'], axis=1)
    labels['Sample'] = '1395_tumor_ffpe_wgs'
    labels['Label'] = labels['Label'].apply(lambda x: 1 if x == 'Real' else 0)

    variants = pd.merge(features, labels, on=['Sample', 'Variant'], how='inner')

    ## Convert hot encodings to tri- and penta-nuc context features (if applicable); 
    ## Scale relevant features 
    if 'hot_encoded_ref_base' in variants.columns:
        variants = convert_hot_encodings(variants)
    variants = scale_features(variants)

    for sample in variants[variants['Cohort'] == 'DeepSomatic']['Sample'].unique():
        testing_set = variants[(variants['Sample'] == sample) & (variants['Variant'].str.startswith('chr1:'))].drop(['Sample', 'Variant', 'Cohort'], axis=1)

        X_test = testing_set.drop(['Label'], axis=1)
        y_test = label_encoder.fit_transform(testing_set['Label'])
        y_pred = ebm.predict(X_test)

        logger.info(f'Metrics in DeepSomatic\'s Held-Out Test Set (chr1 variants)')
        f1 = f1_score(y_test, y_pred)
        logger.info("F1: {:.3f}".format(f1))
        
        accuracy = accuracy_score(y_test, y_pred)
        logger.info("Accuracy: {:.3f}".format(accuracy))
        
        precision = precision_score(y_test, y_pred)
        logger.info("Precision: {:.3f}".format(precision))
        
        recall = recall_score(y_test, y_pred)
        logger.info("Recall: {:.3f}".format(recall))

