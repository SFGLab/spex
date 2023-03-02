import argparse
import numpy as np
import pandas as pd
import logging
import os.path
import hashlib
import xgboost as xgb
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import mean_squared_error
import random
from multiprocessing import Pool
import math, sys

from args import add_args, setup, is_mbp
from configs import load_genes

from filenames import get_tensor_filename

# workaround: https://github.com/dmlc/xgboost/issues/1715#issuecomment-420305786
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'

log_full_predictions = True

def get_no_filename_and_idx(genes, data_folder):
    i, nos_filenames_idxs = 0, list()
    for gene in genes.itertuples():
        tensor_filename = get_tensor_filename(gene, data_folder)
        if os.path.isfile(tensor_filename):
            nos_filenames_idxs.append((i, tensor_filename, gene.Index))
            i += 1
    return nos_filenames_idxs


def read_tensor(no_filename_and_idx):
    i, tensor_filename, index = no_filename_and_idx
    tensor = np.load(tensor_filename)
    return i, tensor, index


def read_tensors(nos_filename_and_idx, data, dataind):
    n = 0
    pool = Pool(8)
    for i, tensor, index in pool.imap(read_tensor, nos_filename_and_idx):
        data[i] = tensor
        n += 1
        dataind.append(index)
    pool.close()
    pool.join()
    return n


def read_train_data(args):
    # load protein coding genes and Chromatin Contact Domains (CCDs), find CCDs containing genes
    genes = load_genes(args.genes_filename, args.gene_type, args.data_folder)
    genes = genes[~genes.seqnames.isin(["chrX", "chrY"])]

    train_filenames_idxs = get_no_filename_and_idx(genes[genes.seqnames != "chr8"], args.data_folder)
    test_filenames_idxs = get_no_filename_and_idx(genes[genes.seqnames == "chr8"], args.data_folder)

    tensor_width = np.load(train_filenames_idxs[0][1]).shape[0]
    train = np.zeros(shape=(len(train_filenames_idxs), tensor_width))
    test = np.zeros(shape=(len(test_filenames_idxs), tensor_width))
    trainind, testind = list(), list()
    train_n, test_n = 0, 0

    # for each gene read tensor and add to appropriate array
    logging.info("Loading train tensors...")
    train_n = read_tensors(train_filenames_idxs, train, trainind)
    logging.info("Loading test tensors...")
    test_n = read_tensors(test_filenames_idxs, test, testind)
    
    # loadig gene expression files
    geneexp = pd.read_csv(args.gene_expression_filename)

    # prepare training data
    train, test = train[:train_n], test[:test_n]
    logging.info(f"Done. Train: {train.shape} Test: {test.shape}")

    dtrain = xgb.DMatrix(train)
    dtest = xgb.DMatrix(test)

    dtrain.set_label(np.asarray(
        np.log(geneexp.iloc[trainind, args.target_index] + args.pseudocount)))
    dtest.set_label(np.asarray(
        np.log(geneexp.iloc[testind, args.target_index] + args.pseudocount)))

    ytest = np.asarray(np.log(geneexp.iloc[testind, args.target_index] + args.pseudocount))

    genes_ids = [x[2] for x in test_filenames_idxs]
    
    return train_n + test_n == len(genes), dtrain, dtest, ytest, genes.loc[testind][["id", "symbol"]]


def train(dtrain, dtest, ytest, params, num_round, early_stopping_round):
    # training
    evallist = [(dtrain, 'train'), (dtest, 'eval')]
    bst = xgb.train(params, dtrain, num_round, evallist, verbose_eval=False,
                    callbacks=[xgb.callback.EarlyStopping(rounds=early_stopping_round, save_best=True)])

    ypred = bst.predict(dtest)
    if np.isinf(ypred).any() or np.isnan(ypred).any():
        return None, None, None, None
    rmse = mean_squared_error(ypred, ytest, squared=False)
    sm = spearmanr(ypred, ytest)
    pr = pearsonr(ypred, ytest)

    return rmse, sm.correlation, sm.pvalue, pr[0], pr[1], bst.best_iteration, bst, ypred, ytest


def nanstr(x):
    return "NaN" if math.isnan(x) else x

def sha256sum(filename):
    h  = hashlib.sha256()
    b  = bytearray(128*1024)
    mv = memoryview(b)
    with open(filename, 'rb', buffering=0) as f:
        while n := f.readinto(mv):
            h.update(mv[:n])
    return h.hexdigest()

def log_training(sm, pvalue, pr, pvalue2, rmse, train_n, test_n, args, xparams):
    if args.google_sheets_credentials_filename:
        from glogger import GSheetLogger
        glogger = GSheetLogger(args.google_sheets_credentials_filename, args.google_spreadsheet_id)
        params = dict()
        for k, v in args.__dict__.items():
            if k.endswith("_filename") and v:
                params[k] = sha256sum(v)[:8] + "|" + v
            else:
                params[k] = v

        with open('last_commit.log', mode='r') as f:
            params["last_commit"] = f.read()

        glogger.log("Log", {"spearman": sm, "p-value": pvalue, "pearson": pr, "p-value2": pvalue2, "rmse": rmse, "train_n": train_n, "test_n": test_n, **params, **os.environ, **xparams})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args()
    setup(args)

    num_round, early_stopping_round = 2000, 50
    
    default_params = {'booster': 'gbtree', 'objective': 'reg:squarederror', 'tree_method': 'hist' if is_mbp else 'gpu_hist',
                      'base_score': 2}
    params = {'alpha': 0.09838211642344241, 'colsample_bylevel': 0.7845488136664186, 'colsample_bynode': 0.3190626081646857,
              'colsample_bytree': 0.5898165366422254, 'eta': 0.051817108848631936, 'gamma': 0.4980984973194085,
              'lambda': 275.1911788945532, 'max_bin': 100, 'max_delta_step': 4.636511934158108, 'max_depth': 6,
              'min_child_weight': 2.9990138316074377, 'random_state': random.randint(-sys.maxsize-1, sys.maxsize)}
    all_params = {**default_params, **params}

    is_all_data, dtrain, dtest, ytest, test_genes = read_train_data(args)

    rmse, sm_correlation, sm_pvalue, pr_correlation, pr_pvalue, best_iteration, bst, ypred, ytest = \
        train(dtrain, dtest, ytest, all_params, num_round, early_stopping_round)

    logging.info(f"rmse: {rmse} spearman: {sm_correlation} pval: {sm_pvalue} pearson: {pr_correlation} pval: {pr_pvalue}")

    if(log_full_predictions):
        test_genes[args.log_note+'_ypred'] = ypred
        test_genes[args.log_note+'_ytest'] = ytest
        
        test_genes.to_csv("predictions/"+args.log_note.replace("/", "_").replace(",", "_") + ".csv")

    log_training(nanstr(sm_correlation), nanstr(sm_pvalue), nanstr(pr_correlation), nanstr(pr_pvalue), nanstr(rmse),
                 dtrain.num_row(), dtest.num_row(), args,
                 {'num_round': num_round, 'best_iteration': best_iteration, "biosample": args.biosample,
                  "Note": args.log_note, **all_params})
