# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 10:26:38 2019
Changed on Mon Oct 19 14:53:58 2020
@author: Li,Zhao
Usage:python3 AS_Finder_Denovo.py test.txt data.txt
Note:First, mcl should be loaded. Then, the test.txt is the result of test data blasting, which inculde 12 lines "gene    label   pident  length  mismatch gapopen qstart  qend    sstart  send    evalue  bitscore". "gene" line combines isoform1 and isoform2 with *, and "label" line is consist of 0 or 1, while 0 presents do not come from the same gene and 1 presents yes. Other lines produced by "blastn -query test.fa -db test -outfmt 6 -out test.txt -evalue 1e-10". At last, the data.txt comes from blast like test.txt, which does not have the "label" line.
"""
import warnings
warnings.filterwarnings('ignore') 
import pandas as pd
import numpy as np
from sklearn import ensemble, tree, neighbors, linear_model, svm, naive_bayes
from sklearn import model_selection
from xgboost import XGBClassifier
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt
from scipy import interp
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import GridSearchCV
import sys
from sklearn import preprocessing
import os
import subprocess
from sklearn.model_selection import train_test_split

def checkNanData(file):
    if file.isnull().any() == True:
        raise ValueError('Null values exist in this file')

def processingData(file):#Z-Score
    scaler = preprocessing.StandardScaler()
    file = scaler.fit_transform(file)
    
def tune_roc_curve(classifier):
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    probas_ = classifier.fit(train_X , train_y ).predict_proba(test_X)
    fpr, tpr, thresholds = roc_curve(test_y, probas_[:, 1])
    tprs.append(interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0
    roc_auc = auc(fpr, tpr)
    aucs.append(roc_auc)
    plt.plot(fpr,
             tpr,
             lw=1,
                 alpha=0.3,
                 label='AUC = %0.2f' % roc_auc)

    plt.plot([0, 1], [0, 1],
             linestyle='--',
             lw=2,
             color='r',
             label='Chance',
             alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr,
             mean_tpr,
             color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2,
             alpha=.8)
    
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Tune Hyperparameter {} ROC'.format(
        classifier.__class__.__name__))
    plt.legend(loc="lower right")
    plt.tight_layout()    
    plt.savefig(os.path.join(
        path, '{}.pdf'.format(classifier.__class__.__name__)),
                format='pdf')
    plt.close()
    
def cleandata():#choose data for mcl
    f1=pd.read_csv('prediction.csv',usecols=['gene','ft-probability'])
    f1['gene'],f1['gene1']=f1['gene'].str.split('*', 1).str
    f2=f1.loc[f1['ft-probability']>=0.6]
    f2=f2[['gene','gene1','ft-probability']]
    f2.to_csv('tmp-mcl.txt',sep='\t',header=None,index=None)
    
def mcl1():#cluster by mcl
    ret=subprocess.run("module load mcl/14-137 \n mcl tmp-mcl.txt --abc -I 2.0 -o pseudo-gene.txt ",shell=True)
    if ret.returncode == 0:
        print("success in mcl")
    else:
        print("error in mcl") 
    
def pred_ft_gene(clf):
    pred = clf.predict(testData.values)
    pred_pro = clf.predict_proba(testData.values)
    predictions = pd.DataFrame(
        {
            'label': pred,
            'ft-probability': pred_pro[:, 1],
            'no-ft-probility': pred_pro[:, 0]
        },
        index=testData.index)
    predictions.sort_values(by=['ft-probability'],
                            ascending=False,
                            inplace=True)
    predictions.to_csv(
        os.path.join(path,
                     'prediction.csv'.format(clf.__class__.__name__)))

if __name__ == '__main__':    
    path = os.getcwd()
    DataSet = pd.read_csv(sys.argv[1], sep='\t', index_col='gene')
    testData = pd.read_csv(sys.argv[2], sep='\t', index_col='gene')
    X = DataSet.drop(['label'], axis=1).values
    y = DataSet.label.values
    processingData(X)
    processingData(testData)
    train_X , test_X, train_y ,test_y = train_test_split(
            X, y, test_size=0.2,random_state=0)
    
    scoring = {
        'accuracy': make_scorer(accuracy_score),
        'precision': make_scorer(precision_score),
        'recall': make_scorer(recall_score),
        'f1_score': make_scorer(f1_score),
        'auc': make_scorer(roc_auc_score)
    }
    
    xgb_param = {
        'learning_rate': [.06],  #default: .3
        'max_depth': [11],  #default 2
        'n_estimators': [300],
        'seed': [0]
    }
    xgb = XGBClassifier()
    xgb_gs = GridSearchCV(estimator=xgb,
                          param_grid=xgb_param,
                          cv=None,
                          scoring='roc_auc',
                          refit=True)
    xgb_gs = xgb_gs.fit(X, y)
    
    Bxgb = xgb_gs.best_estimator_
    Axgb = xgb_gs.best_params_
    tune_roc_curve(Bxgb)
    pred_ft_gene(Bxgb)
    cleandata()
    mcl1()
    
