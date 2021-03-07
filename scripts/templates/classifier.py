#!user/bin/python3

import numpy as np
import pandas as pd
import xgboost as xgb
import sklearn.model_selection as ms
from sklearn import metrics
import shap
import pickle
import sys
import os
from datetime import datetime
import json
from json import JSONEncoder
from ast import literal_eval

with open('configurations.json', 'r') as f:
    config = json.load(f)

N_THREADS = config['N_THREADS']
os.environ['OMP_NUM_THREADS'] = str(N_THREADS)
seed = 42

class NumpyArrayEncoder(JSONEncoder):
    """ JSON encoded for numpy array
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)

def fileparser(path, dlist, sample=0, L=2):
    """ The fileparser to read the events from a csv
        argument:
            path: the path to the file
            dlist: the list of variables to be excluded
            frac: the fraction of sample that will be the test sample when sample is set to 0
            sample: the number of events that will be the train sample.
            L: Luminosity scaling
        returns:
            df_train: the training dataframe
            df_test: the testing dataframe
            weight: the weight (related to crosssection)
    """
    df = pd.read_csv(path)
    df.drop(columns=dlist, inplace=True)    
    n = len(df)
    weight = int(round(np.abs(df['weight'].sum()) * 3. * 1e6 * L))
    df_train = df.sample(n=sample, random_state=seed)
    df_test = df.drop(df_train.index)
    
    return df_train, df_test, weight

def runBDT(df, filename, depth=10, seed=seed):
    """ The BDT/RF runner
        argument:
            df: the dataframe with all the events
            filename: the name of the pickle file to store the model in
            depth: The depth of the trees
            seed: the seed for the random number generator
    """
    X = df.drop(columns=['class', 'weight'])
    y = df['class'].values
    nchannels = np.unique(y).shape[0]

    # Split for training and testing
    x_train, x_test, y_train, y_test = ms.train_test_split(X.values, y, test_size=0.2, random_state=seed)
    eval_set = [(x_train, y_train), (x_test, y_test)]
    
    # Fit the decision tree
    print('running classifier')
    now_r = datetime.now()
    classifier = xgb.XGBClassifier(max_depth=depth, learning_rate=0.01, objective='multi:softprob', num_class=nchannels,
                                         n_jobs=N_THREADS, subsample=0.5, colsample_bytree=1, n_estimators=5000, random_state=seed)
    classifier = classifier.fit(x_train, y_train, early_stopping_rounds=50, eval_set=eval_set,
                                eval_metric=["merror", "mlogloss"], verbose=False)
    print('classification time: {}'.format(datetime.now() - now_r))
    
    # Predictions
    print('making predictions')
    y_pred = classifier.predict(x_test)
    print('Accuracy Score: {:4.2f}% '.format(100*metrics.accuracy_score(y_test, y_pred)))
    print('dumping model')
    pickle.dump(classifier, open(filename, 'wb'))
    
    # Calculate the SHAP scores
    print('calculating SHAP')
    now_s = datetime.now()
    X_shap = pd.DataFrame(x_test, columns=df.drop(columns=['class', 'weight']).columns)
    explainer = shap.TreeExplainer(classifier)
    shap_values = explainer.shap_values(X_shap)
    print('Shapley time: {}'.format(datetime.now() - now_s))
    
    print('dumping Shapley values')
    numpyData = {"shap_values": shap_values}
    with open('shapley_files/shapley_values.json', 'w') as f:
        encodedNumpyData = json.dump(numpyData, f, cls=NumpyArrayEncoder)
    
    # Dump the data used to compute the Shapley values
    X_shap.to_json('shapley_files/shapley_X.json', orient='records')
    
    

def main():
    
    now = datetime.now()
    
    # Set up the file system
    if not os.path.exists('test_files'):
        os.makedirs('test_files')
    if not os.path.exists('shapley_files'):
        os.makedirs('shapley_files')
    
    # filename for getting the training data
    prefix = config['prefix']
    try:
        file_bbxaa = prefix+'bbxaa.csv'
        file_yb2 = prefix+'yb2.csv'
        file_ybyt = prefix+'ybyt.csv'
        file_yt2 = prefix+'yt2.csv'
        file_zh = prefix+'zh.csv'
        file_tth = prefix+'ttH_full.csv'
        file_bkg = prefix+config['bkg_filename']
        file_sig = prefix+config['signal_filename']
    except:
        print('training data files not found')
        exit(1)
    
    # list of kinematic variables to be removed
    dlist = literal_eval(config['dlist'])
    
    # Signal loading
    df_sig, df_sig_test, weight_sig = fileparser(file_sig, dlist, sample=40000)
    df_sig['class'] = 4
    df_sig_test['class'] = 4

    # Primary background loading
    df_bkg, df_bkg_test, weight_bkg = fileparser(file_bkg, dlist, sample=40000)
    df_bkg['class'] = 3
    df_bkg_test['class'] = 3
    
    # Secondary background loading
    df_tth, df_tth_test, weight_tth = fileparser(file_tth, dlist, sample=20000)
    df_yb2, df_yb2_test, weight_yb2 = fileparser(file_yb2, dlist, sample=3890*6)
    df_ybyt, df_ybyt_test, weight_ybyt = fileparser(file_ybyt, dlist, sample=500*6)
    df_yt2, df_yt2_test, weight_yt2 = fileparser(file_yt2, dlist, sample=7360*6)
    df_zh, df_zh_test, weight_zh = fileparser(file_zh, dlist, sample=4960*6)
    
    df_bbh = pd.concat([df_yb2, df_ybyt, df_yt2, df_zh])
    df_bbh_test = pd.concat([df_yb2_test, df_ybyt_test, df_yt2_test, df_zh_test])
    df_bbh_test = df_bbh_test.sample(frac=0.5).reset_index(drop=True)
    df_bbh['class'] = 1
    df_bbh_test['class'] = 1
    weight_bbh = int(weight_yb2*1.5 - weight_ybyt*1.9 + weight_yt2*2.5 + weight_zh*1.3)
    
    df_bbxaa, df_bbxaa_test, weight_bbxaa = fileparser(file_bbxaa, dlist, sample=100000)
    
    weight_dict = {'weight_bbxaa': weight_bbxaa, 
                   'weight_tth': weight_tth, 
                   'weight_bbh': weight_bbh, 
                   'weight_bkg': weight_bkg, 
                   'weight_sig': weight_sig}
    
    with open('test_files/weights.json', 'w') as f:
        json.dump(weight_dict, f)
    
    # dump the test sets
    df_sig_test.to_json('test_files/sig_test.json', orient='records')
    df_bkg_test.to_json('test_files/bkg_test.json', orient='records')
    df_tth_test.to_json('test_files/tth_test.json', orient='records')
    df_bbh_test.to_json('test_files/bbh_test.json', orient='records')
    df_bbxaa_test.to_json('test_files/bbxaa_test.json', orient='records')
    
    # Prepare for training
    channels = [df_sig, df_bkg, df_bbh, df_tth, df_bbxaa]
    nchannels = len(channels)
    df_train = pd.concat(channels, ignore_index=True)
    df_train = df_train.sample(frac=1).reset_index(drop=True)
    
    # Train the BDT and compute Shapley values
    runBDT(df_train, config['classifier_filename'])
    
    print('run time: {}'.format(datetime.now() - now))


if __name__ == "__main__":
    # execute only if run as a script
    main()