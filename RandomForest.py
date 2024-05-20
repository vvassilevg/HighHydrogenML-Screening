#!/usr/bin/env python
# coding: utf-8

# Libraries
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics, model_selection, preprocessing
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import cross_val_score
from sklearn.utils import shuffle
from sklearn.metrics import median_absolute_error
import pickle
import sys

#A random_state of 42 was used to train the models published in the original publication.



# Load datasets
data = sys.argv[1]
df = pd.read_excel(data, header=0)

X_Eads = df.iloc[:,:9]
# Visualization of feature vector space in dataset 
X_Eads.describe()


# Hyperparameter optimization for the model predicting the surface activity (Eads)
# Target
y_Eads = df[["Eads"]]
# Dataset splitting in train and test set
X_train, X_test, y_train, y_test = train_test_split(X_Eads, y_Eads, test_size=0.15, random_state=42)
#Grid
n_est = [1000]
max_feat = np.arange(0,10)
max_d = [100, 300, 600, 900]
min_samples_spl = [2, 3, 4]
#Grid search optimal parameters
param_grid = {'n_estimators': n_est, 'max_features': max_feat, 'max_depth': max_d, 'min_samples_split':min_samples_spl}
extr = ExtraTreesRegressor()
grid_search = GridSearchCV(estimator = extr, param_grid = param_grid, cv = 5, n_jobs = -1, verbose = 2)
grid_search.fit(X_train, y_train.values.ravel())
best_param = grid_search.best_params_

print(best_param)


# Models assessement


# Model predicting the adsorption energy
y_Eads = df[["Eads"]]

X_trainset, X_testset, y_trainset, y_testset = train_test_split(X_Eads, y_Eads, test_size=0.15, random_state=42) 
model = ExtraTreesRegressor(n_estimators=1000, random_state=42, max_features=9, min_samples_leaf=1, min_samples_split=3, max_depth=900) 

rs=30
X_trainset,y_trainset = shuffle(X_trainset,y_trainset, random_state=rs)

model.fit(X_trainset, y_trainset.values.ravel()) 
y_pred = model.predict(X_testset) 
y_predtrain = model.predict(X_trainset)

print(model.score(X_testset, y_testset)) 
print(mean_absolute_error(y_testset, y_pred))  
print(mean_squared_error(y_testset, y_pred, squared=False)) 
print(median_absolute_error(y_testset, y_pred))

#Save model
with open('My_RF_model.pickle', 'wb') as F:
    pickle.dump(model, F)

