#!/usr/bin/python

# Train a random forest model for phage prediction

import os
import datetime
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
import methods as ms


# Train VFM
def VFM_train(model_name, data_csv, threads, data_dir=None, loading_data=False, dims=24):
    if threads == None:
        threads = 1
    print('Start time for training: ' + str(datetime.datetime.now()))
    train_X = None
    train_Y = None
    if loading_data:
        (train_X, train_Y) = ms.load_data(data_csv, True, dims)
    else:
        for index, folder in enumerate(data_dir):
            (ids, temp_X, temp_Y) = ms.train_preprocess(folder, index, threads)
            if train_X is None:
                train_X = temp_X
                train_Y = temp_Y
            else:
                train_X = np.concatenate((train_X, temp_X), axis=0)
                train_Y = np.concatenate((train_Y, temp_Y), axis=0)
        # save data
        data = np.concatenate((train_X, train_Y), axis=1)
        np.savetxt(os.path.join("user_models", data_csv), data, delimiter=",")

    # Param:
    seed = 7
    num_tree = 50

    print("Start training...")
    model = RandomForestClassifier(n_estimators=num_tree, random_state=seed)
    # Deal with nan values
    if np.isnan(train_X).sum() > 0:
        train_X = ms.fill_nan(train_X)
    if train_Y.ndim == 2:
        train_Y = train_Y.flatten()
    model.fit(train_X, train_Y)
    print("Training finished.")
    print("Saving model...")
    joblib.dump(model, "user_models/" + model_name)
    print('End time for training: ' + str(datetime.datetime.now()) + "\n")