import sys
import autosklearn
import autosklearn.classification
import sklearn.model_selection
import sklearn.datasets
import sklearn.metrics
from sklearn.inspection import permutation_importance
from autosklearn.experimental.askl2 import AutoSklearn2Classifier
import pandas as pd
from pprint import pprint
import numpy as np
import pickle

mcc_scorer = autosklearn.metrics.make_scorer(
    "mcc",
    sklearn.metrics.matthews_corrcoef,
    greater_is_better=True
)

def rf_auto_ml_random_drop(file1, random_size):
    df = pd.read_csv(file1, sep="\t", dtype={"name/position": str, "outcome": int, "*": int})

    array = df.to_numpy()
    max_len = array.shape[1] - 1

    mutation_columns = df.columns[1:-1]

    drop_percent = len(mutation_columns) * (10-random_size) / 10

    mutations_selected_random = np.random.choice(mutation_columns, int(drop_percent), replace=False)

    mutations_to_drop = mutations_selected_random.tolist()

    df = df.drop(mutations_to_drop, axis=1)

    array = df.to_numpy()
    max_len = array.shape[1] - 1

    X = array[:, 1:-1].astype(int)
    y = array[:, max_len].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=42, test_size=0.33)

    for e, classifier in enumerate(["random_forest"]):
        cls = autosklearn.classification.AutoSklearnClassifier(
            memory_limit=400000,
            time_left_for_this_task=100000,
            include={'classifier': [classifier]},
            resampling_strategy='holdout',
            resampling_strategy_arguments={'train_size': 0.67},
            delete_tmp_folder_after_terminate=False,
            ensemble_size=1,
            metric=mcc_scorer
        )
        cls.fit(X_train, y_train)
        y_hat = cls.predict(X_test)

        outfile = file1[:-4] + "_RF_AutoML_Result"

        filename = f"{outfile}_model.sav"
        pickle.dump(cls, open(filename, 'wb'))

        with open(outfile, "w") as ofile:
            ofile.write(str(cls.show_models()))
            ofile.write("\n")
            ofile.write("Accuracy score: " + str(sklearn.metrics.accuracy_score(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("Balanced Accuracy score: " + str(sklearn.metrics.balanced_accuracy_score(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("Brier score loss: " + str(sklearn.metrics.brier_score_loss(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("F1 score macro: " + str(sklearn.metrics.f1_score(y_test, y_hat, average='macro')))
            ofile.write("\n")
            ofile.write("F1 score micro: " + str(sklearn.metrics.f1_score(y_test, y_hat, average='micro')))
            ofile.write("\n")
            ofile.write("F1 score weighted: " + str(sklearn.metrics.f1_score(y_test, y_hat, average='weighted')))
            ofile.write("\n")
            ofile.write("F1 score binary: " + str(sklearn.metrics.f1_score(y_test, y_hat, average='binary')))
            ofile.write("\n")
            ofile.write("Precision score: " + str(sklearn.metrics.precision_score(y_test, y_hat, average='binary')))
            ofile.write("\n")
            ofile.write("Recall score: " + str(sklearn.metrics.recall_score(y_test, y_hat, average='binary')))
            ofile.write("\n")
            ofile.write("Confussion matrix: " + str(sklearn.metrics.confusion_matrix(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("ROC Curve: " + str(sklearn.metrics.roc_curve(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("ROC AUC Score: " + str(sklearn.metrics.roc_auc_score(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("Jaccard score: " + str(sklearn.metrics.jaccard_score(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("Hinge loss: " + str(sklearn.metrics.hinge_loss(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("Hamming loss: " + str(sklearn.metrics.hamming_loss(y_test, y_hat)))
            ofile.write("\n")
            ofile.write(
                "Fbeta score macro: " + str(sklearn.metrics.fbeta_score(y_test, y_hat, average='macro', beta=0.5)))
            ofile.write("\n")
            ofile.write(
                "Fbeta score micro: " + str(sklearn.metrics.fbeta_score(y_test, y_hat, average='micro', beta=0.5)))
            ofile.write("\n")
            ofile.write("Fbeta score weighted: " + str(
                sklearn.metrics.fbeta_score(y_test, y_hat, average='weighted', beta=0.5)))
            ofile.write("\n")
            ofile.write("Log loss: " + str(sklearn.metrics.log_loss(y_test, y_hat)))
            ofile.write("\n")
            ofile.write("Matthews correlation coefficient: " + str(sklearn.metrics.matthews_corrcoef(y_test, y_hat)))

