from sklearn.inspection import permutation_importance
import pandas as pd
import sklearn
from sklearn.metrics import classification_report
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import pickle


def svm_ml(file1, abiotic):

    print(abiotic)
    print("\n")

    df = pd.read_csv(file1, sep="\t", dtype={"name/position": str, "outcome": int, "*": int})

    array = df.to_numpy()
    max_len = array.shape[1] - 1

    X = array[:, 1:-1].astype(int)
    y = array[:, max_len].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=42, test_size=0.33)

    loaded_model = pickle.load(open(file1[:-4]+"_RF_AutoML_Result2_model.sav", 'rb'))

    #print("Model loaded")

    #result = loaded_model.predict(X_test)

    #print("Perm importance started")
    r = permutation_importance(loaded_model, X_test, y_test, n_repeats=30, random_state=42)

    #print("Perm importance ended")
    #print("\n")
    for i in r.importances_mean.argsort()[::-1]:
        if r.importances_mean[i] - 2 * r.importances_std[i] > 0:
            print(f"{df.columns[i]:<8}"
                  f"{r.importances_mean[i]:.3f}"
                  f" +/- {r.importances_std[i]:.3f}")


    #print(sklearn.metrics.matthews_corrcoef(y_test, result))

    #result = loaded_model.score(X_test, y_test)

    #print(result)

    print("\n")

antibiotics = ["amikacin", "capreomycin", "ethionamide", "kanamycin", "ofloxacin", "streptomycin"]
#antibiotics = ["ofloxacin"]
for abiotic in antibiotics:
    svm_ml(
        "/scratch/SCRATCH_SAS/alper/Mycobacterium/automl30percent/combined_binary_mutations_non_snp_corrected_0.2_column_corrected_after_pyhlogeny_ordered_30_%s.tsv" % abiotic, abiotic)
