import pandas as pd
import sklearn
from sklearn.metrics import classification_report
from sklearn.svm import SVC
import numpy as np


def add_snippy_info_to_fia(fia_file):
    outfile = fia_file + "_MUT_INFO"

    with open(fia_file, "r") as infile:
        fia_lines = infile.readlines()

    with open("./20032023_sonia_results_MUT_INFO.tsv") as snippy_info:
        snippy_info_lines = snippy_info.readlines()

    with open(outfile, "w") as ofile:
        ofile.write(snippy_info_lines[0])
        for line in fia_lines:
            splitted = line.split("\t")
            mutation = splitted[0].strip()
            for line2 in snippy_info_lines:
                splitted2 = line2.split("\t")
                mutation2 = splitted2[0].strip()
                if mutation == mutation2:
                    ofile.write(line2)


def svm_ml(file1, abiotic):
    df = pd.read_csv(file1, sep="\t", dtype={"name/position": str, "outcome": int, "*": int})

    print(df.shape)

    array = df.to_numpy()
    max_len = array.shape[1] - 1

    mutation_columns = df.columns[1:-1]

    seventy_percent = len(mutation_columns) * 7 / 10

    mutations_selected_random = np.random.choice(mutation_columns, int(seventy_percent), replace=False)

    mutations_to_drop = mutations_selected_random.tolist()

    df = df.drop(mutations_to_drop, axis=1)

    array = df.to_numpy()
    max_len = array.shape[1] - 1

    X = array[:, 1:-1].astype(int)
    y = array[:, max_len].astype(int)

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=42, test_size=0.33)

    best_model_mcc = -1.0
    bm_c = 0

    print(df.shape)

    for c_val in range(1,2):
        svm_cls = SVC(class_weight={0: sum(y_train), 1: len(y_train) - sum(y_train)}, kernel="linear", C=c_val)
        svm_cls.fit(X_train, y_train)

        y_hat = svm_cls.predict(X_test)

        cur_mcc_val = sklearn.metrics.matthews_corrcoef(y_test, y_hat)
        #print(c_val)
        if cur_mcc_val > best_model_mcc:
            best_model_mcc = cur_mcc_val
            best_model = svm_cls
            bm_c = c_val

    outfile = "./randomsvm/" + abiotic + "_SVM"

    idx = (-best_model.coef_[0]).argsort()[:20]

    y_hat = best_model.predict(X_test)

    with open(outfile + "_FIA_TOP20", "w") as fiafile:
        for i in idx:
            fiafile.write(df.columns[int(repr(i))] + "\n")

    add_snippy_info_to_fia(outfile + "_FIA_TOP20")

    with open(outfile + "_results", "w") as ofile:
        ofile.write("C: " + str(bm_c))
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
        ofile.write("F1 score: " + str(sklearn.metrics.f1_score(y_test, y_hat)))
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


if __name__ == "__main__":

    #antibiotics = ["amikacin", "capreomycin", "ethionamide", "kanamycin", "ofloxacin", "streptomycin"]
    antibiotics = ["kanamycin"]
    for abiotic in antibiotics:
        svm_ml(
            "/scratch/SCRATCH_SAS/alper/Mycobacterium/non_dropped/combined_binary_mutations_non_snp_corrected_0.2_%s_name_corrected_no_duplicates_dropped_zero_cols.tsv" % abiotic, abiotic)

