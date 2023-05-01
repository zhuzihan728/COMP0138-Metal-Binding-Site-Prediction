from sklearn.metrics import precision_score, recall_score, roc_auc_score, average_precision_score, matthews_corrcoef, confusion_matrix, precision_recall_curve, PrecisionRecallDisplay
from sklearn.model_selection import KFold
import json
import numpy as np
import torch

np.random.seed(42)


def retrieve_json(path):
    with open(path, 'r') as fp:
        data = json.load(fp)
    return data

class Validation_tool():
    def __init__(self, class_encode_json=None):
        if class_encode_json is None:
            self.class_encode = retrieve_json('label_encode\class_encode.json')
        else:
            self.class_encode = retrieve_json(class_encode_json)
        self.nmetals = len(self.class_encode)


    def dataset_class_summary(self, label_path):
        labels = np.load(label_path)
        dc = {i: 0 for i in range(self.nmetals+1)}  # 18 binding + 1 non-binding
        for key, value in labels.items():
            dc[self.nmetals] += np.where(value.sum(axis=0) == 0)[0].shape[0]  # non-binding
            for i in range(self.nmetals):
                dc[i] += value[i].sum()
        for key, value in self.class_encode.items():
            print(f"{value[1]:25}|{dc[int(key)]}")
        print(f"{'neg':25}|{dc[self.nmetals]}")

    def five_fold_val_split(self, label_path, children = [0, 1, 2, 3, 4, 7, 9, 10,
                        11, 12, 14, 15, 16], parent = [5, 6, 8, 13, 17]):
        labels = np.load(label_path)
        children_metals = children  # rows of children metals
        parent_metals = parent  # rows of parent metals
        assert len(children_metals) + \
            len(parent_metals) == self.nmetals, "wrong number of metals"

        dc = {}  # store the protein id for each metal
        for i in children_metals:
            dc[i] = []
        dc[i+1] = []  # 17 parent metal
        dc[i+2] = []  # 18 non-metal

        for key, value in labels.items():
            is_child = False
            is_parent = False
            if value.sum() == 0:
                dc[self.nmetals].append(key)
                continue
            for i in children_metals:
                if value[i].sum() > 0:
                    is_child = True
                    dc[i].append(key)
                    break
            if not is_child:
                for i in parent_metals:
                    if value[i].sum() > 0:
                        is_parent = True
                        dc[self.nmetals-1].append(key)
                        break
            if not is_child and not is_parent:
                dc[self.nmetals].append(key)

        # just sanity check
        total_len1 = 0
        for i in dc.keys():
            total_len1 += len(dc[i])
        assert total_len1 == len(labels), "wrong label division"

        duplicate_count = 0
        split = [[], [], [], [], []]
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        for key, value in dc.items():
            if len(value) == 0:
                continue
            try:
                for i, (train_index, test_index) in enumerate(kf.split(value)):
                    split[i] += [value[i] for i in test_index]
            except:
                print(
                    f"not enough data to divide {self.class_encode[str(key)]}, add to all split!")
                for i in range(5):
                    split[i] += value
                duplicate_count += len(value)

        five_fold_splits = {}
        for i in range(5):
            temp = {'train': [], 'test': []}
            for j in range(5):
                if j == i:
                    temp['test'] += split[j]
                else:
                    temp['train'] += split[j]
            temp['train'] = list(set(temp['train']))
            temp['test'] = list(set(temp['test']))
            five_fold_splits[f'fold{i}'] = temp

        for i in range(5):
            assert len(five_fold_splits[f'fold{i}']['train']) + len(five_fold_splits[f'fold{i}']
                                                                    ['test']) == len(labels) + duplicate_count, f"wrong number of labels in split {i}"

        return five_fold_splits


    def calculate_2metrics(self, y_true, y_pred, thres, all=False):
        y_bin = (y_pred > thres).astype(np.int8)
        if all:
            MCC = 0
        else:
            MCC = matthews_corrcoef(y_true, y_bin)

        try:
            AUPR = average_precision_score(y_true, y_pred)
        except:
            AUPR = np.nan
        return MCC, AUPR


    def sum_2metrics(self, gts, preds, thres):
        print()
        summary = {}
        all_MCC = []
        all_AUPR = []
        if type(thres) == float:
            thres = [thres for i in range(preds.shape[0])]
        for i in range(preds.shape[0]):
            MCC, AUPR = self.calculate_2metrics(gts[i], preds[i], thres[i])
            all_MCC.append(MCC)
            all_AUPR.append(AUPR)
            metal_name = self.class_encode[str(i)][1]
            summary[f'{metal_name} MCC'] = MCC
            summary[f'{metal_name} AUPR'] = AUPR

        summary['mean MCC'] = np.mean(all_MCC)
        summary['mean AUPR'] = np.mean(all_AUPR)
        print(
            f"{'mean':25}|MCC(-1-1): {summary['mean MCC']:6.3f}|AUPR(0-1): {summary['mean AUPR']:.3f}")
        return summary


    def calculate_full_metrics(self, y_true, y_pred, thres):
        y_bin = (y_pred > thres).astype(np.int8)

        tn, fp, fn, tp = confusion_matrix(y_true, y_bin, labels=[0, 1]).ravel()
        precision = precision_score(y_true, y_bin, labels=[0, 1])
        recall = recall_score(y_true, y_bin, labels=[0, 1])
        specificity = tn/(tn + fp)
        MCC = matthews_corrcoef(y_true, y_bin)
        try:
            AUROC = roc_auc_score(y_true, y_pred)
        except:
            AUROC = np.nan

        try:
            AUPR = average_precision_score(y_true, y_pred)
        except:
            AUPR = np.nan
        return precision, recall, specificity, MCC, AUROC, AUPR


    def sum_full_metrics(self, gts, preds, thres):
        print()
        summary = {}
        all_pre = []
        all_rec = []
        all_spe = []
        all_MCC = []
        all_AUROC = []
        all_AUPR = []
        if type(thres) == float:
            thres = [thres for i in range(preds.shape[0])]
        for i in range(preds.shape[0]):
            precision, recall, specificity, MCC, AUROC, AUPR = self.calculate_full_metrics(
                gts[i], preds[i], thres[i])
            all_pre.append(precision)
            all_rec.append(recall)
            all_spe.append(specificity)
            all_MCC.append(MCC)
            all_AUROC.append(AUROC)
            all_AUPR.append(AUPR)
            metal_name = self.class_encode[str(i)][1]
            summary[f'{metal_name} precision'] = precision
            summary[f'{metal_name} recall'] = recall
            summary[f'{metal_name} specificity'] = specificity
            summary[f'{metal_name} MCC'] = MCC
            summary[f'{metal_name} AUROC'] = AUROC
            summary[f'{metal_name} AUPR'] = AUPR
            print(f"{metal_name:25}|precision: {precision:6.3f}|recall: {recall:6.3f}|specificity: {specificity:6.3f}|MCC: {MCC:6.3f}|AUROC: {AUROC:.3f}|AUPR: {AUPR:.3f}")

        summary['mean precision'] = np.mean(all_pre)
        summary['mean recall'] = np.mean(all_rec)
        summary['mean specificity'] = np.mean(all_spe)
        summary['mean MCC'] = np.mean(all_MCC)
        summary['mean AUROC'] = np.mean(all_AUROC)
        summary['mean AUPR'] = np.mean(all_AUPR)
        print(
            f"{'mean':25}|pre(0-1): {summary['mean precision']:6.3f}|rec(0-1): {summary['mean recall']:6.3f}|spe(0-1): {summary['mean specificity']:6.6f}|MCC(-1-1): {summary['mean MCC']:6.3f}|AUROC(0-1): {summary['mean AUROC']:.3f}|AUPR(0-1): {summary['mean AUPR']:.3f}")
        return summary

    def plot_pr_curve(self, gts, preds):
        for i in range(preds.shape[0]):
            metal_name = self.class_encode[str(i)][1]
            display = PrecisionRecallDisplay.from_predictions(gts[i], preds[i])
            _ = display.ax_.set_title(f"{metal_name} Precision-Recall curve")



    def mcc_on_thres(self, y_true, y_pred, thres):
        dc = {'thres': [], 'MCC': []}
        
        for t in thres:
            y_bin = (y_pred > t).astype(np.int8)
            MCC = matthews_corrcoef(y_true, y_bin)
            dc['thres'].append(t)
            dc['MCC'].append(MCC)
            
        dc1 = {'thres': dc['thres'], 'MCC': np.array(dc['MCC'])}
        return dc1 




    def mcc_sum(self, gts, preds, thres):
        print()
        summary = {}


        for i in range(preds.shape[0]):
            MCC_dc = self.mcc_on_thres(gts[i], preds[i], thres)
            metal_name = self.class_encode[str(i)][1]
            summary[f'{metal_name}'] = MCC_dc
            max_MCC = np.max(MCC_dc['MCC'])
            max_thres = MCC_dc['thres'][np.argmax(MCC_dc['MCC'])]
            print(f"{metal_name:25}|max MCC: {max_MCC:6.3f}|thres: {max_thres:6.3f}")
        return summary

