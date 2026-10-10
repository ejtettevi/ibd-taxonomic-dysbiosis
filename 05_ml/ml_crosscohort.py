import pandas as pd
import sys
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import LabelEncoder

def main():
    if len(sys.argv) < 4:
        print("Usage: python ml_crosscohort.py <train_cohort> <test_cohort> <level>")
        sys.exit(1)
        
    train_cohort = sys.argv[1]
    test_cohort = sys.argv[2]
    level = sys.argv[3] # "species" or "genus"
    
    try:
        # Load data from data/ directory
        taxon_path = f"data/taxon/{train_cohort}/combined_counts.csv" if level == "species" else f"data/taxon/{train_cohort}/genus_counts.csv"
        meta_path = f"data/metadata/{train_cohort}.csv"
        
        if not os.path.exists(taxon_path) or not os.path.exists(meta_path):
            raise FileNotFoundError(f"Missing data for {train_cohort}")
            
        train_tax = pd.read_csv(taxon_path, index_col=0)
        train_meta = pd.read_csv(meta_path, index_col='sample_name')
        common_train = train_tax.columns.intersection(train_meta.index)
        X_train = train_tax[common_train].T
        y_train = train_meta.loc[common_train, 'status'].map({'control': 0, 'PD': 1})
        
        test_taxon_path = f"data/taxon/{test_cohort}/combined_counts.csv" if level == "species" else f"data/taxon/{test_cohort}/genus_counts.csv"
        test_meta_path = f"data/metadata/{test_cohort}.csv"
        
        test_tax = pd.read_csv(test_taxon_path, index_col=0)
        test_meta = pd.read_csv(test_meta_path, index_col='sample_name')
        common_test = test_tax.columns.intersection(test_meta.index)
        X_test = test_tax[common_test].T
        y_test = test_meta.loc[common_test, 'status'].map({'control': 0, 'PD': 1})
        
        common_features = X_train.columns.intersection(X_test.columns)
        if len(common_features) < 5:
            print(f"ERROR: Too few overlapping features ({len(common_features)}) between {train_cohort} and {test_cohort}")
            return

        # Real shotgun profiles are sparse: with no prevalence filter the shared
        # feature set is dominated by near-zero, near-constant columns that carry
        # no signal and inflate variance. Retain features present at >=20% of
        # TRAINING samples, measured on the training cohort only so the test
        # cohort plays no part in feature selection.
        prev = (X_train[common_features] > 0).mean(axis=0)
        common_features = common_features[prev >= 0.20]
        if len(common_features) < 5:
            print(f"ERROR: Too few features at >=20% prevalence ({len(common_features)}) between {train_cohort} and {test_cohort}")
            return

        X_train = X_train[common_features]
        X_test = X_test[common_features]

        # class_weight='balanced': the case/control ratio differs sharply between
        # cohorts (e.g. 103 vs 27 in HMP2), so unweighted fitting would be biased
        # toward predicting the majority class on the test cohort.
        clf = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
        clf.fit(X_train, y_train)
        probs = clf.predict_proba(X_test)[:, 1]

        auc = roc_auc_score(y_test, probs)

        out_dir = 'out/ml_results'
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, 'auc_results.txt'), 'a') as f:
            f.write(f"Train: {train_cohort}, Test: {test_cohort}, Level: {level}, AUC: {auc:.3f}\n")

        # Persist per-sample predictions so ROC curves and CIs can be drawn from
        # observed data rather than reconstructed from the AUC value alone.
        pred = pd.DataFrame({
            'train_cohort': train_cohort,
            'test_cohort': test_cohort,
            'level': level,
            'sample_name': list(y_test.index),
            'y_true': list(y_test.values),
            'prob_ibd': probs,
        })
        pred.to_csv(os.path.join(out_dir, f'predictions_{train_cohort}__{test_cohort}__{level}.csv'), index=False)
        print(f"Train: {train_cohort}, Test: {test_cohort}, Level: {level}, AUC: {auc:.3f} "
              f"(n_test={len(y_test)}, features={len(common_features)})")
            
    except Exception as e:
        print(f"ERROR: {e}")

if __name__ == "__main__":
    main()
