import sys, argparse

import numpy as np
import pandas as pd
from pathlib import Path

from .iforest import SraIsolationForest
#from .plot import outlier_umap, plot_importance

RANDOM_STATE = np.random.RandomState(42)


def detection(features_file):
    """Isolation forest training using features extracted from RNA seq"""
    rnaseq_features = pd.read_parquet(features_file)
    iso = SraIsolationForest(
        rnaseq_features, random_state=RANDOM_STATE, iso_kwargs=dict(n_estimators=100)
    )

    # Check that proportion of outliers in Train and Test is similar.
    assert np.isclose(iso.prop_outliers_test, iso.prop_outliers_train, atol=0.01)

    # Save a list of good rnaseq samples
    inliers = iso.inliers(rnaseq_features).index.tolist()
    outliers = iso.outliers(rnaseq_features).index.tolist()
    # plot_importance(iso, features_file)
    # outlier_umap(features_file, inliers)
    rnaseq_features['labels'] = rnaseq_features.index.to_series().apply(
        lambda x: "Inliers" if x in inliers else "Outlier"
    )
    result_path = Path(features_file).parent.parent / "result.csv"
    rnaseq_features.to_csv(result_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--features_file', required=True, type=str,
                        help="Path to feature dataframe file.")
    args = parser.parse_args()
    detection(args.features_file)

