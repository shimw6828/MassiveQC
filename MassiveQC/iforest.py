from typing import Tuple
import numpy as np
import pandas as pd
import shap
from sklearn.ensemble import IsolationForest
from sklearn.model_selection import train_test_split


class SraIsolationForest:
    def __init__(
        self,
        features: pd.DataFrame,
        random_state=None,
        iso_kwargs=None,
        explain_kwargs=None,
        shap_values_kwargs=None,
    ):
        """Instantiate, fit, and explain an isolation forest.
        Parameters
        ----------
        features : pd.DataFrame
            Set of features to fit. Automatically does X_train/X_test splitting.
        random_state: np.random.RandomState, optional
            Numpy random state.
        iso_kwargs : optional
            Passed directly to `sklearn.ensemble.IsolationForest, by default None
        explain_kwargs : optional
            Passed directly to `shap.TreeExplainer`, by default None
        shap_values_kwargs : optional
            Passed directly to `shap.TreeExplainer.shap_values`, by default None
        Attributes
        ----------
        features : pd.DataFrame
            Stores the features passed at init.
        X_train : pd.DataFrame
            The subset used to X_train the isolation forest model.
        X_test : pd.DataFrame
            The subset used to X_test
        model_ : sklearn.ensemble.IsolationForest
            Fitted isolation forest model.
        explainer_ : shap.TreeExplainer
            Shap model to explain isolation forest results
        expected_value_ : float
            The mean value from `shap.TreeExpliner.expected_value`.
        index : pd.Index
            features.index
        columns : pd.Columns
            features.columns
        Note that many method names are available as suffixed version
        {_test, _train, _all} to return the results for
        {X_test, X_train, or features}
        Methods
        -------
        predict
            Return an array from running isolation forest outlier detection.
        inliers
            Return a dataframe with inliers.
        outliers
            Retrun a dataframe with outliers.
        isinlier
            Return a boolean array where inlier.
        isoutlier
            Return a boolean array where outlier.
        prop_outliers
            Return the proportion of samples that are outliers.
        shap_values
            Return an array with shap_values
        """
        # Prep features
        self.random_state = random_state or np.random.RandomState(42)
        self.features = features
        self.X_train, self.X_test = train_test_split(self.features, random_state=self.random_state)

        # Run Isolation Forest
        iso_kwargs = iso_kwargs or {}
        self.model_ = IsolationForest(**iso_kwargs).fit(self.X_train)

        # Run Tree Explainer
        explain_kwargs = explain_kwargs or {}
        self.explainer_ = shap.TreeExplainer(self.model_, data=self.inliers_test, **explain_kwargs)
        self.expected_value_ = self.explainer_.expected_value

        shap_values_kwargs = shap_values_kwargs or {}
        self.shap_values_ = self.explainer_.shap_values(self.X_test, **shap_values_kwargs)

    @property
    def index(self) -> pd.Index:
        """Return sample index"""
        return self.features.index

    @property
    def columns(self) -> pd.Index:
        """Return feature names."""
        return self.features.columns

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        """Run Isolation Forest outlier detection.
        Parameters
        ----------
        X : pd.DataFrame
            Data frame of features to run prediction
        Returns
        -------
        np.ndarray
            An array where -1 are outliers and 1 are inliers.
        """
        return self.model_.predict(X)

    def isinlier(self, X: pd.DataFrame) -> np.ndarray:
        """Create inlier boolean mask.
        Parameters
        ----------
        X : pd.DataFrame
            Data frame of features to run prediction
        Returns
        -------
        np.ndarray
            An array where True are inliers and False are outliers.
        """
        return np.where(self.model_.predict(X) == 1, True, False)

    @property
    def isinlier_test(self) -> np.ndarray:
        """Inlier boolean mask for X_test data"""
        return self.isinlier(self.X_test)

    @property
    def isinlier_train(self) -> np.ndarray:
        """Inlier boolean mask for X_train data"""
        return self.isinlier(self.X_train)

    @property
    def isinlier_all(self) -> np.ndarray:
        """Inlier boolean mask for all data"""
        return self.isinlier(self.features)

    def isoutlier(self, X: pd.DataFrame) -> np.ndarray:
        """Create outlier boolean mask.
        Parameters
        ----------
        X : pd.DataFrame
            Data frame of features to run prediction
        Returns
        -------
        np.ndarray
            An array where True are outliers and False are inliers.
        """
        return np.where(self.model_.predict(X) == -1, True, False)

    @property
    def isoutlier_test(self) -> np.ndarray:
        """Inlier boolean mask for X_test data"""
        return self.isoutlier(self.X_test)

    @property
    def isoutlier_train(self) -> np.ndarray:
        """Inlier boolean mask for X_train data"""
        return self.isoutlier(self.X_train)

    @property
    def isoutlier_all(self) -> np.ndarray:
        """Inlier boolean mask for all data"""
        return self.isoutlier(self.features)

    def inliers(self, X: pd.DataFrame) -> pd.DataFrame:
        """Return a data frame of inliers"""
        return X[self.isinlier(X)]

    @property
    def inliers_test(self) -> pd.DataFrame:
        """Return a data frame of inliers"""
        return self.inliers(self.X_test)

    @property
    def inliers_train(self) -> pd.DataFrame:
        """Return a data frame of inliers"""
        return self.inliers(self.X_train)

    @property
    def inliers_all(self) -> pd.DataFrame:
        """Return a data frame of inliers"""
        return self.inliers(self.features)

    def outliers(self, X: pd.DataFrame) -> pd.DataFrame:
        """Return a data frame of outliers"""
        return X[self.isoutlier(X)]

    @property
    def outliers_test(self) -> pd.DataFrame:
        """Return a data frame of inliers"""
        return self.outliers(self.X_test)

    @property
    def outliers_train(self) -> pd.DataFrame:
        """Return a data frame of inliers"""
        return self.outliers(self.X_train)

    @property
    def outliers_all(self) -> pd.DataFrame:
        """Return a data frame of inliers"""
        return self.outliers(self.features)

    def prop_outliers(self, X) -> float:
        """Return the proportion of outliers in X."""
        return np.mean(self.isoutlier(X))

    @property
    def prop_outliers_test(self):
        """Return the proportion of outliers in X_test data."""
        return self.prop_outliers(self.X_test)

    @property
    def prop_outliers_train(self):
        """Return the proportion of outliers in X_train data."""
        return self.prop_outliers(self.X_train)

    @property
    def prop_outliers_all(self):
        """Return the proportion of outliers in all data."""
        return self.prop_outliers(self.features)

    @property
    def shap_values(self):
        """Retruns calculated shap values for X_test data."""
        return self.shap_values_

    @property
    def shap_values_inliers(self):
        """Return shap values for inliers for X_test data"""
        return self.shap_values_[self.isinlier_test]

    @property
    def shap_values_outliers(self):
        """Return shap values for outliers for X_test data"""
        return self.shap_values_[self.isoutlier_test]

    def mean_shap_values(self, shap_values) -> Tuple[np.ndarray, pd.Index]:
        """Calculate the mean(abs(shap_values)).
        Returns
        -------
        Tuple[np.ndarray, pd.Index]
            Mean absolute shapely values and feature names ordered by
            importance.
        """
        mean_shap = np.mean(np.abs(shap_values), axis=0)
        idx_order = np.argsort(mean_shap)[::-1]
        return mean_shap[idx_order], self.columns[idx_order]

    @property
    def mean_shap_values_inliers(self) -> Tuple[np.ndarray, pd.Index]:
        """Calculate the mean(abs(shap_values)).
        Returns
        -------
        Tuple[np.ndarray, pd.Index]
            Mean absolute shapely values and feature names ordered by
            importance.
        """
        return self.mean_shap_values(self.shap_values[self.isinlier_test])

    @property
    def mean_shap_values_outliers(self) -> Tuple[np.ndarray, pd.Index]:
        """Calculate the mean(abs(shap_values)).
        Returns
        -------
        Tuple[np.ndarray, pd.Index]
            Mean absolute shapely values and feature names ordered by
            importance.
        """
        return self.mean_shap_values(self.shap_values[self.isoutlier_test])