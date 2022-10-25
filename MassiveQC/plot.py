import os
import numpy as np
import pandas as pd
import umap
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

RANDOM_STATE = np.random.RandomState(42)

CATEGORIES = ["Inliers", "Outlier"]
COLORS = ["C0", "C3"]
ZORDER = [1, 2]
SCATTER_STYLE = dict(s=10, edgecolors="w", linewidths=0.2, rasterized=True)
feature_path = "/home/mwshi/project/MassiveQC/Features/features.parquet"

def outlier_umap(feature_path, inliers):
    feature_df = pd.read_parquet(feature_path).dropna()
    scaled_features = scale_feature(feature_df)
    umap_df = umap_feature(scaled_features)
    svg = Path(feature_path).parent / "umap_outliers.svg"
    plot_umap(umap_df, inliers, svg)



def scale_feature(feature_df):
    scaler = StandardScaler()
    scale_feature_df = pd.DataFrame(scaler.fit_transform(feature_df),
                                    index=feature_df.index, columns=feature_df.columns)
    return scale_feature_df


def umap_feature(scaled_features):
    X = scaled_features.values
    reducer = umap.UMAP(random_state=RANDOM_STATE)
    embeddings = reducer.fit_transform(X)
    umap_df = pd.DataFrame(
        embeddings, columns=["UMAP1", "UMAP2"], index=scaled_features.index
    )
    return umap_df


def plot_umap(umap_df, inliers, svg):
    umap_df['labels'] = umap_df.index.to_series().apply(
        lambda x: "Inliers" if x in inliers else "Outlier"
    )
    for cat, color, zorder in zip(CATEGORIES, COLORS, ZORDER):
        df = umap_df.query(f"labels == '{cat}'")
        plt.scatter(df.UMAP1, df.UMAP2, c=color, label=cat, zorder=zorder, **SCATTER_STYLE)

    ax = plt.gca()
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    plt.legend(loc="upper left")
    plt.savefig(svg)


NAME_MAPPER = {
    "percent_alignment": r"Reads Aligned (%)",
    "percent_duplication": f"Duplicated Reads (%)",
    "percent_genes_on": f"Gene Expressed (%)",
    "number_junction_reads": r"Reads at Exon Junctions (#)",
    "percent_utr_bases": r"UTR Bases (%)",
    "gene_body_three_prime": "3' Gene Body Coverage (avg)",
    "percent_intergenic_bases": r"Intergenic Bases (%)",
    "percent_reverse": r"Reverse Oriented Reads (%)",
    "number_genic_reads": "Genic Reads (#)",
    "gene_body_middle": "Middle Gene Body Coverage (avg)",
    "number_reads": "Total Reads (#)",
    "median_cv_coverage": "Coef Var Coverage (median)",
    "gene_body_five_prime": "5' Gene Body Coverage (avg)",
    "percent_mrna_bases": r"mRNA Bases (%)",
    "reads_MQ0": "Reads 0 Mapping Quality (#)",
    "percent_rrna_reads": r"rRNA Reads (%)",
    "number_junctions_on": "Exon Junctions Touched (#)",
    "percent_intronic_bases": r"Intronic Bases (%)",
    "number_multimapping_reads": "Multimapping Reads (#)",
    "average_quality": "Quality Score (avg)",
    "number_reads_too_short": "Reads Too Short (#)",
}


def plot_importance(iso, feature_path):
    mean_shap, columns = iso.mean_shap_values_outliers
    cols = [NAME_MAPPER[x] for x in columns]
    fig, ax = plt.subplots(figsize=plt.figaspect(1.2))
    sns.barplot(x=mean_shap, y=cols, color="C0", ax=ax)
    ax.set(xlabel="mean(|SHAP Value|)", ylabel="")
    svg = Path(feature_path).parent / "feature_importance.svg"
    plt.savefig(svg, bbox_inches='tight')