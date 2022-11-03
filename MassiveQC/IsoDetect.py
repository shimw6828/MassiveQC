import argparse
from pathlib import Path

from .feature_store import check_done_sample, feature_store
from .build_features import build_features
from .detection import detection



def main():
    """Combine multiple modules to provide commands for cli"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str, help='Input file, containing two columns srx and srr')
    parser.add_argument('-o', '--outdir', required=True, type=str, help="Path to result output directory of main process.")
    args = parser.parse_args()
    done_samples = check_done_sample(args.outdir)
    feature_store(done_samples, args.outdir)
    done_samples = check_done_sample(args.outdir)
    feature_store(done_samples, args.outdir)
    Features = Path(args.outdir) / "Features"
    build_features(args.input, Features)
    features_file = (Features / "features.parquet").as_posix()
    detection(features_file)


if __name__ == "__main__":
    main()