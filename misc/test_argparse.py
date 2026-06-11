
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--iter", help="Flag to use score_genes_iter instead of score_genes",
                    action="store_true")
args = parser.parse_args()

if args.iter:
    print("Scoring iteratively")
else:
    print("Normal score_genes")