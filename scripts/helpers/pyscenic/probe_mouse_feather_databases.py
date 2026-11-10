import sys
from pathlib import Path
from ctxcore.rnkdb import FeatherRankingDatabase

gene_file = Path(sys.argv[1])
db_files = [Path(x) for x in sys.argv[2:]]

with gene_file.open('r', encoding='utf-8') as handle:
    expression_genes = {line.strip().split('\t')[0] for line in handle if line.strip()}

for idx, db_file in enumerate(db_files, start=1):
    try:
        db = FeatherRankingDatabase(fname=str(db_file), name=f'db_{idx}')
        db_genes = set(db.genes)
        overlap = expression_genes.intersection(db_genes)
        fraction = len(overlap) / max(1, len(expression_genes))
        print(
            str(db_file) + '\tPASS\t' + str(db.total_genes) + '\t' +
            str(len(expression_genes)) + '\t' + str(len(overlap)) + '\t' +
            str(fraction) + '\t'
        )
    except Exception as exc:
        text = str(exc).replace('\t', ' ').replace('\n', ' ')
        print(str(db_file) + '\tFAIL\tNA\t' + str(len(expression_genes)) +
              '\tNA\tNA\t' + text)
