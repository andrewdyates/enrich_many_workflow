#!/usr/bin/python
"""Compute gene enrichment curves.


fname_rowlabels:
Tab-delimited file, each line corresponding to row labels.
File produced like in prune_dependency_matrix_workflow

FORMAT :
Original Row Number, Probe ID, Gene Symbol

EXAMPLE:
0	1007_s_at	DDR1
3	121_at	PAX8
11	1487_at	ESRRA
61	1552326_a_at	CCDC11

SAMPLE USE:
python $HOME/enrich_many_workflow/script.py fname_json=$HOME/enrich_many_workflow/sample.json fname_pina=$HOME/Homo-sapiens-20110628.txt fname_rowlabels=$HOME/gse7307/GSE7307_GPL570.symbol_rownums.gt0.25.pina.txt
"""
from gene_enrichment.pina import PINAEnriched
import gene_enrichment
import os
import json
import numpy as np
import sys


def main(fname_json, fname_pina, fname_rowlabels):
  Enrich = PINAEnriched(open(fname_pina))
  gene_syms = [s.split('\t')[2].strip('\n') for s in open(fname_rowlabels) if s.strip()]
  assert len(Enrich.vars) == len(gene_syms)
  assert not (Enrich.vars - set(gene_syms))
  J = json.load(open(fname_json))
  for dep in J["dependencies"]:
    print dep["function"]
    print dep["values_file"]
    M = np.load(os.path.join(J["dir"], dep["values_file"]))

    
if __name__ == "__main__":
  print sys.argv
  main(**dict([s.split('=') for s in sys.argv[1:]]))
             
