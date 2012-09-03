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
from gene_enrichment.py_symmetric_matrix import *
import os
import json
import numpy as np
import sys


def main(fname_json, fname_pina, fname_rowlabels, w=10000):
  Enrich = PINAEnriched(open(fname_pina))
  gene_syms = [s.split('\t')[2].strip('\n') for s in open(fname_rowlabels) if s.strip()]
  assert len(Enrich.vars) >= len(gene_syms)
  for s in gene_syms:
    assert Enrich.is_in(s)
  n = len(gene_syms)
  J = json.load(open(fname_json))
  Results = {}

  for dep in J["dependencies"]:
    print dep["function"]
    print dep["values_file"]
    M = np.load(os.path.join(J["dir"], dep["values_file"]))
    n_expected = sym_idx(n-2,n-1,n) + 1
    assert np.size(M,0) == n_expected
    Q = M.argsort()[::-1]
    g = []
    print M[Q[0]], M[Q[1]], np.sum(M==0)
    for i in xrange(w):
      x, y = inv_sym_idx(Q[i], n)
      if Enrich.exists(gene_syms[x], gene_syms[y]):
        g.append(i)
    print "last value considered:", M[Q[w]]
    print "matches in top (%d) rank: %d" % (w, len(g))
    Results[dep["function"]] = g
    

    # Absolute value
    if "abs" in dep and dep["abs"]:
      print dep["function"], "absolute value"
      Q = np.abs(M).argsort()[::-1]
      g = []
      print M[Q[0]], M[Q[1]], np.sum(M==0)
      for i in xrange(w):
        x, y = inv_sym_idx(Q[i], n)
        if Enrich.exists(gene_syms[x], gene_syms[y]):
          g.append(i)
      print "last value considered:", M[Q[w]]
      print "matches in top (%d) rank: %d" % (w, len(g))
      Results[dep["function"]+"_abs"] = g
      
  # Plot enrichment curves
  gene_enrichment.make_enrichment_curve_figure(\
    title="GSE7307 top %d Enrichment" % w, Ranks=Results, plotpath="enrichment.png")
  

    
if __name__ == "__main__":
  print sys.argv
  main(**dict([s.split('=') for s in sys.argv[1:]]))
             
