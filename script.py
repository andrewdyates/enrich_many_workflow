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
python $HOME/enrich_many_workflow/script.py fname_json=$HOME/gse7307/dependency_matrices_single_gene_qnorm.json fname_pina=$HOME/Homo-sapiens-20110628.txt fname_rowlabels=$HOME/gse7307/GSE7307_GPL570.symbol_rownums.gt0.25.pina.txt outdir=$HOME/gse7307
"""
from gene_enrichment.pina import PINAEnriched
import gene_enrichment
from gene_enrichment.py_symmetric_matrix import *
import os
import json
import numpy as np
import sys
from lab_util import *


def main(fname_json, fname_pina, fname_rowlabels, outdir=""):
  if outdir != "":
    make_dir(outdir)
  
  Enrich = PINAEnriched(open(fname_pina))
  gene_syms = [s.split('\t')[2].strip('\n') for s in open(fname_rowlabels) if s.strip()]
  assert len(Enrich.vars) >= len(gene_syms)
  for s in gene_syms:
    assert Enrich.is_in(s)
  n = len(gene_syms)
  J = json.load(open(fname_json))

  for w in (1000, 10000, 100000, 1000000, 5000000):
    outpath = os.path.join(outdir, "gse7307_enrichment_%d.pdf" % w)
    print "Generating enrichment curves for w=%d..." % w
    Results = make_ranks(Enrich, J, w, n, gene_syms)
    gene_enrichment.make_enrichment_curve_figure(\
      title="GSE7307 top %d Enrichment" % w, Ranks=Results, plotpath=outpath)
    print "Saved enrichment curves as %s." % outpath
    print


def enrich_rank(Enrich, M, w, n, gene_syms):
  Q = M.argsort()[::-1] # Copy argsorted in reverse order
  self.assertEqual(np.size(Q), np.size(M))
  g = []
  print "Rank Check: Top 2:", M[Q[0]], M[Q[1]]
  print "Number of dependencies == 0 (may indicate missing values):", np.sum(M==0)
  for i in xrange(w):
    try:
      x, y = inv_sym_idx(Q[i], n)
    except IndexError:
      print "!", i, n
      raise
    if Enrich.exists(gene_syms[x], gene_syms[y]):
      g.append(i)
  print "last value considered:", M[Q[w]]
  print "matches in top (%d) rank: %d" % (w, len(g))
  print
  return g


def make_ranks(Enrich, J, w, n, gene_syms):
  Results = {}

  # Random.
  M = np.random.random(n*(n-1)/2)
  print "Random..."
  Results["Random"] = enrich_rank(Enrich, M, w, n, gene_syms) 
  
  for dep in J["dependencies"]:
    print dep["function"], dep["values_file"]
    M = np.load(os.path.join(J["dir"], dep["values_file"]))
    n_expected = sym_idx(n-2,n-1,n) + 1
    assert np.size(M,0) == n_expected == n*(n-1)/2
    Results[dep["function"]] = enrich_rank(Enrich, M, w, n, gene_syms)
    
    # Absolute and negative values
    if "abs" in dep and dep["abs"]:
      print dep["function"], "absolute value"
      Results["abs(%s)"%dep["function"]] = enrich_rank(Enrich, np.abs(M), w, n, gene_syms)
      print dep["function"], "negative value"
      Results["neg(%s)"%dep["function"]] = enrich_rank(Enrich, -M), w, n, gene_syms)
      
  return Results


  
if __name__ == "__main__":
  print sys.argv
  main(**dict([s.split('=') for s in sys.argv[1:]]))
             
