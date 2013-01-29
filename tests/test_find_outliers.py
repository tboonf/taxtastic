from Bio import Phylo
from StringIO import StringIO
import unittest
from taxtastic import find_outliers

"""
Here's how some of the data files used in testing this module came to be:

src=/shared/silo_researcher/Matsen_F/MatsenGrp/micro_refset/rdp/10_30/1200bp/named
dest=~/src/taxtastic/testfiles

csvgrep $src/rdp_10_30_named_1200bp.seq_info.csv -c description -m 'Enterococcus faecium' | head -101 > $dest/e_faecium.seq_info.csv
csvcut $dest/e_faecium.seq_info.csv -c 1 | \
seqmagick convert --include-from-file - $src/rdp_10_30_named_1200bp.fasta $dest/e_faecium.fasta
cmalign.py $dest/e_faecium.fasta $dest/e_faecium.sto
seqmagick convert $dest/e_faecium.sto $dest/e_faecium.aln.fasta
FastTree -nt -makematrix $dest/e_faecium.aln.fasta > $dest/e_faecium.distmat

csvgrep $src/rdp_10_30_named_1200bp.seq_info.csv -c description -m 'Enterococcus faecalis' | head -101 > $dest/e_faecalis.seq_info.csv
csvcut $dest/e_faecalis.seq_info.csv -c 1 | \
seqmagick convert --include-from-file - $src/rdp_10_30_named_1200bp.fasta $dest/e_faecalis.fasta
cmalign.py $dest/e_faecalis.fasta $dest/e_faecalis.sto
seqmagick convert $dest/e_faecalis.sto $dest/e_faecalis.aln.fasta
FastTree -nt -makematrix $dest/e_faecalis.aln.fasta > $dest/e_faecalis.distmat
"""

from . import config

class TestReadDists(unittest.TestCase):

    def test_01(self):
        with open(config.data_path('e_faecium.distmat')) as f:
            taxa, mat = find_outliers.read_dists(f)
            self.assertEqual(len(taxa), 100)
            self.assertEqual(mat.shape, (100, 100))

class TestFastTreeDists(unittest.TestCase):

    def test_01(self):
        fa = config.data_path('e_faecium.aln.fasta')
        taxa, mat = find_outliers.fasttree_dists(fa)
        self.assertEqual(len(taxa), 100)
        self.assertEqual(mat.shape, (100, 100))

class TestFindOutliers(unittest.TestCase):

    def test_01(self):
        with open(config.data_path('e_faecium.distmat')) as f:
            taxa, mat = find_outliers.read_dists(f)
            is_outlier = find_outliers.outliers(mat, cutoff = 0.015)
            outliers = {t for t,o in zip(taxa, is_outlier) if o}
            self.asserEqual(len(outliers), 7)
