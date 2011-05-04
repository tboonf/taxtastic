"""Generates a csv file representing the ratio of leaves to edges in a tree."""

from Bio import Phylo

from collections import Counter, defaultdict
from itertools import combinations
import argparse
import csv
import sys

from taxtastic import algotax, refpkg

def build_parser(parser):
    parser.add_argument('refpkg', nargs=1,
        help='the reference package to operate on')
    parser.add_argument('-o', '--outfile',
        type=argparse.FileType('w'), default='-',
        help='the name of the csv file to write out')

def action(args):
    args.outfile = csv.writer(args.outfile)
    args.outfile.writerow(['rank', 'color', 'leaves', 'min_psi', 'total_psi'])
    rp = refpkg.Refpkg(args.refpkg[0])
    rp.load_db()
    with rp.resource('tree_file', 'rU') as fobj:
        tree = Phylo.read(fobj, 'newick')
    clade_map = {c.name: c for c in tree.get_terminals()}

    curs = rp.db.cursor()
    seq_colors = curs.execute("""
        SELECT t.rank,
               seqname,
               t.tax_name
        FROM   parents
               JOIN taxa t
                 ON t.tax_id = parent
               JOIN sequences s
                 ON s.tax_id = child
    """)
    rank_map = defaultdict(list)
    for rank, seq, color in seq_colors:
        rank_map[rank].append((seq, color))

    rank_order = curs.execute("""
        SELECT   rank
        FROM     ranks
        ORDER BY rank_order ASC
    """)

    for rank, in rank_order:
        seq_colors = rank_map[rank]
        if not seq_colors:
            continue
        colors = {}
        for seq, color in seq_colors:
            colors[clade_map[seq]] = color
        metadata = algotax.color_clades(tree, colors)
        psis = defaultdict(Counter)

        for clade in tree.find_elements(order='postorder'):
            if clade is tree: continue
            if clade.is_terminal():
                psis[clade][colors[clade]] += 1
                continue
            for child in clade:
                psis[clade] += psis[child]

        for color in metadata.cut_colors[tree.root]:
            for clade in tree.find_elements(terminal=False):
                if clade is tree: continue
                if any(psis[child][color] == 0 for child in clade):
                    continue
                args.outfile.writerow([
                    rank,
                    color,
                    clade.count_terminals(),
                    min(psis[child][color] for child in clade),
                    sum(psis[child][color] for child in clade),
                ])
