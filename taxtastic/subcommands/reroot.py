"""Reroots a reference package"""
import itertools
import argparse
import tempfile
import logging
import sys
import os

from Bio import Phylo

from taxtastic import algotax, refpkg

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refpkg', nargs=1,
        help='the reference package to operate on')
    parser.add_argument('-p', '--pretend',
        action='store_true', default=False,
        help="don't save the rerooted tree; just attempt the rerooting.")
    parser.add_argument('-f', '--force',
        action='store_true', default=False,
        help="reroot with the first found root, even if there are multiple taxonomic roots.")
    parser.add_argument('--write-taxonomic-roots', metavar='FILE',
        type=argparse.FileType('wb'), default=None,
        help="if specified, write out a file containing all of the taxonomic roots in Newick format.")

def action(args):
    log.info('loading reference package')
    rp = refpkg.Refpkg(args.refpkg[0])
    rp.load_db()
    with rp.resource('tree_file', 'rU') as fobj:
        tree = Phylo.read(fobj, 'newick')
    log.info('rerooting')
    stop_at_first_root = not (args.force or args.write_taxonomic_roots)
    orig_root = tree.root
    root, other_roots = algotax.reroot_from_rp(
        tree.root, rp, stop_at_first_root)
    if other_roots:
        if args.write_taxonomic_roots:
            def new_trees():
                for r in itertools.chain(other_roots, [root]):
                    tree.root_with_outgroup(r)
                    yield tree
            with args.write_taxonomic_roots as fobj:
                Phylo.write(new_trees(), fobj, 'newick')

        log.warn('%d other taxonomic roots found', len(other_roots))
        if stop_at_first_root:
            log.warn('there are also potentially still more roots')

        if not args.force:
            log.error('multiple roots and not a forced rerooting; exiting')
            sys.exit(2)
    tree.root_with_outgroup(root)
    if log.isEnabledFor(logging.DEBUG):
        log.debug('new root is %d steps from the old root',
                  len(tree.get_path(orig_root)))

    if args.pretend:
        return
    log.info('saving reference package')
    fd, name = tempfile.mkstemp(dir=args.refpkg[0])
    with os.fdopen(fd, 'w') as fobj:
        Phylo.write(tree, fobj, 'newick',
                    branchlengths_only=True,
                    format_branch_length='%0.6f')
    os.rename(name, rp.resource_path('tree_file'))
    rp.rehash('tree_file')
    rp.save()
