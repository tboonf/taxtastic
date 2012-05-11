# This file is part of taxtastic.
#
#    taxtastic is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    taxtastic is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with taxtastic.  If not, see <http://www.gnu.org/licenses/>.
"""
Methods and variables specific to the GreenGenes taxonomy.
"""
import bisect
import collections
import logging
import itertools
import os.path
import sqlite3

from . import taxdb
from .errors import IntegrityError

log = logging

ranks = taxdb.ranks

# Map from rank abbreviation to taxdb rank name
_rank_map = {'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o': 'order',
        'f': 'family', 'g': 'genus', 's': 'species'}

download_url = 'http://www.secondgenome.com/go/2011-greengenes-taxonomy/'
tax_map = 'taxonomy_16S_all_gg_2011_1'

db_connect = taxdb.db_connect

def _get_source_id(con):
    """Get the source.id corresponding to greengenes"""
    cursor = con.cursor()
    cursor.execute('SELECT id FROM source where name = ?', ['greengenes'])
    return cursor.fetchone()

def _add_space(species, genus):
    """
    Add a space between genus and species
    """
    l = len(genus)
    if species.startswith(genus) and species[l] != ' ':
        return ' '.join((species[:l], species[l:]))
    return species

def _clean_species(con):
    """
    Try to clean species-level classifications

    For every species that doesn't start with it's genus, try to find the genus
    and insert a space between genus and species.

    :param con: Database connection
    """
    cursor = con.cursor()
    all_genus = dict(cursor.execute("""
    SELECT names.tax_name, names.tax_id
    FROM nodes JOIN names USING (tax_id)
    WHERE nodes.rank = 'genus'"""))

    # Sorted list of all genus in the database, for searching
    keys = sorted(all_genus)

    def find_genus(tax_name):
        """Search for a genus in the tax_name"""
        # We don't know where the genus/species separation is, so search from
        # the first genus starting with tax_name[0]
        start = bisect.bisect_left(keys, tax_name[0])
        poss_genus = itertools.islice(keys, start, None)
        poss_genus = itertools.takewhile(lambda k: k < tax_name, poss_genus)
        for genus in poss_genus:
            if tax_name.startswith(genus):
                return genus

    # All species names that don't start with the genus and a space
    cursor.execute("""SELECT DISTINCT snode.tax_id, sname.tax_name
    FROM nodes snode INNER JOIN nodes gnode ON snode.parent_id = gnode.tax_id
        INNER JOIN names sname ON sname.tax_id = snode.tax_id
        INNER JOIN names gname ON gname.tax_id = gnode.tax_id
    WHERE sname.is_primary = 1 AND
       snode.rank = 'species' AND gnode.rank = 'genus' AND
       SUBSTR(sname.tax_name, 1, LENGTH(gname.tax_name) + 1) <>
           gname.tax_name + ' ';""")

    with con:
        stmt = "UPDATE names SET tax_name = ? WHERE tax_id = ? AND tax_name = ?"
        args = []
        for tax_id, tax_name in cursor:
            genus = find_genus(tax_name)
            if genus:
                new_tax_name = _add_space(tax_name, genus)
                args.append((new_tax_name, tax_id, tax_name))
        cursor.executemany(stmt, args)

def _load_taxonomy(rows, con):
    """
    Load GreenGenes taxonomy into database.

    Sequences are labeled with generated tax_ids starting from GG00000001

    No action is taken if tables are already populated.

    :param rows: Iterable containing lists of (rank_name, rank) tuples by descending rank
    :param con: database connection
    """
    source_id = _get_source_id(con)
    cursor = con.cursor()
    if taxdb.has_row(cursor, 'nodes') or taxdb.has_row(cursor, 'names'):
        log.warning("nodes table has data. not updating.")
        return

    count = itertools.count()
    # map from (name, parent_id) to tax_id
    tax_ids = {}
    # map from tax_id to name
    tax_names = {}
    # Parent tax_ids for a (tax_name, rank) pair
    parents = collections.defaultdict(set)

    def get_or_insert(name, rank, parent_id, source_id, alternate_names=None):
        try:
            return tax_ids[name, parent_id]
        except KeyError:
            pass

        # Generate a new tax_id
        tax_id = 'GG{0:08d}'.format(count.next())
        tax_ids[name, parent_id] = tax_id
        tax_names[tax_id] = name
        parents[name, rank].add(parent_id)

        # If no parent id, node should be it's own parent (root)
        if parent_id is None:
            parent_id = tax_id

        # Add to database
        cursor.execute("""INSERT INTO nodes (tax_id, parent_id, rank, source_id)
                VALUES (?, ?, ?, ?)""", (tax_id, parent_id, rank, source_id))
        cursor.execute("""INSERT INTO names (tax_id, tax_name, is_primary)
                VALUES (?, ?, 1)""", (tax_id, name))

        return tax_id

    def alternate_name(tax_id, name, name_class='greengenes lineage'):
        cursor.execute("SELECT tax_id FROM names WHERE tax_name = ? AND name_class = ?",
                (name, name_class))
        tax_id = cursor.fetchone()
        if tax_id:
            return
        else:
            cursor.execute("""INSERT INTO names (tax_id, tax_name, is_primary, name_class)
VALUES (?, ?, 0, ?)""", (parent, orig, name_class))

    def fix_polyphyly():
        """
        Give a unique name to each polyphyletic group
        """
        polyphyly = ((k, v) for k, v in parents.items() if len(v) > 1)
        for (child, rank), p in polyphyly:
            logging.warn("Renaming items from polyphyletic %s %s", rank, child)
            for parent in p:
                tax_id = tax_ids[child, parent]
                parent_name = tax_names[parent]
                new_name = '{0} [{1}]'.format(child, parent_name)
                logging.warn("  '%s' -> '%s'", child, new_name)
                cursor.execute("""UPDATE names
SET tax_name = ?
WHERE tax_id = ? AND tax_name = ?""", (new_name, tax_id, child))
                assert cursor.rowcount == 1

    with con:
        for classes, orig in rows:
            parent = None
            for rank, name in classes:
                parent = get_or_insert(name, rank, parent, source_id)
            alternate_name(parent, orig)

        # Fixup
        fix_polyphyly()

def db_load(con, fname, maxrows=None):
    """
    Load data from taxonomic map into database identified by con. Data
    is not loaded if target tables already contain data.

    :param con: Database connection
    :param archive: Tar archive path
    :param maxrows: Maximum rows to load (ignored)
    """
    with open(fname) as handle:
        try:
            cursor = con.cursor()
            if not taxdb.has_row(cursor, 'nodes'):
                next(handle) # Skip header
                rows = _parse_greengenes_lineages(handle)
                _load_taxonomy(rows, con)
            _clean_species(con)
        except sqlite3.IntegrityError, err:
            raise IntegrityError(err)

def fetch_data(dest_dir='.', clobber=False, url=None):
    """
    Download data from greengenes required to generate local taxonomy
    database.

    :param dest_dir: directory in which to save output files (created if necessary).
    :param boolean clobber: don't download if False and target of url exists in dest_dir
    :param url: url to archive;

    Returns (fname, downloaded), where fname is the name of the
    downloaded zip archive, and downloaded is True if a new files was
    downloaded, false otherwise.
    """
    log.warn("Downloading is not implemented for GreenGenes.\n"
            "Please download '%s' manually from '%s'", tax_map, download_url)
    dest = os.path.join(dest_dir, tax_map)
    if not os.path.exists(dest):
        raise IOError("Missing file: {0}".format(dest))
    return dest, False

def _parse_classes(classes):
    """
    Parse classes from GreenGenes taxonomy

    classes - GreenGenes taxonomy assignments, as semicolon delimited list of
    [rank_key]__[class] items, e.g.
    'k__Bacteria';p__Proteobacteria'

    If species starts with the genus, a space is inserted after genus

    Returns list of (rank, name) tuples

    :param classes: semicolon delimited lineage from GreenGenes
    :param otu_id: OTU number
    """
    split = [i.strip().split('__') for i in classes.split(';')]
    # Start with a single root node
    result = [[taxdb.root_name, taxdb.root_name]]

    # Add classifications
    result.extend([_rank_map[cls_key], cls or None]
                  for cls_key, cls in split if cls)

    # special handling for species -> genus
    if result[-1][0] == 'species' and result[-2][0] == 'genus':
        species = result[-1][1]
        genus = result[-2][1]
        if species.startswith(genus) and species[len(genus)] != ' ':
            result[-1][1] = _add_space(species, genus)
    return result

def _parse_greengenes_lineages(handle):
    """
    Parse GreenGenes lineages from a file handle.

    Yields a list of parsed lineages for each OTU

    :param handle: File-like object
    """
    for line in handle:
        otu_id, classes = line.rstrip().split('\t')
        parsed = _parse_classes(classes)
        yield parsed, classes
