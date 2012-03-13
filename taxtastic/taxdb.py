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
Methods and variables for taxonomy databases.
"""

import itertools
import logging
import os
import sqlite3
import urllib

log = logging

db_schema = """
-- nodes.dmp specifies additional columns but these are not implemented yet
CREATE TABLE nodes(
tax_id        TEXT UNIQUE PRIMARY KEY NOT NULL,
parent_id     TEXT,
rank          TEXT,
embl_code     TEXT,
division_id   INTEGER,
source_id     INTEGER DEFAULT 1 -- added to support multiple sources
);

CREATE TABLE names(
tax_id        TEXT REFERENCES nodes(tax_id),
tax_name      TEXT,
unique_name   TEXT,
name_class    TEXT,
is_primary    INTEGER -- not defined in names.dmp
);

CREATE TABLE merged(
old_tax_id    TEXT,
new_tax_id    TEXT REFERENCES nodes(tax_id)
);

-- table "source" supports addition of custom taxa (not provided by NCBI)
CREATE TABLE source(
id            INTEGER PRIMARY KEY AUTOINCREMENT,
name          TEXT UNIQUE,
description   TEXT
);

INSERT INTO "source"
  (id, name, description)
VALUES
  (1, "NCBI", "NCBI taxonomy");
INSERT INTO "source"
  (id, name, description)
VALUES
  (2, "GreenGenes", "GreenGenes taxonomy");

-- indices on nodes
CREATE INDEX nodes_tax_id ON nodes(tax_id);
CREATE INDEX nodes_parent_id ON nodes(parent_id);
CREATE INDEX nodes_rank ON nodes(rank);

-- indices on names
CREATE INDEX names_tax_id ON names(tax_id);
CREATE INDEX names_tax_name ON names(tax_name);
CREATE INDEX names_is_primary ON names(is_primary);
CREATE INDEX names_taxid_is_primary ON names(tax_id, is_primary);
CREATE INDEX names_name_is_primary ON names(tax_name, is_primary);

-- CREATE UNIQUE INDEX names_id_name ON names(tax_id, tax_name, is_primary);

"""

# define headers in names.dmp, etc (may not correspond to table columns above)
merged_keys = 'old_tax_id new_tax_id'.split()

undefined_rank = 'no_rank'
root_name = 'root'

_ranks = """
root
superkingdom
kingdom
subkingdom
superphylum
phylum
subphylum
superclass
class
subclass
infraclass
superorder
order
suborder
infraorder
parvorder
superfamily
family
subfamily
tribe
subtribe
genus
subgenus
species group
species subgroup
species
subspecies
varietas
forma
"""

ranks = [k.strip().replace(' ','_') for k in _ranks.splitlines() if k.strip()]

def db_connect(dbname='ncbi_taxonomy.db', schema=db_schema, clobber=False):
    """
    Create a connection object to a database. Attempt to establish a
    schema. If there are existing tables, delete them if clobber is
    True and return otherwise. Returns a connection object.
    """

    if clobber:
        log.info('Creating new database %s' % dbname)
        try:
            os.remove(dbname)
        except OSError:
            pass

    con = sqlite3.connect(dbname)
    cur = con.cursor()

    cmds = [cmd.strip() for cmd in schema.split(';') if cmd.strip()]
    try:
        for cmd in cmds:
            cur.execute(cmd)
            log.debug(cmd)
    except sqlite3.OperationalError as err:
        log.warn(err)

    return con

def has_row(cursor, tablename):
    """
    Check if tablename contains any rows
    """
    cursor.execute("SELECT * FROM \"{0}\" LIMIT 1".format(tablename))
    return cursor.fetchone() is not None

def do_insert(con, tablename, rows, maxrows=None, add=True, colnames=None):
    """
    Insert rows into a table. Do not perform the insert if
    add is False and table already contains data.
    """

    cur = con.cursor()

    cur.execute('select count(*) from "%s" where rowid < 2' % tablename)
    has_data = cur.fetchone()[0]

    if not add and has_data:
        log.info('Table "%s" already contains data; load not performed.' % tablename)
        return False

    # pop first row to determine number of columns
    row = rows.next()
    if colnames:
        cmd = 'INSERT INTO "%s" (%s) VALUES (%s)' % (tablename,
                                         ', '.join(colnames),
                                         ', '.join(['?'] * len(colnames)))
    else:
        cmd = 'INSERT INTO "%s" VALUES (%s)' % (tablename, ', '.join(['?']*len(row)))
    log.info(cmd)

    # put the first row back
    rows = itertools.chain([row], rows)
    if maxrows:
        rows = itertools.islice(rows, maxrows)

    cur.executemany(cmd, rows)
    con.commit()

    return True

def fetch_url(url, dest_dir='.', clobber=False):
    """
    Download data from url to dest_dir

    * url - url to download
    * dest_dir - directory in which to save output files (created if necessary).
    * clobber - don't download if False and target of url exists in dest_dir

    Returns (fname, downloaded), where fname is the name of the downloaded
    file, and downloaded is True if a new file was downloaded, false otherwise.
    """

    dest_dir = os.path.abspath(dest_dir)
    try:
        os.mkdir(dest_dir)
    except OSError:
        pass

    fout = os.path.join(dest_dir, os.path.split(url)[-1])

    if os.access(fout, os.F_OK) and not clobber:
        downloaded = False
        log.warning('%s exists; not downloading' % fout)
    else:
        downloaded = True
        log.warning('downloading %(url)s to %(fout)s' % locals())
        urllib.urlretrieve(url, fout)

    return fout, downloaded


class OnUpdate(object):
    def __init__(self, proxied):
        self.proxied = proxied
        self.setter = None

    def __get__(self, inst, cls):
        if inst is None:
            return self
        return getattr(inst, self.proxied)

    def __set__(self, inst, value):
        if self.setter:
            self.setter(inst, value)
        setattr(inst, self.proxied, value)

    def __call__(self, f):
        self.setter = f
        return self

class _IntermediateTaxon(object):
    def __init__(self, tax_id, parent, rank, tax_name):
        self.children = set()
        self.tax_id = tax_id
        self.parent = parent
        self.rank = rank
        self.tax_name = tax_name

    _parent = _adjacent_to = None

    @OnUpdate('_parent')
    def parent(self, p):
        if self.parent is not None:
            self.parent.children.discard(self)
        if p is not None:
            p.children.add(self)

    @OnUpdate('_adjacent_to')
    def adjacent_to(self, n):
        if n is not None:
            self.parent = n.parent

    def iterate_children(self, on_pop=None, including_self=True):
        search_stack = [(None, set([self]))]
        while search_stack:
            if not search_stack[-1][1]:
                parent, _ = search_stack.pop()
                if on_pop:
                    on_pop(parent)
                continue
            node = search_stack[-1][1].pop()
            if node is not self or including_self:
                yield node
            search_stack.append((node, node.children.copy()))

class Taxdb(object):
    def __init__(self, sqlite_db=None):
        if sqlite_db is None:
            sqlite_db = sqlite3.connect(':memory:')
        self.db = sqlite_db

    def __getattr__(self, attr):
        return getattr(self.db, attr)

    def create_tables(self):
        curs = self.db.cursor()

        curs.execute("""
            CREATE TABLE ranks (
              rank TEXT PRIMARY KEY NOT NULL,
              rank_order INT
            )
        """)

        curs.execute("""
            CREATE TABLE taxa (
              tax_id TEXT PRIMARY KEY NOT NULL,
              tax_name TEXT NOT NULL,
              rank TEXT REFERENCES ranks (rank) NOT NULL
            )
        """)

        curs.execute("""
            CREATE TABLE sequences (
              seqname TEXT PRIMARY KEY NOT NULL,
              tax_id TEXT REFERENCES taxa (tax_id) NOT NULL
            )
        """)

        curs.execute("""
            CREATE TABLE hierarchy (
              tax_id TEXT REFERENCES taxa (tax_id) PRIMARY KEY NOT NULL,
              lft INT NOT NULL UNIQUE,
              rgt INT NOT NULL UNIQUE
            )
        """)

        curs.execute("""
            CREATE VIEW parents AS
            SELECT h1.tax_id AS child,
                   h2.tax_id AS parent
            FROM   hierarchy h1
                   JOIN hierarchy h2
                     ON h1.lft BETWEEN h2.lft AND h2.rgt
        """)

    def insert_from_taxtable(self, fieldnames_cb, table):
        curs = self.db.cursor()

        taxon_map = {}
        for row in table:
            parent = taxon_map.get(row['parent_id'])
            taxon = _IntermediateTaxon(
                row['tax_id'], parent, row['rank'], row['tax_name'])
            taxon_map[taxon.tax_id] = taxon

        root = next(taxon_map.itervalues())
        while root.parent is not None:
            root = root.parent
        counter = itertools.count(1).next
        def on_pop(parent):
            if parent is not None:
                parent.rgt = counter()
        for node in root.iterate_children(on_pop=on_pop):
            node.lft = counter()

        fieldnames = fieldnames_cb()
        curs.executemany("INSERT INTO ranks (rank_order, rank) VALUES (?, ?)",
            enumerate(fieldnames[4:]))
        curs.executemany("INSERT INTO taxa VALUES (?, ?, ?)",
            ((t.tax_id, t.tax_name, t.rank) for t in taxon_map.itervalues()))
        curs.executemany("INSERT INTO hierarchy VALUES (?, ?, ?)",
            ((t.tax_id, t.lft, t.rgt) for t in taxon_map.itervalues()))
        self.db.commit()
