#!/usr/bin/env python

import sys
import os
import unittest
import logging
import shutil

import config
import taxonomy

log = logging

module_name = os.path.split(sys.argv[0])[1].rstrip('.py')
outputdir = os.path.abspath(config.outputdir)
datadir = os.path.abspath(config.datadir)

def newdir(path):
    """
    Create a new directory "path", deleting an existing directory if
    necessary.
    """

    shutil.rmtree(path, ignore_errors = True)
    os.makedirs(path)
    
class TestFetchData(unittest.TestCase):
    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.outdir = os.path.join(outputdir, self.funcname)
        newdir(self.outdir)
        _, self.zfilename = os.path.split(taxonomy.ncbi.ncbi_data_url)
        
    def test01(self):
        zfile = os.path.join(self.outdir, self.zfilename)        
        fout, downloaded = taxonomy.ncbi.fetch_data(dest_dir=self.outdir)

        # file is downloaded the first time
        self.assertTrue(downloaded)
        self.assertTrue(os.path.isfile(fout))
        self.assertTrue(zfile == fout)

        # ... but not the second time
        fout, downloaded = taxonomy.ncbi.fetch_data(dest_dir=self.outdir)
        self.assertFalse(downloaded)

        # ... unless clobber = True
        fout, downloaded = taxonomy.ncbi.fetch_data(dest_dir=self.outdir, clobber=True)
        self.assertTrue(downloaded)

        
class TestCreateSchema(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.dbname = os.path.join(outputdir, self.funcname + '.db')
        log.info(self.dbname)

    def test01(self):
        with taxonomy.ncbi.db_connect(self.dbname, new=True) as con:
            cur = con.cursor()
            cur.execute('select name from sqlite_master where type = "table"')
            tables = set(i for j in cur.fetchall() for i in j) # flattened
            self.assertTrue(set(['nodes','names','merged','source']).issubset(tables))
        
class TestLoadData(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.dbname = os.path.join(outputdir, self.funcname + '.db')
        self.zfile = os.path.join(outputdir, 'taxdmp.zip')

    def test01(self):
        con = taxonomy.ncbi.db_connect(self.dbname, new=True)
        taxonomy.ncbi.db_load(con, self.zfile, maxrows=10)
        con.close()

    def test02(self):
        con = taxonomy.ncbi.db_connect(self.dbname, new=True)
        taxonomy.ncbi.db_load(con, self.zfile, maxrows=10)
        con.close()



