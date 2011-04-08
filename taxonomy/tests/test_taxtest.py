#!/usr/bin/env python

import sys
import os
from os import path
import unittest
import logging
import shutil
import commands
import re

import config

log = logging

_, module = path.split(__file__.replace('.py',''))
outputdir = path.abspath(config.outputdir)
datadir = path.abspath(config.datadir)

class TestBase(unittest.TestCase):

    executable = './taxtest.py'

    def __getitem__(self, i):
        """
        Enables string formatting, eg:
        print 'current function: %(funcname)s' % self
        """
        return getattr(self, i)

    def wrap_cmd(self, args = None, cmd = None):
        if cmd is None:
            cmd = self.executable
        input = (cmd + ' ' + args) % self
        log.warning('--> '+ input)
        status, output = commands.getstatusoutput(input)
        log.info(output)
        return status, output

    def cmd_ok(self, cmd=None, args=None):
        status, output = self.wrap_cmd(cmd, args)
        self.assertTrue(status == 0)

    def cmd_fails(self, cmd=None, args=None):
        status, output = self.wrap_cmd(cmd, args)
        self.assertFalse(status == 0)

    def setUp(self):
        self.funcname = '.'.join(self.id().split('.')[1:])
        self.outputdir = config.outputdir
        self.package = path.join(self.outputdir, self.funcname)

    def tearDown(self):
        pass

class TestHelp(TestBase):

    def test01(self):
        self.cmd_ok('-h')

    def test02(self):
        self.cmd_ok('--help')
        
    def test03(self):
        self.cmd_fails('')

    def test04(self):
        self.cmd_fails('notacommand')
        
    def test05(self):
        self.cmd_ok('help create')

    def test06(self):
        self.cmd_ok('help check')

    def test07(self):
        self.cmd_ok('help help')

    def test08(self):
        self.cmd_fails('help notacommand')

    def test09(self):
        self.cmd_ok('--version')
        
        
class TestCreate(TestBase):

    def test01(self):
        self.cmd_fails('create -P %(package)s')

    def test02(self):
        """
        Create a minimal package.
        """
        
        shutil.rmtree(self.package)
        self.cmd_ok('create -P %(package)s -l 16s')
        self.assertTrue(path.exists(self.package))
        self.assertTrue(path.exists(path.join(self.package,'CONTENTS.json')))
        
        # fails the second time because package already exists
        self.cmd_fails('create -P %(package)s -l 16s')
