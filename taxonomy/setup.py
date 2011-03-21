"""
Create unix package:    python setup.py sdist
"""

from distutils.core import setup
import os
import sys
import glob

from __init__ import __version__

params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'name': 'taxonomy',
          'package_dir': {'taxonomy': '.'},
          'packages': ['taxonomy'],
          'scripts': glob.glob('scripts/*.py'),
          # 'package_data':{'taxonomy': glob.glob('data/*')},
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': __version__,
          'requires': ['Python (>= 2.6)', 'sqlalchemy']}

setup(**params)

