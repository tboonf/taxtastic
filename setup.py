"""
Create unix package:    python setup.py sdist
"""

import os
import subprocess
import shutil
from os.path import join
import glob

try:
    from setuptools import setup, find_packages
except ImportError:
    import distribute_setup
    distribute_setup.use_setuptools()
    from setuptools import setup, find_packages

# The code contained by the with statement below configures and
# updates the insertion of the abbrevated git sha hash in the package
# version number. This assumes that the repository root contains a
# '.gitattributes' file with the contents
# 
# taxtastic/data/sha    filter=sha
#
# and that the (empty) file 'taxtastic/data/sha' exists.
with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
    outputs = {'stdout': stdout, 'stderr': stderr}
    # outputs = {}
    try:
        subprocess.check_call(['git', 'status'], **outputs) 
    except subprocess.CalledProcessError:
        git_available = False
    else:
        git_available = True

    # The first time a repo is cloned, the local configuration needs
    # to be updated to add the attributes 'filter.sha.clean' and
    # 'filter.sha.smudge'. Do this now if necessary.        
    if git_available and subprocess.call(['git', 'config', '--get-regexp', 'filter.sha'], **outputs) != 0:
        print 'updating git config'
        subprocess.call(['git', 'config', '--remove-section', 'filter.sha'], **outputs)    
        gitconf = '.git/config'
        shutil.copyfile(gitconf, 'git-config.bak')
        with open(gitconf, 'a') as fobj:
            fobj.write('\n[filter "sha"]\n')
            fobj.write('clean = cat > /dev/null\n')
            fobj.write(r'smudge = echo "$(git --no-pager log --pretty=format:\"%h\" -1)"')
        subprocess.call(['git', 'config', '--get-regexp', 'filter.sha'])

    # Make sure that the version number is up to date if we're
    # installig from a git repo. Note that the version number is
    # performed even if installation is not (eg, 'python setup.py
    # -h'). The version number is stored in `shafile`; we must make
    # sure that the shafile retains its original attributes since this
    # command may be run as root.
    if git_available:
        shafile = join('taxtastic','data','sha')
        stats = os.stat(shafile)
        os.rename(shafile, shafile+'.bak')
        subprocess.check_call(['git', 'checkout', shafile], **outputs)
        os.chown(shafile, stats.st_uid, stats.st_gid)
        new, old = [open(f).read().strip() for f in [shafile, shafile+'.bak']]
        if old != new:
            print 'updated version sha: %s --> %s' % (old, new)

            
from taxtastic import __version__

# all files with .py extension in top level are assumed to be scripts
scripts = ['taxit'] + list(set(glob.glob('*.py')) - set(['setup.py', 'distribute_setup.py']))

params = {'author': 'Noah Hoffman',
          'author_email': 'ngh2@uw.edu',
          'description': 'Tools for taxonomic naming and annotation',
          'name': 'taxtastic',
          'packages': find_packages(exclude=['tests']),
          'scripts': scripts,
          'url': 'https://github.com/fhcrc/taxtastic',
          'version': __version__,
          'package_data': {'taxtastic': [join('data',f) for f in ['sha']]},
          'install_requires': ['sqlalchemy', 'decorator']}

setup(**params)

