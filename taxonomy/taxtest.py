#!/usr/bin/env python

"""\
taxomatic.py
============

Usage: %prog action <options>

Creation, validation, and modification of reference packages for use
with `pplacer` and related software.

Actions
=======

Use `taxomatic.py -h` or `taxomatic.py --help` to print help text.

help
  Print detailed help for an action below using `taxomatic.py help <action>`
create
  Create a reference package
check
  Verify that reference package is intact and valid.

"""

import argparse
import sys
import logging
from taxonomy.package import manifest_name, package_contents, write_config
from taxonomy import __version__ as version

log = logging
PROG = 'taxtest.py'
DESCRIPTION = 'To Be Named -- Creation, validation, and modification of ' + \
              'reference packages for use with `pplacer` and related software.'

def main():
    """
    """
    action, arguments = parse_arguments()
    verbose = arguments.verbose


def parse_arguments(action_arguments=None):
    """
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description=DESCRIPTION, prog=PROG)

    parser.add_argument('-V', '--version', action='version',
            version='%(prog)s v' + version,
            help='Print the version number and exit')

    ##########################
    # Setup all sub-commands #
    ########################## 

    subparsers = parser.add_subparsers(dest='subparser_name')
    
    # Begin help sub-command 
    parser_help = subparsers.add_parser('help', help='Detailed help for ' + \
        'actions using `help <action>`')
    parser_help.add_argument('action', nargs=1)
    # End help sub-command

    # Begin create sub-command
    parser_create = subparsers.add_parser('create', 
        help='Create a reference package')

    parser_create.add_argument("-f", "--aln-fasta",
        action="store", dest="aln_fasta",
        help='Multiple alignment in fasta format', metavar='FILE')

    parser_create.add_argument("-i", "--seq-info",
        action="store", dest="seq_info",
        help='CSV format file describing the aligned reference ' + \
             'sequences, minimally containing the fields "seqname" ' + \
             'and "tax_id"', metavar='FILE')

    parser_create.add_argument("-m", "--mask",
        action="store", dest="mask",
        help='Text file containing a mask.', metavar='FILE')

    parser_create.add_argument("-p", "--profile",
        action="store", dest="profile",
        help='Alignment profile', metavar='FILE')

    parser_create.add_argument('-P', '--package-name', 
        action='store', dest='package_name',
        default='./taxonomy.refpkg', metavar='PATH',
        help='Name of output directory [default %(default)s]')

    parser_create.add_argument("-r", "--package-version",
        action="store", dest="package_version",
        help='Release version for the reference package', metavar='VERSION')

    parser_create.add_argument("-s", "--tree-stats",
        action="store", dest="tree_stats",
        help='File containing tree statistics (for example ' + \
             'RAxML_info.whatever").', metavar='FILE')

    parser_create.add_argument("-S", "--aln-sto",
        action="store", dest="aln_sto",
        help='Multiple alignment in Stockholm format.', metavar='FILE')

    parser_create.add_argument("-t", "--tree-file",
        action="store", dest="tree_file",
        help='Phylogenetic tree in newick format', 
        metavar='FILE')

    parser_create.add_argument("-T", "--taxonomy",
        action="store", dest="taxonomy",
        help='CSV format file defining the taxonomy. Fields include ' + \
             '"tax_id","parent_id","rank","tax_name" followed by a column ' + \
             'defining tax_id at each rank starting with root', metavar='FILE')
    # End create sub-command

    # Begin check sub-command
    parser_check = subparsers.add_parser('check', 
        help='The check action is not yet implemented')
        #help='Verify that a reference package is intact and valid')
    parser_check.add_argument('-P', '--package-name', 
        action='store', dest='package_name',
        default='./taxonomy.refpkg', metavar='PATH',
        help='Name of output directory. [default %(default)s]')
    # End check sub-command

    # Begin taxtable sub-command
    # TODO: Merged from taxtable.py and move from optparse to argparse.
    # End taxtable sub-command

    # With the exception of 'help', all subcommands can share a  
    # number of arguments, which are added here.
    for subcommand in subparsers.choices.keys():
        if subcommand == 'help': continue
        subparser = subparsers.choices[subcommand]
        subparser.add_argument('-v', '--verbose', action='count', dest='verbose', 
            help='Increase verbosity of screen output (eg, -v is verbose,' + \
                 '-vv more so', default=0)

    # Determine we have called ourself (e.g. "help <action>")
    # Set arguments to display help if parameter is set
    #           *or* 
    # Set arguments to perform an action with any specified options.
    arguments = parser.parse_args(action_arguments) if parser.parse_args(action_arguments) else parser.parse_args()
    # Determine which action is in play.
    action = arguments.subparser_name

    # Support help <action> by simply having this function call itself and 
    # translate the arguments into something that argparse can work with.
    if action == 'help':
        parse_arguments(action_arguments=[str(arguments.action[0]), '-h'])

    print arguments

    return action, arguments



#    loglevel = {
#        0:logging.WARNING,
#        1:logging.INFO,
#        2:logging.DEBUG
#        }.get(options.verbose, logging.DEBUG)
#
#    verbose_format = '%(levelname)s %(module)s %(lineno)s %(message)s'
#
#    logformat = {0:'%(message)s',
#        1:verbose_format,
#        2:verbose_format}.get(options.verbose, verbose_format)
#
#    # set up logging
#    logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)
#
#        try:
#            if hasattr(Taxonomy.package, action):
#                getattr(Taxonomy.package, action)(options)
#            else:
#                log.error('Sorry: the %s action is not yet implemented' % action)
#                sys.exit(1)
#        except OSError:
#            log.error('A package named "%s" already exists' % options.package_name)
#            sys.exit(2)



if __name__ == '__main__':
    sys.exit(main())
