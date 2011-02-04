#!/usr/bin/env python
import sys, os, argparse, re

# Insert one level above project directory to path for testing.

sys.path.insert(0, "../..")
from taxonomy.alignment import Alignment, PROFILE_TYPES, ALIGNMENT_DEFAULTS
#from Taxonomy.alignment import Alignment


def main():
    """
    Entry point for this script.
    """
    # Get command-line arguments.
    action, arguments = parse_arguments()

    # Simplify things by going from arguments object to a dict.
    # argparse puts positional arguments in a list, even when
    # nargs=1.  Clean this up for required positional args.
    # Remove keys that the constructor is not expecting.
    # Finally, set default search options based on selected profile, 
    # if necessary.
    align_args = arguments.__dict__
    align_args['reference_package'] = align_args['refpkg'][0]
    align_args['sequence_file'] = align_args['seqfile'][0]
    del align_args['subparser_name'], align_args['seqfile'], \
        align_args['refpkg']
    if 'alignment_options' in align_args and not align_args['alignment_options']:
        align_args['alignment_options'] = \
        ALIGNMENT_DEFAULTS['alignment_options'][align_args['profile_version']]
    if 'search_options' in align_args and not align_args['search_options']:
        align_args['search_options'] = \
        ALIGNMENT_DEFAULTS['search_options'][align_args['profile_version']]

    if align_args['verbose']: print 'Argument parsing complete.'

    if align_args['debug']: print align_args

    align = Alignment(**align_args)
    if align_args['verbose']: print 'Alignment object instance created.'

    # Only perform hmmer_search for the searchalign action when a hmmer
    # profile is specified.
    if action == 'searchalign' and 'hmmer' in align_args['profile_version']:
        hmmsearch_output_file = align.hmmer_search()
        # Create alignment with hmmer for the file containing 
        # recruited sequences.  hmmsearch output file is in .sto format.
        align.hmmer_align(sequence_file=hmmsearch_output_file,
                          sequence_file_format='stockholm',
                         )

    elif action == 'align' and 'hmmer' in align_args['profile_version']:
        # Create alignment with hmmer. 
        align.hmmer_align()

    elif action == 'align' and 'infernal1' in align_args['profile_version']:
        # Create alignment with hmmer. 
        align.infernal_align()


    else:
        raise Exception, 'Not implemented exception raised. ' + \
                          'action: %s, profile: %s' % \
                          (action, align_args['profile_version'])


def parse_arguments(action_arguments=None):
    """
    Parse command-line arguments.
    """
    # Build up a list of search and align defaults, for help text used below.
    search_defaults = ''
    for profile in ALIGNMENT_DEFAULTS['search_options']:
        options =  ALIGNMENT_DEFAULTS['search_options'][profile]
        search_defaults += '(' + profile + ': "' + \
        options + '")' 

    align_defaults = ''
    for profile in ALIGNMENT_DEFAULTS['alignment_options']:
        options =  ALIGNMENT_DEFAULTS['alignment_options'][profile]
        align_defaults += '(' + profile + ': "' + \
        options + '")' 


    parser = argparse.ArgumentParser(
                        description='pkg_align.py - Wrapper script for ' + \
                                    'Infernal and hmmer alignment ' + \
                                    'binaries.')

    # Setup sub-commands
    subparsers = parser.add_subparsers(dest='subparser_name')
    # Help
    parser_help = subparsers.add_parser('help', help='Help for actions')
    parser_help.add_argument('action', nargs=1)

    # align
    parser_align = subparsers.add_parser('align', help='Create an alignment')

    # searchalign
    parser_searchalign = subparsers.add_parser('searchalign', 
                help='Search sequences and create an alignment')

    parser_searchalign.add_argument('--search-opts', dest='search_options', metavar='OPTS',
                        help='search options, such as "--notextw --noali ' + \
                        '-E 1e-2"' + " Defaults are as follows for the ' + \
                        different profiles: " + \
                        search_defaults)

    # With the exception of 'help', all subcommands share a certain 
    # number of arguments, which are added here.
    for subcommand in subparsers.choices.keys():
        if subcommand == 'help': continue
        subparser = subparsers.choices[subcommand]

        subparser.add_argument('--align-opts', dest='alignment_options', metavar='OPTS',
                            help='Alignment options, such as "--mapali $aln_sto".' + \
                            ' $ characters will need to be escaped if ' + \
                            'using template variables. ' + \
                            'Available template variables are ' + \
                            '$aln_sto, $profile. ' + \
                            " Defaults are as follows for the different ' + \
                            profiles: " + align_defaults)
        subparser.add_argument('--profile-version', dest='profile_version',
                            default=ALIGNMENT_DEFAULTS['profile'],
                            choices=PROFILE_TYPES,
                            help='Profile version to use. Defaults to ' + \
                            ALIGNMENT_DEFAULTS['profile'])                        
        subparser.add_argument('-o', '--outprefix', required=True, dest='out_prefix', 
                                  help='Output file prefix. ' + 'Currently only works ' + \
                                  'with a single sequence file')
        subparser.add_argument('--min-length', dest='min_length', type=int, metavar='N', 
                            default=ALIGNMENT_DEFAULTS['min_length'],
                            help='Minimum sequence length. Defaults to ' + \
                            str(ALIGNMENT_DEFAULTS['min_length']) + '.' + \
                            ' Currently only works with the hmmer profile.')
        subparser.add_argument('--format', dest='sequence_file_format', 
                               default=ALIGNMENT_DEFAULTS['sequence_file_format'],
                               help='Specify format of seqfile.  Defaults ' + \
                               ALIGNMENT_DEFAULTS['sequence_file_format'])
        subparser.add_argument('refpkg', nargs=1, 
                               type=reference_package, help='Reference package directory')
        subparser.add_argument('seqfile', nargs=1,
                               help='A single fasta file')
        subparser.add_argument('--debug', action='store_true', help='Enable debug output')
        subparser.add_argument('--verbose', action='store_true', help='Enable verbose output')

    arguments = parser.parse_args()
    action = arguments.subparser_name
    # Override arguments passed at the command-line to support help <action>
    if action_arguments:
        arguments = parser.parse_args(action_arguments)
    else:
        arguments = parser.parse_args()

    if (action == 'help'):
        parse_arguments(action_arguments=[str(arguments.action[0]), '-h'])

    return action, arguments


def reference_package(reference_package):
    """
    A custom argparse 'type' to make sure the path to the reference package exists.
    """
    if os.path.isdir(reference_package):
        # Remove trailing / from directory if found.
        if reference_package.endswith('/'):
            reference_package = reference_package[:-1]
        return reference_package
    else:
        raise Exception, 'Path to reference package does not exist: ' + reference_package


if __name__ == '__main__':
    sys.exit(main())
