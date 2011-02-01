#!/usr/bin/env python
import json, sys, os, string, subprocess, re
from string import Template
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq, SeqRecord
from itertools import islice

# Supported profiles
PROFILE_TYPES = ( 'hmmer3', 'infernal1', 'infernal1mpi', )

# Default options that can be used by scripts.
# Keys for search_options and alignment options must map to a valid profile.
ALIGNMENT_DEFAULTS = {
                        'min_length' : 1,
                        'profile' : 'hmmer3',
                        'search_options' : { 'hmmer3' : '--notextw --noali',
                                             'infernal1' : '',
                                             'infernal1mpi' : '',
                                           },
                        'alignment_options' : { 'hmmer3' : '--mapali $aln_sto', 
                                                'infernal1' : '',   
                                                'infernal1mpi' : '',   
                                              },
                        'sequence_file_format' : 'fasta',
                     }


class Alignment(object):
    """
    A class to provide alignment-related tools for reference packages.

    Notes:
    A "mask" is just a boolean list, with True meaning include.
    """

    def __init__(self, reference_package, out_prefix, profile_version, 
                 sequence_file, sequence_file_format, alignment_options,
                 search_options=None, min_length=None, debug=False, 
                 verbose=False):
        """
        Constructor - sets up a number of properties when instantiated.
        """
        self.debug = debug
        self.verbose = verbose

        self.reference_package = reference_package        
        # Determine the name of the reference package, excluding any other 
        # elements in its path.
        self.reference_package_name = os.path.split(reference_package)[-1]
        self.reference_package_name_prefix = os.path.splitext(self.reference_package_name)[0]

        # Keep track of whether or not an out_prefix was specifed.
        self.out_prefix_arg = out_prefix
        self.out_prefix = out_prefix
        # Set default prefix if unspecified.
        if out_prefix is None:
            # Sequence file names will be appended within for loops.
            out_prefix = os.path.join(self.reference_package, self.reference_package_name_prefix) + '.'

        if min_length:
            self.min_length = min_length

        self.sequence_file = sequence_file
        self.search_options = search_options
        self.alignment_options = alignment_options
        self.sequence_file_format = sequence_file_format

        # Read in CONTENTS.json and retrieve settings we will be working with.
        json_file = os.path.join(reference_package, 'CONTENTS.json')
        json_contents = self._read_contents_json(json_file)
        aln_sto, aln_fasta, profile = [reference_package] * 3
        self.aln_sto = os.path.join(aln_sto, json_contents['files']['aln_sto'])
        self.aln_fasta = os.path.join(aln_fasta, json_contents['files']['aln_fasta'])
        self.profile = os.path.join(profile, json_contents['files']['profile'])

        self.profile_version = profile_version
        # There is no reason to go beyond this point if the profile is not valid.
        if not self.match_model():
            raise Exception, 'profile: ' + self.profile + ' does not appear to be valid.' + \
                             ' Version match specified: ' + profile_version
 
        # read in the consensus RF line
        self.consensus_list = self._consensus_list_of_sto(self.aln_sto)

        # debugging
        # print ("".join(map(lambda(b): str(int(b)), self.consensus_list)))

        # initialize the masking
        if 'mask' in json_contents['files']:
            self.mask_file = os.path.join(self.reference_package, json_contents['files']['mask'])
   
        sto_len = self._get_sequence_length(self.aln_sto, "stockholm")

        self.trimal_mask = self._mask_of_file(self.mask_file, sto_len)
        # first make sure that the trimal mask only includes consensus columns according to HMMER
        for pos in range(sto_len):
            if self.trimal_mask[pos] & (not self.consensus_list[pos]):
                print("trying to include a non-consensus column %d in the mask" % pos)
                raise Exception, "trying to include a (non-consensus column " + \
                                     str(pos) + " in the mask"

        # Now we make consensus_only_mask, which is the mask after we have taken just the consensus columns.
        self.consensus_only_mask = self._mask_list(mask=self.consensus_list, to_mask=self.trimal_mask)


            
    # Public methods

    def hmmer_search(self):
        """
        Recruit fragments using hmmsearch.  Works with a single sequence file, further 
        work would be required if it is to be expanded to work with multiple sequence files.
        """
        # hmmsearch must be in PATH for this to work.
        if self.verbose: print 'Entering hmmer_search()'
        hmmsearch_output_file = self.out_prefix + '.search_out.sto'
        hmmsearch_command = 'hmmsearch ' + self.search_options + \
                            ' -A ' + hmmsearch_output_file + \
                            ' ' + self.profile + ' ' + self.sequence_file

        if self.debug: print "Command to execute: \n\t" + hmmsearch_command

        child = subprocess.Popen(hmmsearch_command,
                                 stdin=None,
                                 stdout=None,
                                 stderr=None,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()

        # If return code was not 1, hmmsearch completed without errors.
        if not return_code:
            return hmmsearch_output_file
        else:
            raise Exception, "hmmsearch command failed: \n\t" + hmmsearch_command
        if self.verbose: print 'Leaving hmmer_search()'


    def hmmer_align(self, sequence_file=None, sequence_file_format=None):
        # Note to Erick, ref=False was unused, as was frag=False. They have
        # been removed.
        """
        Create an alignment with hmmalign. Then, separate out reference sequences 
        from the fragments into two separate files. Note that mask=True forces 
        a consensus action.
        """
        if self.verbose: print 'Entering hmmer_align()'
        # If sequence_file not specified, we'll just work with
        # self.sequence_file.
        if not sequence_file:
            sequence_file = self.sequence_file
        if not sequence_file_format:
            sequence_file_format = self.sequence_file_format

        # hmmalign must be in PATH for this to work.
        hmmer_template = Template('hmmalign -o $together_aln ' + \
                                   self.alignment_options + ' ' + \
                                   '$profile $sequence_file')

        if self.debug: print "hmmer_template: \n    " + hmmer_template.template
        
        together_aln = self.out_prefix + '.align_out.sto'
        _,sequence_file_name = os.path.split(sequence_file)

        hmmalign_command = hmmer_template.substitute(sequence_file=sequence_file,
                                                     aln_sto=self.aln_sto,
                                                     profile=self.profile,
                                                     together_aln=together_aln,
                                                    )
        if self.debug: print "Command to execute: \n    " + hmmalign_command

        try:
            child = subprocess.Popen(hmmalign_command,
                                     stdin=None,
                                     stdout=None,
                                     stderr=None,
                                     shell=(sys.platform!="win32"))
            return_code = child.wait()

        
            # If return code was not 1, split off the reference sequences from the fragments.
            if not return_code:

                frag_names = self._names(SeqIO.parse(sequence_file, sequence_file_format))

                # pull up the together_aln then mask it
                def make_masked_iterator():
                    full_aln = SeqIO.parse(together_aln, "stockholm")
                    aln_consensus_list = self._consensus_list_of_sto(together_aln)
                    return(self._min_length_filter(
                                        self.min_length,
                                        self._maskerator(
                                            # mask using the consensus_only mask
                                            self.consensus_only_mask, 
                                            # after taking only the alignment columns
                                            self._maskerator(aln_consensus_list, full_aln)
                                                        )
                                                  )
                          )
            
                SeqIO.write(self._id_filter(make_masked_iterator(), lambda(idstr): idstr in frag_names), 
                            self.out_prefix + '.masked.fasta', "fasta")

                SeqIO.write(self._id_filter(make_masked_iterator(), lambda(idstr): idstr not in frag_names), 
                            self.out_prefix + '.refs.fasta', "fasta")

        except:
                raise
            # we may want to tidy things up with an option later, but for now we just leave everything around.
            # finally:
            #     # Always remove the temporary alignment file
            #     os.remove(tmp_file)
        if self.verbose: print 'Leaving hmmer_align()'


    def match_model(self):
        """
        Verify that the header contents look to be from a specific profile
        version.
        """
        header = self._extract_header()

        is_match = False
        if self.profile_version.startswith('infernal1'):
            # Match a couple of things one would expect to find in a .cm file.
            cm_header = re.compile(r"^INFERNAL-1" + r".*MODEL", re.MULTILINE|re.DOTALL)
            is_match = bool(cm_header.search(header))
        if self.profile_version == 'hmmer3':
            # Match a couple of things one would expect to find in a .hmm file.
            hmm_header = re.compile(r"^HMMER3" + r".*CKSUM", re.MULTILINE|re.DOTALL)
            is_match = bool(hmm_header.search(header))

        return is_match


    # Private methods

    def _read_contents_json(self, json_file):
        with open(json_file, 'r') as contents:
            json_contents = json.load(contents)
        return json_contents 


    def _names(self, records):
        """
        Get the sequence names.
        """
        s = set()
        for record in records:
            s.add(record.id)
        return(s)


    def _id_filter(self, in_seqs, f):
        """
        Generator function to filter out sequences by id.
        """
        for record in in_seqs:
            if (f(record.id)):
                yield record

    def _min_length_filter(self, min_length, records):
        """
        Generator function to filter out sequences by sequence length.
        """
        for record in records:
            nongaps = len(filter(lambda(s): s != "-", list(str(record.seq))))
            if nongaps >= min_length:
                yield record
            else:
                print record.id + " is too short"  


    def _extract_header(self):
        """
        Extract the header portion of the file, approximately.
        """
        line_count = 20

        # Only read in the first <line_count> lines to extract the header,
        # approximately.
        with open(self.profile, 'r') as model:
            header = ''.join(list(islice(model, line_count)))
        
        return header


    # Begin mask-related functions
    # Brian: does it make sense to have these be part of the object if none of them refer to self?
    # the answer to me seems "no" but...

    def _mask_list(self, mask, to_mask):
        """
        Mask a list!
	    """
        assert(len(mask) == len(to_mask))
        masked = []
        for i in range(len(mask)):
            if mask[i]:
                masked.append(to_mask[i])

        return(masked)


    def _maskerator(self, mask, records):
        """
	Apply a mask to all sequences in records.
        """
        for record in records:
            sequence = list(str(record.seq))
            yield SeqRecord(Seq(''.join(self._mask_list(mask, sequence))),
                            id=record.id, description=record.description)

    def _int_list_of_file(self, file_name):
        """
        Get a comma-delimited integer list from a file.
        """
        # Regex to remove whitespace from mask file.
        whitespace = re.compile(r'\s|\n', re.MULTILINE)
        with open(file_name, 'r') as handle:
            mask_text = handle.read() 
            mask_text = re.sub(whitespace, '', mask_text)

        # Cast mask positions to integers
        return(map(int, mask_text.split(',')))

    def _mask_of_file(self, mask_file, length):
        """
        Make a mask of a zero-indexed comma-delimited integer list in a file.
        The included indices are turned to True in the mask.
        Will fail if mask is out of range, and that's a good thing.
        """
        mask = [False] * length
        for i in self._int_list_of_file(mask_file):
            mask[i] = True
        return(mask)

    # End mask-related functions


    # Alignment-validation-related functions

    def _get_sequence_length(self, source_file, source_file_type):
        """
        Returns the length of the first sequence in an alignment.  If file is empty 
        or cannot be read in by AlignIO.read, an exception will be thrown.
        """
        alignment = AlignIO.read(source_file, source_file_type)
        return len(alignment.__getitem__(0))


    def _sequence_length_check(self, records, reference_length):
        """
        Generator to validate sequence lengths, as we go through the sequences.
        """
        for record in records:
            if reference_length == len(record):
                yield record
            else:
                raise Exception, 'Sequence "' + record.id + '" has length of ' + str(len(record)) + \
                                 ', expected ' + str(reference_length)

    # End alignment-validation-related functions

    # consensus column related functions

    # I wish biopython did this...
    def _consensus_list_of_sto(self, sto_aln):
        """
        Return a boolean list indicating if the given column is consensus according to the GC RF line.
        """
        rf_rex = re.compile("#=GC RF\s+([x.]*)")
        rf_list = []

        def is_consensus(c):
            if c == 'x':
                return(True)
            if c == '.':
                return(False)
            assert(False)

        with open(sto_aln, 'r') as handle:
            for line in handle.readlines():
                m = rf_rex.match(line)
                if m:
                    rf_list.append(m.group(1))

        return(map(is_consensus, list("".join(rf_list))))


        # End consensus column related functions


