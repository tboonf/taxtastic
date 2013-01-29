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

import collections
import logging
from itertools import combinations
import subprocess
import tempfile

import numpy

log = logging.getLogger(__name__)

def read_dists(fobj):
    """
    Read interleaved phylip distance matrix from file-like object fobj
    into a numpy matrix. Return (taxon_names, matrix).
    """

    # TODO: there has to be a more efficient way to read a numpy matrix
    elements = fobj.read().split()
    N = int(elements.pop(0))
    distmat = numpy.repeat(-1*numpy.inf, N**2)
    distmat.shape = (N, N)

    taxa = []
    for row, i in enumerate(range(0, len(elements), N+1)):
        taxa.append(elements[i])
        distmat[row,:] = [float(x) for x in elements[i+1:i+N+1]]

    return taxa, distmat

def fasttree_dists(fasta):
    """
    Calculate pairwise distances among DNA multiple alignment in
    `fasta` using FastTree and return (taxon_names, matrix).
    """

    # TODO: need a more informative error if FastTree is not installed.

    cmd = ['FastTree','-nt','-makematrix', fasta]

    with tempfile.TemporaryFile('rw') as stdout:
        proc = subprocess.Popen(cmd, stdout = stdout)
        proc.communicate()
        stdout.seek(0)
        taxa, distmat = read_dists(stdout)

    return taxa, distmat

def outliers(mat, cutoff):
    """
    Given pairwise distance matrix `mat`, identify elements with a
    distance to the centrid element of > cutoff. Returns a boolean
    vector corresponding to the margin of mat.
    """

    # index of most central element.
    medoid = numpy.argmin(numpy.median(mat, 0))

    # distance from each element to most central element
    dists = mat[medoid,:]

    return dists > cutoff
