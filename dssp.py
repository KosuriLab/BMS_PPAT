"""Reads output from DSSP program.

DSSP is a program for classifying secondary structures and solvent
accessibility from protein structures. This module contains functions
that process the output.

These functions have been tested on output from the DSSP webserver
(http://swift.cmbi.ru.nl/gv/dssp/) using the version available
as of May-24-2013.

This module was written by Jesse Bloom.


List of functions
----------------------

* *MaxASA* : Maximum accessible surface area for amino acids.

* *ReadDSSP* : reads secondary structures and RSA from DSSP output.


Documentation for individual functions
----------------------------------------
Documentation for individual functions is provided in their
definitions below.

"""

import os


def MaxASA(scale):
    """Returns the maximum accessible surface area for amino acids.

    This function returns the maximum accessible surface area (ASA)
    for the amino acids in units of square angstroms. These ASAs
    are necessary for calculating the relative solvent accessibility
    of a residue.

    There are a large variety of estimates for the exact ASA for each
    amino acid. The calling variable *scale* specifies which scale to
    use. Allowed values are the following strings:

        * 'Tien2013' : The values provided by Tien et al (Maximum
          allowed solvent accessibilities of residues in proteins),
          as defined in Table 1 in the column labeled "Theoretical"
          at http://www.plosone.org/article/info:doi/10.1371/journal.pone.0080635

    The returned variable is a dictionary *asa* keyed by each 
    upper-case one-letter amino-acid code, with the values
    being the ASA for that residue.

    Example:

    >>> asa = MaxASA('Tien2013')
    >>> len(asa) == 20
    True
    >>> asa == {'A':129.0, 'R':274.0, 'N':195.0, 'D':193.0, 'C':167.0, 'E':223.0, 'Q':225.0, 'G':104.0, 'H':224.0, 'I':197.0, 'L':201.0, 'K':236.0, 'M':224.0, 'F':240.0, 'P':159.0, 'S':155.0, 'T':172.0, 'W':285.0, 'Y':263.0, 'V':174.0}
    True

    """
    if scale == 'Tien2013':
        return {'A':129.0, 'R':274.0, 'N':195.0, 'D':193.0, 'C':167.0,
                'E':223.0, 'Q':225.0, 'G':104.0, 'H':224.0, 'I':197.0,
                'L':201.0, 'K':236.0, 'M':224.0, 'F':240.0, 'P':159.0,
                'S':155.0, 'T':172.0, 'W':285.0, 'Y':263.0, 'V':174.0,}
    else:
        raise ValueError("Invalid value for scale")


def ReadDSSP(infile, asa_scale, chain=None):
    """Reads secondary structure and ASA from DSSP output.

    DSSP is a program that calculates secondary structure and
    absolute solvent accessibility (ASA) on a per-residue
    basis from a PDB protein structure file.

    These functions have been tested on output from the DSSP webserver
    (http://swift.cmbi.ru.nl/gv/dssp/) using the version available
    as of May-24-2013.

    **NOTE** that this function currently does not process DSSP output
    from PDB structures where the residue numbers have letter suffixes
    (such as residue 25A). It will raise an error if such residue
    numbers exist in the DSSP output.

    CALLLING VARIABLES:

        * *infile* : name of a text file containing the output of a 
          DSSP analysis.

        * *asa_scale* : calculation of relative solvent accessiblity (RSA)
          from the absolute solvent accessibility (ASA) returned by DSSP
          requires normalization by the maximum solvent accessibility
          of each residue. *asa_scale* specifies those maximum values.
          It should be a string that is passable to *MaxASA* as
          its single calling paramter -- those maximum ASAs are used
          to normalize the DSSP ASAs to RSAs.

        * *chain* is an optional variable. If the PDB structure used to 
          run DSSP has only one chain, then this variable can just have
          its default value of *None*. Otherwise it should provide a
          letter specifying the protein chain from which we read
          the ASA values. If there are multiple chains, this variable
          must be set to a value matching one of those chains.

    RETURN VARIABLE:

    The return variable is a dictionary *dssp*. For each residue
    number in the original PDB analyzed by DSSP, there is an
    integer key in *dssp*. The value for each of these residues
    is in turn a dictionary which has the following string keys:

        * *ASA* : the value is absolute solvent accessibility in
          square angstroms.

        * *RSA* : the value is the relative solvent accessiblity.

        * *SS* : the DSSP secondary structure code, which can be
          the following letters:

            - *G* : 3-10 helix

            - *H* : alpha helix

            - *I* : pi helix

            - *B* : beta bridge

            - *E* : beta bulge

            - *T* : turn

            - *S* : high curvature

            - *_* : loop

        * *SS_CLASS* : the larger secondary structure classes, which 
          are defined by the following strings based on *SS*:

            - *helix* : a *SS* value of *G*, *H*, or *I*

            - *strand* : a *SS* value of *B* or *E*

            - *loop* : any of the other *SS* values.
    """
    maxasa = MaxASA(asa_scale)
    if not os.path.isfile(infile):
        raise IOError("Cannot find infile of %s" % infile)
    lines = open(infile).readlines()
    iline = 0
    while iline < len(lines) and '  #  RESIDUE AA STRUCTURE' not in lines[iline]:
        iline += 1
    if iline == len(lines):
        raise ValueError("Problem parsing DSSP infile %s, never found header" % infile)
    dssp = {}
    chainsfound = {}
    for line in lines[iline + 1 : ]:
        try:
            ires = int(line[5 : 10])
        except ValueError:
            if line[5 : 10].strip() == '':
                # probable missing residue
                if line[13] == '!':
                    continue # definitely missing residue, skip it
            raise ValueError("Cannot convert into a residue number, are they letter suffixes on some of the residue numbers? The problem is with: %s\nFound in line\n%s" % (line[5 : 10], line))
        ichain = line[11]
        if chain:
            if ichain != chain:
                continue # not the chain we are looking at
        else:
            chainsfound[ichain] = True
            if len(chainsfound) > 1:
                raise ValueError("No chain specified but multiple chains found in %s" % infile)
        if ires in dssp:
            raise ValueError("Duplicate residue number of %d" % ires)
        wt = line[13]
        dssp[ires] = {}
        ss = line[16]
        if ss == ' ':
            ss = '_'
        try:
            asa = float(line[35 : 38])
        except TypeError:
            raise ValueError("Cannot process ASA: %s" % line[35 : 38])
        dssp[ires]['ASA'] = asa
        dssp[ires]['RSA'] = asa / float(maxasa[wt])
        dssp[ires]['SS'] = ss
        if ss in ['G', 'H', 'I']:
            dssp[ires]['SS_CLASS'] = 'helix'
        elif ss in ['B', 'E']:
            dssp[ires]['SS_CLASS'] = 'strand'
        elif ss in ['T', 'S', '_']:
            dssp[ires]['SS_CLASS'] = 'loop'
        else:
            raise ValueError("Unrecognized secondary structure of %s" % ss)
    if not dssp:
        raise ValueError("Failed to read any data from %s. Did you specify the right value for chain?" % (infile))
    return dssp



if __name__ == '__main__':
    import doctest
    doctest.testmod()
