# csp.py
#
# S. Ruedisser April 2017
# colors atoms in pymol according to chemical shift differences between two states: e.g. free and bound
#
# compute chemical shift difference from two .prot files and color residues accordingly
#
# usage:
# ------
# within pymol:
# run csp.py
# csp pdbfile.pdb, shifts_apo.prot, shifts_complex.prog, nr_colors, offset, low_limit 
#
##############################################################################################################################################################################
#
# Copyright (c) 2017
#
# ETH Zuerich, Simon Ruedisser
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation 
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without 
# specific prior written permission.
# 
# 
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
################################################################################################################################################################################

cmd.reinitialize()

import re
comment = re.compile('^\s*$|^\s*#')
from operator import add
from operator import sub

# dictionary containing residue_nr, atomtypes and chemical shifts
# {resn_atomtype: chemical shifts}

prot =  {}
prot_lig = {}


# main function 

def csp(pdbfile, shifts_apo, shifts_complex, nr_colors, offset, low_limit=None):
    
    """
    colors residues according to chemical shift difference
    chemical shifts are provided in the files shifts_apo and shifts_complex
    nr_colors defines the number of colors
    low_limit is the lower threshold given as muliples of the standard deviations
    if no low_limit is given, the low_limit is computed automatically
    """
    
    print("pdb file: " + pdbfile + "\n" + "--------")
    cmd.load(pdbfile)

    
    # read chemical shifts for apo protein
    apo_shifts = open(shifts_apo, 'rU')
    for line in apo_shifts:
        
        # igone comment lines
        if not comment.match(line):
            split = line.split()
            resn_atom = split[4] + "_" + split[3]
            shift = split[1]
            
            prot[resn_atom] = shift

    # read chemical shifts for the complex       
    complex_shifts = open(shifts_complex, 'rU')
    for line in complex_shifts:
        
        # igone comment lines
        if not comment.match(line):
            split = line.split()
            resn_atom = split[4] + "_" + split[3]
            shift = split[1]
            
            prot_lig[resn_atom] = shift

    # compute shift difference for NH
    
    # retrieve keys
    prot_resn_atom = list(prot.keys())
    prot_lig_resn_atom = list(prot_lig.keys())

    # print first 10 keys 
    # print(prot_resn_atom[:9])

    # match keys containing _N (amide N)
    N = re.compile(".*_N")
    N_prot = list(filter(N.match, prot_resn_atom))
    N_prot_lig = list(filter(N.match, prot_lig_resn_atom))

    # select N with shift given for the apo and complex state
    N_intersect = set(N_prot).intersection(set(N_prot_lig))

    # helper function to obtain resn from string resn_atom
    def getresnr(resn_atom):
        return resn_atom[0:resn_atom.find("_")]

    # intersection of residue numbers, for intersection with other nuclei
    # N_intersect = set(map(getresnr, N_prot)).intersection(set(map(getresnr, N_prot_lig)))
    # print(N_intersect)

    # match keys containing _H (amide H)
    H = re.compile(".*_N")
    H_prot = list(filter(N.match, prot_resn_atom))
    H_prot_lig = list(filter(N.match, prot_lig_resn_atom))

    # select H with shift given for the apo and complex state
    H_intersect = set(H_prot).intersection(set(H_prot_lig))
     
    NH_csp = {}

    # shift difference as the euclidian distance between 15N, scaled by 0.2, and 1H
    
    for atom in N_intersect:
        resn = getresnr(atom)
        atom_N = resn + "_N"
        atom_H = resn + "_H"
        shift_difference = (
            (0.2 * (float(prot[atom_N]) - float(prot_lig[atom_N]))**2 +
            (float(prot[atom_H]) - float(prot_lig[atom_H]))**2)**(0.5))
        
        NH_csp[atom_N] = shift_difference

    # print(NH_csp)

    # compute average csp
    NH_csps = list(map(float, NH_csp.values()))
    NH_avg = sum(NH_csps)/len(NH_csps)

    # compute standard deviation for NH csp
    NH_stdev = (
        sum([(NH_avg - NH_csps[i])**2 for i in range(len(NH_csps))])**0.5 /
        len(NH_csps))

    # find minimum and maximum of NH_csps
    csp_max = max(map(abs, NH_csps))
    csp_min = min(map(abs, NH_csps))
    
    print(("average CSP: " + "%.4f" % NH_avg +
           "\nstandard deviation: " + "%.4f" % NH_stdev +
           "\nnumber of residues: " + "%i" % len(NH_csps) +
           "\nmin(CSP): " + "%.4f" % csp_min +
           "\nmax(CSP): " + "%.4f" % csp_max))

    # estimate standard deviation of CSP for non-interacting residues by removing residues with a CSP > 3 * stdev
    # recalcuate stdev and repeat procedure until the stdev stays constant
    
    if low_limit == None:

        print("--------\n" + "removing residues with CSP > 3 * stdev")
        stdev = NH_stdev
        csp_current = {}
        stdev_previous = 0

        while abs(stdev - stdev_previous) > 0.02:
            csp_current = {}
            stdev_limit = 3 * stdev
            
            for residue in NH_csp.keys():
                if NH_csp[residue] < stdev_limit:
                    csp_current[residue] =  NH_csp[residue]
                    
            stdev_previous = stdev    
            csps = list(map(float, csp_current.values()))

            if len(csps) != 0:
                avg = sum(csps)/len(csps)
                stdev = (
                    sum([(avg - csps[i])**2 for i in range(len(csps))])**0.5 /
                    len(csps))
            else:
                print("Number of remaining residues is zero!\n Try setting last parameter (low_limit) to 1")
                break
                
            print("stdev: " + str(stdev))
            print("Number of residues: " + str(len(list(csp_current.keys()))))
        low_limit = stdev
        

    # calculate rgb values from shifts
    # normalized to the maximal CSP
        
    def shift2rgb(shifts, low_limit, offset):
        if float(shifts) < (low_limit * NH_stdev):
            return 0.0
        else:
            return (float(shifts) - low_limit * NH_stdev)/(csp_max - csp_min - low_limit * NH_stdev) + offset

    # define n different colors as rgb

    def rgbcolors(nr_colors, offset, low_limit):
        """
        usage: rgbcolors(nr_colors, offset, low_limit)
        nr_colors: number of color-intervals
        offset: offset to zero on rgb channel
        low_limit: lower limit for csp considered in times standard deviation
        """
        
        # create dictionary with nr_colors {color_0:[0, 0, 0], ..., color_n:[1, 0, 0]},
        # these are dummy colors whose rgb values are set later
        
        colors = {}
        for col in range(nr_colors):
            colname = "color_" + str(col)
            
            # use r channel as color indicator
            color = [float(col) / (float(nr_colors) - 1), 0, 0]
            colors[colname] = color
            cmd.set_color("%s"%colname, "%s"%colors[colname])
            
        # list of all r-values
        r_values = sorted([item[0] for item in colors.values()])

        # dictionary with occurences for each color

        occurences = {}
        for col in range(nr_colors):
            occurences[col] = 0
        
        # set color for each residue with csp
        
        for residue in NH_csp.keys():
            residue_nr = residue[0:residue.find("_")]
            r_current = shift2rgb(NH_csp[residue], low_limit, offset)
            
            # list with distances to r_values
            distances = [abs(r_current - r) for r in r_values]

            # set color index to index of color with closest distance to current csp
            index = distances.index(min(distances))
            
            colors_key = list(sorted(colors.keys()))[index]

            # update counts

            occurences.update({index:occurences[index] + 1})
            
            # printing diagnostics
            if False:
                print("for csp > 0.8: min distance: " + str(min(distances)) + " r_current: " + str(r_current))
                print("distances: " + str(distances))
                print(("CSP for residue " + str(residue) + ": " + str(r_current) +
                       " index: " + str(index) + " gives color " + str(colors_key)))
                
            cmd.color("%s"%(colors_key), "resi %s"%residue_nr)

        print("counts: " + str(occurences) + "\n------------")

        # print a histogramm
        for i in list(occurences.keys()):
            bar = "#" * occurences[i]
            print(str(i) + " : " + str(bar))
                  
        # set rgb values for colors
        # rgp = [1, 1, 1] : white
        # rgb = [1, 0, 0] : red
        # rgb = [0, 0, 1] : blue

        # update colors
        for col in range(nr_colors):
            colname = "color_" + str(col)
            fraction = colors[colname][0]
            # colors[colname] = list(map(sub, [1, 1, 1], color))
            # colors[colname] = list(map(add, list(map(lambda x: 1 - x, rgbmap)), color))
            
            colors[colname] = [ 1.0, 1-fraction, 1-fraction]
            cmd.set_color("%s"%colname, "%s"%colors[colname])

        # print("colors: " + str(colors) + "\n------")
        
        return None

    rgbcolors(int(nr_colors), float(offset), float(low_limit))
    # cmd.show("sticks")
    # cmd.hide("sticks", "hydrogens")
    cmd.hide("lines", "all")
    cmd.show("cartoon")
    
# add function to pymol commands
cmd.extend('csp', csp)
