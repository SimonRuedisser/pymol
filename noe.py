#########################################################################
# Script to display upper distance limits (from Cyana final.upl) as lines
# on pdb output from cyana (finalQ.pdb)
# #######################################################################
#
# in pymol type:
# run noe.py
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

# from __future__ import print_function

from pymol import cmd, CmdException
import re

# restraints are placed in dictionary with the following format
# restraints = {R1:[DIS, A1, A2], R2:[...]}

filename = "./final.upl"
pdb_file = "./finalQ.pdb"
strong = 3.0
medium = 4.0


cmd.reinitialize()
cmd.load(pdb_file)

spaces = re.compile("  *")

RID = 0
restraints = {}

for line in open(filename):
    splitline = spaces.split(line)
    print(splitline[3], splitline[6], splitline[7])
    restraints[RID] = [splitline[3], splitline[6], splitline[7]]
    RID+=1

print(restraints)

atom1="/////1HD"
atom2="/////2HD"
pre="/////"

def draw(atom1, atom2, label):
    cmd.distance(label, atom1, atom2, label=0, gap=0)
    return 1

def _color_restraint(object_name, color="yellow"):
    cmd.set("dash_color", color, object_name)

for restraint in restraints.keys():
    atom1 = pre+restraints[restraint][0]
    atom2 = pre+restraints[restraint][1]
    distance = float(restraints[restraint][2])

    if distance < strong:
        color = "red"
        category = "s"
    if ((distance < medium) & (distance >= strong)):
        color = "orange"
        category = "m"
    if distance >= medium:
        color = "yellow"
        category = "w"

    
    NOEID = "NOE_"+category+"_"+str(restraint)+"_"+restraints[restraint][2]
    draw(atom1, atom2, NOEID)

        
    _color_restraint(NOEID, color)

    
cmd.group("NOE_strong", members="NOE_s*", action='auto', quiet=1)
cmd.group("NOE_medium", members="NOE_m*", action='auto', quiet=1)
cmd.group("NOE_weak", members="NOE_w*", action='auto', quiet=1)
# cmd.order("NOE_*") this does nothing

cmd.select("DummyAtoms", selection="name Q*")
