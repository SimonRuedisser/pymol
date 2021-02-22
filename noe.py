#########################################################################
# Script to display upper distance limits (from Cyana final.upl) as lines
# on pdb output from cyana (finalQ.pdb)
# #######################################################################
#
# in pymol type:
# run noe.py
#
# S. Ruedisser, 2018 April, ETH Zürich


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
