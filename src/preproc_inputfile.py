################################################################
################# Preprocessing for CPP_UEL ####################
##    Identify the "auelPart*" part, and declare the UEL and its
##  element properties. Multiple UEL parts are enabled, so that
##  the local node and element label are changed for global
##  consistent numbering.
################################################################
##
##
################################################################
##  Applicable abaqus input file:
##  - multiple parts with name "auelPart1", "auelPart2" ...
##  - part definition can only include *Node and *Element
##  - assembly definition is copied without modification
################################################################
##
##  Record of revision:
##    Date      Programmer            Description of change
##    ========  ====================  ========================================
##    15/01/25  Lu Xin                New code for CPP_UEL
##
################################################################
##
##  Comments:
##    Global node and element counters are defined, and after
##  looping over one part, the node and element index will be
##  shifted by the previously accumulated numbers.
##
###############################################################

import math
# operating system module
import os, sys
# shutil module has functions for different kinds of copying
import shutil

#***************************************************************
#   define parameters
#***************************************************************
# elnode  : no. of nodes defined in this user cohesive element
# ndim    : dimension; currently only support 3D
# nprops  : no. of input material properties used in uel code (min 1)
# nsvars  : no. of sol. dpdnt var. used and ouput by uel code (min 1)
# uelcode : code for the user cohesive element, used in uel input file
# uellinecount  : max no. of entries in a data line in uel input file (16)
# uellinelength : line length in writing uel input file (max. =80 in *Element)
ndim     = 3
nprops   = 10
njprops  = 1
nsvars   = 1
uelcode6 = 306
uelcode8 = 308
uellinecount  = 14
uellinelength = 70

#***************************************************************
#   fetch original input file and open uel files
#***************************************************************
# ask for job name
jobname      = raw_input('abaqus job name:')
# abaqus original input file
abqinputfile = jobname+'.inp'
# uel input file
uelinputfile = 'uel_'+jobname+'.inp'

# define UEL parts identifier
uel_id    = 'auelPart'

# open files
abq_input = open(abqinputfile,'r')
uel_input = open(uelinputfile,'w')


#***************************************************************
#       Read Abaqus input file
#***************************************************************

# Store all lines first
All_lines = abq_input.readlines()
# Find the length of all lines
lenAll = len(All_lines)

#==================================================
# READ PARTS SECTION:
# - Only cohesive parts are extracted and modified
# - Other parts will be directly copied
# ** NOTE **:
# - for now, only ONE type of CE is supported, it
#   is used in parts with "cohesive" in name
#==================================================
parts   = []    # list of all parts in the model

# find the line no. of *Part and store in jparts
jparts = [j for j,line in enumerate(All_lines) if '*Part' in line]
# find the line no. of *End Part and store in eparts
eparts = [j for j,line in enumerate(All_lines) if '*End Part' in line]

# check the consistence
if (len(jparts) != len(eparts)):
    print("ERROR: inconsistent part numbers")
    sys.exit()

node_id = 0
elem_id = 0
node_shift = 0
elem_shift = 0

for ip in range( len(jparts) ):
# loop over all parts

    # get the starting and ending line numbers of this part
    start_line = jparts[ip]
    end_line   = eparts[ip]

    # read Part name
    pname = All_lines[start_line][12:].rstrip()

    # check if uel_id - 'auelPart' in part name
    if(uel_id in pname):
    # this is a cohesive part

        # find the line of *Node and *Element
        jnode = next( j for j in range(start_line,end_line) if '*Node'    in All_lines[j] )
        jelem = next( j for j in range(start_line,end_line) if '*Element' in All_lines[j] )

        #*******************************************#
        # read the nodes inf
        #*******************************************#
        for jn in range(jnode+1,jelem):
            # read the node inf line by line
            nline = All_lines[jn]
            # break out of for loop if end of node section is reached
            if('*' in nline):
                break

            # read the coords of this node into a list of float numbers
            coords = []
            for t in nline.split(','):
                try:
                    coords.append(float(t))
                except ValueError:
                    pass
            # keep tracking the local node id
            node_id = int(coords[0])
            # update the local node id to global considering node shift
            # -- the node_shift will be updated at the end of this part
            #    node_shift: node_shift = node_shift + node_id
            #    always shift the current no. of nodes
            coords[0] = int(coords[0] + node_shift)

            # define the current line string
            nline = str(coords[0]) + ', ' + str(coords[1]) + ', ' +\
            str(coords[2]) + ', ' + str(coords[3])+'\n'
            # update to All_lines variable
            All_lines[jn] = nline

        #*******************************************#
        # Write the UEL declarations
        #*******************************************#
        # get the no of node for this element
        line_lcl  = All_lines[jelem+1]
        el_lcl    = []
        for t in line_lcl.split(','):
            try:
                el_lcl.append(int(t))
            except ValueError:
                pass
        elnode = len(el_lcl) - 1
        if(elnode == 6):
            uelcode = uelcode6
        elif(elnode == 8):
            uelcode = uelcode8
        else:
            print("ERROR: no of element node NOT supported")
            sys.exit()
        # add the uel header inf
        uelline = '*USER ELEMENT, TYPE=U'+str(uelcode)+\
        ', NODES='+str(elnode)+', COORDINATES='+str(ndim)+\
        ', PROPERTIES='+str(nprops)+ ', IPROPERTIES='+str(njprops)+\
        ', VARIABLES='+str(nsvars)+'\n'
        uelline = uelline + '1,2,3\n'
        uelline = uelline + '*ELEMENT, TYPE=U'+str(uelcode)+', ELSET='+\
        uel_id+str(ip+1)+'\n'
        All_lines[jelem] = uelline

        #*******************************************#
        #  read the elem inf
        #*******************************************#
        for je in range(jelem+1,end_line):
            # read the elem inf line by line
            eline = All_lines[je]
            # break out of for loop if end of elem section is reached
            if('*' in eline):
                break
            # read the index and real nodes of this elem into a list of int numbers
            el = []
            for t in eline.split(','):
                try:
                    el.append(int(t))
                except ValueError:
                    pass
            # keep tracking the local element id
            elem_id = el[0]
            # update the local elem id to global considering elem shift
            # -- the elem_shift will be updated at the end of this part
            #    elem_shift: elem_shift = elem_shift + elem_id
            #    always shift the current no. of elems
            el[0] = el[0] + elem_shift

            # the elem node cnc will be shifted by the same amount of node shift
            for ke in range(len(el)-1):
                el[ke+1] = el[ke+1] + node_shift

            eline = ''
            for k in range( len(el) - 1):
                eline = eline + str(el[k]) + ', '
            eline = eline + str(el[k+1]) + '\n'
            # update to All_lines variable
            All_lines[je] = eline

        #*******************************************#
        #   Write the UEL property after elem cnc
        #*******************************************#
        All_lines[je+1] = All_lines[je+1] + \
        '*UEL PROPERTY, ELSET='+uel_id+str(ip+1)+'\n'+\
        '1.0e6, 57., 57., 0.28, 0.28, 1.0, 0.0, 0.0\n'+\
        '1.0, 1.0, '+str(ip+1)+'\n' +\
        '** K,  tnc, tsc, Gnc,  Gsc, alpha, neta, cof\n'+\
        '** dmax, power_fg, part_id\n'

        #*******************************************#
        #   Update shift place
        #*******************************************#
        node_shift = node_shift + node_id
        elem_shift = elem_shift + elem_id

    # if 'auelPart' is not in the part name, end of loop and
    # write markers for the end of UEL part
    if(not (uel_id in pname)):
        All_lines[start_line] = '** end of auelPart ** \n' + All_lines[start_line]
        break

#==================================================
# WRITE UEL INP AND FORTRAN FILES:
#==================================================
# write uel file
for i in range(lenAll):
    uel_input.write(All_lines[i])
uel_input.close()
