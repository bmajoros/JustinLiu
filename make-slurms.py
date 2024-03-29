#!/usr/bin/env python
#=========================================================================
# Copyright (C)William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import os
import sys
import ProgramName
from SlurmWriter import SlurmWriter
from Rex import Rex
rex=Rex()

BASE="/hpc/group/majoroslab/BIRD/justin"
MCMC_SAMPLES=1000
JOB_NAME="BIRD"
MAX_PARALLEL=300
additional_SBATCH_lines=""

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <slurms-dir> <chunks-dir> <model> <out-dir>\n")
(slurmsDir,chunksDir,model,outDir)=sys.argv[1:]

slurm=SlurmWriter()
files=os.listdir(chunksDir)
for filename in files:
    if(not rex.find("^chunk(\d+).txt",filename)): continue
    chunkNum=rex[1]
    cmd="cd "+BASE+"\ngit/justin.py "+model+" "+chunksDir+\
        "/chunk"+chunkNum+".txt"+\
        " . "+str(MCMC_SAMPLES)+\
        " > "+outDir+"/out"+chunkNum+".txt"
    slurm.addCommand(cmd)
slurm.mem(1500)
slurm.setQueue("majoroslab-gpu")
slurm.writeArrayScript(slurmsDir,JOB_NAME,MAX_PARALLEL,
                       additional_SBATCH_lines)

