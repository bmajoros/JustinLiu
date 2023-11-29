#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2018 William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import math
import ProgramName
from Rex import Rex
rex=Rex()
import TempFilename
import getopt
from scipy.optimize import minimize
from scipy.stats import beta as scipybeta
import numpy as np

DEBUG=False
WARMUP=1000
ALPHA=0.05
STDERR=TempFilename.generate(".stderr")
INPUT_FILE=TempFilename.generate(".staninputs")
INIT_FILE=TempFilename.generate(".staninit")
OUTPUT_TEMP=TempFilename.generate(".stanoutputs")

def fit(samples):
    def f(x):
        logLik=sum(betaModeConc(samples,x[0],x[1]))
        return -logLik
    x0 = np.array([0.5,30]) # mode and concentration
    rec = minimize(f,x0,
                   method="Powell",
                   bounds=[[0,1],[2,1000000]],
                   options={'maxiter': 50, 'disp': False})
    (mode,conc)=rec.x
    return (mode,conc)
    
def betaModeConc(parm,m,c):
    logPDF=scipybeta.logpdf(parm, m*(c-2)+1, (1-m)*(c-2)+1)
    return logPDF

def printFields(fields,hFile):
    numFields=len(fields)
    for i in range(7,numFields):
        print(i-6,"=",fields[i],sep="",end="",file=hFile)
        if(i<numFields-1): print("\t",end="",file=hFile)
    print(file=hFile)

def getFieldIndex(label,fields):
    numFields=len(fields)
    index=None
    for i in range(7,numFields):
        if(fields[i]==label): index=i
    return index

def writeToFile(fields,OUT):
    numFields=len(fields)
    for i in range(7,numFields):
        print(fields[i],end="",file=OUT)
        if(i<numFields-1): print("\t",end="",file=OUT)
    print(file=OUT)

def writeReadCounts(counts,index,varName,OUT):
    print(varName,"<- c(",file=OUT,end="")
    NUM_RNA=len(counts)
    for rep in range(NUM_RNA):
        rec=counts[rep]
        print(rec[index],file=OUT,end="")
        if(rep+1<NUM_RNA): print(",",file=OUT,end="")
    print(")",file=OUT)

def writeInitializationFile(counts,filename):
    NUM_RNA=len(counts)
    totalRef=0; totalAlt=0
    for i in range(NUM_RNA):
        (alt,ref)=counts[i]
        totalRef+=ref
        totalAlt+=alt
    v=float(totalAlt+1)/float(totalAlt+totalRef+2)
    if(v==0): v=0.01
    OUT=open(filename,"wt")
    print("q <-",v,file=OUT)
    print("qi <- c(",file=OUT,end="")
    for i in range(NUM_RNA-1):
        print(v,",",sep="",end="",file=OUT)
    print(v,")",sep="",file=OUT)
    OUT.close()

def writeInputsFile(fields,filename):
    NUM_RNA=len(counts)
    OUT=open(filename,"wt")
    print("N_RNA <-",str(NUM_RNA),file=OUT)
    writeReadCounts(counts,0,"k",OUT) # alt
    writeReadCounts(counts,1,"m",OUT) # ref
    OUT.close()

def getMedian(thetas):
    # Precondition: thetas is already sorted
    n=len(thetas)
    mid=int(n/2)
    if(n%2==0): 
        return (thetas[mid-1]+thetas[mid])/2.0
    return thetas[mid]

def getCredibleInterval(thetas,alpha):
    halfAlpha=alpha/2.0
    n=len(thetas)
    leftIndex=int(halfAlpha*n)
    rightIndex=n-leftIndex
    left=thetas[leftIndex+1]
    right=thetas[rightIndex-1]
    return (left,right)

def runVariant(model,counts,numSamples,outfile):
    # Write inputs file for STAN
    NUM_RNA=len(counts)
    writeInputsFile(counts,INPUT_FILE)
    writeInitializationFile(counts,INIT_FILE)

    # Run STAN model
    init=" init="+INIT_FILE
    cmd=model+" sample thin=1"+\
        " num_samples="+numSamples+\
        " num_warmup="+str(WARMUP)+\
        " data file="+INPUT_FILE+\
        init+\
        " output file="+OUTPUT_TEMP+" refresh=0 > "+STDERR
    if(DEBUG):
        print(cmd)
    os.system(cmd)

    # Parse MCMC output
    samples=[]; qIndex=None
    OUT=None if outfile=="." else open(outfile,"wt")
    with open(OUTPUT_TEMP,"rt") as IN:
        for line in IN:
            if(len(line)==0 or line[0]=="#"): continue
            fields=line.rstrip().split(",")
            numFields=len(fields)
            if(numFields>0 and fields[0]=="lp__"):
                if(OUT is not None): printFields(fields,OUT)
                qIndex=getFieldIndex("q",fields)
            else:
                if(OUT is not None): writeToFile(fields,OUT)
                q=float(fields[qIndex])
                samples.append(q)
    if(OUT is not None): OUT.close()
    samples.sort(key=lambda x: x)
    return samples


#=========================================================================
# main()
#=========================================================================
(options,args)=getopt.getopt(sys.argv[1:],"s:")
if(len(args)!=4):
    exit(ProgramName.get()+" [-s stanfile] <model> <input.txt> <output.txt> <#MCMC-samples>\n   -s = save raw STAN file\n")
(model,inFile,outfile,numSamples)=args
stanFile=None
for pair in options:
    (key,value)=pair
    if(key=="-s"): stanFile=value
    
# Process all input lines, each line = one variant (one MCMC run)
with open(inFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        numFields=len(fields)
        if(numFields<11): raise Exception("Wrong number of fields in input")
        q=fields[0]
        NUM_RNA=int(fields[1])
        counts=[]
        nextField=2
        for i in range(NUM_RNA):
            (k,m,modeI,concI,lnA,lnB,lnGamAB,lnGamA,lnGamB)=\
                fields[nextField:(nextField+9)]
            nextField+=9
            counts.append([int(k),int(m)])
        samples=runVariant(model,counts,numSamples,outfile)
        (mode,conc)=fit(samples)
        alpha=mode*(conc-2)+1
        beta=(1-mode)*(conc-2)+1
        print(round(mode,3),round(conc,1),round(alpha,3),round(beta,3),
              sep="\t",end="")
        for i in range(NUM_RNA):
            print("\t",end="")
            print(k,m,modeI,concI,lnA,lnB,lnGamAB,lnGamA,lnGamB,
                  sep="\t",end="")
        print()
os.remove(STDERR)
os.remove(INPUT_FILE)
if(stanFile is None):
    os.remove(OUTPUT_TEMP)
else:
    os.system("cp "+OUTPUT_TEMP+" "+stanFile)

