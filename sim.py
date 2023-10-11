#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2021 William H. Majoros <bmajoros@alumni.duke.edu>
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
import numpy as np
import random
from random import randint

LAMBDA=30

def betaModeConc(m,c):
    if(m==0 or m==1): raise Exception("Invalid mode")
    if(c<=2): raise Exception("Invalid concentration")
    alpha=m*(c-2)+1
    beta=(1-m)*(c-2)+1
    return np.random.beta(alpha,beta)

def sampleFreq():
    p=0
    while(p==0 or p==1):
        p=round(np.random.uniform(0,1),3)
    return p
    
def sim():
    p=sampleFreq()
    theta=0; q=0
    while(q==0 or q==1):
        theta=np.random.lognormal(0,1)
        q=round(theta*p/(1-p+theta*p),3)
    if(p==0 or p==0 or q==0 or q==1): return 0
    N1=random.randint(1,10); N2=random.randint(1,10)
    p=round(p,3); theta=round(theta,3); q=round(q,3)
    print(p,q,theta,N1,N2,sep="\t",end="")
    conc=0
    while(conc<=2): conc=np.random.gamma(1.1, 1/0.0005)
    for i in range(N1):
        #print("p=",p,conc,flush=True);print(flush=True)
        pi=betaModeConc(p,conc)
        n=np.random.poisson(LAMBDA)
        k=np.random.binomial(n,pi)
        m=n-k
        pi=round(pi,3)
        print("\t",pi,"\t",k,"\t",m,sep="",end="")
    for i in range(N2):
        #print("q=",q,conc,flush=True);print(flush=True)
        qi=betaModeConc(q,conc)
        n=np.random.poisson(LAMBDA)
        k=np.random.binomial(n,qi)
        m=n-k
        qi=round(qi,3)
        print("\t",qi,"\t",k,"\t",m,sep="",end="")
    print()
    return 1

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <#cases>\n")
(numCases,)=sys.argv[1:]
numCases=int(numCases)

i=0
while(i<numCases):
    if(sim()>0): i+=1


