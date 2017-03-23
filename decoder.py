#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 17:00:37 2017

@author: valeriya.potapova
"""

def initializeDataStructures():
    global inactivePathIndices;
    inactivePathIndices=[[]]*L;
    global activePath;
    activePath=[0]*L;
    global arrayPointer_P;
    arrayPointer_P=[[]]*(m+1);
    global arrayPointer_C;
    arrayPointer_C=[[]]*(m+1);
    global pathIndexToArrayIndex;
    pathIndexToArrayIndex=[[]]*(m+1);
    global inactiveArrayIndices;
    inactiveArrayIndices=[[]]*(m+1);
    global arrayReferenceCount;
    arrayReferenceCount=[0]*(m+1);

    for i in range(m+1):
        arrayReferenceCount[i]=[0]*L;
        pathIndexToArrayIndex[i]=[0]*L;

    for lambd in range(m+1):
        arrayPointer_P[lambd]=[[]]*L;
        arrayPointer_C[lambd]=[[]]*L;

    for lambd in range(m+1):
        for s in range(L):
            arrayPointer_P[lambd][s]=[[]]*(2**(m-lambd));
            arrayPointer_C[lambd][s]=[[]]*(2**(m-lambd));

    for lambd in range(m+1):
        for s in range(L):
            for d in range(2**(m-lambd)):
                arrayPointer_P[lambd][s][d]=[0]*2;
                arrayPointer_C[lambd][s][d]=[0]*2;

    for lambd in range(m+1):
        for s in range(L):
            inactiveArrayIndices[lambd].append(s);

    for l in range(L):
        activePath[l]=False;
        inactivePathIndices.append(l);

def assignInitialPath():
    l=inactivePathIndices.pop();
    activePath[l]=True;
    for lambd in range(m):
        s=inactiveArrayIndices[lambd].pop();
        pathIndexToArrayIndex[lambd][l]=s;
        arrayReferenceCount[lambd][s]=1;
    return l

def clonePath(l):
    l_t=inactivePathIndices.pop();
    activePath[l_t]=True;
    for lambd in range(m):
        s=pathIndexToArrayIndex[lambd][l];
        pathIndexToArrayIndex[lambd][l_t]=s;
        arrayReferenceCount[lambd][s]=arrayReferenceCount[lambd][s]+1;
    return l_t;

def killPath(l):
    activePath[l]=False;
    inactivePathIndices.append(l);
    for lambd in range(m):
        s=pathIndexToArrayIndex[lambd][l];
        arrayReferenceCount[lambd][s]=arrayReferenceCount[lambd][s]-1;
        if arrayReferenceCount[lambd][s]==0:
            inactiveArrayIndices[lambd].append(s);

def getArrayPointer_P(lambd,l):
    s=pathIndexToArrayIndex[lambd][l];
    if arrayReferenceCount[lambd][s]==1:
        s_t=s;
    else:
        s_t=inactiveArrayIndices[lambd].pop();
        for p in range(2**(m-lambd)):
            for q in range(2):
                arrayPointer_P[lambd][s_t][p][q]=arrayPointer_P[lambd][s][p][q];
        arrayReferenceCount[lambd][s]=arrayReferenceCount[lambd][s]-1;
        arrayReferenceCount[lambd][s_t]=1;
        pathIndexToArrayIndex[lambd][l]=s_t;

    return arrayPointer_P[lambd][s_t]

def getArrayPointer_C(lambd,l):
    s=pathIndexToArrayIndex[lambd][l];
    if arrayReferenceCount[lambd][s]==1:
        s_t=s;
    else:
        s_t=inactiveArrayIndices[lambd].pop();
        for p in range(2**(m-lambd)):
            for q in range(2):
                arrayPointer_C[lambd][s_t][p][q]=arrayPointer_C[lambd][s][p][q];
        arrayReferenceCount[lambd][s]=arrayReferenceCount[lambd][s]-1;
        arrayReferenceCount[lambd][s_t]=1;
        pathIndexToArrayIndex[lambd][l]=s_t;

    return arrayPointer_C[lambd][s_t]

def pathIndexInactive(l):
    if activePath[l]==True:
        return False
    else:
        return True

def stop(p,q):
    res=np.exp((-((p-q)**2))/(2*var));
    res=res/np.sqrt(2*np.pi*var);
    return res

def recursivelyCalcP(lambd,phi):
    if lambd==0:
        return;
    psi=phi//2;
    if phi%2==0:
        recursivelyCalcP(lambd-1,psi);
    sigma=0;
    for l in range(L):
        if pathIndexInactive(l):
            continue;
        P_lambd=getArrayPointer_P(lambd,l);
        P_lambd_1=getArrayPointer_P(lambd-1,l);
        C_lambd=getArrayPointer_C(lambd,l);
        for beta in range(2**(m-lambd)):
            if phi%2==0:
                P_lambd[beta][0]=P_lambd_1[2*beta][0]*P_lambd_1[2*beta+1][0]/2+P_lambd_1[2*beta][1]*P_lambd_1[2*beta+1][1]/2;
                P_lambd[beta][1]=P_lambd_1[2*beta][1]*P_lambd_1[2*beta+1][0]/2+P_lambd_1[2*beta][0]*P_lambd_1[2*beta+1][1]/2;
                sigma=max(sigma,P_lambd[beta][0],P_lambd[beta][1]);
            else:
                u_t=int(C_lambd[beta][0]);
                P_lambd[beta][0]=P_lambd_1[2*beta][u_t]*P_lambd_1[2*beta+1][0]/2;
                P_lambd[beta][1]=P_lambd_1[2*beta][(u_t+1)%2]*P_lambd_1[2*beta+1][1]/2;
                sigma=max(sigma,P_lambd[beta][0],P_lambd[beta][1]);
    for l in range(L):
        if pathIndexInactive(l):
            continue;
        P_lambd=getArrayPointer_P(lambd,l);
        for beta in range(2**(m-lambd)):
            for p in {0,1}:
                P_lambd[beta][p]=P_lambd[beta][p]/sigma;

def recursivelyUpdateC(lambd,phi):
    psi=phi//2;
    for l in range(L):
        if pathIndexInactive(l):
            continue;
        C_lambd_1=getArrayPointer_C(lambd-1,l);
        C_lambd=getArrayPointer_C(lambd,l);
        for beta in range(2**(m-lambd)):
            C_lambd_1[2*beta][psi%2]=(C_lambd[beta][0]+C_lambd[beta][1])%2;
            C_lambd_1[2*beta+1][psi%2]=C_lambd[beta][1];
    if psi%2==1:
        recursivelyUpdateC(lambd-1,psi);

def continuePaths_FrozenBit(phi):
    for l in range(L):
        if pathIndexInactive(l):
            continue;
        C_m=getArrayPointer_C(m,l);
        C_m[0][phi%2]=u[phi];

def continuePaths_UnfrozenBit(phi):
    probForks=[[]]*L;
    for l in range(L):
        probForks[l]=[0]*2;

    i=0;
    for l in range(L):
        if pathIndexInactive(l):
            probForks[l][0]=-1;
            probForks[l][1]=-1;
        else:
            P_m=getArrayPointer_P(m,l);
            probForks[l][0]=P_m[0][0];
            probForks[l][1]=P_m[0][1];
            i=i+1;
    rho=min(2*i,L);
    contForks=[[]]*L;
    for l in range(L):
        contForks[l]=[0]*2;
    k=0;
    while k<rho:
        maximum=-1;
        ind_x=-1;
        ind_y=-1;
        for l in range(L):
            for d in range(2):
                if probForks[l][d]>maximum and contForks[l][d]==False:
                    maximum=probForks[l][d];
                    ind_x=l;
                    ind_y=d;
        contForks[ind_x][ind_y]=True;
        k=k+1;

    for l in range(L):
        if pathIndexInactive(l):
            continue;
        if contForks[l][0]==False and contForks[l][1]==False:
            killPath(l);

    for l in range(L):
        if contForks[l][0]==False and contForks[l][1]==False:
            continue;
        C_m=getArrayPointer_C(m,l);
        if contForks[l][0]==True and contForks[l][1]==True:
            C_m[0][phi%2]=0;
            l_t=clonePath(l);
            C_m=getArrayPointer_C(m,l_t);
            C_m[0][phi%2]=1;
        else:
            if contForks[l][0]==True:
                C_m[0][phi%2]=0;
            else:
                C_m[0][phi%2]=1;

def findMostProbablePath():
    l_t=0;
    p_t=0;
    for l in range(L):
        if pathIndexInactive(l):
            continue;
        C_m=getArrayPointer_C(m,l);
        P_m=getArrayPointer_P(m,l);
        if p_t<P_m[0][C_m[0][1]]:
            l_t=l;
            p_t=P_m[0][C_m[0][1]];
    return l_t

"----------------------------------"
"Here we start the main loop"
import numpy as np
global c;
global u;
global n;
F=[0,1,2,4];
UnF=[7,3,5,6]; "We have to decode the vector"
u=[0,0,0,-1,0,-1,-1,-1];"the set of frozen (known) bits."
s=[0,0,0,0]; "The source symbols"
k=4; "the dimension of the code"
c=[1,1,1,1,0,0,0,0]; "codeword to test"
n=len(c);
global L; "the length of the list"
L=1;
global m;
m=int(np.log2(n)); "the log of the length of the code"
global snr;
snr=10;
global var;
var=0.001;
y=[];
"Now I have to add some noise";
for i in range(n):
    y.append(round((c[i]+np.random.randn()*var),5));
"This vector y I have to decode"
initializeDataStructures();
l=assignInitialPath();
P_0=getArrayPointer_P(0,l);
for beta in range(n):
    P_0[beta][0]=stop(y[beta],0);
    P_0[beta][1]=stop(y[beta],1);
for phi in range(n):
    recursivelyCalcP(m,phi);
    if phi in F:
        continuePaths_FrozenBit(phi);
    else:
        continuePaths_UnfrozenBit(phi);
    if phi%2==1:
        recursivelyUpdateC(m,phi);

l=findMostProbablePath();
C_0=getArrayPointer_C(0,l);
codeword=[];
for beta in range(n):
    codeword.append(C_0[beta][0]);
print(codeword)