# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 12:15:54 2020

@author: Runda
"""

import json
import string

def readTXT(fileName):
    f = open(f'{fileName}','r')
    s = f.read()
    f.close()
    return string.Template(s)

def readJSON(fileName):
    f = open(f'{fileName}','r')
    data = json.load(f)
    return data

def writeC(fileName, s):
    f = open(f'{fileName}','w')
    f.write(s)
    f.close()
    return 0

def main():
    global scheme
    scheme = [None]*3
    for i in range(0, len(scheme)):
        scheme[i] = readJSON(f'json/AF{2*(i+1)+1}.json')
    
    global findBasisTemplate
    findBasisTemplate = readTXT('txt/findBasisTemplate.txt')
    findDiffBasisTemplate = readTXT('txt/findDiffBasisTemplate.txt')
    findDiff2BasisTemplate = readTXT('txt/findDiff2BasisTemplate.txt')
    findIntBasisTemplate = readTXT('txt/findIntBasisTemplate.txt')
    
    for i in range(0,len(scheme)):
        phi      = scheme[i]['phi']
        diffphi  = scheme[i]['diffphi']
        diff2phi = scheme[i]['diff2phi']
        intphi   = scheme[i]['intphi']
        for j in range(0, len(phi)):
            findBasis              = findBasisTemplate.safe_substitute({f'AF{2*(i+1)+1}phi{j}': phi[j]})
            findBasisTemplate      = string.Template(findBasis)
            findDiffBasis          = findDiffBasisTemplate.safe_substitute({f'AF{2*(i+1)+1}phi{j}': diffphi[j]})
            findDiffBasisTemplate  = string.Template(findDiffBasis)
            findDiff2Basis         = findDiff2BasisTemplate.safe_substitute({f'AF{2*(i+1)+1}phi{j}': diff2phi[j]})
            findDiff2BasisTemplate = string.Template(findDiff2Basis)
            findIntBasis           = findIntBasisTemplate.safe_substitute({f'AF{2*(i+1)+1}phi{j}': intphi[j]})
            findIntBasisTemplate   = string.Template(findIntBasis)
            
    findBasis      = findBasis.replace('**', '')
    findDiffBasis  = findDiffBasis.replace('**', '')
    findDiff2Basis = findDiff2Basis.replace('**', '')
    findIntBasis   = findIntBasis.replace('**', '')
    s = findBasis + '\n' \
      + findDiffBasis + '\n' \
      + findDiff2Basis + '\n' \
      + findIntBasis + '\n'
    writeC('c/BasisFunc.c', s)
    return 0

if __name__ == "__main__":
    main()