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
    global DG
    DG = [None]*4
    for i in range(1, len(DG)):
        DG[i] = readJSON(f'json/DG{i}.json')
    
    global findBasisTemplate
    findBasisTemplate = readTXT('txt/findBasisTemplate.txt')
    findTestTemplate = readTXT('txt/findTestTemplate.txt')
    findDiffTestTemplate = readTXT('txt/findDiffTestTemplate.txt')
    
    for i in range(1,len(DG)):
        phi = DG[i]['phi']
        psi = DG[i]['psi']
        diffpsi = DG[i]['diffpsi']
        for j in range(0, len(phi)):
            findBasis             = findBasisTemplate.safe_substitute({f'DG{i}phi{j}': phi[j]})
            findBasisTemplate     = string.Template(findBasis)
        for j in range(0, len(psi)):
            findTest              = findTestTemplate.safe_substitute({f'DG{i}psi{j}': psi[j]})
            findTestTemplate      = string.Template(findTest)
            findDiffTest          = findDiffTestTemplate.safe_substitute({f'DG{i}diffpsi{j}': diffpsi[j]})
            findDiffTestTemplate  = string.Template(findDiffTest)
            
    findBasis     = findBasis.replace('**', '')
    findTest      = findTest.replace('**', '')
    findDiffTest  = findDiffTest.replace('**', '')
    
    writeC('c/findBasis.c', findBasis)
    writeC('c/findTest.c', findTest)
    writeC('c/findDiffTest.c', findDiffTest)
    return 0

if __name__ == "__main__":
    main()