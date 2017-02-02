#!/usr/bin/env python

from sys import argv,exit
from re import sub

suffixes = { 'float':'F', 
             'int':'I',
             'uint':'i',
             'uint64':'Ulong64_t',
           }
ctypes = {  
            'uint':'unsigned int',
            'float':'float',
            'int':'int',
            'uint64':'ULong64_t'
         }

class Branch:
  def __init__(self,name,dtype,counter='NMAX'):
    self.name = name
    self.dtype = dtype
    self.counter = counter
    self.suffix = ''
    try:
      self.suffix = '/'+suffixes[sub('Array','',dtype)]
      if 'Array' in dtype:
        self.suffix = '[%s]%s'%(counter,self.suffix)
    except KeyError:
      # must be a TObject
      self.suffix = dtype

class TreeClass:
  def __init__(self,name):
    self.name=name
    self.members = []
  def writeDef(self):
    s = '''
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "genericTree.h"
class %s : public genericTree {
  public:
    %s();
    ~%s();
    void ReadTree(TTree *t);
    void WriteTree(TTree *t);
    void Reset();
    \n'''%(self.name,self.name,self.name)
    for member in self.members:
      s += createDefString(member)
    s += '//ENDDEF\n};\n'
    return s
  def writeConst(self):
    s = '''
%s::%s() {
    '''%(self.name,self.name)
    for member in self.members:
      s+= createInitString(member)
    s += '//ENDCONST\n}\n'
    return s
  def writeDest(self):
    s = '''
%s::~%s() {
    '''%(self.name,self.name)
    for member in self.members:
      if ' new ' in createInitString(member):
        s += '\tdelete %s; %s = 0;\n'%(member.name,member.name)
    s += '//ENDDEST\n}\n'
    return s
  def writeResetTree(self):
    s = '''
void %s::Reset() {
    '''%(self.name)
    for member in self.members:
      s+= createResetString(member) 
    s += '//ENDRESET\n}\n'
    return s
  def writeReadTree(self):
    s = '''
void %s::ReadTree(TTree *t) {
      treePtr = t;
      treePtr->SetBranchStatus("*",0);
    '''%(self.name)
    for member in self.members:
      s+= createInString(member,'treePtr')
    s += '//ENDREAD\n}\n'
    return s
  def writeWriteTree(self):
    s = '''
void %s::WriteTree(TTree *t) {
      treePtr = t;
    '''%(self.name)
    for member in self.members:
      s+= createOutString(member,'treePtr')
    s += '//ENDWRITE\n}\n'
    return s
    

def createDefString(branch):
  sInit = ''
  # initialize
  if branch.dtype=='TClonesArray':
    sInit += '\tTClonesArray *%s=0;\n'%(branch.name)
  elif branch.dtype=='TLorentzVector':
    sInit += '\tTLorentzVector *%s=0;\n'%(branch.name)
  else:
    if 'Array' in branch.dtype:
      basetype = ctypes[sub('Array','',branch.dtype)]
      sInit += '\t%s *%s=0;\n'%(basetype,branch.name)
    else:
      sInit += '\t%s %s=0;\n'%(ctypes[branch.dtype],branch.name)

  return sInit

def createInitString(branch):
  sInit = ''
  # initialize
  if branch.dtype=='TClonesArray':
    sInit += '\t%s = new TClonesArray("TLorentzVector",%s);\n'%(branch.name,str(branch.counter))
  elif branch.dtype=='TLorentzVector':
    sInit += '\t%s = new TLorentzVector();\n'%(branch.name)
  else:
    if 'Array' in branch.dtype:
      basetype = ctypes[sub('Array','',branch.dtype)]
      counter = 'NGENMAX' if 'mc' in branch.name else 'NMAX'
      sInit += '\t%s = new %s[%s];\n'%(branch.name,basetype,counter)
    else:
      sInit += '\t%s=0;\n'%(branch.name)
  
  return sInit

def createResetString(branch):
  sReset = ''
  if 'Array' in branch.dtype  or branch.dtype=='TLorentzVector':
    pass
  elif 'chi' in branch.name:
    sReset += '\t%s = -99;\n'%(branch.name)
  elif branch.dtype=='float':
    sReset += '\t%s = -1;\n'%(branch.name)
  else: #int-like
    sReset += '\t%s = 0;\n'%(branch.name)

  return sReset


def createInString(branch,treeName='t'):
  # set branch address
  sSet = ''
  sSet += '\t%s->SetBranchStatus("%s",1);\n'%(treeName,branch.name)
  if 'Array' in branch.dtype and not(branch.dtype[0]=='T'):
    sSet += '\t%s->SetBranchAddress("%s",%s);\n'%(treeName,branch.name,branch.name)
  else:
    sSet += '\t%s->SetBranchAddress("%s",&%s);\n'%(treeName,branch.name,branch.name)

  return sSet

def createOutString(branch,treeName='tOut'):
  sCreate = ''
  if branch.dtype[0]=='T':
    sCreate += '\t%s->Branch("%s","%s",&%s,128000,0);\n'%(treeName,branch.name,branch.dtype,branch.name)
  else:
    sCreate += '\t%s->Branch("%s",&%s,"%s%s");\n'%(treeName,branch.name,branch.name,branch.name,branch.suffix)
  
  return sCreate

def parseCfg(cfg,className):
  cfg = list(cfg)
  tree = TreeClass(className)
  for line in cfg:
    ll = line.split()
    if ll[0]=='#':
      continue
    varname = ll[0]
    dtype = ll[1]
    counter='NMAX' if len(ll)==2 else ll[2]
    b = Branch(varname,dtype,counter)
    tree.members.append(b)
  s = tree.writeDef()
  s += tree.writeConst()
  s += tree.writeDest()
  s += tree.writeResetTree()
  s += tree.writeReadTree()
  s += tree.writeWriteTree()
  return s


if len(argv)<3:
  print 'Usage %s className cfgFile'%(argv[0])
  exit(1)

className = argv[1]
cfg = open(argv[2])
s = parseCfg(cfg,className)
with open(sub('config','interface',sub('cfg','h',argv[2])),'w') as outFile:
  outFile.write(s)

