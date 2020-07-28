'''
python score.py -i input [-o output]

for scoring peptides not comparing 13mer peptide array labels
Matrix is Sites(5) x Residues(20) x Terms(6)
  Terms = ['vdw','elec','harm','solvfe','cp','dunfe']
  Sites = [-2,-1,0,+1,+2] (index to [0<->4])
  Residues = {'E':0,'D':1,'K':2,'R':3,'Q':4,'N':5,'P':6,'H':7,'T':8,'S':9,'G':10,'A':11,'V':12,'M':13,'C':14,'I':15,'L':16,'Y':17,'F':18,'W':19}


'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
pd.set_option("display.precision", 3)

site_index   = [0,1,2,3,4]
term_index   = [1,2,3,4,5,6,7]
termsDict    = {1:'vdw',2:'elec',3:'harm',4:'cp',5:'dunfe',
                6:'solvFE-resN',7:'solvFE-site'}

residues     = ['E','D','K','R','Q','N','P','H','T','S',
                'G','A','V','M','C','I','L','Y','F','W']
residuesDict = {'E':0,'D':1,'K':2,'R':3,'Q':4,'N':5,'P':6,
                'H':7,'T':8,'S':9,'G':10,'A':11,'V':12,'M':13,
                'C':14,'I':15,'L':16,'Y':17,'F':18,'W':19}

wsite=[0.2, 0.6, 1.0, 0.4, 0.2]
wterm0=[0.20, 0.60, 0.0, 0.0, 0.0, 1.0, 0.40]
wterm_out=wterm0

####
#load params from numpy file rq
params=np.load('params.npy',allow_pickle=True)
np.set_printoptions(precision=1,linewidth=80)

# model form:
# E(XXXXX) = Sum_site( w_site Sum_term( w_term E_{aa,term,site}  )  )

tmpparam = []
def score(fivemer):
  try:
    score = 0
    for ws,siteIdx,residue in zip(wsite,site_index,fivemer):
      score_site = 0
      resIdx=residuesDict[residue] #returns matrix index of residue I
      tmpparam = params[siteIdx][resIdx]
      if siteIdx is 2: wterms=wterm0
      else: wterms=wterm_out
      for wt,termIdx in zip(wterms,term_index):
        score_site += wt * tmpparam[termIdx]
      score += ws * score_site
    return score
  except KeyError:
    print(fivemer,'has an unnatural amino acid')

def score_peptide(seq):
  seq_score = []
  for i in range(0,len(seq)-4):
    pep = seq[i:int(i+5)]
    s = score(pep)
    seq_score.append([pep,s])
  df = pd.DataFrame(seq_score)
  df.columns = ['seq','score']
  return df 

def read_fasta(ifile):
  Pepdict={}
  with open(ifile) as f: content = [i.strip() for i in f.readlines()]
  for line in content:
    if line[0] is '>':
      key = line[1:]
    else:
      Pepdict[key] = line
  return Pepdict

def write_scores(sequences,ofile):
  out={}
  for name,seq in sequences.items():
    scores = score_peptide(seq)  
    print(scores)
    out[name] = scores
  #print(out)
  out = pd.DataFrame.from_dict(out,orient='index',columns=['names'])
  out.to_csv(ofile)



def main():
  # cl arguments
  ifile='test.fasta'
  ofile=ifile+'.out'
  parser = argparse.ArgumentParser(description='calculate scores for peptides')
  parser.add_argument('-i','--infile',default=ifile,dest='IFile',help='source of peptides to be scored')
  parser.add_argument('-o','--outfile',default=ofile,dest='OFile',help='destination of scores')
  parser.add_argument('-d','--detail',action='store_true',dest='detail',help='flag to give detailed score breakdown')
  args=parser.parse_args()
  
  # read fasta-format input file
  # fasta format:
  #>name
  #SEQUENCE
  sequences=read_fasta(ifile)
  write_scores(sequences,ofile)

if __name__ == "__main__":
  main()
