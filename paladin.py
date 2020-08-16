'''
python score.py [-i input] [-o output] [-d] [-h]

for scoring peptides not comparing 13mer peptide array labels
Matrix is Sites(5) x Residues(20) x Terms(6)
  Terms = ['vdw','elec','harm','solvfe','cp','dunfe']
  Sites = [-2,-1,0,+1,+2] (index to [0<->4])
  Residues = {'E':0,'D':1,'K':2,'R':3,'Q':4,'N':5,'P':6,'H':7,'T':8,'S':9,'G':10,'A':11,'V':12,'M':13,'C':14,'I':15,'L':16,'Y':17,'F':18,'W':19}


'''

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


# default weights
wsite=[0.5, 0.5, 1.0, 0.2, 0.1]
wterm=[0.20, 0.60, 0.0, 0.0, 0.0, 1.0, 0.40]

# penalty for reverse-binding peptides
reverse_penalty = 1

####
#load params from numpy file rq
params=np.load('params.npy',allow_pickle=True)
#np.set_printoptions(precision=1,linewidth=80)

# model form:
# E(XXXXX) = Sum_site( w_site Sum_term( w_term E_{aa,term,site}  )  )

def score_fivemer(fivemer):
  tmpparam = []
  try:
    score = 0
    for ws,siteIdx,residue in zip(wsite,site_index,fivemer):
      score_site = 0
      resIdx=residuesDict[residue] #returns matrix index of residue I
      tmpparam = params[siteIdx][resIdx]
      for wt,termIdx in zip(wterm,term_index):
        termij = wt * tmpparam[termIdx]
        score_site += termij
      score += ws * score_site
    return score
  except KeyError:
    print(fivemer,'has some weird input')

def reverse(seq):
  return seq[::-1]

def score_peptide(seq,detail_set):
  seq_score = []
  for i in range(0,len(seq)-4):
    fivemer= seq[i:int(i+5)]
    if detail_set:
      columns=['seq','score','for-detail','reverse','rev-detail']
      score, fdetails = detail(fivemer)
      rscore, rdetails = detail(fivemer)
      seq_score.append([fivemer,score,fdetails,rscore,rdetails])
    else:
      columns=['seq','score','reverse']
      score = score_fivemer(fivemer)
      rscore = score_fivemer(reverse(fivemer)) + reverse_penalty
      details = ''
      seq_score.append([fivemer,score,rscore])
  df = pd.DataFrame(seq_score)
  df.columns=columns
  return df 

def detail(fivemer):
  details = {}
  try:
    score = 0
    for ws,siteIdx,residue in zip(wsite,site_index,fivemer):
      terms = {}
      score_site = 0
      resIdx=residuesDict[residue] #returns matrix index of residue I
      tmpparam = params[siteIdx][resIdx]
      for wt,termIdx in zip(wterm,term_index):
        t = ws * wt * tmpparam[termIdx]
        terms[termsDict[termIdx]] = t
        score_site += t
      details["site "+str(siteIdx-2)] = terms
      score += score_site
    return score, pd.DataFrame(details)

  except KeyError:
    print(fivemer,'has some weird input')


def read_fasta(ifile):
  Pepdict={}
  with open(ifile) as f: content = [i.strip() for i in f.readlines()]
  for line in content:
    if line[0] is '>':
      seq = line[1:]
    else:
      Pepdict[seq] = line
  return Pepdict

def get_scores(fasta,detail):
  out={}
  for name,seq in fasta.items(): out[name] = score_peptide(seq,detail)

  scores = pd.concat(out)
  scores.index.rename(['label','n'],inplace=True)
  scores.set_index(['seq'],inplace=True,append=True)
  return scores
  
def print_summary(scores):
  print('SUMMARY')
  print(scores)

def write_scores(scores,ofile):
  out = pd.DataFrame.from_dict(scores)
  out.to_csv(ofile,float_format='%.2f')

def main():
  # cl arguments
  parser = argparse.ArgumentParser(description='calculate scores for peptides saved in fasta format')
  parser.add_argument('-i','--infile', required=True,       dest='IFile',  help='fasta format, source of peptides to be scored')
  parser.add_argument('-o','--outfile',                     dest='OFile',  help='destination of scores')
  parser.add_argument('-d','--detail', action='store_true', dest='Detail', help='flag to give detailed score breakdown')
  parser.add_argument('-q','--quiet',  action='store_true', dest='Quiet',  help='flag to suppress summary ouput printed to screen')
  args=parser.parse_args()
  
  # read fasta-format input file
  # fasta format:
  #>name
  #SEQUENCE

  fasta=read_fasta(args.IFile)
  scores=get_scores(fasta,args.Detail)
  if not args.Quiet: print_summary(scores)
  if args.OFile is not None: write_scores(scores,args.OFile)
  if args.Quiet and args.OFile is None: print('use -o outputfile to save output to text file')

if __name__ == "__main__":
  main()
