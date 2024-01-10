
#RECOMBINATION FREQUENCY:

def recom_freq():

   print("You have to enter the frequency of 2 non recombinant and 2 recombinant genotypes :")

   nr1 = int(input("Frequency of first non recombinant genotype: "))
   nr2 = int(input("Frequency of second non recombinant genotype: "))
   r1 = int(input("Frequency of first recombinant genotype: "))
   r2 = int(input("Frequency of second recombinant genotype: "))

   rf = (r1+r2)/(nr1 + nr2 + r1 + r2)

   print(f"Recombination frequency is : {rf}")

recom_freq()

#CALCULATE ALLELE FREQUENCY
def allele_freq():

 pp = int(input("Enter number of individuals with allele PP"))
 pq = int(input("Enter number of individuals with allele PQ"))
 qq = int(input("Enter number of individuals with allele QQ"))

 pfreq = ((2*pp)+pq)/(2*(pp+pq+qq))
 qfreq = ((2*qq)+pq)/(2*(pp+pq+qq))

 print (f"Frequency of allele P is {pfreq}")
 print (f"Frequency of allele Q is {qfreq}")

allele_freq()

#FIND POSITION OF START AND STOP CODON

code = input("please enter the sequence of codons :")

start = code.find('AUG') #returns index of the substring
termination = code.find('UGA')

print(f"the start codon is at position {start} and the stop codon is at position {termination}")

orf = code[start:(termination+3)] #we have to print stop codon as part of orf

#TRANSCRIPTION

code = input("please enter the sequence of codons :")
newcode = code.upper()
splitcode = newcode[4:10]
rna = splitcode.replace("t","u")
print (rna)

#DNA OR RNA

code = input("please enter the sequence :  ")

if "U" in code:
  print ("This is a RNA sequence")

elif "T" in code:
  print ("This is a DNA sequence")

else:
  print ("This is a wrong sequence")

#POLYPEPTIDE

code = input("please enter the sequence :  ")

if "AUG" in code:
  print ("This sequence contains methionine")
  if len(code)>=29:
    print ("This sequence has over 10 AAs")
  else:
    print ("This sequence has less than 10 AAs")

if "AUG" not in code:
  print ("This sequence does not contain methionine")
  if len(code)>=29:
    print ("This sequence has over 10 AAs")
  else:
    print ("This sequence has less than 10 AAs")

#PAIRING
while True:
  code = input("please enter the sequence :  ")
  if code in "ATGC":
    if code == "A":
      print("A PAIRS WITH T")
      break
    elif code == "T":
      print("T PAIRS WITH A")
      break
    elif code == "G":
      print("G PAIRS WITH C")
      break
    elif code == "C":
      print("C PAIRS WITH G")
      break
  else:
    print("incorrect base,please retry")
    continue

#HAMMING DISTANCE:

code1= input("enter first sequence")
code2= input("enter second sequence")
count = 0
length = len(code1)

for num in range(0,length):
  if code1[num]!= code2[num]:
    count+=1

print(f"hamming distance is : {count}")

per = (count/length)*100

print(f"hamming distance in percentage is : {per}")

#TRANSCRIPTION WITH TEMPLATE STRAND

code = input("enter DNA sequence : ")
rna = ""

for item in code:
  if item == "T":
    rna+="U"
  else:
    rna+=item

print(rna)

#TRANSCRIPTION WITH CODING STRAND

code = input("enter DNA sequence : ")
trans = {"A":"U","G":"C","C":"G","T":"A"}
rna = ""

for item in code:
  rna+= trans[item]

print (rna)

#POLYPEPTIDE CHAIN

translation = {'UUU':'Phe','UUC':'Phe',
               'UUA':'Leu','UUG':'Leu','CUU':'Leu','CUC':'Leu','CUA':'Leu','CUG':'Leu',
               'AUU':'Ile','AUC':'Ile','AUA':'Ile',
               'AUG':'Met',
               'GUU':'Val','GUC':'Val','GUA':'Val','GUG':'Val',
               'UCU':'Ser','UCC':'Ser','UCA':'Ser','UCG':'Ser',
               'CCU':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
               'ACU':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
               'GCU':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
               'UAU':'Tyr','UAC':'Tyr',
               'CAU':'His','CAC':'His',
               'CAA':'Gln','CAG':'Gln',
               'AAU':'Asn','AAC':'Asn',
               'AAA':'Lys','AAG':'Lys',
               'GAU':'Asp','GAC':'Asp',
               'GAA':'Glu','GAG':'Glu',
               'UGU':'Cys','UGC':'Cys',
               'UGG':'Trp',
               'CGU':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
               'AGU':'Ser','AGC':'Ser',
               'AGA':'Arg','AGG':'Arg',
               'GGU':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly',
               'UAA':'XXX','UAG':'XXX','UGA':'XXX'}

RNA = input("enter rna sequence")

peptide = ""

start = 0
end = 3

while start < len(RNA)-3:
  peptide = peptide + translation[RNA[start:end]] + '-'
  start += 3
  end += 3

#CHALLENGE - ORFS
#test code : GCAAUGCCUUUCUUACCUCCUAACCCUUGACCUUAUACAUGAUUUAUGUUUAUGCCUCCUCCCGGGAACUAACACGCAGCGUAG


code = input("enter rna sequence : ")
start = 0
end = 3
aug = ""
while start<len(code)-3:

  aug = code[start:end]
  start +=3
  end +=3

  if aug == "AUG":
    orf = "AUG"

    start1 = start
    end1 = end
    stop = ""
    while start1<len(code)-3:

      orf = orf + code[start1:end1]
      stop = code[start1:end1]
      start1 +=3
      end1 +=3


      if stop == "UAA" or stop == "UGA" or stop == "UAG":
        print (f"ORF: {orf}")
        orf = ""
        stop = ""
        start=start1
        end=end1
        break
      else:
        continue
  else:
     continue

#CALCULATION OF % OF A BASE IN DNA

def percentage(dna):
  A = 0
  T=0
  G=0
  C=0
  len(dna)=length

  for item in dna:
    if item == "A":
      A+=1
    elif item == "T":
      T+=1
    elif item == "G":
      G+=1
    elif item == "C":
      C+=1



  print("percentage of A is : " + (A/length)*100)
  print("percentage of T is : " + (T/length)*100)
  print("percentage of G is : " + (G/length)*100)
  print("percentage of C is : " + (C/length)*100)

#RECOMB FREQ

def recom(nonrec1,nonrec2,rec1,rec2):

  sum = nonrec1+nonrec2+rec1+rec2

  freq = (rec1+rec2)/sum

  return freq

#HARDY WEIRBERG

def population (PP,PQ,QQ,sum):

  Pfreq = (2*PP + PQ)/(2*sum)
  Qfreq = (2*QQ + PQ)/(2*sum)

  print (f"Frequency of P is : {Pfreq}")
  print (f"Frequency of Q is : {Qfreq}")

  expPP = (Pfreq**2)*sum
  expQQ = (Qfreq**2)*sum
  expPQ = (Pfreq*Qfreq)*sum*2

  print (f"expected PP is : {expPP}")
  print (f"expected PQ is : {expPQ}")
  print (f"expected QQ is : {expQQ}")

  t1 = ((PP - expPP)**2)/expPP
  t2 = ((QQ - expQQ)**2)/expQQ
  t3 = ((PQ - expPQ)**2)/expPQ
  chisq= t1+t2+t3

  print(f"chi square is equal to :  {chisq}")

population(1200,745,555,2500)

#TRY EXP FINALLY

def validInt(question):
  while True:
    try:
      n = int(input(question))
    except:
      print('Invalid data type! Type an integer number!')
    else:
      break
  return n


nonrecgen1 = validInt('How many individuals with non-recombinant genotype 1: ')
nonrecgen2 = validInt('How many individuals with non-recombinant genotype 2: ')

recgen1 = validInt('How many individuals with recombinant genotype 1: ')
recgen2 = validInt('How many individuals with recombinant genotype 2: ')

#MATH MODULE

import math
math.sqrt(x)
math.log(number,base)
math.e
math.pi
math.sin(x)
math.cos(x)
math.tan(x)
math.factorial(x)
math.floor(x) #roundoff lower
math.ceil(x)#roundoff higher

#RANDOM MODULE

import random
random.random()
random.randint(lowerlimit,upperlimit)
random.randint(lowerlim,upperlimit,step)
random.choice(sometuple)

#TIME MODULE

import time as tm

before = tm.time()#to calculate time before execution

#some code here

after = tm.time()#to calculate time after execution

interval = before - after

tm.sleep(somedelayincode)

#GENERATING A RANDOM DNA

length = int(input("enter length of the sequence : "))
dnapair = ("A","T","G","C")
dna = ""
num = 0
import random
while num<length:
  dna+=random.choice(dnapair)
  num+=1

print(dna)

#TIMING A RANDOM DNA

length = int(input("enter length of the sequence : "))
dnapair = ("A","T","G","C")
dna = ""
num = 0
import time
import random

before = time.time()
while num<length:
  dna+=random.choice(dnapair)
  num+=1
after = time.time()

interval =  after - before

print(interval)

#CROSSING GENOTYPES :

parent1 = input("enter genotype of the first parent")
parent2 = input("enter genotype of the second parent")
offspring = int(input("enter the number of offsprings"))

par1=[]
par2 =[]
num=0
off={}
offlist=[]
import random

for item in parent1:
  par1.append(item)

for item in parent2:
  par2.append(item)

while num < offspring:
  off = random.choice(par1)+random.choice(par2)
  offlist.append(off)
  print (off)
  off= ""
  num+=1

#READING AND WRITING ON TEXT FILES

with open('/content/python_text.txt','r') as python_text: #upload file on drive and then input file location
  for line in python_text:
    print(line) #print the file

with open('/content/python_text.txt','r') as python_text:
  lines = python_text.readlines()
  print(lines) #prints every line with a /n as a list

with open('python_phrase.txt','w') as python_phrase:
  python_phrase.write('''Welcome to my Python course
  I hope you are en joying it!
  ''') #overwrites a file

#READING AND WRITING ON FASTA FILES

with open('/content/sequence (14).fasta') as viral_seq:
  sequence = []
  count = 0 #to remove the first header line
  for line in viral_seq:
    if count > 0:
      sequence.append(line)

sequence2 = ''.join(sequence)

with open('sequence.txt','w') as sequence: #CREATING A TEXT FILE
  sequence.write(sequence2)

#SARS COV 2 TRANSCRIBE

with open('/content/sequence.fasta copy','r') as viralseq:
  sequence = []
  count = 0
  for line in viralseq:
    if count > 0:
      sequence.append(line)
    count+=1

sarscov = ''.join(sequence)

dna_to_rna = {'A':'U','T':'A','C':'G','G':'C'}
rnaSeq = ''
for character in sarscov:
  if character in dna_to_rna.keys():
    rnaSeq += dna_to_rna[character]
  else:
    rnaSeq += character

print(rnaSeq)

with open('rnaSeq.txt','w') as rna:
  rna.write(rnaSeq)

#FINAL CHALLENGE : USER CHOICE WITH SARS COV


def transcribe():

  dna_to_rna = {'A':'U','T':'A','C':'G','G':'C'}
  rnaSeq = ''
  for character in sarscov:
     if character in dna_to_rna.keys():
        rnaSeq += dna_to_rna[character]
     else:
         rnaSeq += character

  print(rnaSeq)

  with open('rnaSeq.txt','w') as rna:
     rna.write(rnaSeq)


def dnacomp():

  dna_to_rna = {'A':'T','T':'A','C':'G','G':'C'}
  comp = ''
  for character in sarscov:
     if character in dna_to_rna.keys():
        comp += dna_to_rna[character]
     else:
         comp += character

  print(comp)

  with open('comp.txt','w') as compdna:
     compdna.write(comp)


def ATcont():


  at=0
  for character in sarscov:
    if character == "A" or character == "T":
      at+=1
  gc=0
  for character in sarscov:
    if character == "C" or character == "G":
      gc+=1


  atcont = (at/(at+gc))*100

  print(f"The A-T content of the sequence is : {atcont}")


def GCcont():

  length = len(sarscov)

  at=0
  for character in sarscov:
    if character == "A" or character == "T":
      at+=1
  gc=0
  for character in sarscov:
    if character == "C" or character == "G":
      gc+=1

  gccont = (gc/(gc+at))*100

  print(f"The A-T content of the sequence is : {gccont}")

sequence = input('Path of the file: ')

with open(sequence,'r') as viralseq:
  seq = []
  count = 0
  for line in viralseq:
    if count > 0:
      seq.append(line)
    count+=1

sarscov = ''.join(seq)

print ("Please choose one of the following options :")
print ("To transcribe the sequence press 1 ")
print ("To get the complementary dna sequence press 2")
print ("To find the A-T content press 3")
print("To find the G-C content press 4")

while True :
    try:
      choice = int(input("Please enter your choice : "))
    except :
      print("Please choose an option : 1,2,3 or 4 ")
    else:
      if choice == 1 :
        transcribe()
        break
      if choice ==2:
        dnacomp()
        break
      if choice == 3:
        ATcont()
        break
      if choice==4:
        GCcont()
        break
