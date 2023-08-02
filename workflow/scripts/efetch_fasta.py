#!/usr/bin/env python3
import sys
from Bio import Entrez

#define email for entrez login
db           = "nuccore"
Entrez.email = "some_email@somedomain.com"

#load accessions from arguments
if len(sys.argv[1:]) > 1:
  accs = sys.argv[1:]
else: #load accesions from stdin
  accs = [ l.strip() for l in sys.stdin if l.strip() ]
#fetch
sys.stderr.write( "Fetching %s entries from GenBank: %s\n" % (len(accs), ", ".join(accs[:10])))
for i,acc in enumerate(accs):
  try:
    sys.stderr.write( " %9i %s          \r" % (i+1,acc))
    handle = Entrez.efetch(db=db, rettype="fasta", id=acc)
    #print output to stdout
    sys.stdout.write(handle.read())
  except:
    sys.stderr.write( "Error! Cannot fetch: %s        \n" % acc)