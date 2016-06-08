import sys, os
import argparse
from subprocess import check_call
import numpy as np

USEARCH = "usearch"


def readFasta(fastafile):

  with open(fastafile, 'r') as fp:
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            line=line[1:]
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

class readable_file(argparse.Action):
  def __call__(self,parser, namespace, values, option_string=None):
      prospective_file=values
      if not os.path.isfile(prospective_file):
          raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid file path".format(prospective_file))
      if os.access(prospective_file, os.R_OK):
          setattr(namespace,self.dest,prospective_file)
      else:
          raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable file".format(prospective_file))

class writeable_dir(argparse.Action):
  def __call__(self,parser, namespace, values, option_string=None):
      prospective_dir=values
      if not os.path.isdir(prospective_dir):
          raise argparse.ArgumentTypeError("writeable_dir:{0} is not a valid path".format(prospective_dir))
      if os.access(prospective_dir, os.W_OK):
          setattr(namespace,self.dest,prospective_dir)
      else:
          raise argparse.ArgumentTypeError("writeable_dir:{0} is not a writeable dir".format(prospective_dir))

class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def checkInputAndClean(readfile, outdir, verbose):
  #check file is in fasta format and that each read has a sample identifier.

  if verbose:
    print "setting up fasta file..."

  read_set = set()

  outputfile = outdir + os.path.splitext(os.path.basename(readfile))[0] + "_renamed.fasta"
  with open(outputfile, 'w') as outfile:
    for h,s in readFasta(readfile):

      r = h.split(";")[0]
      if r in read_set:
        MyError("Reads may be duplicated!")
      read_set.add(r)

      if "sample=" in h:
        #is already in usearch format.
        if "size=" in h:
          #want to remove the original size information
          h = h.split(";")
          h = [t for t in h if "size" not in t]
          h = ";".join(h)
      elif "." in h:
        #Assume the string prior to the first period (.) is the sample identifier (Rask format)
        sample = h.split(".")[0]
        h = h + ";sample=" + sample
      else:
        MyError("Fastafile is neither usearch or rask like.")

      outfile.write(">"+h+"\n"+s+"\n")

  return outputfile

def cluster(readfile, per_id, cpu, verbose):

  file_preface = os.path.splitext(readfile)[0]

  if verbose:
    print "clustering with usearch..."

  usearch_cmd = (USEARCH
    + " -derep_prefix " + readfile
    + " -fastaout " + file_preface + "_unqiues.fasta"
    + " -sizeout"
    + " -threads " + str(cpu))

  if verbose:
    print "running... ", usearch_cmd

  check_call(usearch_cmd, shell=True)

  usearch_cmd = (USEARCH
      + " -cluster_fast"
      + " " + file_preface + "_unqiues.fasta"
      + " -centroids " + file_preface + "_centroids.fasta"
      + " -sort size"
      + " -id " + str(per_id)
      + " -threads " + str(cpu))

  if verbose:
    print "running... ", usearch_cmd

  check_call(usearch_cmd, shell=True)

  usearch_cmd = (USEARCH
    + " -usearch_global " + readfile
    + " -db " + file_preface + "_centroids.fasta"
    + " -strand both"
    + " -id " + str(per_id)
    + " -otutabout " + file_preface + "_otuTable.txt"
    + " -blast6out " + file_preface + "_blast6out.txt"
    + " -threads " + str(cpu))

  if verbose:
    print "running... ", usearch_cmd

  check_call(usearch_cmd, shell=True)

  return file_preface + "_otuTable.txt"

def convertToBinaryMatrix(otufile, verbose):
  file_preface = os.path.splitext(otufile)[0]

  if verbose:
    print "converting to binary matrix..."

  with open(file_preface + "_binary.txt", 'w') as outfile:
    with open(otufile, 'r') as infile:
      outfile.write(infile.next())
      for line in infile:
        line=line.split()
        array = np.array(map(int, line[1:]))
        array[array>0] = 1
        outfile.write("\t".join([line[0]]
          + list(array.astype('str'))) + "\n")

  return file_preface + "_binary.txt"




def main():

  parser = argparse.ArgumentParser(description='Cluster cleaned DBLalpha sequence tags.')

  parser.add_argument('-o', '--outputDir', action=writeable_dir
    , dest='outputdir'
    , help="location of output directory. Will be created if it doesn't exist"
    , required=True)

  parser.add_argument('-r', '--read', dest='read', action=readable_file
    , help="location of fasta file containing sequences."
    , required=True)

  parser.add_argument('--perID', dest='perID', type=float, default=0.96
    , help="percentage ID threshold. (default=0.96)")

  parser.add_argument('--cpu', dest='cpu', type=int, default=1
    , help="number of cpus to use. (default=1)")

  parser.add_argument('--verbose', dest='verbose', action='store_true'
    , default=False
    , help='print verbose output (default=False)')

  args = parser.parse_args()

  #get full path names
  args.outputdir = os.path.abspath(args.outputdir) + "/"
  args.read = os.path.abspath(args.read)

  renamed_reads = checkInputAndClean(args.read, args.outputdir, args.verbose)

  otu_table = cluster(renamed_reads, args.perID, args.cpu, args.verbose)

  convertToBinaryMatrix(otu_table, args.verbose)

  return


if __name__ == '__main__':
  main()
