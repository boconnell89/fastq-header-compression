#!/home/groups/oroaklab/src/conda/miniconda3/bin/python3

#from permute import permute ### Helper library, copied below for simplicity's sake.
from gzip import open as gopen
import argparse, sys
from time import time
'''
Nearly all of our fastq data has a custom read header format, where most of the read name is replaced with a DNA sequence composed of multiple barcodes (e.g, sci-ATAC or 10X barcodes).  Since a typical read appears as follows:

@TCGACAGGACAATAGAATCAATAGACATCATAATGTG:num=7:tag=NS500556:522:H3HTGAFX3:1:11101:9799:1006#0/1
TCGTCGGCAGCGTCAGATGTGTAT
+
EEEEEEAEEEEEEEEEEE/EEAEE

It should be pretty obvious that there are some major space saving opportunities to be realized by encoding/compressing the header.  Furthermore, the SAM/BAM file format uses unicode characters for read names.  Since most of our datasets are composed of multiple NovaSeq runs and easily range in the 0.5TB+ just for the gzipped fastqs per experiment, any savings are worth the effort.

I considered a number of different methods, but settled on the this one for several reasons, including general compatibility with our pipelines (all the NGS analysis software I've tested this with is fine with the 256 bit unicode characters. I specifically wanted to avoid using an external key file or similar, despite the potential speed benefit, since I wanted this to be more broadly useable.

The method is as follows:
1. Split the read header on ':', keeping the barcode header, the 'num' interger, and the '#0/1' or '#0/2', since that is recognized by BWA as the R1/R2 identifier.  Ditch the original read tag, since we never use it, and we have a unique identifier anyway.
2. Create a dict mapping all 4-base kmers to unicode characters.  This includes 'N' as a padding base.  As part of demultiplexing, we speciafically correct 'N' bases to the original barcode sequence, so there is no risk of conflict.
3. Encode and compress the header sequence. In the above example, the end result would be a header like:

@{o&Gc4ÇÒRü:7#0/1
TCGTCGGCAGCGTCAGATGTGTAT
+
EEEEEEAEEEEEEEEEEE/EEAEE

In other words, we just turned a 93 character header into a 17 character header (including the initial '@'). That's 52% space saving.
In testing, this is born out for uncompressed data-types, while you see a more modest ~20% space saving on .gzip'd files and BAM formatted files.

4. This also allows for header decompression (with the 'D' flag), just in case one wants to return to a more human-readable version.

5. Future steps: write an extension for SAM/BAM format if people request it.  That's pretty straightforward with the pysam API.
'''


def encodeHeader(header, encodeDict, l):
	"""
	Encodes a fastq header 
	"""
	D=header.split(':')
	BC=D[0][1:]
	BC=BC+'N'*(l-(len(BC)%l))
	comp=''
	for i in range(int(len(BC)/l)):
		comp+=encodeDict[BC[i*l:(i*l)+l]]
	return('@'+comp+':'+D[1].strip('num=')+'#0/'+D[-1][-1])	

def decodeHeader(compHeader, decodeDict):
	"""
	Decodes a header based on the same sceme as encodeHeader
	"""
	deComp=''
	barcode=compHeader.split(':')[0][1:]
	count=compHeader.split(':')[1]
	for entry in compHeader:
		deComp+=decodeDict[entry]
	return('@'+deComp.strip('N')+':num='+count[:-1]+'#0/'+count[-1])
	
def compress_fastq(r1, r2, encodeDict):
	"""
	Compresses paired fastq file headers
	"""
	if r1.endswith('.gz'):
		h1 = gopen(r1, 'rt')
	else:
		h1 = open(r1, 'rU')
	if r2.endswith('.gz'):
		h2 = gopen(r2, 'rt')
	else:
		h2 = open(r2, 'rU')
	line_data = ['','','','','','']
	for line in h1:
		# print (h2.readline().split())
		r2L = h2.readline().strip('\n')
		if line_data[0] == '' and line.startswith('@'):
			line_data[0] = encodeHeader(line.strip('\n'), encodeDict, 4)
			line_data[3] = encodeHeader(r2L, encodeDict, 4)
		elif line_data[1] == '': ## Sequence empty
			line_data[1] = line.strip('\n')
			line_data[4] = r2L
		elif line_data[2] == '' and not line.startswith('+'):
			line_data[2] = line.strip('\n')
			line_data[5] = r2L
			ld = line_data
			line_data = ['','','','','','']
			yield ld
	h1.close()
	h2.close()
	
def decompress_fastq(r1, r2, decodeDict):
	"""
	Decompresses paired fastq headers
	"""
	if r1.endswith('.gz'):
		h1 = gopen(r1, 'rt')
	else:
		h1 = open(r1, 'rU')
	if r2.endswith('.gz'):
		h2 = gopen(r2, 'rt')
	else:
		h2 = open(r2, 'rU')
	line_data = ['','','','','','']
	for line in h1:
		# print (h2.readline().split())
		r2L = h2.readline().strip('\n')
		if line_data[0] == '' and line.startswith('@'):
			line_data[0] = decodeHeader(line.strip('\n'), decodeDict)
			line_data[3] = decodeHeader(r2L, decodeDict)
		elif line_data[1] == '': ## Sequence empty
			line_data[1] = line.strip('\n')
			line_data[4] = r2L
		elif line_data[2] == '' and not line.startswith('+'):
			line_data[2] = line.strip('\n')
			line_data[5] = r2L
			ld = line_data
			line_data = ['','','','','','']
			yield ld
	h1.close()
	h2.close()

def make_dicts():
"""
Generates all possible 4-mers, then uses the index as the basis for the character encode (index + 32).
"""
	tmp=sorted(permute(['GGGG'], 4)) # This serves to generate all possible 4-mers

	### Do some logic to remove the sequences where N comes before non N:
	### While this is slightly hack-y, it is repeatable, robust, and takes about a millisecond.

	j=1
	while j>0:
		j= remove_incorrect_entries(tmp)

	for i in tmp:
		if 'NA' in i or 'NT' in i or 'NC' in i or 'NG' in i:
			tmp.remove(i)
	### Make the dict
	encodeDict={}
	decodeDict={}
	for i,j in enumerate(tmp):
		encodeDict[j]=chr(i+32)
		decodeDict[chr(i+32)]=j
		
	return(encodeDict, decodeDict)
	
	
### Note: the following is part of a set of helper libraries I've been writing. Copy pasted here for the sake of convenience

def permute(seqs, h):
	"""
	Just a module to do hamming distance permutations of sequences. Takes a list of sequences (if one sequence only, use a length 1 list). Returns a set. Best used in a for loop, e.g.:
		for permutation in permute(["SOME LIST"], 2):
			<do something>
	"""
	basedict = {'A': ['G','T','C','N'], 'T':['G','A','C','N'] ,'C':['G','A','T','N'], 'G':['A','T','C','N'], 'N': ['G','A','T','C' ]}
	h-=1
	# print h, seqs
	perms = set()
	for seq in seqs:
		# print seq
		for i,b in enumerate(seq):
			# print i, b
			tmp = list(seq)
			for p in basedict[b]:
				tmp[i]=p
				perms.add( "".join(tmp))
	# print perms, h
	if h!=0:
		return permute(perms.union(seqs), h)
	else:
		return perms.union(seqs)
		
		
def remove_incorrect_entries(tmp):
	"""
	Takes the all 4-mers and removes cases where N comes before another base.  Since N's are only used as padding, these are unnecessary.
	"""
	ct = 0
	for i in tmp:
		if 'NA' in i or 'NT' in i or 'NC' in i or 'NG' in i:
			tmp.remove(i)
			ct+=1
	return ct
	
	
def main(argv):
	encodeDict, decodeDict=make_dicts()
	options = parse_args(sys.argv)
	
	ct=0
	if options.E is True:
		sys.stderr.write('Endoding headers... ')
		files = [gopen('{}/{}.c.R1.fq.gz'.format(options.directory, options.output), 'wt'),gopen('{}/{}.c.R2.fq.gz'.format(options.directory,options.output), 'wt')]
		for r in compress_fastq(options.r1, options.r2, encodeDict):
			ct+=1
			files[0].write('{}\n{}\n+\n{}\n'.format(r[0], r[1],r[2]))
			files[1].write('{}\n{}\n+\n{}\n'.format(r[3], r[4],r[5]))
	elif options.D is True:
		sys.stderr.write('Decoding headers... ')
		files = [gopen('{}/{}.d.R1.fq.gz'.format(options.directory, options.output), 'wt'),gopen('{}/{}.d.R2.fq.gz'.format(options.directory,options.output), 'wt')]
		for r in decompress_fastq(options.r1, options.r2, decodeDict):
			ct+=1
			files[0].write('{}\n{}\n+\n{}\n'.format(r[0], r[1],r[2]))
			files[1].write('{}\n{}\n+\n{}\n'.format(r[3], r[4],r[5]))
"""
Boilerplate command line handling
"""
def parse_args(args):
	parser = argparse.ArgumentParser(description='''Encodes or decodes paired fastq file readnames''', epilog = "")
	parser.add_argument('-r1', '-R1', type=str, help='Read1 fastq file')
	parser.add_argument('-r2', '-R2', type=str, help='Read2 fastq file ')
	parser.add_argument('-output' , nargs='?', type=str,help='output prefix.')
	parser.add_argument('-directory' , nargs='?', type=str, default='.', help='path to output')
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('-E',action="store_true", help='Encode Headers')
	group.add_argument('-D',action="store_true", help='Decode Headers')
	# parser.add_argument('-v', action="store_true", help='suppress verbose output [default off]')
	return parser.parse_args()			
start = time()


if __name__ == "__main__":
	main(sys.argv);
	stop = time()
	m, s = divmod(stop-start, 60)
	h, m = divmod(m, 60)
	sys.stderr.write('Processing completed in {} hours, {} minutes, {} seconds\n'.format(int(h),int(m),int(s)))
	raise SystemExit