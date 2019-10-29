import os, re, argparse, sys

#input for the desired Nxx/Lxx
Nxx = argparse.ArgumentParser(prog='Contig_Analysis_Nxx', usage='python3 Nxx_calc.py --xx [xx] --input [input file/dir] '
				'--output [output dir] --GC [True/False]',
                              description='Contig Nxx/Lxx and GC-content calculator. Simply fill in your desired xx '
                                          'for your test (so xx = 50 for N50) and let the script take care of the rest')
Nxx.add_argument('--xx', '-x',
		action='store',
		default=50,
		help='Specify the Nxx to calculate (only the number is required, default = 50 for N50)')
Nxx.add_argument('--input', '-i',
		action='store',
		help='Specify the input file or directory. Only accepts .fasta or .fa files.',
		required=True)
Nxx.add_argument('--version', '-v',
		action='version',
		version='Nxx_calc 1.3')
Nxx.add_argument('--output', "-o",
                 action ='store',
                 help ='Specify the output directory',
                 required=True)
Nxx.add_argument('--GC', '-gc',
		action='store',
		help='Calculate the GC value or not. Default = False',
		default=False)
namespace = Nxx.parse_args()
xx = int(namespace.xx)
input = str(namespace.input)
outdir = str(namespace.output)
gc = bool(namespace.GC)
if os.path.isfile(input):
	if input.endswith('.gz'):
		fname = input.split('/')
		print('Input file: ' + str(fname[-1]))
		print('Gunzipping file')
		os.system('gunzip ' + input)
		input = input.strip('.gz')
		fname2 = input.split('/')
		print('Gunzipping of ' + str(fname2[-1]) + ' complete')
	else:
		fname = input.split('/')
		print('Input file: ' + str(fname[-1]))
elif os.path.isdir(input):
	for files in os.walk(input):
		for file in files:
			file = str(file).strip('[]')
			if file.endswith('.gz'):
				print('Gunzipping files')
				os.system('gunzip ' + str(input) + str(file))
	print('Input folder: ' + str(os.path.realpath(input)))
else:
	print('Error: User provided a file/nonexistent directory. Please specify an (existing) directory')
if os.path.isdir(outdir):
	print('Output folder: ' + str(os.path.realpath(outdir)))
else:
	print('Error: User provided a file/nonexistent directory. Please specify an (existing) directory')
	sys.exit()
print('Chosen Nxx: N' + str(xx))

def count_contig(file):
	if os.path.isfile(input):
		doc = open(str(input), 'r')
		isfile = True
	elif os.path.isdir(input):
		doc = open(str(file), 'r')
		isfile = False
	else:
		print('Error: Specified path is no file or directory or does not exist')
	fasta = list(doc.readlines())
	#selecting the contigs
	con = list(filter(lambda x: not '>' in x, fasta))
	#counting the amount of elements within the list of contig annotations
	contigs = len(con)
	if isfile == True:
		fname = input.split('/')
		print('You have ' + str(contigs) + ' contigs in ' + str(fname[-1]))
	else:
		fname = file.split('/')
		print('You have ' + str(contigs) + ' contigs in ' + str(fname[-1]))

def Nxx_calc(file,xx):
	if os.path.isfile(input):
		doc = open(str(input), 'r')
		isfile = True
	elif os.path.isdir(input):
		doc = open(str(file), 'r')
		isfile = False
	else:
		print('Error: Specified path is no file or directory')
	fasta = list(doc.readlines())
	#filtering out the contig annotations
	seq = list(filter(lambda x: not '>' in x, fasta))
	#determining the total nucleotides
	seq.sort(key=len, reverse=True)
	seq_len = int()
	for s in seq:
		seq_len += len(str(s).strip('[]'))
	if isfile:
		fname = input.split('/')
		print(str(fname[-1]) + ' has a total of ' + str(seq_len) + ' nucleotides')
	else:
		fname = file.split('/')
		print(str(fname[-1]) + ' has a total of ' + str(seq_len) + ' nucleotides')
	print('The longest contig is ' + str(len(seq[0])) + ' bp and the shortest is ' + str(len(seq[-1])) + ' bp')
	#calculating the GC content
	if gc:
		seq_str = str(seq).strip('[]')
		G = seq_str.count('G')
		C = seq_str.count('C')
		total_len = len(seq_str)
		GC = round((G+C)/total_len*100,2)
		if isfile == True:
			fname = input.split('/')
			print('The GC content of ' + str(fname[-1]) + ' is ' + str(GC) + '%')
		else:
			fname = file.split('/')
			print('The GC content of ' + str(fname[-1]) + ' is ' + str(GC) + '%')
	else:
		GC = 'X'
	#calculating the Nxx
	print('Calculating the N' + str(xx) + ' and L' + str(xx) + ':')
	threshold = int(seq_len)*(xx/100)
	contigs = int()
	var1 = int()
	#here the actual calculations are called and the output is stored in a list
	step1 = list(plus_100(contigs,seq,threshold,var1))
	contigs1 = step1[0]
	var1_1 = step1[1]
	step2 = list(min_10(contigs1,seq,threshold,var1_1))
	contigs2 = step2[0]
	var1_2 = step2[1]
	step3 = list(plus_1(contigs2,seq,threshold,var1_2))
	Lxx = int(step3[1])
	Nxx = len(str(seq[Lxx]))
	print('The N' + str(xx) + ' is ' + str(Nxx) + ' and the L' + str(xx) + ' is ' + str(Lxx))
	if isfile == True:
		fname = input.split('/')
		f.write(str(fname[-1]) + '\t' + str(Nxx) + '\t' + str(Lxx) + '\t' + str(GC) + '\n')
	else:
		fname = file.split('/')
		f.write(str(fname[-1]) + '\t' + str(Nxx) + '\t' + str(Lxx) + '\t' + str(GC) + '\n')

def plus_100(contigs,seq,threshold,var1):
	var0 = 0
	var1 += 100
	while contigs < threshold:
		for i in range(var0,var1):
			contigs += len(str(seq[i]).strip('[]'))
		if contigs > threshold:
			return contigs, var1
		var0 += 100
		var1 += 100

def min_10(contigs,seq,threshold,var1):
	while contigs > threshold:
		for i in range((var1  - 10),var1):
			contigs -= len(str(seq[i]))
		var1 -= 10
		if contigs < threshold:
			return contigs,var1

def plus_1(contigs,seq,threshold,var1):
	while contigs < threshold:
		var1 += 1
		contigs += len(str(seq[var1]))
		if contigs > threshold:
			return contigs,var1

#finding your files and calling the calculating functions
print('Creating report file')
if os.path.isdir(input):
    # making lists out of the input files
	files = sorted(os.listdir(input))
	new_files = [] # with the full path attached
	for file in files:
		if file.endswith('.fasta') or file.endswith('.fa'):
			file = os.path.join(input, file)
			new_files.append(file)
		elif file.endswith('.gz'):
			file = os.path.join(input, file)
			print('Gunzipping ' + file)
			os.system('gunzip ' + file)
			nfile = file.strip('.gz')
			new_files.append(nfile)
    # writing stuff if correct
	if new_files:
		with open(str(outdir) + '/' + 'contig_report.txt','w') as f:
			f.write('Filename\tN' + str(xx) + '\tL' + str(xx) + '\t' + 'GC content' + '\n')
			for file in new_files:
				fname = str(file).split('/')
				print('File detected: ' + str(fname[-1]))
			for file in new_files:
				fname = str(file).split('/')
				print('Processing ' + str(fname[-1]))
				count_contig(file)
				Nxx_calc(file, xx)
	else: # exit when not correct
		print("No fasta files detected, please make sure your input files are in fasta format")
		sys.exit()

# or is the input a file?
else:
	if re.search(".+\.fasta", input) or re.search(".+\.fa", input):
		with open(str(outdir) + '/' + 'contig_report.txt','w') as f:
			f.write('Filename\tN' + str(xx) + '\tL' + str(xx) + '\t' + 'GC content' + '\n')
			fname = str(input).split('/')
			print("Processing " + str(fname[-1]))
			count_contig(input)
			Nxx_calc(input, xx)
	else:
		print("Please make sure your input file is in fasta format")
		sys.exit()
f.close()
print('A report called \'contig_report.txt\' can be found in your output directory')







