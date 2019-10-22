import os, re, argparse, sys

#input for the desired Nxx/Lxx
Nxx = argparse.ArgumentParser(prog = 'Contig_Analysis_Nxx',
                              usage = 'python3 Nxx_calc.py --xx [xx] --input [input dir] ',
                              description = 'Contig Nxx/Lxx and GC-content calculator. Simply fill in your desired xx '
                                          'for your test (so xx = 50 for N50) and let the script take care of the rest')
Nxx.add_argument('--xx', "-x",
                 action = 'store',
                 help = 'specify the Nxx to calculate (only the number is required)',
                 default = 50,
                 required = False)
Nxx.add_argument('--input', "-i",
                 action = 'store',
                 help = 'give the path to the input',
                 required = True)
Nxx.add_argument('--output', "-o",
                 action = 'store',
                 help = 'specify the output directory',
                 required = True)
Nxx.add_argument('--version',
                 action = 'version',
                 version = 'Nxx_calc 1.2')
args = Nxx.parse_args()
xx = int(args.xx)
input = args.input
print('Input folder: ' + input)
print('Chosen Nxx: N' + str(xx))

def count_contig(file):
    doc = open(str(input + '/' + file), 'r')
    fasta = list(doc.readlines())
    #selecting the contigs
    con = list(filter(lambda x: not '>' in x, fasta))
    #counting the amount of elements within the list of contig annotations
    contigs = len(con)
    print('You have ' + str(contigs) + ' contigs in ' + file)

def Nxx_calc(file, xx):
    doc = open(str(input + '/' + file), 'r')
    fasta = list(doc.readlines())
    #filtering out the contig annotations
    seq = list(filter(lambda x: not '>' in x, fasta))
    #determining the total nucleotides
    seq.sort(key=len, reverse=True)
    seq_str = str(seq).strip('[]')
    seq_len = len(seq_str)
    print(str(file) + ' has a total of '+ str(seq_len) + ' nucleotides')
    print('The longest contig is ' + str(len(seq[0])) + ' bp and the shortest is ' + str(len(seq[-1])) + ' bp')
    #calculating the GC content
    G = seq_str.count('G')
    C = seq_str.count('C')
    total_len = len(seq_str)
    GC = round((G+C)/total_len*100,2)
    print('The GC content of ' + file + ' is ' + str(GC) + '%')
    #calculating the N50
    print('Calculating the N' + str(xx) + ' and L' + str(xx) + '; this may take a while, depending on your hardware, '
    'file size and specified xx')
    threshold = int(seq_len)*(xx/100)
    contigs = []
    var1 = int()
    #here the actual calculations are called and the output is stored in a list
    step1 = list(plus_100(contigs,seq,threshold,var1))
    contigs1 = step1[0]
    var1_1 = step1[1]
    step2 = list(min_10(contigs1,threshold,var1_1))
    contigs2 = step2[0]
    var1_2 = step2[1]
    step3 = list(plus_1(contigs2,seq,threshold,var1_2))
    contigs3 = step3[0]
    Lxx = step3[1]
    Nxx = len(str(contigs3[-1]))
    print('The N' + str(xx) + ' is ' + str(Nxx) + ' and the L' + str(xx) + ' is ' + str(Lxx))
    f.write(str(file) + '\t' + str(Nxx) + '\t' + str(Lxx) + '\t' + str(GC) + '\n')

def plus_100(contigs,seq,threshold,var1):
    var0 = 0
    var1 += 100
    while int(len(str(contigs).strip('[]'))) < threshold:
        for i in range(var0,var1):
            contigs.append(str(seq[i]).strip('[]'))
        if int(len(str(contigs).strip('[]'))) > threshold:
            return contigs, var1
        var0 += 100
        var1 += 100

def min_10(contigs,threshold,var1):
    while int(len(str(contigs).strip('[]'))) > threshold:
        del contigs[-10:]
        var1 -= 10
        if int(len(str(contigs).strip('[]'))) < threshold:
            return contigs,var1

def plus_1(contigs,seq,threshold,var1):
    while int(len(str(contigs).strip('[]'))) < threshold:
        contigs.append(seq[var1])
        var1 += 1
        if int(len(str(contigs).strip('[]'))) > threshold:
            return contigs,var1

#finding your files and calling the calculating functions
#print('This Nxx/Lxx calculator is versatile. It can calculate the N1 to N100. Please specify what you want to calculate'
#' below.')
#xx = int(input('Enter a number between 1 and 100 for xx (only one!) and press enter \n'))
print('Creating report file')

# is the input a directory?
if os.path.isdir(input):
    # making lists out of the input files
    files = sorted(os.listdir(input))
    new_files = [] # with the full path attached
    for file in files:
        if re.search(".+\.fasta", file) or re.search(".+\.fa", file):
            file = os.path.join(input, file)
            new_files.append(file)
    # writing stuff if correct
    if new_files:
        with open('contig_report.txt','w') as f:
            f.write('Filename\tN' + args.xx + '\tL' + args.xx + '\t' + 'GC content' + '\n')
            print('Fasta file(s) detected:' + ''.join(new_files))
            for file in new_files:
                print('Processing ' + str(file))
                count_contig(file)
                Nxx_calc(file, xx)
    else: # exit when not correct
        print("No fasta files detected, please make sure your input files are in fasta format")
        sys.exit()

# or is the input a file?
else:
    if re.search(".+\.fasta", input) or re.search(".+\.fa", input):
        with open('contig_report.txt','w') as f:
            f.write('Filename\tN' + args.xx + '\tL' + args.xx + '\t' + 'GC content' + '\n')
            print("Processing" + input)
            count_contig(input)
            Nxx_calc(input, xx)
    else:
        print("Please make sure your input file is in fasta format")
        sys.exit()

print('A report can be found in your running directory')