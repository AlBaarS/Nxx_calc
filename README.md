# Nxx_calc
Nxx calculator (N50, N80, N90, etc.) which also gives the respective Lxx and GC content of contigs.

### Info
Nxx_calc is a python 3 based tool to calculate the Nxx, Lxx and GC content of your contigs. Simply fill in your desired xx for your test (so xx = 50 for N50) and let the script take care of the rest. Currently, the script is built for linux command line, but I may upload a version to use on a windows desktop, with a very simple command interface.

### Prerequisites
Python 3

### Getting started
Install the script as follows: 
'''
git clone https://github.com/AlBaarS/Nxx_calc
'''
That's all you need to do. You can use the script as described below

### Usage
The script is used as follows: 'python3 Nxx_calc [args]'

#### Required args:
Use '''--input''' or '''-i''' [path] to specify an input directory or path to a file.
Use '''--output''' or '''-o''' [path] to specify the output directory, where your report 'contig_report.txt' will appear.

#### Optional args:
Use '''--xx''' or '''-x''' [int] you can specify your desired xx (--xx 50 will give N50/L50).
Use '''--GC''' or '''-gc''' to specify whether you want the GC content to be calculated. Since this is now the slowest feature, I decided to make it toggleable.
Use '''--array''' or '''-a''' to run in array mode. This will create an array per file with the N10-N90, L10-L90.
 * Can be combined with the -gc flag.
 * Can NOT be combined with -x. Or do. I have no idea what will happen. Maybe a glitch will occur on the matrix. Who knows. I don't recommend it. But we live in a free world, so do what you desire. Just don't blame me if the FBI is at your doorstep. Thanks in advance.
Use '''--length_distribution''' or '''-ld''' to generate a length distribution plot.
Use '''--array_plot''' or '''-ap''' to generate an array plot, which makes a bar graph with the several Nxx's and their respective length.
The flags '''--help''' and '''--version''' will give you the list of possible flags and the version, respectively.
  
### Disclaimer
This tool is a highly experimental python script made by a biology/bioinformatics student that wanted to refresh their python after a summer break. I cannot guarantee high speed or efficient use of computational power. Feel free to suggest improvements and submit feedback. Feel free to use the code in any way you want. 

### To-do list (things that MAY appear in future versions)
- Support for multiple files without the need for a folder
- More professional installer