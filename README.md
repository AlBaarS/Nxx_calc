# Nxx_calc
Nxx calculator (N50, N80, N90, etc.) which also gives the respective Lxx and GC content of contigs.

### Info
Nxx_calc is a python 3 based tool to calculate the Nxx, Lxx and GC content of your contigs. Simply fill in your desired xx for your test (so xx = 50 for N50) and let the script take care of the rest. Currently, the script is built for linux command line, but I may upload a version to use on a windows desktop, with a very simple interface.

### Usage
The script is used as follows: 'python3 Nxx_calc --xx <xx> --input <input directory>'
Use '--xx' or '-x' [int] you can specify your desired xx (--xx 50 will give N50/L50).
Use '--input' or '-i' [path] to specify an input directory or path to a file.
Use '--output' or '-o' [path] to specify the output directory, where your report 'contig_report.txt' will appear.
Use '--GC' or '-gc' [True/False] to specify whether you want the GC content to be calculated. Since this is now the slowest feature, I decided to make it toggleable.
The flags '--help' and '--version' will give you the list of possible flags and the version, respectively.
  
### Disclaimer
This tool is a highly experimental python script made by a biology/bioinformatics student that wanted to refresh their python after a summer break. I cannot guarantee high speed or efficient use of computational power. Feel free to suggest improvements and submit feedback. Feel free to use the code in any way you want. 

### To-do list (things that MAY appear in future versions)
- support for gzipped files
- support for multiple files without the need for a folder
- Nxx-array: create a tab-delimited array of Nxx per file (N10-N20-N30-etc.)
- Graphs depicting length distribution
