# Nxx_calc
Nxx calculator (N50, N80, N90, etc.) which also gives the respective Lxx and GC content of contigs.

### Info
Nxx_calc is a python 3 based tool to calculate the Nxx, Lxx and GC content of your contigs. Simply fill in your desired xx for your test (so xx = 50 for N50) and let the script take care of the rest. Currently, the script is built for linux/mac command line, but I may upload a version to use on a windows desktop, with a very simple interface.

### Usage
The script is used as follows: 'python3 Nxx_calc --xx <xx> --input <input directory>'
Use '--xx' or '-x' you can specify your desired xx (--xx 50 will give N50/L50).
Use '--input' or '-i' to specify an input directory or path to a file.
As of now, the '--output' or '-o' function does nothing, but it will be added somewhere in the future, if I feel like it.
The flags '--help' and '--version' will give you the list of possible flags and the version, respectively.
  
### Disclaimer
This tool is a highly experimental python script made by a biology/bioinformatics student that wanted to refresh their python after a summer break. I cannot guarantee high speed or efficient usage of your computational power, but it does do the job. Feel free to use the code in any way you want. 
