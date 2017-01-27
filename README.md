ThIEF: Finding Genome-wide Trajectories of Epigenetics Marks

INSTALL:
- you need to have python
- you need to have gpsol (part of the GLPK, for mac os install via homebrew). 

INPUTS:
Input format is a set of comma-separated-value text files using the "," delimiter (if you need to, you can add other formats in function ReadNucTable of the process.py script). Each line corresponds to a nucleosome (or other genetic feature). The first column is the genomic location. Any other column is ignored at this point (but could be used in the future). All files should have the same number of columns. 

OUTPUT:
The output format for the THiEF is table, where each line is the resulting track of features (note that if you have multi column input files, i.e. < x1 a1 b1> and <y2 p2 q2>, then in the output the line corresponding to corresponding track will look like < x1 a1 b1 y2 p2 q2>). In case the track has "empty spots" (gaps), they are marked as "-1" (in multicolumn case it will be < -1 -1 -1>).

USAGE:
As said, the input files (named *.track in the examples provided) should contain a list of nucleosome positions as the first column. It may contain other information such as the start and end positions of the nucleosomes but those information are not used for building the trajectories. See input1.track, input2.track and input3.track for examples

Given the input files (*.track) inside the THiEF directory, the trajectory is build using the following command

python process.py < output_file_name> < input1.track> < input2.track> < input3.track>

OR

python process.py < output_file_name> *.track

The output file should contain trajectories build by THiEF. A “-1.00" in the output file represents a gap in a trajectory. See output_trajectories for example.

If you have questions please contact Abid mhasa006@ucr.edu, Anton polishka@cs.ucr.edu, or Stefano stelo@cs.ucr.edu
