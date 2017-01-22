In order to have THiEF results you need to have gpsol intalled on your machine (part of the GLPK, for mac os intall via homebrew). If you don't need THiEF results comment lines #99-143 in "process.py"


To generate a track file simply run
	>python process.py <filename> [<linker size> <STD of the spread>]
linker size - is the distance between locations to track
STD of the spread - is standard deviation of the jump step for location from one layer to the other

To generate a set of tracks with different parameters check Start.sh

The track file "<filename>_3track" would be put in OUTPUT folder
	if you would like to have tracks consisting of 4 or 5 layers uncomment lines #87-97. In this case corresponding files "<filename>_4[5]track" would be placed in OUTPUT folder

The format of the track file is 
	each line is a track "x1 y1 x2 y2 x3 y3 [x4 y4] [x5 y5]"
	every odd column x<i> is a location coordinate
	every even column y<i> is a boolen flag whether that location is part of the track or not (relevant for THiEF tool - only locations with TRUE flag are considered to be present)

Output format for the THiEF is the same as above, the only difference is "-1 -1" marks an empty spot on the track. This means that true track "... X 0 ..." will correspond to  "... -1 -1 ..." recovered by THiEF.

any questions send to: polishka@cs.ucr.edu

