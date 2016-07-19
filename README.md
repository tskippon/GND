# GND
Matlab code for calculation of geometrically necessary dislocation densities from EBSD data

This code requires Matlab and MTEX

Here is a brief description of the files included in this project:


GND_auto.m 
	The main file containing the matlab function that should be called to calculate GND values.

doBurgersBCC.m
	A function used to set up the Burgers vector information for Body Centred Cubic materials.
	This gets called by GND_auto.
	
doBurgersFCC.m
	A function used to set up the Burgers vector information for Face Centred Cubic materials.
	This gets called by GND_auto.
	
doBurgersHCP.m
	A function used to set up the Burgers vector information for Hexagonal Close Packed materials.
	This gets called by GND_auto.
	
cpp_curvatures.cpp 
	A c++ source file containing the code that calculates the lattice curvatures for an EBSD dataset.
	This is called by GND_auto.  C++ is used because matlab is much, much slower at doing these calculations.

cpp_curvatures.mexw64
	This is a "mex" file, which, allows Matlab to execute c++ code.  If you are experiencing errors, you can
	have matlab rebuild this file by issuing the following command in Matlab: mex cpp_curvatures.cpp
	
sum_dislocations
	An optional file included for convenience.  This function can take the results from GND_auto and package
	them into a structure where similar dislocation types are grouped together for easy plotting.  For example
	with HCP materials, the dislocation densities will be arranged into prism<a>,screw<a>,basal<a>,pyramidal<c+a>,
	screw<c+a>, and total GND density.


*************************************	
Running The Code:
*************************************
First, make sure you have installed Matlab and MTEX, the latter of which can be found here: (http://mtex-toolbox.github.io/download.html).

Import your EBSD data into Matlab using MTEX.  If you are unfamiliar with this, the MTEX website contains tutorials that you may find
helpful.  Once you have imported your data, you must perform grain reconstruction on the dataset in order for the GND code to be used.
It is also highly reccomended that you smooth out the orientation data in your EBSD map, using the smooth command in MTEX.  Smoothing
is not strictly necessary to obtain results, but doing so drastically improves the quality of the results you get.  Smoothing can be 
done with this command:

ebsd_smoothed=smooth(ebsd,splineFilter);

Where ebsd is your dataset.  The splineFilter is the simplest filter to use, but there are others if you wish to experiment.

Once your data is smoothed, simply run the following command:

[dislocations,systems]=GND_auto(ebsd,nthreads,poisson,cubictype);

The GND data will be stored in "dislocations".

"systems" will be a structure containing info on the slip systemds for each phase (burgers vectors, planes, etc). 

ebsd is your ebsd dataset that you would like to analyze

nthreads is the number of cores you want to run the code in parallel (requires that matlab has parallel computing toolbox, otherwise the code will run in serial)

poisson is an array containing Poisson's ratio for each indexed phase in the dataset, in order

cubictype is a cell array denoting the types of any cubic phases, in the order they appear (for example if there are two cubic phases in your dataset, and the first one is FCC and the second BCC then cubictype={'FCC' 'BCC'}

NOTE: If you do not include the poisson or cubictype arguments when calling the function, you will be prompted by a pop-up window to enter that information.


****************************************
Notes for HCP systems:
****************************************

HCP metals are not as clear-cut as cubic ones in which slip systems they exhibit.  This code by default uses slip systems typical of Zr and its alloys for any HCP metal.  Specifically, the systems below:

prism<a>: {10-10}/<11-20>

basal<a>: {0001}/<11-20>

pyramidal<c+a>: {10-11}/<11-23>

I know that some HCP metals (eg. Mg) can accomplish pyramidal slip through a slightly different slip system.  If you wish to change which slip system the code uses, you need only open the doBurgersHCP.m file and make
a small edit.  If you look through the file you will find a comment designating where a slip system setup begins (example: %%Pyramidal Slip System (Edge)).  Immediately following this will be a pair of lines that assigns
values to a variable b (the slip direction) and n (the slip plane).  Just change the numbers inside the call to the "Miller" functions in these lines to the slip direction and plane you wish to use.  You don't need to worry
about symmetries, the next few lines of code do that automatically, just get the right families of planes/directions and everything should be fine.