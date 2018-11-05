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
	
cpp_curvatures.cpp (deprecated)*
	A c++ source file containing the code that calculates the lattice curvatures for an EBSD dataset.
	This is called by GND_auto.  C++ is used because matlab is much, much slower at doing these calculations.

cpp_curvatures.mexw64 (deprecated)*
	This is a "mex" file, which, allows Matlab to execute c++ code.  If you are experiencing errors, you can
	have matlab rebuild this file by issuing the following command in Matlab: mex cpp_curvatures.cpp
	
sum_dislocations
	An optional file included for convenience.  This function can take the results from GND_auto and package
	them into a structure where similar dislocation types are grouped together for easy plotting.  For example
	with HCP materials, the dislocation densities will be arranged into prism<a>,screw<a>,basal<a>,pyramidal<c+a>,
	screw<c+a>, and total GND density.

*These files are no longer included in new releases, as curvature calculations are now done entirely within Matlab/MTEX

*************************************	
Running The Code:
*************************************
Some thorough tutorials with worked examples of the code are available here: https://cochranec.github.io/



First, make sure you have installed Matlab and MTEX, the latter of which can be found here: (http://mtex-toolbox.github.io/download.html).

Import your EBSD data into Matlab using MTEX.  If you are unfamiliar with this, the MTEX website contains tutorials that you may find
helpful.  Once you have imported your data, you must perform grain reconstruction on the dataset in order for the GND code to be used, which is 
done with the following command:

[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',5*degree);

(You may substitute any number in for the 5, and this will be the threshold in degrees for defining a grain boundary)

It is also highly reccomended that you smooth out the orientation data in your EBSD map, using the smooth command in MTEX.  Smoothing
is not strictly necessary to obtain results, but doing so drastically improves the quality of the results you get.  Smoothing can be 
done with this command:

ebsd_smoothed=smooth(ebsd,splineFilter);

Where ebsd is your dataset.  The splineFilter is the simplest filter to use, but there are others if you wish to experiment. (I personally reccomend the infimalConvolutionFilter)

I have found that it is best to re-run the grain reconstruction command after smoothing, as the smoothing process can somethings cause small issues with grain boundaries.

Once your data is smoothed, simply run the following command:

[dislocations,systems]=GND_auto(ebsd,nthreads,poisson,cubictype);

The GND data will be stored in "dislocations".

"systems" will be a structure containing info on the slip systemds for each phase (burgers vectors, planes, etc). 

To get the data into a more managable form, run the following command:

GND=sum_dislocations(dislocations,systems,ebsd);

This will create a structure called GND which sorts the dislocation densities according to slip system.

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

pyramidal2<c+a>: {10-22}/<11-23>

If you wish to change which slip system(s) the code uses, you need only open the doBurgersHCP.m file and make
a small edit.  If you look through the file you will find a comment designating where a slip system setup begins (example: %%Pyramidal Slip System (Edge)).  Immediately following this will be a pair of lines that assigns
values to a variable b (the slip direction) and n (the slip plane).  Just change the numbers inside the call to the "Miller" functions in these lines to the slip direction and plane you wish to use.  You don't need to worry
about symmetries, the next few lines of code do that automatically, just get the right families of planes/directions and everything should be fine.  Using the doBurgersHPC.m file as a template, you should be able to create similar files for any set of crystal symmetries and slip systems you care to use.  I may add more crystal symmetries in future releases, but FCC, HCP, and BCC are all currently supported.

****************************************
Notes for Older Versions of MATLAB:
****************************************
For versions of MATLAB earlier than 2013:
- line 94 of GND_auto.m should be changed from:
options=optimoptions('linprog','Algorithm','dual-simplex','Display','off');
to
options=optimset('Display','off');

- add a line to GND_auto.m before line 101 ('%loop through all points') with the commend
matlabpool(nthreads);

****************************************
Further Reference
****************************************

This code is based on the approach detailed by W. Pantleon in the paper "Resolving the geometrically necessary dislocation content by conventional
electron backscattering diffraction", Scripta Materialia, 2008, 58, p994-997, DOI: 10.1016/j.scriptamat.2008.01.050
