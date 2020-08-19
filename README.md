# Large-Ensemble_seasonal_extremes

### Problem statement
This is code that I wrote in response to an interview assignment. The goal was to create an end-to-end climate data analysis that meets the following requirements: 
1. Is no longer than 500 lines of code(+comments)
2. Illustrates objective-oriented programming (complete with functions or classes)
3. Possibly exhibits knowledge of bash scripting, climate model ensembles, downscaling methods, and regional climate analysis.

### Objective
Because I am most familiar with model ensembles, I chose to extract and plot the change in daily extremes during summer and winter seasons from 1920 through 2100 in the CESM1 Large Ensemble. 

### Running the code yourself
That said, this code is most easily utilized from an account on the supercomputer hosted by the National Center for Atmospheric Research, because that's where the CESM1 Large Ensemble data is stored. If you have an account, you can simply utilize the NCAR_pylib module to access the relevant packages, and use the bash script to submit the python script to a Casper node.

If you are working from your local computer, you need to (1) [download at least some the large ensemble data](http://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html) to a local directory (and then change the `rootdir` on line 191 and the `cesmle_20802100_pth` on line 271), and (2) install the environment.yml file in order to use the relevant packages.

### Further action is needed to: 
1. Employ Dask in order to speed up the process of analyzing ensemble members
2. Debug the map plotting function - it currently does not work as intended.