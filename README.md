# Large-Ensemble_seasonal_extremes

### Problem statement
This is code that I wrote in response to a job interview assignment. The assignment was to create an end-to-end climate data analysis that meets the following requirements: 
1. Is no longer than 500 lines of code(+comments)
2. Illustrates objective-oriented programming (complete with functions or classes)
3. (optional) Possibly exhibits knowledge of bash scripting, climate model ensembles, downscaling methods, and regional climate analysis.

### Objective
Among the three topics listed in no. 3 above, I am most familiar with model ensembles. Thus, I chose to extract and plot the change in the number of daily extremes during summer and winter seasons from 1920 through 2100 in the CESM1 Large Ensemble. I visualize the change in frequency of four different types of extremes: cold summers, hot summers, cold winters, and hot winters.

### Running the code yourself
This code is most easily utilized from an account on the supercomputer hosted by the National Center for Atmospheric Research, because that's where the CESM1 Large Ensemble data is stored. If you have an account, you can simply load the `NCAR_pylib` environment to access the relevant packages, and use the provided bash script to submit the python script to a Casper node.

If you are working from your local computer, you need to (1) [download at least some the large ensemble data](http://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html) to a local directory (and then change the `rootdir` on line 191 and the `cesmle_20802100_pth` on line 271), and (2) install the environment.yml file in order to use the relevant packages.

### Further updates to the code is required to: 
1. Employ Dask in order to speed up the process of analyzing ensemble members
2. Debug the map plotting function - it currently does not work as intended.

### Analysis pipeline functions
#### Calculation functions
`extract_seasons` - returns only the days that fall within JJA or DJF.<br>
`calculate_baseline` - extracts the JJA & DJF days from 1920 through 1949 in order to use them as a baseline for defining an 'extreme'.<br>
`find_extremes` - calls `calculate_baseline` and returns days that exceed 2 standard deviations of the baseline period (<2&sigma; & >2&sigma; for a total of 4 different types of extremes: hot summers, cold summers, hot winters, cold winters).<br>
`count_extremedays_peryear` - count the number of times each of the four extremes occur each year.<br>
`count_extreme_freq_in_time` - a function that combines all of the above functions.<br>
`calculation_driver` - this code performs filename I/O in order to execute the `count_extreme_freq_in_time` over every ensemble member. The calculate outputs three small files for each ensemble member: (1) a baseline file for the summer days, (2) a baseline file for the winter days, (3) a file with a timeseries of the number of extremes per year for each of the four types of extremes.<br><br>
#### Plotting functions - run after the completion of calculation functions over all ensemble members
`return_daily_temp_dists` - given the ensemble number and season, opens the baseline & extreme files and returns the distribution of seasons during the end of the 22nd century.<br>
`dist_plot` - given the ensemble number, season, and latitude & longitudinal coordinates, plots the baseline distribution overlaid the end-of-century distribution calculated with `return_daily_temp_dists`. <br>
`read_all_ensembles_extremes` - extracts the frequency of extremes for all ensemble members from the files produced in the calculation pipeline.<br>
`roundup` - a funtion to round numbers to the nearest eight.<br>
`make_map_plot` - a function to specify characteristics of the map plot created in `plot_extremedays_map`<br>
`plot_extremedays_map` - This function needs some work. The goal is to create a plot of 4 subplots that show the regions that result in the largest changes in extreme over the course of a century.<br>
`plot_extremes_in_time` - creates four timeseries plots that show the number of daily extremes for each of the four types of extremes from 1920-2100.<br>
`plot_driver` - a function that calls all the plotting functions.<br>
