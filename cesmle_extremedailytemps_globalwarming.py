#************************************************************************
# Written by Eleanor Middlemas
# July 27, 2020

import xarray as xr
import numpy as np
import time
from scipy import stats
import math
from glob import glob
from itertools import chain

# plotting tools
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
import cartopy
import cartopy.crs as ccrs
import cmocean

#************************************************************************
# Calculations

#------------------------------------------------------------------------
# Find baseline daily climatology distribution for summer & winter 
# seasons for each ensemble member
#------------------------------------------------------------------------

def extract_seasons(infilename):
    """
    This function takes in a filename corresponding to daily-resolution timeseries
    data and returns only the values over the U.S. during JJA and DJF in two separate
    datasets.
    
    INPUTS:
    infilename: string indicating full path of daily timeseries netcdf filename
    
    RETURNS:
    jja: dataset with only daily data during June, July, and August months
    djf: dataset with only daily data during December, January, and February months
    """
    ds = xr.open_dataset(infilename).sel(lat=slice(24,72), lon=slice(190,310))
    
    def is_jja(month):
        return (month >= 6) & (month <= 8)

    def is_djf(month):
        return (month == 1) | (month == 2) | (month==12)
    
    djf = ds.sel(time=is_djf(ds['time.month']))
    jja = ds.sel(time=is_jja(ds['time.month']))
    
    return jja, djf

def calculate_baseline(infilename, jja_outfilename, djf_outfilename):
    """
    Subselects DJF & JJA daily values from years 1920 to the end of 1949 
    for the purpose of serving as a baseline period for global warming.
    A 30-year period was selected to avoid effects of large ENSO events. 
    
    The resulting daily values are saved as a NetCDF file in the baseline 
    folder (see calculation_driver() function below for more information.)
    
    INPUTS:
    1. infilename: string indicating full path of daily timeseries netcdf filename
    2. jja_outfilename: string indicating full path of output daily data during JJA months from years 1920 - 1949.
    3. djf_outfilename: "  " " during DJF months from years 1920-1949.
    
    RETURNS:
    Nothing, but saves two NetCDF files per input file: one for JJA and one for DJF.
    """
    jja, djf = extract_seasons(infilename)

    djf_baseline = djf.sel(time=(djf['time.year']>= 1920) & (djf['time.year']< 1950))
    jja_baseline = jja.sel(time=(jja['time.year']>= 1920) & (jja['time.year']< 1950))
    
    # CESMLE output data have all these extra coordinates that we don't need.
    # Let's remove them.
    ds_ndcoords = {'P0', 'ch4vmr', 'co2vmr', 'date', 'date_written', 'datesec', 'f11vmr', 'f12vmr', 'gw', 'hyai', 'hyam', 'hybi', 'hybm', 'mdt', 'n2ovmr', 'nbdate', 'nbsec', 'ndbase', 'ndcur', 'nlon', 'nsbase', 'nscur','nsteph','ntrk','ntrm','ntrn','sol_tsi','time_bnds','time_written','w_stag','wnummax'}
    
    # Could merge these into the same dataset for faster
    # computation, but requires defining an additional coordinate
    # because the number(JJA days)!=number(DJF days)
    jja_baseline = jja_baseline.drop(ds_ndcoords)
    jja_baseline.to_netcdf(jja_outfilename,mode="w")
    
    djf_baseline = djf_baseline.drop(ds_ndcoords)
    djf_baseline.to_netcdf(djf_outfilename, mode="w")

#------------------------------------------------------------------------
# Find extremes & save to netcdf
#------------------------------------------------------------------------

def find_extremes(season, baseline):
    """
    This function is quite simple: it simply returns temperature "extremes"
    based on the 1920-1950 JJA/DJF baselines calculated above.
    
    INPUTS: 
    - season: string indicating season ('jja' or 'djf')
    - baseline: XArray dataset calculated from the calculate_baseline() 
    function above.
    
    RETURNS:
    - extremes: an XArray dataset that contains the days that correspond
    to the cold extremes and warm extremes of either a summer or winter
    seasons, i.e., those with temperature values that exceed 2 standard 
    deviations of the baseline distribution.
    """
    avg, std = np.mean(baseline.TS, axis=0), np.std(baseline.TS,axis=0)
    
    hot = np.greater(season.TS, avg + 2*std)
    hot.name = 'hot'
    cold = np.less(season.TS, avg - 2*std)
    cold.name = 'cold'
    extremes = xr.merge([hot,cold],
                        compat='broadcast_equals',join='inner')
    return extremes

def count_extremedays_peryear(jja_extremes, djf_extremes):
    """
    This function is also very simple: it counds the number of days
    that are designated an "extreme" from the find_extreme() function above.
    
    The resulting dataset includes the four extremes (hot/cold jja/djf), each with
    a number of extreme days per year and at every gridpoint.
    
    INPUTS: 
    - jja_extremes/djf_extremes: both datasets that were output as a result
    of the find_extremes() function above.
    
    RETURNS:
    - extreme_days_count: a dataset that contains counts for each extreme type
    for each gridpoint and for every year in the input dataset.
    """
    jjahot = jja_extremes.hot.groupby('time.year').sum()
    jjahot.name = 'jjahot'
    jjacold = jja_extremes.cold.groupby('time.year').sum()
    jjacold.name = 'jjacold'
 
    djfhot = djf_extremes.hot.groupby('time.year').sum()
    djfhot.name = 'djfhot'
    djfcold = djf_extremes.cold.groupby('time.year').sum()
    djfcold.name = 'djfcold'

    extreme_days_count = xr.merge([jjahot,jjacold,djfhot,djfcold],
                        compat='broadcast_equals',join='inner')
    return extreme_days_count

def count_extreme_freq_in_time(infilename,jja_baseline_path,
                               djf_baseline_path,outfilename):
    """
    This combines the last two functions into one function. Users must first 
    compute seasonal baselines (see calculate_baselines() above) before 
    running this function.
    
    INPUT: 
     - infilename: a string of the full path to any daily-timeseries netcdf 
    dataset
    - jja_baseline_path: a string to the full path of the pre-calculated 
    jja baselines
    - djf_baseline_path: " " but for djf
    - outfilename: a string indicating the preferred output filename that
    will hold the dataset of the number of extreme JJA & DJF days each year.
    """
    jja, djf = extract_seasons(infilename)
    
    jja_baselines = xr.open_dataset(jja_baseline_path)
    djf_baselines = xr.open_dataset(djf_baseline_path)
    
    jja_extremes = find_extremes(jja, jja_baselines)
    djf_extremes = find_extremes(djf, djf_baselines)
    
    extreme_count = count_extremedays_peryear(jja_extremes,djf_extremes)
    
    extreme_count.to_netcdf(outfilename,mode="w")

#------------------------------------------------------------------------
# Calculation code driver
#------------------------------------------------------------------------

def calculation_driver(baselines_path,extremects_path):
    """
    Hard-coded path to CESM-LE daily data on NCAR's supercomputer Cheyenne, because... well, that's the easiest way for most people I know to access it.
    
    Also, because I didn't have time, users must create the following directory hierarchy somewhere:
    1. folder holding two folders that are referred to as "baselines_path" & "extremects_path" below. 
    2. two folders in the "baselines_path" folder: one titled "jja" and one titled "djf"
    """
    # Set up string manipulation for handling large ensemble data
    rootdir = '/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE'
    T20C = 'b.e11.B20TRC5CNBDRD.f09_g16'
    TR85 = 'b.e11.BRCP85C5CNBDRD.f09_g16'    

    # Define which ensemble members should be analyzed
    member_ids = [i for i in chain(range(2,36),range(101,106))]
    component = 'atm'
    freq = 'daily'
    variable = 'TS'
    pth = f'{rootdir}/{component}/proc/tseries/{freq}/{variable}'
    
    # create dictionaries that define input filenames, baseline filenames,
    # and extreme count filenames
    baseline_fs = {}
    extremects_fs = {}
    filedict = {}
    for ens in member_ids:
        filedict[ens] = []
        for sc in [T20C, TR85]:
            glob_res = sorted(glob(f'{pth}/{sc}.{ens:03d}*.nc'))
            filedict[ens].extend(glob_res)
        baseline_fs[ens] = [f'{baselines_path}jja/ens.{filedict[ens][0].split(".")[4]}.jja.{variable}.19200101-19491231.nc',
                               f'{baselines_path}djf/ens{filedict[ens][0].split(".")[4]}.djf.{variable}.19200101-19491231.nc']
        extremects_fs[ens] = [f'{extremects_path}ens.{mem.split(".")[4]}.extreme_counts.{".".join(mem.split(".")[7::])}' for mem in filedict[ens]]
        
    # Step through each ensemble member
    start = time.time()
    for mem in member_ids:
        print("--- ensemble member "+str(mem)+" ---")
        
        # First, calculate the baseline DJF/JJA for each ensemble member.
        calculate_baseline(filedict[mem][0],baseline_fs[mem][0],baseline_fs[mem][1])
        print("Baselines computed.")
        
        # Then, calculate the number of JJA/DJF extremes over time. Some 
        # ensembles have 2 files spanning 1920-2100, and others have 3 files.
        if mem>=34:
            for i in range(0,2):
                count_extreme_freq_in_time(filedict[mem][i],
                                  baseline_fs[mem][0],baseline_fs[mem][1],
                                  extremects_fs[mem][i])
        else:
             for i in range(0,3):
                count_extreme_freq_in_time(filedict[mem][i],
                                  baseline_fs[mem][0],baseline_fs[mem][1],
                                  extremects_fs[mem][i])
        print("Extremes found from 1920-2100.")
    
    end = time.time()
    print(f"Time to produce data for all {len(member_ids)} ensemble members: ")
    print(str(end-start)+" seconds.")
    
#************************************************************************
#------------------------------------------------------------------------
# *** Plot everything! ****
# 
# Now that the calculations/heavy lifting is done (should take a couple
# of hours if no parallel processing is used), we can start to make plots
# to visualize what the changes in frequency of winter/summer extremes
# look like with climate change.
#------------------------------------------------------------------------

def return_daily_temp_dists(ensno, season):
    """ 
    This simply extracts distributions of temperature from the baseline
    file as well as from the file that contains the time period 2080-2100.
    
    INPUTS:
    - ensno: an integer representing the ensemble number
    - season: a string, either "jja" or "djf"
    
    RETURNS:
    - Two data arrays, one for the baseline daily temperatures and another
    for daily temperatures during 2080-2100 after RCP8.5 forcing
    """
    baseline_pth = "$somepath/baselines/"
    baselinedatafile = glob(f"{baseline_pth}jja/*.{ensno:03d}*.nc")[0]
    bsds = xr.open_dataset(baselinedatafile)
    basets = bsds.TS-273.15
    
    cesmle_20802100_pth = "/glade/collections/cdg/data/cesmLE/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/TS/b.e11.BRCP85C5CNBDRD.f09_g16"
    extremesdatafile = glob(f"{cesmle_20802100_pth}*.{ensno:03d}*2100*nc")[0]
    if season=="jja":
        ex,__ = extract_seasons(extremesdatafile)
    else:
        __,ex = extract_seasons(extremesdatafile)
    extremets = ex.TS-273.15

    return basets, extremets

def dist_plot(ensno, season, lat, lon):
    """
    Overlays the distribution of temperature from the baseline and the
    distribution of temperature from years 2080-2100 after RCP8.5 forcing
    during a provided season and at a given lat & lon point. Saves in .eps 
    format.
    
    INPUTS:
    - ensno: an integer indicating the model ensemble number
    - season: a string, either 'djf' or 'jja'
    - lat: a float indicating the latitudinal coordinate
    - lon: a float indicating the longitudinal coordinate
    """
    fname = 'daily_'+season+'dist_with_warming_lat{0:2.2f}'.format(mean_or_coord[0])+'_lon{0:2.2f}'.format(mean_or_coord[1])+".eps"
    basets, extremets = return_daily_temp_dists(ensno, season)
    basets_local = basets.sel(lat=lat,lon=lon,method='nearest')
    extremets_local = extremets.sel(lat=lat,lon=lon,method='nearest')
    
    # plot baseline daily temperatures (1920-1949)
    ax = sns.distplot(basets_local,color='lightblue',label="1920-1949");
    bottom, top = ax.get_ylim()
    
    # Indicate 2 standard deviations of baseline temps
    avg, stddev = np.mean(basets_local, axis=0), np.std(basets_local,axis=0)
    plt.plot([-2*stddev+avg, -2*stddev+avg], 
             [bottom,top],linestyle='--',color='black')
    plt.plot([2*stddev+avg, 2*stddev+avg], 
             [bottom,top],linestyle='--',color='black')
    
    # overlay 2080-2100 daily temperature distribution    
    sns.distplot(extremets_local,color="orangered",label="2080-2100")
    
    ax.legend()
    plt.title(f"Distribution of {season.upper()} daily surface temperatures \n in ens. member {str(ensno)} of CESM-LE",fontsize=16)
    plt.xlabel("TS (Celsius)",fontsize=14)
    plt.ylabel("Occurrence",fontsize=14)
    plt.ylim(bottom,top)
   
    plt.savefig(fname)

def read_all_ensembles_extremes():
    """
    This function reads all the extremes into a dataset, concatenating in
    time and along a new ensemble dimension.
    The resulting dimensions of the dataset should be
    (ensemble member, year, lat, lon).
    
    Returns a dataset, with four arrays according to each type of extreme
    (hot/cold djf/jja).
    """
    extremes_pth = "$somepath/extreme_counts/"
    member_ids = [i for i in chain(range(2,36),range(101,106))]
    
    dsets = []
    for mem in member_ids:
        glob_res = sorted(glob(f'{extremes_pth}/*.{mem:03d}*.nc'))
        dsets.append(xr.open_mfdataset(glob_res,concat_dim='year',combine='nested'))

    dsets_aligned = xr.align(*dsets,join='inner')
    first = dsets_aligned[0]
    rest = [ds.reset_coords(drop=True) for ds in dsets_aligned[1:]]
    objs_to_concat = [first] + rest

    # concatenate
    ensemble_dim = xr.DataArray(member_ids, dims='member_id', name='member_id')
    ds = xr.concat(objs_to_concat, dim=ensemble_dim, coords='minimal')

    # restore non_dim_coords to variables
    non_dim_coords_reset = set(ds.coords) - set(ds.dims)
    ds = ds.reset_coords(non_dim_coords_reset)

    return ds

def roundup(x):
    return int(math.ceil(x / 8.0)) * 8

def make_map_plot(nplot_rows, nplot_cols, plot_index, data, cmap):
    """ Create a single map subplot. Used in function below."""
    ax = plt.subplot(nplot_rows, nplot_cols, plot_index, projection = ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cartopy.feature.OCEAN, zorder=100, edgecolor='k')
    mxlevel = data.max().values
    levels = np.linspace(0,roundup(0.9*mxlevel),9)
    cplot = plt.contourf(data.lon,data.lat,
                         data, 
                         levels = levels,
                         cmap = cmap,
                         extend="both",
                         transform=ccrs.PlateCarree())
    cbar = plt.colorbar(cplot, ax=ax, orientation='horizontal', pad = 0.1) #shrink = 0.7
    cbar.set_ticks(levels)
    
    return cplot, ax

def plot_extremedays_map(season_extremes):
    """
    Displays a map of a count of extremes in four different time periods:
    2000-2010, 2030-2040, 2060-2070, 2090-2100. Averaged across ensemble members.
    Saves in .eps format.
    
    season_extremes is a string indicating one of the following types of extremes:
    'jjahot','jjacold','djfhot','djfcold'
    """
    ds = read_all_ensembles_extremes()
    
    fname = season_extremes+"_spatialmap_counts.eps"
    ends = [2010,2040,2070,2100]

    numPlotCols = 2
    numPlotRows = 2    
    fig, axs = plt.subplots(numPlotRows, numPlotCols, figsize=(10,6))

    if season_extremes=='jjahot' or season_extremes=='djfhot':
        cmap = cmocean.cm.amp
    else:
        cmap = cmocean.cm.deep
    
    for i,v in enumerate(ends):
        plotRow = i // numPlotCols
        plotCol = i % numPlotCols
        plot_index = i + 1
        plot_data = ds[season_extremes].sel(year=slice(v-10,v)).mean('member_id').sum('year')

        # Using functino above to make each subplot.
        cplot, axs[plotRow, plotCol] = make_map_plot(numPlotRows, numPlotCols, plot_index, plot_data, cmap)
        plt.title(f'{v-10}-{v}',fontsize=15);

    plt.suptitle('Average number of extreme days across ensembles', fontsize = 20)
    plt.tight_layout()
    plt.savefig(fname)

def plot_extremes_in_time(mean_or_coord):
    """
    Finally, the most exciting plot of all.
    Plotting the number of extremes per year in time from 1920-2100.
    Each ensemble is plotted in light blue, with the ensemble average plotted over.
    
    The user has the choice of plotting either the average over the whole U.S., 
    in which case, mean_or_coord = "areamean", or the user may put in a single
    coordinate array, i.e., [latvalue,lonvalue], to see how the extremes change
    at a single point.
    
    Figure outputs in .eps format.
    """
    ds = read_all_ensembles_extremes()
        
    if type(mean_or_coord)==str and mean_or_coord=="areamean":
        area_weights = np.cos(np.deg2rad(ds.lat))
        area_weights.name = "weights"
        area_weighted = ds.weighted(area_weights)
        plot_data = area_weighted.mean(("lon", "lat"))
        fname_suff = "area_mean"
    else:
        plot_data = ds.sel(lat=mean_or_coord[0],lon=mean_or_coord[1],method='nearest')
        fname_suff = 'lat{0:2.2f}'.format(mean_or_coord[0])+'_lon{0:2.2f}'.format(mean_or_coord[1])
    fname = "./extreme_count_timeseries_"+fname_suff+".eps"

    # set some colors
    c = {'cesmlea': '#1f78b4', 
         'cesmle': '#a6cee3'}

    # define ensemble member ids:
    member_ids = [i for i in chain(range(2,36),range(101,106))]
    
    # specify the xtick locations
    xtick = [1920] + list(np.arange(1950,2125,25))
    xticklbl = [f'{y:0d}' for y in xtick]

    # create figure
    fig = plt.figure(figsize=(5*1.2,7*1.2))

    # loop over hemisphere and plot sea-ice extent for each
    for i, v in enumerate(['jjahot','djfhot','jjacold','djfcold']):
        ax = fig.add_subplot(4,1,i+1)
        # plot each ensemble member
        for ens in member_ids:
            pe = ax.plot(ds.year, plot_data[v].sel(member_id=ens),
                         color=c['cesmle'], linewidth=0.5, label='CESM-LE')        
            
        maxvl = plot_data[v].max(['year','member_id']).values
        print(f"Max frequency for {v} is {maxvl}")     
        # plot the ensemble mean
        pea = ax.plot(ds.year, plot_data[v].mean('member_id'),
                      color=c['cesmlea'], linewidth=2., label='<CESM-LE>')
        # ticks, labels and legend
        ax.set_xticks(xtick)
        ax.set_xticklabels(xticklbl)
        if i == 0:
            ax.set_xticklabels([])
            plt.legend([pe[0], pea[0]] ,['CESM-LE', '<CESM-LE>'])
            plt.title("Number of extreme days across the U.S.")
        else:
            plt.title("")

        ax.set_ylabel(f'{v.upper()}')
        ylm = ax.get_ylim()
        fuzz = np.diff(ylm)*0.05
        ytick = [y for y in ax.get_yticks() if ylm[0] < y and y < ylm[1]]
    
    plt.savefig(fname)

def plot_driver():
    """
    A function which should be edited by the user to take advantage of
    previous functions to output three plots.
    """    
    dist_plot(13, 'jjahot', lat, lon)

    plot_extremes_in_time("emean")

    # Try changing the coordinates - don't forget to add 180  
    # to the longitude.
    [lat, lon] = [25.44174, 80.09543+180.] # Miami, FL
    plot_extremes_in_time([lat,lon])
    
    # plot maps for each extreme type
    for season in ['jjahot','djfhot','jjacold','djfcold']:
        plot_extremedays_map(season)
    