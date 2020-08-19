#!/bin/bash -l
#SBATCH --job-name=cesmle_extreme_daily_analysis
#SBATCH --account=ACCOUNTNO
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --time=04:00:00
#SBATCH --partition=dav
#SBATCH --output=us_extremes_cesmle_analysis.out.%j

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

### Run program
module load python/3.7.5
module load ncarenv
ncar_pylib
python -c "from cesmle_extremedailytemps_globalwarming import calculation_driver; calculation_driver('/glade/work/eleanorm/testdata/baselines/','/glade/work/eleanorm/testdata/extreme_counts/')" >> output_extremes_timing.txt

python -c "from cesmle_extremedailytemps_globalwarming import plot_driver; plot_driver()" >> output_plot_errors.txt
