# Driver to make the terrace profile plots
# FJC 01/11/17

# import modules
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import sys
import os

import TerracePlotter

#=============================================================================
# This is the main function that runs the whole thing
#=============================================================================
def main(argv):

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()

    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the terrace analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # Some filtering info for terrace pixels
    parser.add_argument("-min_size", "--min_size", type=int, help="The minimum size (in pixels) of a terrace patch. Default = 5", default=5)
    parser.add_argument("-min_elev", "--min_elev", type=int, help="The minimum elevation above the channel of a terrace pixel. Default = 0", default=0)
    parser.add_argument("-max_elev", "--max_elev", type=int, help="The maximum elevation above the channel of a terrace pixel. Default = large", default=10000000)

    # What sort of analyses you want to do
    parser.add_argument("-LP", "--long_profiler", type=bool, default=False, help="If this is true, I'll make plots of the terrace long profiles (Default = true)")
    parser.add_argument("-compiled", "--compiled", type=bool, default=False, help="If this is true, I'll combine all reaches to make a super-plot for the whole river")
    parser.add_argument("-PR", "--plot_rasters", type=bool, default=False, help="If this is true, I'll make raster plots of the terrace locations (Default=false)")
    parser.add_argument("-HM", "--heat_map", type=bool, default=False, help="if true I'll make a heat map of terrace locations along the river long profile")
    parser.add_argument("-dips", "--dips", type=bool,default=False, help="If this is true, I'll calculate the dip and dip direction of each terrace.")
    parser.add_argument("-3d", "--plot_3d", type=bool,default=False, help="If this is true, I'll make a 3d plot of each terrace surface.")

    # These control the format of your figures
    parser.add_argument("-fmt", "--FigFormat", type=str, default='png', help="Set the figure format for the plots. Default is png")
    parser.add_argument("-size", "--size_format", type=str, default='ESURF', help="Set the size format for the figure. Can be 'big' (16 inches wide), 'geomorphology' (6.25 inches wide), or 'ESURF' (4.92 inches wide) (defualt esurf).")

    args = parser.parse_args()

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        # check if you remembered a / at the end of your path_name
        # if not this_dir.endswith("/"):
        #     print("You forgot the '/' at the end of the directory, appending...")
        #     this_dir = this_dir+"/"
    else:
        this_dir = os.getcwd()

    # get the path separator for this os
    path_sep = os.path.sep

    # check if you supplied the DEM prefix
    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()
    # print the arguments that you used to an output file for reproducibility
    with open(this_dir+args.fname_prefix+'_report.csv', 'w') as output:
        for arg in vars(args):
            output.write(str(arg)+','+str(getattr(args, arg))+'\n')
        output.close()

    if not args.compiled:
        # read in the baseline channel csv
        lp = pd.read_csv(this_dir+args.fname_prefix+'_baseline_channel_info.csv')
        lp = lp[lp['Elevation'] != -9999]

        # read in the terrace csv
        terraces = pd.DataFrame()
        dist_file = this_dir+args.fname_prefix+'_terrace_info_filtered_dist.csv'
        # check if you have already calculated the distance along the baseline for each point
        if os.path.isfile(dist_file):
            terraces = pd.read_csv(this_dir+args.fname_prefix+'_terrace_info_filtered_dist.csv')
        else:
            terraces = pd.read_csv(this_dir+args.fname_prefix+'_terrace_info_filtered.csv')
            # find the nearest point along the baseline for each terrace ID
            terraces = TerracePlotter.get_distance_along_baseline_points(terraces, lp)
            terraces.to_csv(this_dir+args.fname_prefix+'_terrace_info_filtered_dist.csv', index=False)


        if args.long_profiler:
            TerracePlotter.long_profiler_all_terraces(this_dir, args.fname_prefix, terraces, lp)
            #TerracePlotter.long_profiler(this_dir, args.fname_prefix, terraces, lp)

        if args.plot_3d:
            TerracePlotter.PlotTerraceSurfaces(this_dir, args.fname_prefix, terraces)

        # if args.plot_rasters:
        #     TerracePlotter.MakeRasterPlotTerraceIDs(this_dir, args.fname_prefix, args.FigFormat, args.size_format)
        #     TerracePlotter.MakeRasterPlotTerraceElev(this_dir, args.fname_prefix, args.FigFormat, args.size_format)
        if args.dips:
            TerracePlotter.write_dip_and_dipdir_to_csv(this_dir,args.fname_prefix, args.digitised_terraces, args.shapefile_name)
            # TerracePlotter.MakeRasterPlotTerraceDips(this_dir,args.fname_prefix,FigFormat=args.FigFormat,size_format=args.size_format)
        if args.heat_map:
            TerracePlotter.MakeTerraceHeatMap(this_dir, args.fname_prefix, prec=100, bw_method=0.03, FigFormat=args.FigFormat, ages="")

    else: # Compile all the reaches to make a terrace plot for the whole river.
        lp = this_dir+'UMV_combined'+path_sep+args.fname_prefix+'_points.csv'
        dist_file = this_dir+'UMV_combined'+path_sep+args.fname_prefix+'_terrace_info_filtered_dist.csv'

        # read in the long profile csv
        lp_file = this_dir+'UMV_combined'+path_sep+args.fname_prefix+'_baseline_channel_info.csv'
        lp_df = pd.DataFrame()
        if os.path.isfile(lp_file):
            lp_df = pd.read_csv(lp_file)
        else:
            lp_df = TerracePlotter.merge_baselines(lp, lp_file)
        lp_df = lp_df[lp_df['Elevation'] != -9999]

        # check if you have already calculated the distance along the baseline for each point
        if os.path.isfile(dist_file):
            terraces = pd.read_csv(dist_file)
            print(terraces['DistAlongBaseline_new'])
        else:
            # find each sub-directory and get the distance from the shapefile
            subdirs = next(os.walk(this_dir))[1]
            master_df = pd.DataFrame()
            for dir in subdirs:
                if 'UMV_DEM5m_' in dir:
                    terraces = pd.read_csv(this_dir+dir+path_sep+dir+'_final_terrace_info_filtered.csv')
                    # find the nearest point along the baseline for each terrace ID
                    terraces = TerracePlotter.get_distance_along_baseline_points(terraces, lp_df)
                    terraces['reach'] = dir
                    master_df = master_df.append(terraces)
            master_df.to_csv(dist_file, index=False)
        #
        # make the long profile plot
        TerracePlotter.long_profiler_all_reaches(this_dir+'UMV_combined'+path_sep, args.fname_prefix, terraces, lp_df)



#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! Welcome to the terrace long profiler tool.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("For help type:")
    print("   python terrace_profile_plots.py -h\n")
    print("=======================================================================\n\n ")

#=============================================================================
if __name__ == "__main__":
    main(sys.argv[1:])
