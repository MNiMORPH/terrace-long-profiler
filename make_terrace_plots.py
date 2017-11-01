# Driver to make the terrace profile plots
# FJC 01/11/17

# import modules
import sys
import os

import LongProfilePlotter

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
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error.")

    # What sort of analyses you want to do
    parser.add_argument("-LP", "--long_profiler", type=bool, default=True, help="If this is true, I'll make plots of the terrace long profiles (Default = true)")

    args = parser.parse_args()

    # get the base directory
    if args.base_directory:
        this_dir = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not this_dir.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            this_dir = this_dir+"/"
    else:
        this_dir = os.getcwd()

    # check if you supplied the DEM prefix
    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()

    if args.long_profiler:
        LongProfilePlotter.long_profiler(this_dir, args.fname_prefix)

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
