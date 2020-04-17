# terrace-long-profiler
Plots profiles of river long profiles and their terraces; may eventually fit ancient long profiles to the terraces. Uses outputs from the terrace detection code in [LSDTopoTools](https://github.com/LSDtopotools).

For more details of how the terrace code works then see [Clubb et al. (2017)](https://www.earth-surf-dynam.net/5/369/2017/esurf-5-369-2017.html):

Clubb, F.J., Mudd, S.M., Milodowski, D.T., Valters, D.A., Slater, L.J., Hurst, M.D., and Limaye, A.B. (2017) Geomorphometric delineation of floodplains and terraces from objectively defined topographic thresholds, Earth Surface Dynamics, doi:10.5194/esurf-2017-21

## Install the python packages
I recommend that you run this code using an Anaconda/Miniconda distribution. You need to create an environment which will have
all the necessary packages to run the code. I have provided an `environment.yml` file which you can use to create the environment.

If you are using a Linux system you can do this natively. I will assume you have Anaconda/Miniconda installed with Python 3.7.

First, clone this repository:
```
git clone https://github.com/umn-earth-surface/terrace-long-profiler.git
```
Then navigate to the terrace long profiler directory and create the environment:
```
cd terrace-long-profiler
conda env create -f environment.yml
```
This will create an environment called `terraces` with all the necessary packages installed. You can activate and deactivate this environment with:
```
conda activate terraces
conda deactivate terraces
```

### Run the code
