## Name
Using Correlated GBM, BAT, and XRT GRBs to Create a Progenitor Classifier

## Description
Code and data that will reproduce the analysis done in the paper of the same name, plus a (regularly-updated) python code that will return an initial classification of new GRBs based on their peak energy, fluence, and duration. The general gamma-ray burst data was taken from the [Fermi-GBM catalog](https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dfermigbrst&Action=More+Options) through 22 May 2023, [the Swift-Bat catalog](https://swift.gsfc.nasa.gov/results/batgrbcat/summary_cflux/summary_general_info/) and [Swift preliminary catalog](https://swift.gsfc.nasa.gov/archive/grb_table/) through 14 April 2023,  a [table of Swift-Bat's best fit spectral parameters over t100](https://swift.gsfc.nasa.gov/results/batgrbcat/index_tables.html) through 14 April 2023, both [J. Griener's](https://www.mpe.mpg.de/~jcg/grbgen.html) and the [official Swift-BAT](https://swift.gsfc.nasa.gov/results/batgrbcat/summary_cflux/summary_general_info/GRBlist_redshift_BAT.txt) redshift tables through 14 April 2023, [the Swift-XRT catalog](https://www.swift.ac.uk/xrt_live_cat/) (with all columns selected) through 22 May 2023, and the burst analyzers from the [swifttools.grb submodule](https://www.swift.ac.uk/API/ukssdc/data/GRB.md) provided by the University of Leicester. The creation of the known progenitors is detailed in the associated paper. All the other data files were derived from these data files through the included python notebooks. 

To reproduce the analysis provided in the paper, it is sufficient to run the Ep_S_t90_Data_Creation.py, energy-energy_data_creation.py, and Flu_X-ray_Flux_Data_Creation.py, followed by the side_analysis.rmd. This will output some of the statistics used in the main paper, while also creating some of the files needed for the Main_Analysis.ipynb, which outputs the other statistics and the figures. If one only wants to recreate part of the analysis, it should also be sufficient to run just that part of the data creation, followed by just that "chunk" of the R and main python codes (the analysis of the classifier spans many chunks). 

swift_data_downloading.py can also be combined with xrt_lcs/find_afterglow_start_and_end_times.py to find when each burst's afterglow starts and ends and how much fluence it contained. This step is not needed to run the main analysis, but rather a side analysis performed in the Main_Analysis.ipynb (the data for which is included). This step may take up to an hour.

Finally, the classifer python code can be run from the command line by simply naviagating to the installed folder, and running "python3 single_burst_predictor.py -Ep [] -S [] -t90 []", where peak energy is in keV, fluence is in ergs cm^-2, and duration is in seconds. It will print the initial classification of your GRB and create a file called "your_GRB_in_context.pdf", which displays a figure containing your GRB, a few famous GRBs for context, our known progenitors, and the scale on which the classifier sorts GRBs. This code will be refined regualrly as more known progenitors occur, so the version number may change.

## Installation
R can be downloaded directly from [Comprehensive R Archive Network (CRAN)](https://cran.rstudio.com/). Then [RStudio can be installed](https://posit.co/download/rstudio-desktop/) as appropriate for the operating system, allowing one to edit the included R markdown files. Both the R and python (pip) packages can then be installed from the commands in installed_packages.txt. We recommend running the python notebooks in an enviornment that allows for one to run then in pieces, such as jupyter or Visual Studio. This allows the user to run on the relevant pieces if they so choose.

## Usage
This code exists to serve as both a reproduction and extension of the analysis done in the attached paper. Through the progenitor-prediction python notebook, we hope to refine our machine learning model while simultaneously learning more about why certain long-merger and short-collapsar GRBs exist. The rest simply expands on the statstics and graphics doen in the paper itself.

## Support
This is the main help document, though some of the assistnace is spread through the files themselves.
Most of the active assistance will be provided by Pi Nuessle (n.nuessle@nasa.gov or nuessle167@gwu.edu) but the #software-help channel in the BurstCube Slack may also have some hints.

## Contributing
We are open to contributions.

## Authors and acknowledgment
We thank Israel Martinez, Joe Ascersion, and David Friedlander for various contributions to this project. 

## License
This code is written using open source packages and software and available under an MIT License. This merely states that the liscense be preserved with the code, and that commercial use can only be pursued with express permission. We dedicate the code itself to open use, and we accept modification through this Git repository.

## Project status
This will likely be the final update before the paper's release, though we will update the predictor every 3-6 months with new progenitors.
