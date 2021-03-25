# code_for_EFT_analysis

Creat the following directories :

/data/

/data/madgraph/

/data/madgraph/output/

/data/heppy/

/data/heppy/output/

/results/

/results/heppy/

/results/heppy/ratio_madgraph/

/results/heppy/top_reco/

/results/heppy/top_reco/ATCGRoot/

Put the rootfiles in the data directory : madgraph files in the madgraph directory and heppy files in the heppy directory.

Add the files names in runSingleTopLHEAnalyzer.cpp if you want to run using madgraph files.
Add the files names in runHeppyRootFile.cpp if you want to run using heppy files.

Use the command "Make" from root directory to compile.
Run ./bin/runSingleTopLHEAnalyzer or ./bin/runHeppyRootFile depending of the files you edited.

Then edit the edit_Plot.cpp or Plot.cpp files to get the plot or files to run combine/ATGC.

In heppy_Plot.cpp :

  - normalizeWithChi() compute Chi square
  - Stacked_histo_Fit() takes all MC files as input and compute a fit letting QCD free.
  - Stacked_histo() plot the DATA/MC agreement taking the fraction 0.X returned by TFractionFitter that is promped when running Stacked_histo_Fit().
  
  - Stacked_histo_reversed() and Stacked histo_reversed_fit() are similar to Stacked_histo() and Stacked_histo_Fit() but take the estimated QCD background as input.
  - ATGCRoot() makes the files ATGC files with data and background.
  
  
In Plot.cpp :

  - Ratio_EFT_SM() makes the TF1 file used by ATGC as signal_proc_"..".root
  
  
 After modifying a file, run "Make" before using ./bin/... commands.






# Code_EFT_Analysis
