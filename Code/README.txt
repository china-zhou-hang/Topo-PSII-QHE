MATLAB CODE

This package contains MATLAB scripts for Figures 2-4 and Supplementary Figures 1-6. Numerical results are archived separately in the Source Data package. In particular, some plotting scripts read precomputed .mat files and cannot be run directly without these files.

SETUP

Extract this package and retain its subfolder structure. Set the MATLAB current folder to the directory containing the figure scripts. Download the accompanying Source Data package and place the following files at the indicated relative paths:

  fig2_data.mat       -> fig2data/fig2_data.mat
  fig4_data.mat       -> fig4data/fig4_data.mat
  supp_fig4_data.mat  -> supp_fig4data/supp_fig4_data.mat

The plotting scripts also use fig2data/cmap_RdYlBu.mat, which contains the colormap. This auxiliary file is not part of the numerical Source Data and must be supplied with the code for the existing scripts to run unchanged.

EXECUTION

After the files are in place, run fig2.m, fig4.m, or supp_fig4.m to plot their archived results. The other figure scripts perform their calculations directly and may save additional .mat files when executed. Run scripts from the MATLAB current folder specified above so that relative file paths resolve correctly.

Figure 2 data generation (data_fig2.m) uses a thread-based parallel pool of 50 MATLAB workers (parpool('threads',50)) and requires the Parallel Computing Toolbox; the worker count may be reduced to match available computing resources.

Supplementary Figure 4 uses a consolidated dataset. Its plotting script must load supp_fig4data/supp_fig4_data.mat. An older version that loads individual tlist.mat and p*.mat files is not compatible with the consolidated dataset without updating its load statements.

The figure scripts and data files serve different purposes: the Source Data package permits reuse of the reported numerical arrays, while the MATLAB scripts specify the calculations and plotting procedures. Physical quantities, parameter definitions, and units are described in the article and Supplementary Information.
