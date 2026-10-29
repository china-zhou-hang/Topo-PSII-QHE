SOURCE DATA

This package contains the numerical results underlying Figures 2-4 and Supplementary Figures 1-6. Source_Data.xlsx provides a spreadsheet representation for inspection and reuse. The nine accompanying .mat files preserve the MATLAB arrays used to prepare the figures.

FILE CONTENTS

  fig2_data.mat       Figure 2
  fig3_data.mat       Figure 3
  fig4_data.mat       Figure 4
  supp_fig1_data.mat  Supplementary Figure 1
  supp_fig2_data.mat  Supplementary Figure 2
  supp_fig3_data.mat  Supplementary Figure 3
  supp_fig4_data.mat  Supplementary Figure 4
  supp_fig5_data.mat  Supplementary Figure 5
  supp_fig6_data.mat  Supplementary Figure 6

The .mat files retain the native numerical arrays and are the preferred inputs to the accompanying MATLAB plotting scripts. The spreadsheet is intended for access without MATLAB. Variable names correspond to those used in the scripts. Some plot axes are constructed in the scripts rather than stored in the .mat files; for example, the Figure 2 scan coordinates are reconstructed by fig2.m.

USING THE DATA WITH THE CODE

The MATLAB code is provided as a separate package. To reproduce the plots that load archived results, copy fig2_data.mat into the code's fig2data/ folder, fig4_data.mat into fig4data/, and supp_fig4_data.mat into supp_fig4data/. Keep the filenames unchanged. The code additionally requires the plotting colormap fig2data/cmap_RdYlBu.mat; this auxiliary file is not included among the nine numerical datasets.

For the remaining figures, consult the corresponding MATLAB scripts for the calculation and plotting procedures. The source data document the reported results but do not replace the required scripts or dependencies. See the article and Supplementary Information for the definitions, parameters, and units of the plotted quantities.
