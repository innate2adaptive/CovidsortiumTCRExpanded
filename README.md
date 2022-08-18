# CovidsortiumTCRExpanded
Scripts that produce data and figures for the COVIDsortium TCR paper (Aug 2022)

### Figures – data and scripts

Current version of MS is Word file TCR_COV2_V5.docx; Current main figures : Figures_V4.pptx; Current supplementary figures : Figures_S_V2.pptx 

Data can be found at: https://www.dropbox.com/sh/neg3lpofw85l9ng/AAAw-HckR7LYxKz7TMNQ90zca?dl=0

Scripts in: scripts/

Output plots in: output_figs/

#### Main figs
- [x] **1A** script_1.R, this makes the dot plots for each patient, between parwise timepoints. The saveImage() function is currently commented out. 

Data input is in the dropbox, inside the data folder. 

Need to run twice, once for alpha, and once for beta. Produces plots; and a list of up and downregulated TCRs, which is saved as "TCR_change_all_alpha.RData" and "TCR_change_all_beta.RData"

*I am missing a script here that goes from TCR_change_all_[alpha/beta].RData to exp_AB_wide1.RData (?)*

On exp_AB_wide1.RData, you run Fig_1_remove_invariant.R to remove MAIT and IkT cells from the expanded list -> obtain exp_AB_wide3.RData

*I didn't have exp_AB_wide1.RData from Benny, so I used the version I had in my old repo*

*The new axes here are shifted by a power of 2. When I plot with Benny's plotting without renaming axes, I get the same axes as on my new plots. Was there a reason Benny had shifted the axes by a power of 2?*

- [x] **1B** bulk_analysis_expanded.R calculates richness, Shannon diversity, Renyi etc from exp_AB_wide3.R, then saves it to a result file. Then plot with Fig1b.R

- [x] **1C, D** script_2.R, Section B

- [ ] **1E, F**

- [ ] **1D** script_3.R

- [ ] **2A, B** Numbers in matches_summary.xlsx

- [ ] **2C** Stacked Histogram plotted in Numbers in matches_summary.xlsx

- [ ] **2D** script_4.R Section C

- [ ] **2E** script_4.R Section A

- [ ] **2F** Scripts from Martina to add.

- [ ] **Fig 3**

- [ ] **Fig 4** Scripts from Cedric to add.

#### Additional scripts

Fig_1_remove_invariant.R removes MAIT and IkT cells from the expanded list.

#### Supplementary

- [ ] **S1** FACS Plots from Michal

- [x] **S2** script_2.R, Section A - 2 separate heatmaps for PCR+ and PCR-

- [ ] **S3** script_S1.R

- [x] **S4** bulk_analysis_all.R, this takes all_A_long.Rdata and all_B_long.Rdata as input, and calculates some of the standard bulk metrics. The script then saves the results. This is equivalent to bulk_analysis_expanded.R. Then plot with FigS4.R.

- [x] **S5** script_1.R

- [ ] **S6** script_2.R Section D

- [ ] **S7** ???

- [ ] **S8** script_3.R

- [ ] **S9** Plots from Yanchun Peng

- [ ] **S10** ???

- [ ] **S11** script_4.R Section B

- [ ] **S12** Script from Martina to add

- [ ] **S13** Script from Martina to add
