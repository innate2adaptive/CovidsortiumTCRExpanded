# CovidsortiumTCRExpanded
Scripts that produce data and figures for the COVIDsortium TCR paper (Aug 2022)

### Figures – data and scripts

Current version of MS is Word file TCR_COV2_V5.docx; Current main figures : Figures_V4.pptx; Current supplementary figures : Figures_S_V2.pptx 

Data can be found at: https://www.dropbox.com/sh/neg3lpofw85l9ng/AAAw-HckR7LYxKz7TMNQ90zca?dl=0

Note: the scripts automatically load the data from DropBox. This works for me (although sometimes I need to change the timeout), but does not seem to always work outside of the UK. You might need to download the data locally and tweak the scripts to load them.

Scripts in: scripts/

Output plots in: output_figs/

#### Main figs
- [x] **1A** script_1.R, this makes the dot plots for each patient, between parwise timepoints. 

Need to run twice, once for alpha, and once for beta. Produces plots; and a list of up and downregulated TCRs, which is saved as "TCR_change_all_alpha.RData" and "TCR_change_all_beta.RData"

**to_wide_data.R** converts the TCR_change_all_a/b into the format of exp_AB_wide.RData (i.e. has the abundances for all time points for the expanding TCRs) - *adapted, but I am missing a PCR_positive.csv file so cannot run*

**clean_wide_data.R** takes exp_AB_wide.RData and removes a few TCRs which are very high at -3 or -4 and then go down, because I think they are irrelevant. This provides  exp_AB_wide1.RData. The pGEN values were added manually to this. The version with manually added pGen is found in data/ - *adapted, but cannot run because I need to run to_wide_data.R first*

On exp_AB_wide1.RData (from data/output_data/ and not from data/, so using the version without manually annotated pGen), you run **Fig_1_remove_invariant.R** to remove MAIT and IkT cells from the expanded list -> obtain exp_AB_wide3.RData - *adapted, but cannot run because I need to run to_wide_data.R and clean_wide_data.R first*

*The new axes here are shifted by a power of 2. When I plot with Benny's plotting without renaming axes, I get the same axes as on my new plots. Was there a reason Benny had shifted the axes by a power of 2?*

- [x] **1B** bulk_analysis_expanded.R calculates richness, Shannon diversity, Renyi etc from exp_AB_wide3.R, then saves it to a result file. Then plot with Fig1b.R

- [x] **1C, D** script_2.R, Section B

- [x] **1E, F** script_3.R - plots dynamics for each patient separately.

- [x] **2A, B, C** Numbers and stacked histogram in matches_summary.xlsx 

*I don't fully understand matches_summary.xlsx, so I left them as they were, just changed the font and the styling*

- [x] **2D** script_4.R Section C. 

script_4.R uses exp_AB_merge2.txt, which is the list of annotated TCRs. There is no “script” for generating this. It is constructed semi-manually from the annotation data from Tao as well as what was in VDJdb and the Franics et al. paper.

*script_4.R also uses summary_with_HLA.txt, which I did not have, so got from my other repo*

- [x] **2E** script_4.R, Section A. Save heatmap and legend separately.

- [x] **2F** script_5.R - calculates and plots clustering

- [ ] **Fig 3**

- [ ] **Fig 4** Scripts from Cedric to add.

- [ ] **Fig 5** Scripts? Data?

#### Additional scripts

Fig_1_remove_invariant.R removes MAIT and IkT cells from the expanded list.

#### Supplementary

- [ ] **S1** FACS Plots from Michal - *they are a little pixelated, is there another version?*

- [x] **S2** script_2.R, Section A - 2 separate heatmaps for PCR+ and PCR-

- [x] **S3** script_S3.R *this script used to be called script_S1.R, I changed it for consistency*

- [x] **S4** bulk_analysis_all.R, this takes all_A_long.Rdata and all_B_long.Rdata as input, and calculates some of the standard bulk metrics. The script then saves the results. This is equivalent to bulk_analysis_expanded.R. Then plot with FigS4.R.

- [x] **S5** script_1.R

- [x] **S6** script_2.R Section D. Script 2 also produces a summary file - summary_data.csv with patient information and some key metrics

- [ ] **S7** ???

- [x] **S8A, B** script_3.R

- [ ] **S9** Plots from Yanchun Peng

- [ ] **S10** ???

- [x] **S11** script_4.R Section B

- [ ] **S12** script_6.R - calculates clustering at various threshold and compares expanded vs non-expanded controls

*(The script is ok, but crashes my laptop. Need to run it on big computer or cluster)*

- [x] **S13** script_5.R - calculates and plots clustering

**This script still uses a file where there are a CD4 and CD8 epitope which are the same, which is why some pies have 2 colours. Can we fix this? We were waiting for someone to let us know? - check if Benny's newest version is different (I don't think it will be because I am using the _clean.txt file, and Benny's AB_merge_2.txt is the same as I had before*

- [x] **S14** Generated by TP, Pymol?

### Notes:

**script_4.R**: missing summary_with_HLA.txt, so does not run fully

- Script to select controls?