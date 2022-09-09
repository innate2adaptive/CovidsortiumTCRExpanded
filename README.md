# CovidsortiumTCRExpanded
Scripts that produce data and figures for the COVIDsortium TCR paper (Aug 2022)

### Figures – data and scripts

Current version of MS is Word file TCR_COV2_V5.docx; Current main figures : Figures_V5.pptx; Current supplementary figures : Figures_S_V4.pptx 

*Figures_V4 and S_V2 are before restyling, V5 and S_V3 has restyling (mostly) done, S_V4 is without figure S10.  Note that figures still save with same numbering as in S_V3*

Data can be found at: https://www.dropbox.com/sh/neg3lpofw85l9ng/AAAw-HckR7LYxKz7TMNQ90zca?dl=0

Note: the scripts automatically load the data from DropBox. This works for me (although sometimes I need to change the timeout), but does not seem to always work outside of the UK. You might need to download the data locally and tweak the scripts to load them.

Scripts in: scripts/

Output plots in: output_figs/

#### Main figs
- [x] **1A** script_1.R, this makes the dot plots for each patient, between parwise timepoints. 

Need to run twice, once for alpha, and once for beta. Produces plots; and a list of up and downregulated TCRs, which is saved as "TCR_change_all_alpha.RData" and "TCR_change_all_beta.RData"

**to_wide_data.R** converts the TCR_change_all_a/b into the format of exp_AB_wide.RData (i.e. has the abundances for all time points for the expanding TCRs)

**clean_wide_data.R** takes exp_AB_wide.RData and removes a few TCRs which are very high at -3 or -4 and then go down, because I think they are irrelevant. This provides  exp_AB_wide1.RData. The pGEN values were added manually to this. The version with manually added pGen is found in data/

On exp_AB_wide1.RData (from data/output_data/ and not from data/, so using the version without manually annotated pGen), you run **remove_invariant.R** to remove MAIT and IkT cells from the expanded list -> obtain exp_AB_wide3.RData

*The new axes here are shifted by a power of 2. When I plot with Benny's plotting without renaming axes, I get the same axes as on my new plots. Was there a reason Benny had shifted the axes by a power of 2?*

- [x] **1B** bulk_analysis_expanded.R calculates richness, Shannon diversity, Renyi etc from exp_AB_wide3.R, then saves it to a result file. Then plot with Fig1b.R

- [x] **1C, D** script_2.R, Section B

- [x] **1E, F** script_3.R - plots dynamics for each patient separately.

- [x] **2A, B, C** Numbers and stacked histogram in matches_summary.xlsx 

*I don't fully understand matches_summary.xlsx, so I left them as they were, just changed the font and the styling*

- [x] **2D** script_4.R Section C. 

script_4.R uses exp_AB_merge2.txt, which is the list of annotated TCRs. There is no “script” for generating this. It is constructed semi-manually from the annotation data from Tao as well as what was in VDJdb and the Franics et al. paper.

*script_4.R also uses summary_with_HLA.txt - just got most updated version from Benny, need to re-run because small difference.*

- [x] **2E** script_4.R, Section A. Save heatmap and legend separately.

- [x] **2F** script_5.R - calculates and plots clustering. Note: this clustering keeps duplicates, so that you see multiple nodes when the same cdr3 is present in multiple patients.

- [x] **Fig 3** - emerson1.R: uses the whole (unique) non-expanded set and the whole (unique) expanded set to compare sharing levels and precursor frequencies.

*Note: Emerson data does not fit on my DropBox nor on git, so I don't know where to save it*

- [x] **Fig 4** Scripts from Cedric to add.

- [ ] **Fig 5** Scripts? Data? *Benny to send*

*Found LCMV data on DropBox. 5B is plotted and looks correct, I still need to figure out 5C*

#### Supplementary

- [x] **S1** FACS Plots from Michal - *they are a little pixelated, is there another version?*

- [x] **S2** script_2.R, Section A - 2 separate heatmaps for PCR+ and PCR-

- [x] **S3** script_S3.R *this script used to be called script_S1.R, I changed it for consistency*

- [x] **S4** bulk_analysis_all.R, this takes all_A_long.Rdata and all_B_long.Rdata as input, and calculates some of the standard bulk metrics. The script then saves the results. This is equivalent to bulk_analysis_expanded.R. Then plot with FigS4.R.

- [x] **S5** script_1.R

- [x] **S6** script_2.R Section D. Script 2 also produces a summary file - summary_data.csv with patient information and some key metrics

- [x] **S7** script_3.R *the values are all 0 at week -2? Also looks slightly different from the figure Benny had before, but I don't know from what data the other one had been generated - looks like there were more lines*

- [x] **S8A, B** script_3.R

- [x] **S9** Plots from Yanchun Peng

- [x] **S10** script_4.R Section B

- [x] **S11** script_6.R - calculates clustering at various threshold and compares expanded vs non-expanded controls. This uses controls that are taken from unique CDR3s (controls_long.R) and compares them to the list of unique cdr3s that expand.

**controls_long.R** generates the controls from unique CDR3s. These are the controls used in the current figures - *adapted*

**contronls_long_1.R** generates the controls slightly differently, by taking unique dcb_id/HCW ID combinations, rather than unique cdr3s. 

- [x] **S12** script_5.R - calculates and plots clustering

- [x] **S13** Generated by TP

#### Important things to figure out

- **HCW195** - is it week -1 and week -2?

- Re-make controls to the way I make them for emerson and re-run clustering