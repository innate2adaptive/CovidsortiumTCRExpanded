# CovidsortiumTCRExpanded
Scripts that produce data and figures for the COVIDsortium TCR paper (Aug 2022)

### Figures – data and scripts

Current version of MS is Word file TCR_COV2_V6.docx; Current main figures : Figures_V6.pptx; Current supplementary figures : Figures_S_V4.pptx 

*Figures_V6 and Figures_S_V5 are after correction of HCW 184 and 195 metadata. MS V6 is the one that incorporates comments from all authors.*

Data can be found at: https://www.dropbox.com/sh/neg3lpofw85l9ng/AAAw-HckR7LYxKz7TMNQ90zca?dl=0

Note: the scripts automatically load the data from DropBox. This works for me (although sometimes I need to change the timeout), but does not seem to always work outside of the UK. You might need to download the data locally and tweak the scripts to load them.

Scripts in: scripts/

Output plots in: output_figs/

#### Main figs
- [x] **1A** script_1.R, this script finds the TCRs that are changing between any two timepoints. Takes the original .tsv files as input, and outputs TCR_change_all_a/b. Also makes the dot plots for each patient, between parwise timepoints. 

Need to run twice, once for alpha, and once for beta. Produces plots; and a list of up and downregulated TCRs, which is saved as "TCR_change_all_alpha.RData" and "TCR_change_all_beta.RData"

- [x] **1B** bulk_analysis_expanded.R calculates richness, Shannon diversity, Renyi etc from exp_AB_wide3.R, then saves it to a result file. Then plot and calculate p-value with Fig1b.R. *Note: if bulk_analysis_expanded.R is crashing, restart the session*

- [x] **1C, D** script_2.R, Section B

- [x] **1E, F** script_3.R - plots dynamics for each patient separately.

- [x] **2A, B, C** annotated_overlap.R - looks at the overlap with expanded and with controls (generated with controls_long.R), plots the venn diagram and the bar plot.

*controls numbers generated with old controls to figure out code. Need to re-run with controls from myriad when done*

- [x] **2D** script_4.R Section C. 

script_4.R uses exp_AB_merge2.txt, which is the list of annotated TCRs. There is no “script” for generating this. It is constructed semi-manually from the annotation data from Tao as well as what was in VDJdb and the Franics et al. paper.

- [x] **2E** script_4.R, Section A. Save heatmap and legend separately.

- [x] **2F** script_5.R - calculates and plots clustering. Note: this clustering keeps duplicates, so that you see multiple nodes when the same cdr3 is present in multiple patients.

- [x] **Fig 3** - emerson1.R: uses the whole (unique) non-expanded set and the whole (unique) expanded set to compare sharing levels and precursor frequencies. This should be run after annotate_timepoint.R, which generates exp_AB_wide4.R. The emerson data is collated by merge_emerson_data.ipynb

*Note: Emerson data does not fit on my DropBox nor on git, so I don't know where to save it*

- [x] **Fig 4** Added newest version that CT provided 30/06/2022.

- [x] **Fig 5B** LCMV_analysis.R

- [x] **Fig 5C** LCMV_analysis_final.R

#### Supplementary

- [x] **S1** FACS Plots from Michal - *they are a little pixelated, is there another version?*

- [x] **S2** script_2.R, Section A - 2 separate heatmaps for PCR+ and PCR-

- [x] **S3** script_S3.R

- [x] **S4** bulk_analysis_all.R, this takes all_A_long.Rdata and all_B_long.Rdata as input, and calculates some of the standard bulk metrics. The script then saves the results. This is equivalent to bulk_analysis_expanded.R. Then plot with FigS4.R.

- [x] **S5** script_1.R

- [x] **S6** script_2.R Section D. Script 2 also produces a summary file - summary_data.csv with patient information and some key metrics

- [x] **S7** script_3.R

- [x] **S8A, B** script_3.R

- [x] **S9** Plots from Yanchun Peng

- [x] **S10** script_4.R Section B

- [x] **S11** script_6.R - calculates clustering at various threshold and compares expanded vs non-expanded controls. This uses controls that are taken from unique CDR3s (controls_long.R) and compares them to the list of unique cdr3s that expand.

- [x] **S12** script_5.R - calculates and plots clustering

- [] **S13** annotate_timepoint.R - calculates timepoint at which decombinator id is max in HCW data and then plots it

- [] **S14** calculate_pgen.ipynb - applies OLGA to calculate the pGen of expanded sequences; plot_pgen.R plots the results

- [x] **S13** Generated by TP

- [x] **S14** LCMV_analysis_final.R

#### Other scripts

**to_wide_data.R** converts the TCR_change_all_a/b into the format of exp_AB_wide.RData (i.e. has the abundances for all time points for the expanding TCRs), using PCR_positive.csv

**clean_wide_data.R** takes exp_AB_wide.RData and removes a few TCRs which are very high at -3 or -4 and then go down, because I think they are irrelevant. This provides  exp_AB_wide1.RData. 

On exp_AB_wide1.RData (from data/output_data/ and not from data/, so using the version without manually annotated pGen), you run **remove_invariant.R** to remove MAIT and IkT cells from the expanded list -> obtain exp_AB_wide3.RData.

On exp_AB_wide3.RData you run **annotate_timepoint.R**. This adds a column to annotated TCRs that come up early vs late, then saves exp_AB_wide4.RData.

On exp_AB_wide4.RData you run **calculate_pgen.ipynb**, which generates a separate files for alphas and betas with pgens. This same script also runs on the LCMV data and generates LCMV_pgen.csv.

**controls_long.R** generates the controls from unique CDR3s. These are the controls used in the current figures - does not run locally on my laptop, needs a large computer or a cluster.

**contronls_long_1.R** generates the controls slightly differently, by taking unique dcb_id/HCW ID combinations, rather than unique cdr3s. 