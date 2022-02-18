Data and figure generation files for manuscript: Spatially structured eco-evolutionary dynamics in a host-pathogen interaction render isolated populations vulnerable to disease

Notes: 
All files generated and tested on Matlab 2020b on MacOS 12.0.1.
".c" files will compile automatically in Matlab provided a suitable compiller is installed (if not, run mex -setup)
Typical install time < 1 min.
Expected runtime for typical output < 10 min

Instructions:
To reproduce simulation results from scratch, delete ".mat" files and run "theory_fig.m".

Files:
array_values.m - Sets values for constant arrays used in the simulations
fitness_cost_fig.m - Generates the fitness cost figure for the supplementary material
fitness_cost_fig.pdf - Fitness cost figure for the supplementary material
inf_function.c - Calculates the infection rates for the whole metapopulation
metapop_entry.m - Entry function for the metapopulation simulation
metapop_loop.m - Carries out the main loop for the metapopulation simulation
readme.txt - readme file
theory_fig.m - Generates the theory figure for the main text
theory_fig.pdf - Theory figure for the main text
theory_fig_analysis.m - Analyses simulation data and returns disease prevalence, resistance and infectivity
theory_fig_data.m - Carries out the main parameter sweep
theory_fig_dynamics_data.m - Analyses simulation data for the theory figure dynamics plots
theory_fig_dynamics_data_analysed.mat - Analysed simulation data
theory_fig_dynamics_data_raw.m - Renamed data from the parameter sweep for the main theory figure
theory_fig_dynamics_data_raw.mat - Raw simulation dynamics data for main theory figure
theory_fig_heatmaps_analysis.m - Analyses simulation data for the theory figure heatmaps (host and parasite costs)
theory_fig_heatmaps_analysis.mat - Analysed simulation data for the theory figure heatmaps 
theory_fig_snapshot_data.m - Generates simulation data for the metapopulation snapshot figure
theory_fig_snapshot_data.mat - Simulation data for the metapopulation snapshot figure
theory_supp_table_3_analysis.m - Analyses simulation data for all parameter sets to generate data for Supplementary Table 3
theory_supp_table_3_analysis.mat - Analysed simulation data for Supplementary Table 3