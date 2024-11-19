This repository is a companion to the manuscript

N. Ducci, L. Grilli, M. Pittavino (2024)
"Comparing flexible modelling approaches: the varying-thresholds model versus quantile regression" 

The repository is divided in three main folder: 1) 00_functions, 2) 01_simulations, 3) 02_charts.

The folder "OO_functions" contains three R scripts with the codes for generating the related functions named in the script:
a) interval_score.R
b) parametric_splits_quantiles.R
c) varying_thresholds_model.R

The folder "01_simulations" contains three subfolders for the simulations of the CRPR and the models 1 and 8, with the two different coverages at 80 and 95%. 
In particular, there are:
1) CRPR_SIMULATIONS, with two R scripts: model2_crps_simulation.R and model8_crps_simulation.R
2) MODEL_1_8_COVERAGE_80, with nine R scripts for the different models' simulations, named respectively:
2a) model1(norm_error)_model2(chisq_error)_cloglogLink.R
2b) model1(norm_error)_model2(chisq_error)_probitLink.R
2c) model3_model4(heterosc.)_probitLink.R
2d) model5(quadratic_effect)_probitLink.R
2e) model6(cubic_effect)_probitLink.R
2f) model6(cubic_effect)_probitLink_sampleSize100.R
2g) model6(cubic_effect)_probitLink_sampleSize500.R
2h) model7(cont_norm_error)_probitLink.R
2i) model8(student_error)_probitLink_robitLink.R
3) MODEL_1_8_COVERAGE_90, with nine R scripts for the different models' simulations for the new coverage, named respectively:
3a) model1(norm_error)_model2(chisq_error)_cloglogLink.R
3b) model1(norm_error)_model2(chisq_error)_probitLink.R
3c) model3_model4(heterosc.)_probitLink.R
3d) model5(quadratic_effect)_probitLink.R
3e) model6(cubic_effect)_probitLink.R
3f) model6(cubic_effect)_probitLink_sampleSize100.R
3g) model6(cubic_effect)_probitLink_sampleSize500.R
3h) model7(cont_norm_error)_probitLink.R
3i) model8(student_error)_probitLink_robitLink.R

The folder "02_charts" contains one subfolder "charts_in_the_paper" with the related graphs, named afterwards, and eight R scripts for the graphs generations:
1) The subfolder "charts_in_the_paper" contains five PNG files with updated graphs: 
1a) comparison_chisq_scatter.png
1b) comparison_normal_beta.png
1c) comparison_normal_scatter.png
1d) comparison_scatter_cubic.png
1e) scatter_plot_het2png
2) The eight R scripts for the graphs generations are included here, with the following names:
2a) CHARTS_model1(normal_error).R
2b) CHARTS_model2(chisq_error).R
2c) CHARTS_model3(heterosc.).R
2d) CHARTS_model4(heterosc.).R
2e) CHARTS_model5(quadratic).R
2f) CHARTS_model6(cubic).R
2g) CHARTS_model7(outlier).R
2h) CHARTS_model8(tstudent_error).R