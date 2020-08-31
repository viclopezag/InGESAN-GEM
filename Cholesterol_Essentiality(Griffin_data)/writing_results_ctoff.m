%% Function to write xlsx files of essentialgenes_mtb
% 18 December 2018 - Víctor López University of Antioquia -
% valonso.lopez@udea.edu.co.
function writing_results_ctoff(genecomparison_table,metrics_table, metabolic_model,number_of_genes)
system('taskkill /F /IM EXCEL.EXE');
range_genes = sprintf('A1:D%0.0f', number_of_genes+1);
system('taskkill /F /IM EXCEL.EXE');
filename_genecomp = sprintf('results/%s_genecompctoff.xlsx', metabolic_model);
writetable(genecomparison_table,filename_genecomp,'Sheet',1,'Range',range_genes);
system('taskkill /F /IM EXCEL.EXE');
filename_metrics = sprintf('results/%s_metricsctoff.xlsx', metabolic_model);
writetable(metrics_table,filename_metrics,'Sheet',1,'Range','A1:B19');
system('taskkill /F /IM EXCEL.EXE');
end 