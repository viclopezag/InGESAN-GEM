%% Function to write xlsx files of ROCessentialgenes_mtb
% 18 December 2018 - Víctor López University of Antioquia -
% valonso.lopez@udea.edu.co.
function writing_results_roc(metrics_roc, fva, metabolic_model,number_of_rxns,AUC)

system('taskkill /F /IM EXCEL.EXE');
range_fva = sprintf('A1:C%0.0f',number_of_rxns+1);
system('taskkill /F /IM EXCEL.EXE');
filename_fva = sprintf('results/%s_FVA.xlsx', metabolic_model);
writetable(fva,filename_fva,'Sheet',1,'Range',range_fva);
system('taskkill /F /IM EXCEL.EXE');
filename_metrics_roc = sprintf('results/%s_metrics_roc.xlsx', metabolic_model);
writetable(metrics_roc,filename_metrics_roc,'Sheet',1,'Range','A1:N201')
system('taskkill /F /IM EXCEL.EXE');
AUC_filename = sprintf('results/%s_AUC.mat', metabolic_model);
save(AUC_filename,'AUC')

end 