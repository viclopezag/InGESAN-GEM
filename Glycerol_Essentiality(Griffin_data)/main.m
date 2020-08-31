%% Main Gene Essentiality Analysis to Mtb genome scale models 
% December 18 2018 
% Víctor López - University of Antioquia - valonso.lopez@udea.edu.co
% Griffin Cholesterol essentiality.
clc;
clear;

%% Non Growth Associated Maintenance

NGAM = 1.0 ;%mmol/gDW/h


%% Define models to be evaluated.

% metabolic_models = {'sMtb_Griffin_Glycerol', 'iEK1011_Griffin_Glycerol',...
%         'iCG760_Griffin_Glycerol', 'iSM810_Griffin_Glycerol',...
%         'GSMN_TB_1.1_Griffin_Glycerol','sMtb2018_Griffin_Glycerol',...
%         'iOSDD890_Griffin_Glycerol', 'iNJ661v_modified_Griffin_Glycerol',...
%          'iEK1011_2.0_Griffin_Glycerol', 'sMtb2.0};
metabolic_models = {'sMtb_2.0_Griffin_Glycerol'};
%     
%% Call essentialgenes_mtb function with 0.05 cutoff in growthRates to 
% classify essential or non essential genes.
cutoff = 0.05;
for i = 1:length(metabolic_models)
    essentialgenes_mtb(metabolic_models{i},cutoff,NGAM)
end

%% Call ROCessentialgenes_mtb function with 
clear all;
increments = 0.005;  % growth rate threshold increments
maxcutoff = 0.995; % percentage of growth rate to determine essentiality
metabolic_models = {'sMtb_Griffin_Glycerol', 'iEK1011_Griffin_Glycerol',...
        'iCG760_Griffin_Glycerol', 'iSM810_Griffin_Glycerol',...
        'GSMN_TB_1.1_Griffin_Glycerol','sMtb2018_Griffin_Glycerol',...
        'iOSDD890_Griffin_Glycerol', 'iNJ661v_modified_Griffin_Glycerol'};

for i = 1:length(metabolic_models)
    ROCessentialgenes_mtb(metabolic_models{i},increments,maxcutoff,NGAM)
end



