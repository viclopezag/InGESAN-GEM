%% Main Gene Essentiality Analysis to Mtb genome scale models 
% October 28 2019 
% Víctor López - University of Antioquia - valonso.lopez@udea.edu.co
% Griffin Cholesterol essentiality.
clc;
clear;

%% non-growth associated maintenance to 0.1 or 1.0 or 3.8 mmol gdw-1 h-1 as suggested by Rienksma et al., 2014

NGAM = 1.0; 
%% Define models to be evaluated.
 

metabolic_models = {'sMtb_Griffin_Cholesterol', 'iEK1011_Griffin_Cholesterol',...
        'iCG760_Griffin_cholesterol', 'iSM810_Griffin_Cholesterol',...
        'GSMN_TB_1.1_Griffin_Cholesterol','sMtb2018_Griffin_Cholesterol','iEK1011_2.0_Griffin_Cholesterol',...
        'sMtb2.0_Griffin_Cholesterol'};

%% Call essentialgenes_mtb function with 0.05 cutoff in growthRates to 
% classify essential or non essential genes.
cutoff = 0.05;
for i = 1:length(metabolic_models)
    essentialgenes_mtb(metabolic_models{i},cutoff, NGAM)
end

%% Call ROCessentialgenes_mtb function with 
clear;
NGAM = 1.0;
increments = 0.005;  % growth rate threshold increments
maxcutoff = 0.995; % percentage of growth rate to determine essentiality
metabolic_models = {'sMtb_Griffin_Cholesterol', 'iEK1011_Griffin_Cholesterol',...
        'iCG760_Griffin_cholesterol', 'iSM810_Griffin_Cholesterol',...
        'GSMN_TB_1.1_Griffin_Cholesterol','sMtb2018_Griffin_Cholesterol','iEK1011_2.0_Griffin_Cholesterol',...
        'sMtb2.0_Griffin_Cholesterol'};
for i = 1:length(metabolic_models)
    ROCessentialgenes_mtb(metabolic_models{i},increments,maxcutoff, NGAM)
end



