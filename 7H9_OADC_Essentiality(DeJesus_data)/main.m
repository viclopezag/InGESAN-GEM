%% Main Gene Essentiality Analysis to Mtb genome scale models 
% December 18 2018 
% Víctor López - University of Antioquia - valonso.lopez@udea.edu.co
% Griffin Cholesterol essentiality.
clc;
clear;

%% Non Growth Associated Maintenance

NGAM = 1.0; %mmol/gDW/h

%% Define models to be evaluated.

 metabolic_models = {'sMtb_7H9_10_OADC', 'iEK1011_7H9_10_OADC',...
         'iCG760_7H9_10_OADC', 'iSM810_7H9_10_OADC',...
         'GSMN_TB_1.1_7H9_10_OADC','sMtb2018_7H9_10_OADC',...
         'iOSDD890_7H9_10_OADC', 'iNJ661v_modified_7H9_10_OADC',...
         'iEK1011_2.0_7H9_10_OADC', 'sMtb2.0_TICS_7H9_10_OADC'};

%metabolic_models = {'sMtb2018_7H9_10_OADC'};
    
%% Call essentialgenes_mtb function with 0.05 cutoff in growthRates to 
% classify essential or non essential genes.
cutoff = 0.05;
for i = 1:length(metabolic_models)
    DJessentialgenes_mtb(metabolic_models{i},cutoff, NGAM)
end

%% Call ROCessentialgenes_mtb function with 
clear;
NGAM = 1.0
increments = 0.005;  % growth rate threshold increments
maxcutoff = 0.995; % percentage of growth rate to determine essentiality
metabolic_models = {'sMtb_7H9_10_OADC', 'iEK1011_7H9_10_OADC',...
        'iCG760_7H9_10_OADC', 'iSM810_7H9_10_OADC',...
        'GSMN_TB_1.1_7H9_10_OADC','sMtb2018_7H9_10_OADC',...
        'iOSDD890_7H9_10_OADC', 'iNJ661v_modified_7H9_10_OADC'};
    
for i = 1:length(metabolic_models)
    DJROCessentialgenes_mtb(metabolic_models{i},increments,maxcutoff)
end



