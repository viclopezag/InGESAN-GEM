%% Main Gene Essentiality Analysis to Mtb genome scale models 
% July 09 2019 
% Víctor López - University of Antioquia - valonso.lopez@udea.edu.co
% Minato MtbYM essentiality.
clc;
clear;

%% Non Growth Associated Maintenance

NGAM = 1.0 ;%mmol/gDW/h


%% Define models to be evaluated.

% metabolic_models = {'sMtb_Minato_MtbYM', 'iEK1011_Minato_MtbYM',...
%         'iCG760_Minato_MtbYM', 'iSM810_Minato_MtbYM',...
%         'GSMN_TB_1.1_Minato_MtbYM','sMtb2018_Minato_MtbYM',...
%         'iOSDD890_Minato_MtbYM', 'iNJ661v_modified_Minato_MtbYM',...
%          'iEK1011_2.0_Minato_MtbYM', 'sMtb_2.0_Minato_MtbYM'};
 metabolic_models = {'sMtb_2.0_Minato_MtbYM'} ;   

%% Call essentialgenes_mtb function with 0.05 cutoff in growthRates to 
% classify essential or non essential genes.
cutoff = 0.05;
for i = 1:length(metabolic_models)
    essentialgenes_mtb(metabolic_models{i},cutoff,NGAM)
end

%% Call ROCessentialgenes_mtb function
clear all;
increments = 0.005;  % growth rate threshold increments
maxcutoff = 0.995; % percentage of growth rate to determine essentiality
metabolic_models = {'GSMN_TB_1.1_Minato_MtbYM','sMtb2018_Minato_MtbYM',...
                    'iOSDD890_Minato_MtbYM', 'iNJ661v_modified_Minato_MtbYM'};

for i = 1:length(metabolic_models)
    ROCessentialgenes_mtb(metabolic_models{i},increments,maxcutoff,NGAM)
end



