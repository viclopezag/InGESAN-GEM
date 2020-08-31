function DJROCessentialgenes_mtb(metabolic_model, increments, maxcutoff, NGAM)

%% SINGLE GENE DELETION STUDY - GRIFFIN CHOLESTEROL DATA ANALYZED BY THE GUMBEL METHOD  
%%21 Septiembre 2017 - VICTOR A.LOPEZ-AGUDELO University of Antioquia
%% UNIVERSITY OF SURREY - DANY BESTE 

clc; 
%clear all;

addpath('model');
addpath('dataset');
addpath('medium');


dispstr = sprintf('ROC curves Gene Essentiality for %s',metabolic_model);
disp(dispstr)
time0 = cputime;

%% SOLVER SET UP

% metabolic_model: 'iCG760_7H9_10_OADC',
%                  'iSM810_7H9_10_OADC', 'sMtb_7H9_10_OADC', 'iOSDD890_7H9_10_OADC',
%                  'iNJ661v_modified_7H9_10_OADC';
%                  'GSMN_TB_1.1_7H9_10_OADC';
%                  'sMtb1.9_7H9_10_OADC'
%                  'iNJ661_7H9_10_OADC'
%                  'iNJ661v_7H9_10_OADC'
%                  'sMtb2.0_7H9_10_OADC'
%                  'sMtb2.0_TICS_7H9_10_OADC'
%                  'iEK1011__7H9_10_OADC'

%METABOLIC_MODEL = 'sMtb2.0_TICS_7H9_10_OADC';

solverOK = changeCobraSolver('gurobi7','LP');
%Increments = 0.01;  % growth rate threshold increments
%GR_Threshold = 0.99; % percentage of growth rate to determine essentiality
ROC_points = maxcutoff/increments;  % Number of Point for the ROC Curve per model.
GR_Threshold_c = zeros(ROC_points,1);

%% VECTORS OF RESULTS FOR ROC CURVES

 SENSITIVITY = zeros(ROC_points,1);
 SPECIFICITY = zeros(ROC_points,1);
 PPV = zeros(ROC_points,1);
 NPV = zeros(ROC_points,1);
 FALL_OUT = zeros(ROC_points,1);
 FDR = zeros(ROC_points,1);
 FNR = zeros(ROC_points,1);
 ACCURACY = zeros(ROC_points,1);
 MCC = zeros(ROC_points,1);
 PRESICION = zeros(ROC_points,1);
 RECALL = zeros(ROC_points,1);
%% READ GENOME SCALE MODELS
dispstr = sprintf('%5.1f second: reading network model with media constraints...',cputime-time0);
disp(dispstr)

model = load_model(metabolic_model, NGAM);
WT_FBA_solution = optimizeCbModel(model,'max');
[minFVA, maxFVA] = fluxVariability(model,100);
Rxns = model.rxns;
FVAmin = minFVA;
FVAmax = maxFVA;
fva = table(Rxns, FVAmin, FVAmax);
%% SINGLE GENE DELETION

dispstr = sprintf('%5.1f second: In silico single gene deletion analysis...',cputime-time0);
disp(dispstr)
Fluxes = optimizeCbModel(model,'max');
[grRatio,grRateKO,grRateWT,delRxns,hasEffect] = singleGeneDeletion(model,'FBA');

%% GENE ESSENTIAL DATABASE

dispstr = sprintf('%5.1f second: Loading Griffins Gene Essential Database...',cputime-time0);
disp(dispstr)

GEdatabase='2017_Glycerol_Essentiality_deJesus.xlsx'; % Esssentiality categorization of Cholesterol data (Griffin_2011) by the gumbel method of DeJesus 2013..

[number_data, text_data, mGEdatabase] = xlsread(GEdatabase,1);

Vector_LocusNames = cell(length(model.genes),1) ;
ExperimentalEssential= cell(length(model.genes),1)  ;
Vector_gene_names = cell(length(model.genes),1) ;

for i = 1:length(model.genes)
    
    for j = 1:length(mGEdatabase)
        
            tf = isequal(model.genes{i,1}, text_data{j,1});
        
        if tf == 1
            
            Vector_LocusNames{i} = text_data{j,1}; % Locus
            Vector_gene_names{i} = text_data{j,2}; % Gene Names GUMBEL METHOD DOESN'T GIVE GENE NAMES
            ExperimentalEssential{i} = mGEdatabase{j,4};   % Call
        
        end 
        
              
    end
end

%% FINDING EMPTY MATCHES OF GENES

dispstr = sprintf('%5.1f second: Identify gene name mistakes in the model...',cputime-time0);
disp(dispstr)

emptyCells = cellfun(@isempty,Vector_LocusNames); % Set all empty elements as EMPTY

NOTFOUND_GENES_POSITION = find(emptyCells);
NOTFOUND_GENE = model.genes(NOTFOUND_GENES_POSITION);

%% ERASE UNCERTAIN 'U' GENES
Uindex = find(strcmp(ExperimentalEssential,'U'));% Finding Uncertain Genes Positions
Uncertain_genes = Vector_LocusNames(Uindex);

%% ELIMINATE UNCERTAIN AND NOT FOUND GENES
index_genes_not_compared = union(Uindex,NOTFOUND_GENES_POSITION);
ExperimentalEssential(index_genes_not_compared) = []; % Erase Position of genes will not be compared
Vector_LocusNames(index_genes_not_compared) = []; % Erase position of genes will not be compared.
Vector_gene_names(index_genes_not_compared) = []; % Erase positions of genes will not be compared
grRateKO(index_genes_not_compared) = []; 

N_Genes2Compare = length(NOTFOUND_GENE)+length(Uncertain_genes); 
%% IN SILICO ESSENTIALITY CATEGORIZATION

dispstr = sprintf('%5.1f second: Identifyng In silico Essentiality...',cputime-time0);
disp(dispstr)

InsilicoEssential = cell(length(grRateKO),1);
 %% ARRAYS OF CALCULATIONS

N_TP = zeros(ROC_points,1);           % Number array
N_TN = zeros(ROC_points,1);           % Number array
N_FP = zeros(ROC_points,1);           % Number array
N_FN = zeros(ROC_points,1);           % Number array

for z = 1 : ROC_points
 
    GR_Threshold_c(z)= increments*z;
    
    
for i = 1: length(grRateKO)
    
    
    if grRateKO(i) <= GR_Threshold_c(z)*grRateWT
        
        InsilicoEssential{i} = 'ES';  % Essential
        
    else 
        
        InsilicoEssential{i} = 'NE';  % Non-Essential
        
    end

end
       
    

%% ERASE NOT MATCHES OF GENES

dispstr = sprintf('%5.1f second: Deleting not matched genes...',cputime-time0);
disp(dispstr)

%% COMPARING IN SILICO VERSUS EXPERIMENTAL DATA

dispstr = sprintf('%5.1f second: Comparing In Silico VS Experimental data...',cputime-time0);
disp(dispstr)

%% COMPARING IN SILICO VERSUS EXPERIMENTAL DATA

dispstr = sprintf('%5.1f second: Comparing In Silico VS Experimental data...',cputime-time0);
disp(dispstr)

 % Number of TOtal Genes to compare essentiality

Size = length(model.genes)-N_Genes2Compare;

Confusion_Matrix = cell(length(Size),1);


for i = 1:Size
    
    %all(size('Word1') == size('Word2')) && all('Word1' == 'Word2')
	tf = isequal(InsilicoEssential{i},ExperimentalEssential{i});

	idiS_NE = strfind(InsilicoEssential(i),'NE'); % cell type 
	idiS_NE = logical(cell2mat(idiS_NE)); % logical type

	idiS_E = strfind(InsilicoEssential(i),'ES'); %cell type
    idiS_E = logical(cell2mat(idiS_E)); % logical type

	idexp_NE = strfind(ExperimentalEssential(i),'NE');
	idexp_NE = logical(cell2mat(idexp_NE)); % logical type

	idexp_E = strfind(ExperimentalEssential(i),'ES');
	idexp_E = logical(cell2mat(idexp_E)); % logical type
    
    idexp_nd=strfind(ExperimentalEssential(i),'no-data'); % Datos TraSH
	idexp_nd=logical(cell2mat(idexp_nd)); % logical type


    if  tf==1 

    	if idiS_NE==1
		   Confusion_Matrix{i,1}='TN';
		else 
		   Confusion_Matrix{i,1}='TP';
		end
    
    elseif tf==0

        if idexp_E==1
		  Confusion_Matrix{i,1}='FN';
		end

	    if idexp_NE==1
		  Confusion_Matrix{i,1}='FP';
		end

		if idexp_nd==1 
		  Confusion_Matrix{i,1}='no-data';
	    end

	end


end 

%% CONFUSION MATRIX
dispstr = sprintf('%5.1f second: Confusion matrix analysis...',cputime-time0);
disp(dispstr)


TP = cell(length(Confusion_Matrix),ROC_points);
TP (:,z) = strfind(Confusion_Matrix, 'TP');
TP = logical(cell2mat(TP)); % logical type
N_TP(z) = length(TP);            % Number of true positive genes

TN = cell(length(Confusion_Matrix),ROC_points);
TN(:,z) = strfind(Confusion_Matrix, 'TN');
TN = logical(cell2mat(TN)); % logical type
N_TN(z) = length(TN);            % Number of true negative genes

FP = cell(length(Confusion_Matrix), ROC_points);
FP(:,z) = strfind(Confusion_Matrix, 'FP');
FP = logical(cell2mat(FP)); % logical type
N_FP(z) = length(FP);            % Number of false positive genes

FN = cell(length(Confusion_Matrix),ROC_points);
FN(:,z) = strfind(Confusion_Matrix, 'FN');
FN = logical(cell2mat(FN)); % logical type
N_FN(z) = length(FN);    

	
%% DERIVATIONS FROM THE CONFUSION MATRIX

 % Defining Vectors of results

SENSITIVITY(z) = N_TP(z)/(N_TP(z)+N_FN(z));
SPECIFICITY(z) = N_TN(z)/(N_FP(z)+N_TN(z));
PPV(z) = N_TP(z)/(N_TP(z)+N_FP(z));
NPV(z) = N_TN(z)/(N_TN(z)+N_FN(z));
FALL_OUT(z) = N_FP(z)/(N_FP(z)+N_TN(z));
FDR(z) = 1-PPV(z); 
FNR(z) = N_FN(z)/(N_FN(z)+N_TP(z));
ACCURACY(z) = (N_TP(z)+N_TN(z))/(N_TP(z)+N_TN(z)+N_FP(z)+N_FN(z));
MCC(z) = ((N_TP(z)*N_TN(z))-(N_FP(z)*N_FN(z)))/(sqrt((N_TP(z)+N_FP(z))*(N_TP(z)+N_FN(z))*(N_TN(z)+N_FP(z))*(N_TN(z)+N_FN(z))));
PRESICION(z) = N_TP(z)/(N_TP(z)+N_FP(z));
RECALL(z) = N_TP(z)/(N_TP(z)+N_FN(z));

end

%%  METRICS VALUES
metrics_roc = table(GR_Threshold_c, N_TN, N_TP, N_FN, N_FP,...
                 SENSITIVITY, SPECIFICITY, PPV, NPV,...
                 FALL_OUT, FDR, FNR, ACCURACY, MCC);

%% Area Under the Curve
SENSITIVITY = [0; SENSITIVITY;1];
SPECIFICITY = [1; SPECIFICITY;0];

one = ones(length(SPECIFICITY),1);
AUC = trapz((one-SPECIFICITY),SENSITIVITY);

%% Writing Files
system('taskkill /F /IM EXCEL.EXE');

writing_results_roc(metrics_roc,fva, metabolic_model,length(model.rxns),AUC)

end 