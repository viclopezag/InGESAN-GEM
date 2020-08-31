function ROCessentialgenes_mtb(metabolic_model, increments, maxcutoff, NGAM)

%Increments = 0.005;  % growth rate threshold increments
%maxcutoff = 0.995; % percentage of growth rate to determine essentiality

clc; 
% clear all;
addpath('model');
addpath('dataset');
addpath('medium');
addpath('results');


dispstr = sprintf('ROC curves Gene Essentiality for %s',metabolic_model);
disp(dispstr)
time0 = cputime;



%% SOLVER SET UP

% Options: 'UnifiedMTB_Griffin_cholesterol_nov2016', 'UnifiedMTB_Griffin_cholesterol_feb2017', 'iCG760_Griffin_cholesterol',
%          'iSM810_Griffin_Cholesterol', 'sMtb_Griffin_Cholesterol', 'iVL2017_iOSDD890_Griffin_Cholesterol',
%          'GSMN_TB_1.1_Griffin_Cholesterol';
%
solverOK = changeCobraSolver('gurobi7','LP');
ZBAR_ES_Threshold = 0.9925; % ZBAR > 0.9925 ES genes  GUMBEL METHOD DE JESUS et al., 2013
ZBAR_NE_Threshold = 0.0493; % 0 < ZBAR < 0.0493 NE genes
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

%% READ GENOME SCALE MODELS
dispstr = sprintf('%5.1f second: reading network model with media constraints...',cputime-time0);
disp(dispstr)

model = load_model(metabolic_model,NGAM);

%% SINGLE GENE DELETION

dispstr = sprintf('%5.1f second: In silico single gene deletion analysis...',cputime-time0);
disp(dispstr)
Fluxes = optimizeCbModel(model,'max');
[minFVA, maxFVA] = fluxVariability(model,100);
Rxns = model.rxns;
FVAmin = minFVA;
FVAmax = maxFVA;
fva = table(Rxns, FVAmin, FVAmax);

[grRatio,grRateKO,grRateWT,delRxns,hasEffect] = singleGeneDeletion(model,'FBA');

%% GENE ESSENTIAL DATABASE

dispstr = sprintf('%5.1f second: Loading Griffins Gene Essential Database...',cputime-time0);
disp(dispstr)

GEdatabase='H37Rv_cholesterol_griffin_GUMBEL_sum.xlsx'; % Esssentiality categorization of Cholesterol data (Griffin_2011) by the gumbel method of DeJesus 2013..

[number_data, text_data, mGEdatabase] = xlsread(GEdatabase,1);

Vector_LocusNames = cell(length(model.genes),1) ;
Vector_zbar_values = zeros(length(model.genes),1)  ;
Vector_gene_names = cell(length(model.genes),1) ;

for i = 1:length(model.genes)
    
    for j = 1:length(mGEdatabase)
        
            tf = isequal(model.genes{i,1}, text_data{j,1});
        
        if tf == 1
            
            Vector_LocusNames{i} = text_data{j,1}; % Locus
            Vector_gene_names{i} = text_data{j,2}; % Gene Names GUMBEL METHOD DOESN'T GIVE GENE NAMES
            Vector_zbar_values(i) = mGEdatabase{j,6};   % zbar value
        end 
        
              
    end
end

%% FINDING EMPTY MATCHES OF GENES

dispstr = sprintf('%5.1f second: Identify gene name mistakes in the model...',cputime-time0);
disp(dispstr)

emptyCells = cellfun(@isempty,Vector_LocusNames); % Set all empty elements as EMPTY

NOTFOUND_GENES_POSITION = find(emptyCells);
NOTFOUND_GENE = model.genes(NOTFOUND_GENES_POSITION);

%% IN SILICO ESSENTIALITY CATEGORIZATION

dispstr = sprintf('%5.1f second: Identifyng In silico Essentiality...',cputime-time0);
disp(dispstr)

InsilicoEssential = cell(length(grRateKO),1);

%% ARRAYS OF CALCULATIONS OF TP, TN, FP, FN

N_TP = zeros(ROC_points,1);           % Number array
N_TN = zeros(ROC_points,1);           % Number array
N_FP = zeros(ROC_points,1);           % Number array
N_FN = zeros(ROC_points,1);           % Number array

%% EXPERIMENTAL ESSENTIALITY CATEGORIZATION

dispstr = sprintf('%5.1f second: Identifyng Experimental Essentiality...',cputime-time0);
disp(dispstr)

ExperimentalEssential = cell(length(Vector_zbar_values),1);

for i = 1: length(Vector_zbar_values)
    
    if (Vector_zbar_values(i) > ZBAR_ES_Threshold ) 

    ExperimentalEssential{i} = 'E'; % gene is essential

    elseif (Vector_zbar_values(i) >= ZBAR_NE_Threshold) && (Vector_zbar_values(i) <= ZBAR_ES_Threshold)

    ExperimentalEssential{i} = 'U'; % gene is uncertain

    elseif (Vector_zbar_values(i) >= 0) && (Vector_zbar_values(i) < ZBAR_NE_Threshold)

    ExperimentalEssential{i} = 'NE'; % gene is non-essential
    
    else
    
    ExperimentalEssential{i} = 'S'; % S too few TA sites in the gene
    
    end
    
end 

%% ERASE NOT MATCHES OF GENES

dispstr = sprintf('%5.1f second: Deleting not matched genes...',cputime-time0);
disp(dispstr)

Vector_zbar_values(emptyCells) = [] ;
Vector_gene_names(emptyCells) = [] ;
% ExperimentalEssential(emptyCells) = [];
% InsilicoEssential(emptyCells) =  [];

%% ERASE UNCERTAIN 'U' AND 'S' GENES and erase Genes not to compare
Uindex = find(strcmp(ExperimentalEssential,'U'));% Finding Uncertain Genes Positions
Uncertain_genes = Vector_LocusNames(Uindex);
Sindex = find(strcmp(ExperimentalEssential,'S')); % Finding S Genes Positions
S_genes = Vector_LocusNames(Sindex);

index_genes_not_compared = union([Sindex;Uindex],NOTFOUND_GENES_POSITION);
ExperimentalEssential(index_genes_not_compared) = []; % Erase Position of genes will not be compared
Vector_LocusNames(index_genes_not_compared) = []; % Erase position of genes will not be compared.

N_Genes2Compare = length(NOTFOUND_GENE)+length(Uncertain_genes)+length(S_genes);  % Number of TOtal Genes to compare essentiality

Size = length(model.genes)-N_Genes2Compare;

for z = 1 : ROC_points
 
    GR_Threshold_c(z)= increments*z;
    
    
 
for i = 1: length(grRateKO)
    
    
    if grRateKO(i) <= GR_Threshold_c(z)*grRateWT
        
        InsilicoEssential{i} = 'E';  % Essential
        
    else 
        
        InsilicoEssential{i} = 'NE';  % Non-Essential
        
    end

end
       
InsilicoEssential(index_genes_not_compared) = []; % Erase position of genes will not be compared



%% COMPARING IN SILICO VERSUS EXPERIMENTAL DATA

dispstr = sprintf('%5.1f second: Comparing In Silico VS Experimental data...',cputime-time0);
disp(dispstr)


Confusion_Matrix = cell(length(Size),1);


for i = 1:Size
    
    %all(size('Word1') == size('Word2')) && all('Word1' == 'Word2')
	tf = isequal(InsilicoEssential{i},ExperimentalEssential{i});

	idiS_NE = strfind(InsilicoEssential(i),'NE'); % cell type 
	idiS_NE = logical(cell2mat(idiS_NE)); % logical type

	idiS_E = strfind(InsilicoEssential(i),'E'); %cell type
    idiS_E = logical(cell2mat(idiS_E)); % logical type

	idexp_NE = strfind(ExperimentalEssential(i),'NE');
	idexp_NE = logical(cell2mat(idexp_NE)); % logical type

	idexp_E = strfind(ExperimentalEssential(i),'E');
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
TP(:,z) = strfind(Confusion_Matrix, 'TP');
%TP(:,z) = logical(cell2mat(TP)); % logical type
TP = logical(cell2mat(TP)); % logical type
N_TP(z) = length(TP);       % Number of true positive genes

TN = cell(length(Confusion_Matrix),ROC_points);
TN(:,z) = strfind(Confusion_Matrix, 'TN');
%TN(:,z) = logical(cell2mat(TN(:,z))); % logical type
TN = logical(cell2mat(TN));
N_TN(z) = length(TN);            % Number of true negative genes

FP = cell(length(Confusion_Matrix),ROC_points);
FP(:,z) = strfind(Confusion_Matrix, 'FP');
%FP(:,z) = logical(cell2mat(FP(:,z))); % logical type
FP = logical(cell2mat(FP));
N_FP(z) = length(FP);           % Number of false positive genes

FN = cell(length(Confusion_Matrix),ROC_points);
FN(:,z) = strfind(Confusion_Matrix, 'FN');
%FN(:,z) = logical(cell2mat(FN(:,z))); % logical type
FN = logical(cell2mat(FN));
N_FN(z) = length(FN);           % Number of false negative genes.

	
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

