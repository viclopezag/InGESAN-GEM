function essentialgenes_mtb(metabolic_model,cutoff,NGAM)

clc; 
%clear all;

addpath('model');
addpath('dataset');
addpath('medium');


dispstr = sprintf('In silico Gene Essentiality for %s',metabolic_model);
disp(dispstr)
time0 = cputime;



solverOK = changeCobraSolver('gurobi7','LP');
ZBAR_ES_Threshold = 0.9925; % ZBAR > 0.9925 ES genes  GUMBEL METHOD DE JESUS et al., 2013
ZBAR_NE_Threshold = 0.0493; % 0 < ZBAR < 0.0493 NE genes




%% READ GENOME SCALE MODELS
dispstr = sprintf('%5.1f second: reading network model with media constraints...',cputime-time0);
disp(dispstr)

model = load_model(metabolic_model,NGAM);
WT_FBA_solution = optimizeCbModel(model,'max');

%% SINGLE GENE DELETION

dispstr = sprintf('%5.1f second: In silico single gene deletion analysis...',cputime-time0);
disp(dispstr)
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
            Vector_gene_names{i} = text_data{j,1}; % Gene Names GUMBEL METHOD DOESN'T GIVE GENE NAMES
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


for i = 1: length(grRateKO)
    
    
    if grRateKO(i) <= cutoff*grRateWT
        
        InsilicoEssential{i} = 'E';  % Essential
        
    else 
        
        InsilicoEssential{i} = 'NE';  % Non-Essential
        
    end

end
       
    
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

%Vector_LocusNames(emptyCells)= [];
Vector_zbar_values(emptyCells) = [] ;
Vector_gene_names(emptyCells) = [] ;
%ExperimentalEssential(emptyCells) = [];
%InsilicoEssential(emptyCells) =  [];

%% ERASE UNCERTAIN 'U' AND 'S' GENES
Uindex = find(strcmp(ExperimentalEssential,'U'));% Finding Uncertain Genes Positions
Uncertain_genes = Vector_LocusNames(Uindex);
Sindex = find(strcmp(ExperimentalEssential,'S')); % Finding S Genes Positions
S_genes = Vector_LocusNames(Sindex);

%% COMPARING IN SILICO VERSUS EXPERIMENTAL DATA

dispstr = sprintf('%5.1f second: Comparing In Silico VS Experimental data...',cputime-time0);
disp(dispstr)

index_genes_not_compared = union([Sindex;Uindex],NOTFOUND_GENES_POSITION);
InsilicoEssential(index_genes_not_compared) = []; % Erase position of genes will not be compared
ExperimentalEssential(index_genes_not_compared) = []; % Erase Position of genes will not be compared
Vector_LocusNames(index_genes_not_compared) = []; % Erase position of genes will not be compared.

N_Genes2Compare = length(NOTFOUND_GENE)+length(Uncertain_genes)+length(S_genes);  % Number Genes that will not be compared

Size = length(model.genes)-N_Genes2Compare;

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


TP = strfind(Confusion_Matrix, 'TP');
TP = logical(cell2mat(TP)); % logical type
N_TP = length(TP);            % Number of true positive genes

TN = strfind(Confusion_Matrix, 'TN');
TN = logical(cell2mat(TN)); % logical type
N_TN = length(TN);            % Number of true negative genes

FP = strfind(Confusion_Matrix, 'FP');
FP = logical(cell2mat(FP)); % logical type
N_FP = length(FP);           % Number of false positive genes

FN = strfind(Confusion_Matrix, 'FN');
FN = logical(cell2mat(FN)); % logical type
N_FN = length(FN);           % Number of false negative genes.

	
%% DERIVATIONS FROM THE CONFUSION MATRIX

 % Defining Vectors of results

SENSITIVITY = N_TP/(N_TP+N_FN);
SPECIFICITY = N_TN/(N_FP+N_TN);
PPV = N_TP/(N_TP+N_FP);
NPV = N_TN/(N_TN+N_FN);
FALL_OUT = N_FP/(N_FP+N_TN);
FDR = 1-PPV; 
FNR = N_FN/(N_FN+N_TP);
ACCURACY = (N_TP+N_TN)/(N_TP+N_TN+N_FP+N_FN);
MCC = ((N_TP*N_TN)-(N_FP*N_FN))/(sqrt((N_TP+N_FP)*(N_TP+N_FN)*(N_TN+N_FP)*(N_TN+N_FN)));


%% GROWTH RATE
GROWTH_RATE = {'GrowthRate' WT_FBA_solution.f};
%%  METRICS VALUES

METRICS_VALUES = {'Model Genes' length(model.genes);'Evaluated Genes' Size; 'True Positives' N_TP; 'True Negatives' N_TN; 'False Negatives' N_FN;... 
                  'False Positives' N_FP;'Sensitivity' SENSITIVITY; 'Specificity' SPECIFICITY; 'Positive Predictive Value' PPV; 'Negative Predictive Value' NPV;...
                  'Fall out' FALL_OUT; 'False Discovery Rate' FDR; 'False Negative Rate' FNR; 'Accuracy' ACCURACY; 'Mathews Correlation Coefficient' MCC;...
                  'Uncertain Genes' length(Uncertain_genes); 'Too Small Genes' length(S_genes); 'Not found genes' length(NOTFOUND_GENE)};
%% GENE COMPARISON TABLE.
GENE_COMPARISON = cell(length(Vector_LocusNames),4);
GENE_COMPARISON(:,1) = Vector_LocusNames; 
GENE_COMPARISON(:,2) = InsilicoEssential;              
GENE_COMPARISON(:,3) = ExperimentalEssential;  
GENE_COMPARISON(:,4) = Confusion_Matrix;

Genecomparison_table = table (Vector_LocusNames, InsilicoEssential,...
                              ExperimentalEssential,Confusion_Matrix);


%% METRICS TABLE

LastNames = {'Model Genes';'Evaluated Genes'; 'True Positives'; 'True Negatives'; 'False Negatives';... 
                  'False Positives'; 'Sensitivity'; 'Specificity'; 'Positive Predictive Value'; 'Negative Predictive Value';...
                  'Fall out'; 'False Discovery Rate'; 'False Negative Rate'; 'Accuracy'; 'Mathews Correlation Coefficient';...
                  'Uncertain Genes'; 'Too Small Genes'; 'Not found genes'};
Metr_values = [length(model.genes); Size; N_TP; N_TN; N_FN;... 
            N_FP; SENSITIVITY; SPECIFICITY; PPV; NPV;...
            FALL_OUT; FDR; FNR; ACCURACY; MCC;...
            length(Uncertain_genes); length(S_genes); length(NOTFOUND_GENE)];

Metrics_table = table(LastNames, Metr_values);

%% Writting Files
system('taskkill /F /IM EXCEL.EXE');

writing_results_ctoff(Genecomparison_table,Metrics_table,metabolic_model,Size)
% Range = sprintf('A1:D%0.0f', length(Vector_LocusNames)+1);
% system('taskkill /F /IM EXCEL.EXE');
% filename_genecomp = sprintf('results/%s_genecompctoff.xlsx', METABOLIC_MODEL);
% writetable(Genecomparison_table,filename_genecomp,'Sheet',1,'Range',Range);
% system('taskkill /F /IM EXCEL.EXE');
% filename_metrics = sprintf('results/%s_metricsctoff.xlsx', METABOLIC_MODEL);
% writetable(Metrics_table,filename_metrics,'Sheet',1,'Range','A1:B19');
% system('taskkill /F /IM EXCEL.EXE');

end

        

