function model = load_model(metabolic_model,NGAM)
% Load a genome scale model for the selected organism.
% Options:
% metabolic_model: 'sMtb2.0_Griffin_Cholesterol', 'iCG760_Griffin_cholesterol', 'iSM810_Griffin_Cholesterol',
%                  'GSMN_TB_1.1_Griffin_Cholesterol', 'sMtb2018_Griffin_Cholesterol', 'iEK1011_Griffin_Cholesterol',
%                  'sMtb_Griffin_Cholesterol', 'iEK1011_2.0_Griffin_Cholesterol';
%
% Author: Victor Lopez, September 18/12/2018.

%% non-growth associated maintenance to 0.1 mmol gdw-1 h-1 as suggested by Rienksma et al., 2014
%NGAM = 1;

%% MODEL PATH
    
    GARAY_MODEL = 'model/iCG760_grRuleCorrected.mat'; 
    MA_MODEL = 'model/iSM810_grRuleCorrected.mat';
    LOFTHOUSE_MODEL = 'model/GSMN_TB1.1_grRuleCorrected.mat';
    SMTB_2_0 = 'model/sMtb2.0.mat';
    RIENKSMA_2014 = 'model/sMtb_grRuleCorrected.mat';
    RIENKSMA_2018 = 'model/sMtb2018_EX.mat';
    KAVVAS_MODEL = 'model/iEK1011_xlsx.mat';
    KAVVAS_2_0 = 'model/iEK1011_2.0.mat';
  
%% MEDIUM PATH

    MA_MEDIUM = 'medium/ExchangeRxns_MA2015.mat';
    LOFTHOUSE_MEDIUM = 'medium/ExchangeRxns_MA2015.mat';   

%% CHOICES    
    
switch metabolic_model
    
            case 'sMtb2.0_Griffin_Cholesterol'
            
                    load(SMTB_2_0);
                    model = sMtb2_0;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                     Medium_Griffin_Cholesterol={'EX_ASN[e]','EX_CL[e]','EX_PI[e]','EX_H[e]','EX_K[e]','EX_NA[e]','EX_ZN[e]','EX_MG[e]', 'EX_CA[e]','EX_FE2[e]','EX_FE3[e]','EX_H2O[e]','EX_O2[e]','EX_NH3[e]','EX_SLF[e]','EX_CIT[e]','EX_CHOLESTEROL[e]', 'EX_ETH[e]'};  
                    Griffin_Cholesterol_data = [ 1, 1000, 1, 1000, 1000, 1000, 1000, 1000, 1000, 5, 5, 1000, 20, 10, 1000, 1, 1, 1];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol,-Griffin_Cholesterol_data,'l');
                    %model = changeRxnBounds(model,'EX_CHOLESTEROL[e]',-1,'b'); % Cannot consume cholesterol when LB is assigned.
                    model = changeRxnBounds(model,'Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.

                
             case 'iCG760_Griffin_cholesterol'
                    
                    load(GARAY_MODEL);
                    model = iCG760;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Cholesterol={'EX_ASN(e)','EX_PI(e)','EX_H(e)','EX_FE3(e)','EX_O2(e)','EX_NH3(e)','EX_SLF(e)','EX_CIT(e)','EX_CHOLESTEROL(e)', 'EX_ETH(e)'};
                    Griffin_Cholesterol_data = [ 1, 1, 1000, 5, 20, 10, 1000, 1, 1, 1];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol,-Griffin_Cholesterol_data,'l');
                    %model = changeRxnBounds(model,'EX_CHOLESTEROL(e)',-1,'b');
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % Mainteinence was fixt to 1 mmol/gDW/h
                    
    
             case 'iSM810_Griffin_Cholesterol'
                               
                    load (MA_MODEL);
                    model = iSM810;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000; 
                    load (MA_MEDIUM);
                    model = changeRxnBounds(model,ExchangeRxns,0,'u');
                    Medium_Griffin_Cholesterol = { 'R800', 'R804', 'R932', 'R822', 'R838', 'R841', 'R882', 'R924', 'R858' };                   
                    Griffin_Cholesterol_data = [ 10, 20, 1, 1, 1, 1000, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol, Griffin_Cholesterol_data,'u');
                    %model = changeRxnBounds(model,'R932',1,'b'); % Cannot consume cholesterol with the UB constraint, so I needed to fix it with 'b'.
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % Maintainance was fix to 1 mmol/gDWH
                   

              case 'GSMN_TB_1.1_Griffin_Cholesterol'
                    
                    load(LOFTHOUSE_MODEL);
                    model = GSMN_TB1_1;
                    model.ub(model.ub > 1) = 1000;
                    model.lb(model.lb < 0) = -1000;
                    load(LOFTHOUSE_MEDIUM);
                    model = changeRxnBounds(model,ExchangeRxns,0,'u');
                    Medium_Griffin_Cholesterol = {'R800', 'R804', 'R932', 'R822', 'R838', 'R841', 'R882', 'R924', 'R858'};                   
                    Griffin_Cholesterol_data = [ 10, 20, 1, 1, 1, 1000, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol, Griffin_Cholesterol_data,'u');
                   % model = changeRxnBounds(model,'R932',1,'u'); % not feasible solution when Cholesterol uptake is fix 1 'b' value
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h
                    
                
              case 'sMtb2018_Griffin_Cholesterol'
                  
                    load(RIENKSMA_2018);
                    model = sMtb2018;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Cholesterol={'R_EX_ASN[e]','R_EX_CL[e]','R_EX_PI[e]','R_EX_H[e]','R_EX_K[e]','R_EX_NA[e]','R_EX_ZN[e]','R_EX_MG[e]', 'R_EX_CA[e]','R_EX_FE2[e]','R_EX_FE3[e]','R_EX_H2O[e]','R_EX_O2[e]','R_EX_NH3[e]','R_EX_SLF[e]','R_EX_CIT[e]','R_EX_CHOLESTEROL[e]', 'R_EX_ETH[e]'};  
                    Griffin_Cholesterol_data = [ 1, 1000, 1, 1000, 1000, 1000, 1000, 1000, 1000, 5, 5, 1000, 20, 10, 1000, 1, 1, 1];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol,-Griffin_Cholesterol_data,'l');
                    %model = changeRxnBounds(model,'EX_CHOLESTEROL[e]',-1,'b'); % Cannot consume cholesterol when LB is assigned.
                    model = changeRxnBounds(model,'R_Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
               
              case 'iEK1011_Griffin_Cholesterol'
                 
                    load(KAVVAS_MODEL);
                    model = iEK1011;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Cholesterol={'EX_h_e', 'EX_h2o_e', 'EX_o2_e', 'EX_asn__L_e', 'EX_nh4_e' ,...
                                                'EX_cit_e', 'EX_etoh_e', 'EX_ca2_e', 'EX_cl_e', 'EX_mg2_e',...
                                                'EX_so4_e', 'EX_fe3_e', 'EX_pi_e', 'EX_chsterol_e'};              
                    Griffin_Cholesterol_data = [ 1000, 1000, 20, 1, 10, 1, 1, 1000, 1000, 1000, 1000, 5, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol,-Griffin_Cholesterol_data,'l');
                   % model = changeRxnBounds(model,'EX_CHOLESTEROL[e]',-1,'b'); % Cannot consume cholesterol when LB is assigned.
                    model = changeRxnBounds(model,'ATPM',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.

                
               case 'sMtb_Griffin_Cholesterol'
        
                   load(RIENKSMA_2014);
                   model = sMtb;
                   model.lb(model.lb < 0) = -1000;
                   model.ub(model.ub > 1) = 1000;
                   exchangeRxns  = model.rxns(findExcRxns(model)); 
                   model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                   Medium_Griffin_Cholesterol={'EX_ASN[e]','EX_CL[e]','EX_PI[e]','EX_H[e]','EX_K[e]','EX_NA[e]','EX_ZN[e]','EX_MG[e]', 'EX_CA[e]','EX_FE2[e]','EX_FE3[e]','EX_H2O[e]','EX_O2[e]','EX_NH3[e]','EX_SLF[e]','EX_CIT[e]','EX_CHOLESTEROL[e]', 'EX_ETH[e]'};  
                   Griffin_Cholesterol_data = [ 1, 1000, 1, 1000, 1000, 1000, 1000, 1000, 1000, 5, 5, 1000, 20, 10, 1000, 1, 1, 1];
                   model = changeRxnBounds(model,Medium_Griffin_Cholesterol,-Griffin_Cholesterol_data,'l');
                   %model = changeRxnBounds(model,'EX_CHOLESTEROL[e]',-1,'b'); % Cannot consume cholesterol when LB is assigned.
                   model = changeRxnBounds(model,'Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
                   
               case 'iEK1011_2.0_Griffin_Cholesterol'
                 
                    load(KAVVAS_2_0);
                    model = iEK1011_2_0;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Cholesterol={'EX_h_e', 'EX_h2o_e', 'EX_o2_e', 'EX_asn__L_e', 'EX_nh4_e' ,...
                                                'EX_cit_e', 'EX_etoh_e', 'EX_ca2_e', 'EX_cl_e', 'EX_mg2_e',...
                                                'EX_so4_e', 'EX_fe3_e', 'EX_pi_e', 'EX_chsterol_e'};              
                    Griffin_Cholesterol_data = [ 1000, 1000, 20, 1, 10, 1, 1, 1000, 1000, 1000, 1000, 5, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol,-Griffin_Cholesterol_data,'l');
                   % model = changeRxnBounds(model,'EX_CHOLESTEROL[e]',-1,'b'); % Cannot consume cholesterol when LB is assigned.
                    model = changeRxnBounds(model,'ATPM',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.

                            
end