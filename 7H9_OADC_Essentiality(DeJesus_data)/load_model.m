%% LOAD MODEL IN RICH MEDIUM: 
% 7H9 or 7H10 [7H9/10] plus oleic acid-albumin-dextrose-catalase [OADC]) with a variety of potential carbon sources available 
% (glycerol, dextrose, oleate, citrate, and glutamate)
% Author: Victor Lopez, Agosto 09, 2017.
function model = load_model(metabolic_model, NGAM)
% Load a genome scale model for the selected organism.
% Options:
% metabolic_model: 'iCG760_7H9_10_OADC',
%                  'iSM810_7H9_10_OADC', 'sMtb_7H9_10_OADC', 
%                   'iOSDD890_7H9_10_OADC',
%                  'GSMN_TB_1.1_7H9_10_OADC';
%                  'sMtb1.9_7H9_10_OADC'


%% MODEL PATH
    JAMSHIDI_MODEL = 'model/iNJ661.mat';
    FANG_MODEL = 'model/iNJ661v.mat';
    VASHISHT_MODEL = 'model/iOSDD890.mat';
    GARAY_MODEL = 'model/iCG760_grRuleCorrected.mat'; 
    MA_MODEL = 'model/iSM810_grRuleCorrected.mat';
    LOFTHOUSE_MODEL = 'model/GSMN_TB1.1_grRuleCorrected.mat';
    RIENKSMA_MODEL = 'model/sMtb_grRuleCorrected.mat';
    XAVIER_MODEL = 'model/iNJ661v_modified.mat';
    SMTB2_0_TICS_MODEL = 'model/sMtb2.0_tics_v13_1.mat';
    KAVVAS_MODEL = 'model/iEK1011_xlsx.mat';
    RIENKSMA_2018 = 'model/sMtb2018_EX.mat';
    KAVVAS_2_0 = 'model/iEK1011_2.0.mat';
    
% MEDIUM
    MA_MEDIUM = 'medium/ExchangeRxns_MA2015.mat';
    LOFTHOUSE_MEDIUM = 'medium/ExchangeRxns_MA2015.mat';   

%% CHOICES    
    
switch metabolic_model

	         case 'iEK1011_7H9_10_OADC'
                    load(KAVVAS_MODEL);
                    model = iEK1011;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_glu__L_e','EX_cu2_e','EX_btn_e',...
                        'EX_pydxn_e','EX_ca2_e', 'EX_mg2_e', 'EX_h_e', 'EX_k_e',...
                        'EX_nh4_e', 'EX_h2o_e','EX_pi_e', 'EX_cl_e', 'EX_o2_e',...
                        'EX_na1_e', 'EX_so4_e','EX_cit_e','EX_fe3_e', 'EX_glyc_e',...
                        'EX_glc__D_e','EX_ocdca_e'};
                    Medium_7H9_10_data = [ 1, 1000, 1, 1, 1000, 1000, 1000, 1000, 10,...
                                        1000, 1, 1000, 20, 1000, 1000, 1, 5, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-Medium_7H9_10_data,'l');
                    model = changeRxnBounds(model,'ATPM',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h
                    

             case 'sMtb2.0_TICS_7H9_10_OADC'
                    
                    load(SMTB2_0_TICS_MODEL);
                    model = sMtb2_0_tics_v13_1;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_PYRI[e]', 'EX_BIOTIN[e]', 'EX_GLU[e]','EX_CL[e]','EX_PI[e]','EX_H[e]','EX_K[e]','EX_NA[e]','EX_ZN[e]','EX_MG[e]', 'EX_CA[e]','EX_FE2[e]','EX_FE3[e]','EX_H2O[e]','EX_O2[e]','EX_NH3[e]','EX_SLF[e]','EX_CIT[e]','EX_GL[e]', 'EX_GLC[e]' ,'EX_9OCTADECENOATE[e]'};  
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-10,'l');
                    %model = changeRxnBounds(model,'EX_COII[e]',-0.00001,'l'); % cannot growth if there is not cobalt in the medium.
                    model = changeRxnBounds(model,'Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.                 


             case 'iCG760_7H9_10_OADC'
                    
                    load(GARAY_MODEL);
                    model = iCG760;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_GLU(e)','EX_H(e)','EX_PI(e)',...
                       'EX_O2(e)','EX_NH3(e)','EX_SLF(e)','EX_CIT(e)',...
                       'EX_FE3(e)', 'EX_GL(e)', 'EX_GLC(e)', 'EX_BIOTIN(e)',...
                       'EX_9_OCTADECENOATE(e)'};
                    Medium_7H910_data = [ 1, 1000, 1, 20, 10, 1000, 1, 5, 1, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-Medium_7H910_data,'l'); 
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % Mainteinence was fixt to 1 mmol/gDW/h

    
            case 'iSM810_7H9_10_OADC'
                               
                    load (MA_MODEL);
                    model = iSM810;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000; 
                    load (MA_MEDIUM);
                    model = changeRxnBounds(model,ExchangeRxns,0,'u');
                    Medium_7H9_10_OADC = {'R800', 'R804', 'R812','R863',...
                        'R838', 'R841', 'R851', 'R882', 'R918', 'R924',...
                        'R923', 'R925' ,'R864'};
                    Medium_7H910_data = [ 10, 20, 1, 1, 0, 1000, 1, 1, 1,...
                                          1, 0, 1, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC, Medium_7H910_data,'u');
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % Mainteinence was fixt to 1 mmol/gDW/h

                   
              case 'sMtb_7H9_10_OADC'
            
                   load(RIENKSMA_MODEL);
                    model = sMtb;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_PYRI[e]', 'EX_BIOTIN[e]',...
                        'EX_GLU[e]','EX_CL[e]','EX_PI[e]','EX_H[e]','EX_K[e]',...
                        'EX_NA[e]','EX_ZN[e]','EX_MG[e]', 'EX_CA[e]','EX_FE2[e]',...
                        'EX_FE3[e]','EX_H2O[e]','EX_O2[e]','EX_NH3[e]','EX_SLF[e]',...
                        'EX_CIT[e]','EX_GL[e]', 'EX_GLC[e]' ,'EX_9OCTADECENOATE[e]'};  
                    Medium_7H910_data = [ 1, 1, 1, 1000, 1, 1000, 1000,...
                        1000, 1000, 1000, 1000, 5, 5, 1000, 20, 10,...
                        1000, 1, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-Medium_7H910_data,'l'); 
                    model = changeRxnBounds(model,'Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.

                     
              case 'GSMN_TB_1.1_7H9_10_OADC'
                    
                    load(LOFTHOUSE_MODEL);
                    load(LOFTHOUSE_MEDIUM);
                    model = GSMN_TB1_1;
                    model.ub(model.ub > 1) = 1000;
                    model = changeRxnBounds(model,ExchangeRxns,0,'u');
                    Medium_7H9_10_OADC = {'R800', 'R804', 'R812','R863',...
                        'R838', 'R841', 'R851', 'R882', 'R918', 'R924',...
                        'R923', 'R925' ,'R864'};
                    Medium_7H910_data= [ 10, 20, 1, 1, 0, 1000, 1, 1, 1,...
                                         1, 0, 1, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC, Medium_7H910_data,'u');
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % Mainteinence was fixt to 1 mmol/gDW/h

                    
              case 'iOSDD890_7H9_10_OADC'
                    
                    load(VASHISHT_MODEL);
                    model = iOSDD890;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_glu_L(e)','EX_glc(e)' ,'EX_ca2(e)', 'EX_mg2(e)', 'EX_h(e)', 'EX_k(e)', 'EX_nh4(e)', 'EX_h2o(e)','EX_pi(e)', 'EX_cl(e)', 'EX_o2(e)', 'EX_na1(e)', 'EX_so4(e)','EX_cit(e)','EX_fe3(e)', 'EX_fe2(e)' , 'EX_glyc(e)'};
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-10,'l');
                    model = changeRxnBounds(model,'ATPS1',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h

              case 'iNJ661_7H9_10_OADC'
            
                    load (JAMSHIDI_MODEL);
                    model = iNJ661;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_glu_L(e)','EX_glc(e)',...
                        'EX_ca2(e)', 'EX_mg2(e)', 'EX_h(e)', 'EX_k(e)',...
                        'EX_nh4(e)', 'EX_h2o(e)','EX_pi(e)', 'EX_cl(e)',...
                        'EX_o2(e)', 'EX_na1(e)', 'EX_so4(e)','EX_cit(e)',...
                        'EX_fe3(e)', 'EX_fe2(e)' , 'EX_glyc(e)'};  
                    Medium_7H910_data = [ 1, 1, 1000, 1000, 1000, 1000,...
                        10, 1000, 1, 1000, 20, 1000, 1000, 1, 5, 5, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-Medium_7H910_data,'l'); 
                    model = changeRxnBounds(model,'ATPS1',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h

              case 'iNJ661v_7H9_10_OADC'
            
                    load(FANG_MODEL);
                    model = iNJ661v;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_glu_L(e)','EX_ocdcea(e)',' EX_glc(e)','EX_ca2(e)', 'EX_mg2(e)', 'EX_h(e)', 'EX_k(e)', 'EX_nh4(e)', 'EX_h2o(e)','EX_pi(e)', 'EX_cl(e)', 'EX_o2(e)', 'EX_na1(e)', 'EX_so4(e)','EX_cit(e)','EX_fe3(e)', 'EX_fe2(e)' , 'EX_glyc(e)'};  
                    Medium_7H910_data = [ 1, 1, 1000, 1000, 1000, 1000,...
                        10, 1000, 1, 1000, 20, 1000, 1000, 1, 5, 5, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-Medium_7H910_data,'l');  
                    ADD_iNJ661v = {'EX_xyl_D(e)'};
                    model = changeRxnBounds(model,ADD_iNJ661v,-0.01, 'l');
                    model = changeRxnBounds(model,'ATPS1',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h

              case 'iNJ661v_modified_7H9_10_OADC'
            
                    load(XAVIER_MODEL);
                    model = iNJ661v_modified;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_pydxn(e)','EX_btn(e)','EX_glu_L(e)','EX_ocdcea(e)','EX_glc(e)','EX_ca2(e)', 'EX_mg2(e)', 'EX_h(e)', 'EX_k(e)', 'EX_nh4(e)', 'EX_h2o(e)','EX_pi(e)', 'EX_cl(e)', 'EX_o2(e)', 'EX_na1(e)', 'EX_so4(e)','EX_cit(e)','EX_fe3(e)', 'EX_fe2(e)' , 'EX_glyc(e)'};  
                    Medium_7H910_data= [ 1, 1, 1, 1, 1, 1000, 1000, 1000, 1000, 10, 1000, 1, 1000, 20, 1000, 1000, 1, 5, 5, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-Medium_7H910_data,'l');  
                    ADD_iNJ661v = {'EX_xyl_D(e)'};
                    model = changeRxnBounds(model,ADD_iNJ661v,-0.01, 'l');
                    model = changeRxnBounds(model,'ATPS1',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h

              case 'sMtb2018_7H9_10_OADC'
            
                   load(RIENKSMA_2018);
                    model = sMtb2018;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                     Medium_7H9_10_OADC={'R_EX_PYRI[e]', 'R_EX_BIOTIN[e]',...
                        'R_EX_GLU[e]','R_EX_CL[e]','R_EX_PI[e]','R_EX_H[e]','R_EX_K[e]',...
                        'R_EX_NA[e]','R_EX_ZN[e]','R_EX_MG[e]', 'R_EX_CA[e]','R_EX_FE2[e]',...
                        'R_EX_FE3[e]','R_EX_H2O[e]','R_EX_O2[e]','R_EX_NH3[e]','R_EX_SLF[e]',...
                        'R_EX_CIT[e]','R_EX_GL[e]', 'R_EX_GLC[e]' ,'R_EX_9OCTADECENOATE[e]'};    
                    Medium_7H910_data = [ 1, 1, 1, 1000, 1, 1000, 1000,...
                        1000, 1000, 1000, 1000, 5, 5, 1000, 20, 10,...
                        1000, 1, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-Medium_7H910_data,'l'); 
                    model = changeRxnBounds(model,'R_Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.

              case 'iEK1011_2.0_7H9_10_OADC'
                    load(KAVVAS_2_0);
                    model = iEK1011_2_0;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_7H9_10_OADC={'EX_glu__L_e','EX_cu2_e','EX_btn_e',...
                        'EX_pydxn_e','EX_ca2_e', 'EX_mg2_e', 'EX_h_e', 'EX_k_e',...
                        'EX_nh4_e', 'EX_h2o_e','EX_pi_e', 'EX_cl_e', 'EX_o2_e',...
                        'EX_na1_e', 'EX_so4_e','EX_cit_e','EX_fe3_e', 'EX_glyc_e',...
                        'EX_glc__D_e','EX_ocdca_e'};
                    Medium_7H9_10_data = [ 1, 1000, 1, 1, 1000, 1000, 1000, 1000, 10,...
                                        1000, 1, 1000, 20, 1000, 1000, 1, 5, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_7H9_10_OADC,-Medium_7H9_10_data,'l');
                    model = changeRxnBounds(model,'ATPM',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h
                                     
end
              
              
end