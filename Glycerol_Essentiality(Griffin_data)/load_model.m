function model = load_model(metabolic_model, NGAM)
% Load a genome scale model for the selected organism.
% Options:
% metabolic_model: 'iCG760_Griffin_glycerol',
%                  'iSM810_Griffin_Glycerol', 'sMtb_Griffin_Glycerol', 'iOSDD890_Griffin_Glycerol',
%                  'GSMN_TB_1.1_Griffin_Glycerol', 'sMtb_2.0_Griffin_Glycerol', 'iEK1011_2.0_Griffin_Glycerol';
%
% Author: Victor Lopez, April 12, 2017.

%% non-growth associated maintenance to 0.1 mmol gdw-1 h-1 as suggested by Rienksma et al., 2014
%NGAM = 1;

%% MODEL PATH

    VASHISHT_MODEL = 'model/iOSDD890.mat';
    GARAY_MODEL = 'model/iCG760_grRuleCorrected.mat'; 
    MA_MODEL = 'model/iSM810_grRuleCorrected.mat';
    LOFTHOUSE_MODEL = 'model/GSMN_TB1.1_grRuleCorrected.mat';
    XAVIER_MODEL = 'model/iNJ661v_modified.mat';
    RIENKSMA_2014 = 'model/sMtb_grRuleCorrected.mat';
    RIENKSMA_2018 = 'model/sMtb2018_EX.mat';
    KAVVAS_MODEL = 'model/iEK1011_xlsx.mat';
    KAVVAS_2_0 = 'model/iEK1011_2.0.mat';
    sMTB_2_0 = 'model/sMtb2.0.mat';
    
%% MEDIUM PATH

    MA_MEDIUM = 'medium/ExchangeRxns_MA2015.mat';
    LOFTHOUSE_MEDIUM = 'medium/ExchangeRxns_MA2015.mat';   

%% CHOICES    
    
switch metabolic_model
    
             case 'iCG760_Griffin_Glycerol'
                    
                    load(GARAY_MODEL);
                    model = iCG760;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % Mainteinence was fixt to 1 mmol/gDW/h
                    Medium_Griffin_Glycerol={'EX_ASN(e)','EX_PI(e)','EX_H(e)','EX_FE3(e)',...
                                              'EX_O2(e)','EX_NH3(e)','EX_SLF(e)','EX_CIT(e)','EX_GL(e)','EX_ETH(e)'};
                    Griffin_Glycerol_data = [ 1, 1, 1000, 5, 20, 10, 1000, 1, 1, 1];
                    model = changeRxnBounds(model,Medium_Griffin_Glycerol,-Griffin_Glycerol_data,'l');
                    
    
            case 'iSM810_Griffin_Glycerol'
                               
                    load (MA_MODEL);
                    model = iSM810;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000; 
                    load (MA_MEDIUM);
                    model = changeRxnBounds(model,ExchangeRxns,0,'u');
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % Maintainance was fix to 1 mmol/gDWH
                    Medium_Griffin_Glycerol = { 'R800', 'R804', 'R866', 'R822', 'R838', 'R841', 'R882', 'R924', 'R858' };                   
                    Griffin_Glycerol_data = [ 10, 20, 1, 1, 1, 1000, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Glycerol, Griffin_Glycerol_data,'u');
                    
                   
              case 'sMtb_Griffin_Glycerol'
            
                   load(RIENKSMA_2014);
                   model = sMtb;
                   model.lb(model.lb < 0) = -1000;
                   model.ub(model.ub > 1) = 1000;
                   exchangeRxns  = model.rxns(findExcRxns(model)); 
                   model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                   Medium_Griffin_Glycerol={'EX_ASN[e]','EX_CL[e]','EX_PI[e]','EX_H[e]','EX_K[e]',...
                                            'EX_NA[e]','EX_ZN[e]','EX_MG[e]', 'EX_CA[e]','EX_FE2[e]',...
                                            'EX_FE3[e]','EX_H2O[e]','EX_O2[e]','EX_NH3[e]','EX_SLF[e]',...
                                            'EX_CIT[e]','EX_GL[e]', 'EX_ETH[e]'};  
                   Griffin_Glycerol_data = [ 1, 1000, 1, 1000, 1000, 1000, 1000, 1000, 1000, 5, 5, 1000, 20, 10, 1000, 1, 1, 1];
                   model = changeRxnBounds(model,Medium_Griffin_Glycerol,-Griffin_Glycerol_data,'l');
                   model = changeRxnBounds(model,'Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
                     
                     
              case 'GSMN_TB_1.1_Griffin_Glycerol'
                    
                    load(LOFTHOUSE_MODEL);
                    load(LOFTHOUSE_MEDIUM);
                    model = GSMN_TB1_1;
                    model.ub(model.ub > 1) = 1000;
                    model = changeRxnBounds(model,ExchangeRxns,0,'u');
                    model = changeRxnBounds(model,'R129',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h
                    Medium_Griffin_Glycerol = { 'R800', 'R804', 'R866', 'R822', 'R838', 'R841', 'R882', 'R924', 'R858' };                   
                    Griffin_Glycerol_data = [ 10, 20, 1, 1, 1, 1000, 1, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Glycerol, Griffin_Glycerol_data,'u');
                    
              case 'iOSDD890_Griffin_Glycerol' %EX_etoh not in model.
                    
                    load(VASHISHT_MODEL);
                    model = iOSDD890;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Glycerol={'EX_h(e)', 'EX_h2o(e)', 'EX_o2(e)', 'EX_asn_L(e)', 'EX_nh4(e)' ,...
                                                'EX_cit(e)', 'EX_ca2(e)', 'EX_cl(e)', 'EX_mg2(e)',...
                                                'EX_so4(e)', 'EX_fe3(e)', 'EX_pi(e)', 'EX_glyc(e)'};              
                    Griffin_Glycerol_data = [ 1000, 1000, 20, 1, 10, 1, 1000, 1000, 1000, 1000, 5, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Glycerol,-Griffin_Glycerol_data,'l');
                    model = changeRxnBounds(model,'ATPS1',1,'b'); % maintaince reaction was fix to 1 mmol/gDW/h

               case 'iNJ661v_modified_Griffin_Glycerol' % Etoh not in model.
            
                    load(XAVIER_MODEL);
                    model = iNJ661v_modified;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Glycerol={'EX_h(e)', 'EX_h2o(e)', 'EX_o2(e)', 'EX_asn_L(e)', 'EX_nh4(e)' ,...
                                                'EX_cit(e)', 'EX_ca2(e)', 'EX_cl(e)', 'EX_mg2(e)',...
                                                'EX_so4(e)', 'EX_fe3(e)','EX_fe2(e)', 'EX_pi(e)', 'EX_glyc(e)',...
                                                'EX_k(e)','EX_na1(e)'};              
                    Griffin_Glycerol_data = [ 1000, 1000, 20, 1, 10, 1, 1000, 1000, 1000, 1000, 5, 5, 1, 1, 1000, 1000];
                    model = changeRxnBounds(model,Medium_Griffin_Glycerol,-Griffin_Glycerol_data,'l'); 
                    ADD_iNJ661v = {'EX_xyl_D(e)'};
                    model = changeRxnBounds(model,ADD_iNJ661v,-0.01, 'l');
                    model = changeRxnBounds(model,'ATPS1',NGAM,'b'); % maintaince reaction was fix to 1 mmol/gDW/h
                
                 case 'sMtb2018_Griffin_Glycerol'
                  
                    load(RIENKSMA_2018);
                    model = sMtb2018;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Glycerol={'R_EX_ASN[e]','R_EX_CL[e]','R_EX_PI[e]','R_EX_H[e]','R_EX_K[e]','R_EX_NA[e]','R_EX_ZN[e]',...
                                             'R_EX_MG[e]', 'R_EX_CA[e]','R_EX_FE2[e]','R_EX_FE3[e]','R_EX_H2O[e]','R_EX_O2[e]',...
                                             'R_EX_NH3[e]','R_EX_SLF[e]','R_EX_CIT[e]','R_EX_GL[e]', 'R_EX_ETH[e]'};  
                    Griffin_Glycerol_data = [ 1, 1000, 1, 1000, 1000, 1000, 1000, 1000, 1000, 5, 5, 1000, 20, 10, 1000, 1, 1, 1];
                    model = changeRxnBounds(model,Medium_Griffin_Glycerol,-Griffin_Glycerol_data,'l');
                    model = changeRxnBounds(model,'R_Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
               
                case 'iEK1011_Griffin_Glycerol'
                 
                    load(KAVVAS_MODEL);
                    model = iEK1011;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Cholesterol={'EX_h_e', 'EX_h2o_e', 'EX_o2_e', 'EX_asn__L_e', 'EX_nh4_e' ,...
                                                'EX_cit_e', 'EX_etoh_e', 'EX_ca2_e', 'EX_cl_e', 'EX_mg2_e',...
                                                'EX_so4_e', 'EX_fe3_e', 'EX_pi_e', 'EX_glyc_e'};              
                    Griffin_Cholesterol_data = [ 1000, 1000, 20, 1, 10, 1, 1, 1000, 1000, 1000, 1000, 5, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol,-Griffin_Cholesterol_data,'l');
                    model = changeRxnBounds(model,'ATPM',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.

                 case 'iEK1011_2.0_Griffin_Glycerol'
                 
                    load(KAVVAS_2_0);
                    model = iEK1011_2_0;
                    model.lb(model.lb < 0) = -1000;
                    model.ub(model.ub > 1) = 1000;
                    exchangeRxns  = model.rxns(findExcRxns(model)); 
                    model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                    Medium_Griffin_Cholesterol={'EX_h_e', 'EX_h2o_e', 'EX_o2_e', 'EX_asn__L_e', 'EX_nh4_e' ,...
                                                'EX_cit_e', 'EX_etoh_e', 'EX_ca2_e', 'EX_cl_e', 'EX_mg2_e',...
                                                'EX_so4_e', 'EX_fe3_e', 'EX_pi_e', 'EX_glyc_e'};              
                    Griffin_Cholesterol_data = [ 1000, 1000, 20, 1, 10, 1, 1, 1000, 1000, 1000, 1000, 5, 1, 1 ];
                    model = changeRxnBounds(model,Medium_Griffin_Cholesterol,-Griffin_Cholesterol_data,'l');
                    model = changeRxnBounds(model,'ATPM',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
    
                 case 'sMtb_2.0_Griffin_Glycerol'
            
                   load(sMTB_2_0);
                   model = sMtb2_0;
                   model.lb(model.lb < 0) = -1000;
                   model.ub(model.ub > 1) = 1000;
                   exchangeRxns  = model.rxns(findExcRxns(model)); 
                   model = changeRxnBounds(model,exchangeRxns,0,'l'); 
                   Medium_Griffin_Glycerol={'EX_ASN[e]','EX_CL[e]','EX_PI[e]','EX_H[e]','EX_K[e]',...
                                            'EX_NA[e]','EX_ZN[e]','EX_MG[e]', 'EX_CA[e]','EX_FE2[e]',...
                                            'EX_FE3[e]','EX_H2O[e]','EX_O2[e]','EX_NH3[e]','EX_SLF[e]',...
                                            'EX_CIT[e]','EX_GL[e]', 'EX_ETH[e]'};  
                   Griffin_Glycerol_data = [ 1, 1000, 1, 1000, 1000, 1000, 1000, 1000, 1000, 5, 5, 1000, 20, 10, 1000, 1, 1, 1];
                   model = changeRxnBounds(model,Medium_Griffin_Glycerol,-Griffin_Glycerol_data,'l');
                   model = changeRxnBounds(model,'Maintenance',NGAM,'b'); % Maintenance reaction was fix to 1 mmol/gDW/h.
                     
                    
                    
end
              
              
end