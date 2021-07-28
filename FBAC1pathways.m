%%
initCobraToolbox(false);
changeCobraSolver('gurobi','all');

%% Create a base model E. coli

model=readCbModel('iML1515.mat');

% Allow passive formate import
model=changeRxnBounds(model,'FORtppi', -1000, 'l');
    
% No glucose
model=changeRxnBounds(model,'EX_glc__D_e', 0, 'l');   

% Add C1 interconversion reactions
model=addReaction(model, 'FDH',{'co2_c','nadh_c','h_c','for_c','nad_c'},[-1 -1 -1 1 1]);
model=addReaction(model, 'formylKinase', {'for_c','atp_c','forp_c','adp_c'},[-1 -1 1 1]);
model=addReaction(model, 'formylTransferase', {'forp_c', 'coa_c', 'forcoa_c', 'pi_c'}, [-1 -1 1 1]);
model=addReaction(model, 'acylAldRed', {'fald_c','nad_c','coa_c','forcoa_c','nadh_c','h_c'},[-1 -1 -1 1 1 1]);
model=addReaction(model,'MeOHDH',{'MeOH_c','nad_c','fald_c','nadh_c','h_c'},[-1 -1 1 1 1],1,-1000,1000);
model=addReaction(model,'MMO',{'methane_c', 'nadh_c', 'o2_c', 'h_c', 'MeOH_c', 'nad_c', 'h2o_c'},[-1 -1 -1 -1 1 1 1],1,-1000,1000);
model=addReaction(model, 'hydrogenase',{'nad_c','h2_c','nadh_c','h_c'},[-1 -1 1 1]);


%% Implement pathways

% Alpha reduction to glycolate
ARPGly=addReaction(model, 'HACL', {'forcoa_c','fald_c','glyclcoa_c'},[-1 -1 1]);
ARPGly=addReaction(ARPGly,'glycltcoaTes',{'glyclcoa_c','h2o_c','coa_c','glyclt_c'},[-1 -1 1 1]);
ARPGly=changeRxnBounds(ARPGly, 'HACL',0,'l');

% Alpha reduction to acetyl-CoA
ARPACoA=addReaction(model, 'HACL', {'forcoa_c','fald_c','glyclcoa_c'},[-1 -1 1]);
ARPACoA=addReaction(ARPACoA,'glycltcoaTes',{'glyclcoa_c','h2o_c','coa_c','glyclt_c'},[-1 -1 1 1]);
ARPACoA=addReaction(ARPACoA, 'ACR', {'glyclcoa_c','nadh_c','h_c','gcald_c','nad_c','coa_c'},[-1 -1 -1 1 1 1]);
ARPACoA=addReaction(ARPACoA, 'DOR', {'gcald_c','nadh_c','h_c','ethgly_c','nad_c'},[-1 -1 -1 1 1]);
ARPACoA=addReaction(ARPACoA, 'DDR', {'ethgly_c','acald_c','h2o_c'},[-1 1 1],1,0,1000);
ARPACoA=changeRxnBounds(ARPACoA, 'HACL',0,'l');

% Alpha reduction to glyceraldehyde
ARPGlycer=addReaction(model, 'HACL', {'forcoa_c','fald_c','glyclcoa_c'},[-1 -1 1]);
ARPGlycer=addReaction(ARPGlycer, 'ACR', {'glyclcoa_c','nadh_c','h_c','gcald_c','nad_c','coa_c'},[-1 -1 -1 1 1 1]);
ARPGlycer=addReaction(ARPGlycer, 'HACLC3', {'forcoa_c','gcald_c','glycercoa_c'},[-1 -1 1]);
ARPGlycer=addReaction(ARPGlycer,'glycercoaTes',{'glycercoa_c','h2o_c','coa_c','glyc__R_c'},[-1 -1 1 1]);
ARPGlycer=addReaction(ARPGlycer,'glycercoaRed',{'glycercoa_c','nadh_c','h_c','glyald_c','nad_c','coa_c'},[-1 -1 -1 1 1 1]);
ARPGlycer=changeRxnBounds(ARPGlycer, 'HACL',0,'l');

% RuMP
RuMP=addReaction(model, 'HPS', {'ru5p__D_c','fald_c','h6p_c'},[-1 -1 1]);
RuMP=addReaction(RuMP, 'PHI', {'h6p_c','f6p_c'},[-1 1]);

% Serine
Serine = addReaction(model, 'THFlig', {'fald_c','thf_c','mlthf_c','h2o_c'},[-1 -1 1 1]);
Serine = addReaction(Serine, 'SGA', {'ser__L_c','glx_c','hpyr_c','gly_c'},[-1 -1 1 1]);
Serine = addReaction(Serine, 'MTK', {'mal__L_c','atp_c','coa_c','malylcoa_c','adp_c','pi_c'},[-1 -1 -1 1 1 1],1,0,1000);
Serine = addReaction(Serine, 'MCL', {'malylcoa_c','glx_c','accoa_c'},[-1 1 1]);

% Formolase
FLS = addReaction(model, 'FLS', {'fald_c','dha_c'},[-3 1]);

% SACA
SACA = addReaction(model, 'GALS', {'fald_c','gcald_c'},[-2 1]);
SACA = addReaction(SACA, 'ACPS', {'gcald_c','pi_c','actp_c','h2o_c'},[-1 -1 1 1]);

% Reductive Glycine
RGLS = changeRxnBounds(model, 'GLYCL', -1000,'l');

models = [
    ARPGly
    ARPACoA
    ARPGlycer
    RuMP
    FLS
    SACA
    Serine
    RGLS
    ];

%% FBA
biomassEYield = zeros(length(models),4);
biomassCYield = zeros(length(models),4);
mu = zeros(length(models),4);
osenseStr = 'max';

% Solve the models for different C1 substrates
for i = 1:length(models)
    modelTemp = models(i);
    % Formate
    modelTempFA(i)=changeRxnBounds(modelTemp, 'FDH', 0,'u'); % Do not allow reuse of oxidized carbon when designing pathways for reduced C1
    modelTempFA(i) = changeRxnBounds(modelTempFA(i),'EX_for_e',-10,'l');
    modelTempFA(i) = changeRxnBounds(modelTempFA(i),'EX_for_e',0,'u');
    minNorm = 1e-1;
    solved = false;
    while ~solved
        try
            solTempFA(i)=optimizeCbModel(modelTempFA(i), osenseStr, minNorm);
            if solTempFA(i).stat == 1
                solved = true;
            else
                minNorm = minNorm/10;
            end
        catch
            minNorm = minNorm/10;
        end 
    end
    biomassEYield(i,1) = solTempFA(i).f/10*1000;
    biomassCYield(i,1) = solTempFA(i).f/10*1000;
    mu(i,1) = solTempFA(i).f;
    substrateRxn = {'EX_for_e'};
    solTempFA(i).x = fixFlux(modelTempFA(i),solTempFA(i),substrateRxn);
    
    % Formaldehyde
    modelTempFALD(i)=changeRxnBounds(modelTemp, 'FDH', 0,'u'); % Do not allow reuse of oxidized carbon when designing pathways for reduced C1
    modelTempFALD(i)=changeRxnBounds(modelTempFALD(i),{'formylKinase','formylTransferase'},0,'b');
    modelTempFALD(i)=changeRxnBounds(modelTempFALD(i),'EX_fald_e',-10,'l');
    modelTempFALD(i)=changeRxnBounds(modelTempFALD(i),'EX_fald_e',0,'u');
    minNorm = 1e-1;
    solved = false;
    while ~solved
        try
            solTempFALD(i)=optimizeCbModel(modelTempFALD(i), osenseStr, minNorm);
            if solTempFALD(i).stat == 1
                solved = true;
            else
                minNorm = minNorm/10;
            end
        catch
            minNorm = minNorm/10;
        end 
    end
    biomassEYield(i,2) = solTempFALD(i).f/20*1000;
    biomassCYield(i,2) = solTempFALD(i).f/10*1000;
    mu(i,2) = solTempFALD(i).f;
    substrateRxn = {'EX_fald_e'};
    solTempFALD(i).x = fixFlux(modelTempFALD(i),solTempFALD(i),substrateRxn);
    
    % Methanol
    modelTempMeOH(i)=changeRxnBounds(modelTemp, 'FDH', 0,'u');
    modelTempMeOH(i)=changeRxnBounds(modelTempMeOH(i),{'formylKinase','formylTransferase'},0,'b');
    modelTempMeOH(i)=addReaction(modelTempMeOH(i),'EX_MeOH_e', {'MeOH_e'},[-1],1,-10,0);
    modelTempMeOH(i)=addReaction(modelTempMeOH(i),'MeOHIn', {'MeOH_e', 'MeOH_c'},[-1 1],1,-1000,1000);
    minNorm = 1e-1;
    solved = false;
    while ~solved
        try
            solTempMeOH(i)=optimizeCbModel(modelTempMeOH(i), osenseStr, minNorm);
            if solTempMeOH(i).stat == 1
                solved = true;
            else
                minNorm = minNorm/10;
            end
        catch
            minNorm = minNorm/10;
        end 
    end
    biomassEYield(i,3) = solTempMeOH(i).f/30*1000;
    biomassCYield(i,3) = solTempMeOH(i).f/10*1000;
    mu(i,3) = solTempMeOH(i).f;
    substrateRxn = {'EX_MeOH_e'};
    solTempMeOH(i).x = fixFlux(modelTempMeOH(i),solTempMeOH(i),substrateRxn);    
    
    % Methane
    modelTempCH4(i)=changeRxnBounds(modelTemp, 'FDH', 0,'u');
    modelTempCH4(i)=changeRxnBounds(modelTempCH4(i),{'formylKinase','formylTransferase'},0,'b');    
    modelTempCH4(i)=addReaction(modelTempCH4(i),'EX_methane_e', {'methane_e'},[-1],1,-10,0);
    modelTempCH4(i)=addReaction(modelTempCH4(i),'methaneIn', {'methane_e', 'methane_c'},[-1 1],1,-1000,1000);
    minNorm = 1e-1;
    solved = false;
    while ~solved
        try
            solTempCH4(i)=optimizeCbModel(modelTempCH4(i), osenseStr, minNorm);
            if solTempCH4(i).stat == 1
                solved = true;
            else
                minNorm = minNorm/10;
            end
        catch
            minNorm = minNorm/10;
        end 
    end
    biomassEYield(i,4) = solTempCH4(i).f/40*1000;
    biomassCYield(i,4) = solTempCH4(i).f/10*1000;
    mu(i,4) = solTempCH4(i).f;
    substrateRxn = {'EX_methane_e'};
    solTempCH4(i).x = fixFlux(modelTempCH4(i),solTempCH4(i),substrateRxn);  
end
