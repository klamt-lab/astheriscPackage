function [] = loadAndSavePublicationCommunityModelsAsMat()

%% General settings
% Relative path (seen from CNA's main folder) to the SBML and dG0 JSON
% folder
basePath = "./0Astherisc/publication_runs/ecoli_models/models_and_dG0_data/";
%% Add metabolite concentratons
defaultMinConcentration = 1e-6;
defaultMaxConcentration = 0.02;
exchangeMinConcentration = 1e-6;
exchangeMaxConcentration = 10;
nondefaultConcentrationsMat = [
    "M_co2_c_ecoli1" 1e-6 1e-3;
    "M_co2_c_ecoli2" 1e-6 1e-3;
    "M_co2_p_ecoli1" 1e-6 1e-3;
    "M_co2_p_ecoli2" 1e-6 1e-3;
    "M_h2o_c_ecoli1" 1 1;
    "M_h2o_c_ecoli2" 1 1;
    "M_h2o_p_ecoli1" 1 1;
    "M_h2o_p_ecoli2" 1 1;
    "M_h2o_e_ecoli1" 1 1;
    "M_h2o_e_ecoli2" 1 1;
    "M_h2o_exchg" 1 1;
    "M_h_c_ecoli1" 1 1;
    "M_h_c_ecoli2" 1 1;
    "M_h_p_ecoli1" 1 1;
    "M_h_p_ecoli2" 1 1;
    "M_h_e_ecoli1" 1 1;
    "M_h_e_ecoli2" 1 1;
    "M_h_exchg" 1 1;
];

%% ecolicore2double
% Load ecolicore2 community model with exchanges for all periplasmic and
% cytosolic metabolites
baseName = "ecolicore2double";
sbmlFilePath = basePath + baseName + "_model.xml";
dG0FilePath = basePath + baseName + "_dG0.json";
savePath = basePath + baseName + ".mat";

sbmlFilePath = convertStringsToChars(sbmlFilePath);
dG0FilePath = convertStringsToChars(dG0FilePath);
savePath = convertStringsToChars(savePath);

loadAndSaveCommModelPyCommunityModel(sbmlFilePath, dG0FilePath, savePath,...
    defaultMinConcentration, defaultMaxConcentration, nondefaultConcentrationsMat,...
    exchangeMinConcentration, exchangeMaxConcentration);

%% iML1515double
% Load ecolicore2 community model with exchanges for all periplasmic and
% cytosolic metabolites
baseName = "iML1515double";
sbmlFilePath = basePath + baseName + "_model.xml";
dG0FilePath = basePath + baseName + "_dG0.json";
savePath = basePath + baseName + ".mat";

sbmlFilePath = convertStringsToChars(sbmlFilePath);
dG0FilePath = convertStringsToChars(dG0FilePath);
savePath = convertStringsToChars(savePath);

loadAndSaveCommModelPyCommunityModel(sbmlFilePath, dG0FilePath, savePath,...
    defaultMinConcentration, defaultMaxConcentration, nondefaultConcentrationsMat,...
    exchangeMinConcentration, exchangeMaxConcentration);

%% ecolicore2triple
nondefaultConcentrationsMat = [
    "M_co2_c_ecoli1" 1e-6 1e-3;
    "M_co2_c_ecoli2" 1e-6 1e-3;
    "M_co2_c_ecoli3" 1e-6 1e-3;
    "M_co2_p_ecoli1" 1e-6 1e-3;
    "M_co2_p_ecoli2" 1e-6 1e-3;
    "M_co2_p_ecoli3" 1e-6 1e-3;
    "M_h2o_c_ecoli1" 1 1;
    "M_h2o_c_ecoli2" 1 1;
    "M_h2o_c_ecoli3" 1 1;
    "M_h2o_p_ecoli1" 1 1;
    "M_h2o_p_ecoli2" 1 1;
    "M_h2o_p_ecoli3" 1 1;
    "M_h2o_e_ecoli1" 1 1;
    "M_h2o_e_ecoli2" 1 1;
    "M_h2o_e_ecoli3" 1 1;
    "M_h2o_exchg" 1 1;
    "M_h_c_ecoli1" 1 1;
    "M_h_c_ecoli2" 1 1;
    "M_h_c_ecoli3" 1 1;
    "M_h_p_ecoli1" 1 1;
    "M_h_p_ecoli2" 1 1;
    "M_h_p_ecoli3" 1 1;
    "M_h_e_ecoli1" 1 1;
    "M_h_e_ecoli2" 1 1;
    "M_h_e_ecoli3" 1 1;
    "M_h_exchg" 1 1;
];

% Load ecolicore2 community model with exchanges for all periplasmic and
% cytosolic metabolites
baseName = "ecolicore2triple";
sbmlFilePath = basePath + baseName + "_model.xml";
dG0FilePath = basePath + baseName + "_dG0.json";
savePath = basePath + baseName + ".mat";

sbmlFilePath = convertStringsToChars(sbmlFilePath);
dG0FilePath = convertStringsToChars(dG0FilePath);
savePath = convertStringsToChars(savePath);

loadAndSaveCommModelPyCommunityModel(sbmlFilePath, dG0FilePath, savePath,...
    defaultMinConcentration, defaultMaxConcentration, nondefaultConcentrationsMat,...
    exchangeMinConcentration, exchangeMaxConcentration);
end
