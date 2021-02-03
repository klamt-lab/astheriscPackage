defaultMinConcentration = 1;
defaultMaxConcentration = 10;
exchangeMinConcentration = 1;
exchangeMaxConcentration = 10;
nondefaultConcentrationsMat = [
    "M_P_exchg" 5 10;
];
basePath = "./0Astherisc/publication_runs/toy_model/";
baseName = "toymodelDouble";
sbmlFilePath = basePath + baseName + ".xml";
dG0FilePath = basePath + "dG0_" + baseName + ".json";
savePath = basePath + baseName + ".mat";

sbmlFilePath = convertStringsToChars(sbmlFilePath);
dG0FilePath = convertStringsToChars(dG0FilePath);
savePath = convertStringsToChars(savePath);

loadAndSaveCommModelPyCommunityModel(sbmlFilePath, dG0FilePath, savePath,...
    defaultMinConcentration, defaultMaxConcentration, nondefaultConcentrationsMat,...
    exchangeMinConcentration, exchangeMaxConcentration);


% Global setting variables which are used all across
% the scripts with the publication's shown scenarios
targetSubstrateReaction = "R_EX_C_S_exchg";
targetMaxSubstrateUptake = 1;
minAbsoluteYield = 1e-6;
minimalMdfAdvantage = .1;
minimalMdfWithCommunity = -inf;
temperature = 298.15;
vEpsilon = 1e-6;
reportHeader = "ASTHERISC community MDF advantage search\nWith maximal MILP time 1000 s\n";
showCommunityModelBottlenecks = false;
showSingleSpeciesBottlenecks = true;
calculateIndirectBottlenecks = false;
exchangeReactionSelection = ["R_EX_C_P_exchg"];
modelPath = savePath;
ignoredMetabolites = [];
fixedRatiosNameMat = [];
minimalYieldFactor = 0.1;
numMaxExchanges = inf;
reportPath = "./0Astherisc/publication_runs/toy_model/toymodelDouble_report.txt";
reactionsWithSetMinimumZero = [];
maximalMilpRunTime = 10;
verboseSettings = true;
printValFiles = false;
getDetailedOriginalSingleStrainSolution = false;

finalReport = astherisc(targetSubstrateReaction, targetMaxSubstrateUptake,...
minimalYieldFactor, minAbsoluteYield,...
minimalMdfAdvantage, minimalMdfWithCommunity,...
numMaxExchanges, temperature, vEpsilon,...
modelPath, reportPath, reportHeader,...
showCommunityModelBottlenecks, showSingleSpeciesBottlenecks, calculateIndirectBottlenecks,...
reactionsWithSetMinimumZero, ignoredMetabolites, fixedRatiosNameMat, exchangeReactionSelection,...
maximalMilpRunTime, verboseSettings, printValFiles, getDetailedOriginalSingleStrainSolution);
