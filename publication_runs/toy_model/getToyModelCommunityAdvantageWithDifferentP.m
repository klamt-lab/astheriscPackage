defaultMinConcentration = 1;
defaultMaxConcentration = 10;
exchangeMinConcentration = 1;
exchangeMaxConcentration = 10;
basePath = "./0Astherisc/toy_model/";
baseName = "toymodelDouble";
sbmlFilePath = basePath + baseName + ".xml";
dG0FilePath = basePath + "dG0_" + baseName + ".json";
savePath = basePath + baseName + ".mat";

sbmlFilePath = convertStringsToChars(sbmlFilePath);
dG0FilePath = convertStringsToChars(dG0FilePath);
savePath = convertStringsToChars(savePath);


% Global setting variables which are used all across
% the scripts with the publication's shown scenarios
targetSubstrateReaction = "R_EX_C_A_exchg";
targetMaxSubstrateUptake = 1;
minAbsoluteYield = 1e-6;
minimalMdfAdvantage = 0;
minimalMdfWithCommunity = -inf;
temperature = 298.15;
vEpsilon = 1e-6;
reportHeader = "ASTHERISC community MDF advantage search\nWith maximal MILP time 1000 s\n";
showCommunityModelBottlenecks = false;
showSingleSpeciesBottlenecks = false;
calculateIndirectBottlenecks = false;
exchangeReactionSelection = ["R_EX_C_P_exchg"];
modelPath = savePath;
ignoredMetabolites = [];
fixedRatiosNameMat = [
];
minimalYieldFactor = 0.1;
numMaxExchanges = inf;
reportPath = "";
reactionsWithSetMinimumZero = [];
maximalMilpRunTime = 10;
verboseSettings = false;
printValFiles = false;

minConcentration = .1;
maxConcentration = 25;
concentrationStep = .1;

currentMinConcentration = minConcentration;
distances = [];
output = "";
while currentMinConcentration <= maxConcentration
    if currentMinConcentration <= 10
        nondefaultConcentrationsMat = [
            "M_P_c_species1" currentMinConcentration 10;
            "M_P_c_species2" currentMinConcentration 10;
        ];
    else
        nondefaultConcentrationsMat = [
            "M_P_c_species1" currentMinConcentration currentMinConcentration;
            "M_P_c_species2" currentMinConcentration currentMinConcentration;
        ];
    end
    loadAndSaveCommodepyCommunityModel(sbmlFilePath, dG0FilePath, savePath,...
        defaultMinConcentration, defaultMaxConcentration, nondefaultConcentrationsMat,...
        exchangeMinConcentration, exchangeMaxConcentration);

    [finalReport, mdfsWithCommunity, mdfsWithoutCommunity] = astheriscTest(targetSubstrateReaction, targetMaxSubstrateUptake,...
    minimalYieldFactor, minAbsoluteYield,...
    minimalMdfAdvantage, minimalMdfWithCommunity,...
    numMaxExchanges, temperature, vEpsilon,...
    modelPath, reportPath, reportHeader,...
    showCommunityModelBottlenecks, showSingleSpeciesBottlenecks, calculateIndirectBottlenecks,...
    reactionsWithSetMinimumZero, ignoredMetabolites, fixedRatiosNameMat, exchangeReactionSelection,...
    maximalMilpRunTime, verboseSettings, printValFiles);

    
    if ~isempty(mdfsWithCommunity)
        distance = abs(mdfsWithCommunity - mdfsWithoutCommunity);
        distances = [distances distance];
        output = output + currentMinConcentration + ";" + mdfsWithCommunity + ";" + mdfsWithoutCommunity + "\n";
    end
    
    currentMinConcentration = currentMinConcentration + concentrationStep;
end
fprintf(output)
