generalPublicationRunSettingsDoubleModels;
minimalYieldFactor = 0.98;
numMaxExchanges = inf;
generaliML1515doubleModelPath;
reportPath = getPublicationRunReportPath(minimalYieldFactor, numMaxExchanges, modelPath);
%reportPath = "./0Astherisc/UptakeTestReport.txt";
%exchangeReactionSelection = ["R_EX_C_r1p_exchg"]; % Used in tests
%targetMaxSubstrateUptake = 1;
%maximalMilpRunTime = 50;

astherisc(targetSubstrateReaction, targetMaxSubstrateUptake,...
minimalYieldFactor, minAbsoluteYield,...
minimalMdfAdvantage, minimalMdfWithCommunity,...
numMaxExchanges, temperature, vEpsilon,...
modelPath, reportPath, reportHeader,...
showCommunityModelBottlenecks, showSingleSpeciesBottlenecks, calculateIndirectBottlenecks,...
reactionsWithSetMinimumZero, ignoredMetabolites, fixedRatiosNameMat, exchangeReactionSelection,...
maximalMilpRunTime, verboseSettings, printValFiles, getDetailedOriginalSingleStrainSolution)
