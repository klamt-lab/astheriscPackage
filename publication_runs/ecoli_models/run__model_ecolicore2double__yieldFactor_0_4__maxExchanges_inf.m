generalPublicationRunSettingsDoubleModels;
minimalYieldFactor = 0.4;
numMaxExchanges = inf;
generalEcolicore2doubleModelPath;
reportPath = getPublicationRunReportPath(minimalYieldFactor, numMaxExchanges, modelPath);
% exchangeReactionSelection = ["R_EX_C_2pg_exchg"];

astherisc(targetSubstrateReaction, targetMaxSubstrateUptake,...
minimalYieldFactor, minAbsoluteYield,...
minimalMdfAdvantage, minimalMdfWithCommunity,...
numMaxExchanges, temperature, vEpsilon,...
modelPath, reportPath, reportHeader,...
showCommunityModelBottlenecks, showSingleSpeciesBottlenecks, calculateIndirectBottlenecks,...
reactionsWithSetMinimumZero, ignoredMetabolites, fixedRatiosNameMat, exchangeReactionSelection,...
maximalMilpRunTime, verboseSettings, printValFiles)
