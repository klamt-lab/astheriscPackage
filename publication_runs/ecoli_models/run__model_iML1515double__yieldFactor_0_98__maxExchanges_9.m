generalPublicationRunSettingsDoubleModels;
minimalYieldFactor = 0.98;
numMaxExchanges = 9;
generaliML1515doubleModelPath;
reportPath = getPublicationRunReportPath(minimalYieldFactor, numMaxExchanges, modelPath);

astherisc(targetSubstrateReaction, targetMaxSubstrateUptake,...
minimalYieldFactor, minAbsoluteYield,...
minimalMdfAdvantage, minimalMdfWithCommunity,...
numMaxExchanges, temperature, vEpsilon,...
modelPath, reportPath, reportHeader,...
showCommunityModelBottlenecks, showSingleSpeciesBottlenecks, calculateIndirectBottlenecks,...
reactionsWithSetMinimumZero, ignoredMetabolites, fixedRatiosNameMat, exchangeReactionSelection,...
maximalMilpRunTime, verboseSettings, printValFiles, getDetailedOriginalSingleStrainSolution)
