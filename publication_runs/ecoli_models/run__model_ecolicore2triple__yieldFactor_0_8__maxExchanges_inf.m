generalPublicationRunSettingsTripleModels;
minimalYieldFactor = 0.8;
numMaxExchanges = inf;
generalEcolicore2tripleModelPath;
reportPath = getPublicationRunReportPath(minimalYieldFactor, numMaxExchanges, modelPath);

astherisc(targetSubstrateReaction, targetMaxSubstrateUptake,...
minimalYieldFactor, minAbsoluteYield,...
minimalMdfAdvantage, minimalMdfWithCommunity,...
numMaxExchanges, temperature, vEpsilon,...
modelPath, reportPath, reportHeader,...
showCommunityModelBottlenecks, showSingleSpeciesBottlenecks, calculateIndirectBottlenecks,...
reactionsWithSetMinimumZero, ignoredMetabolites, fixedRatiosNameMat, exchangeReactionSelection,...
maximalMilpRunTime, verboseSettings, printValFiles, getDetailedOriginalSingleStrainSolution)
