% Global setting variables which are used all across
% the scripts with the publication's shown scenarios
targetSubstrateReaction = "R_EX_C_glc__D_exchg";
targetMaxSubstrateUptake = 1;
minAbsoluteYield = 1e-6;
minimalMdfAdvantage = .2;
minimalMdfWithCommunity = 0;
temperature = 298.15;
vEpsilon = 1e-6;
reportHeader = "ASTHERISC community MDF advantage search\n\n";
showCommunityModelBottlenecks = false;
showSingleSpeciesBottlenecks = true;
calculateIndirectBottlenecks = false;
exchangeReactionSelection = [];
maximalMilpRunTime = 1000;
verboseSettings = false;
printValFiles = false;
