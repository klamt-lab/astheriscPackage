generalPublicationRunSettingsDoubleModels;
minimalYieldFactor = 0.98;
numMaxExchanges = 9;
generalEcolicore2doubleModelPath;
reportPath = "./0Astherisc/publication_runs/ecoli_models/run_results/kdo8p_example.txt";
exchangeReactionSelection = ["R_EX_C_kdo8p_exchg"];
maximalMilpRunTime = 50;

% Create temporary model in which only the example
% reactions are allowed.
allowedReactionsNames = [
    "R_A5PISO_ecoli1"
    "R_ADK1_ecoli1"
    "R_ATPS4rpp_ecoli1"
    "R_CO2tpp_ecoli1"
    "R_CYTBO3_4pp_ecoli1"
    "R_EDA_ecoli1"
    "R_EDD_ecoli1"
    "R_ENO_ecoli1"
    "R_F6PA_ecoli1"
    "R_GAPD_ecoli1"
    "R_GLCt2pp_ecoli1"
    "R_GND_ecoli1"
    "R_H2Otpp_ecoli1"
    "R_HEX1_ecoli1"
    "R_KDOPS_ecoli1"
    "R_NADH16pp_ecoli1"
    "R_NADTRHD_ecoli1"
    "R_O2tpp_ecoli1"
    "R_PGI_ecoli1"
    "R_PGK_ecoli1"
    "R_PGM_ecoli1"
    "R_PIt2rpp_ecoli1"
    "R_PPS_ecoli1"
    "R_RPE_ecoli1"
    "R_RPI_ecoli1"
    "R_TALA_ecoli1"
    "R_TKT1_ecoli1"
    "R_TKT2_ecoli1"
    "R_EX_C_kdo8p_exchg"
    "R_EX_C_h_exchg"
    "R_EX_C_co2_exchg"
    "R_EX_C_glc__D_exchg"
    "R_EX_C_h2o_exchg"
    "R_EX_C_o2_exchg"
    "R_EX_C_pi_exchg"
    "R_EXCHG_ecoli1_h_p_to_h"
    "R_EXCHG_ecoli1_co2_p_to_co2"
    "R_EXCHG_ecoli1_kdo8p_c_to_kdo8p"
    "R_EXCHG_ecoli1_glc__D_p_to_glc__D"
    "R_EXCHG_ecoli1_h2o_p_to_h2o"
    "R_EXCHG_ecoli1_6pgc_c_to_6pgc"
    "R_EXCHG_ecoli1_pi_p_to_pi"
    "R_EXCHG_ecoli1_o2_p_to_o2"
    "R_EXCHG_ecoli1_dha_c_to_dha"
    "R_EXCHG_ecoli2_h_p_to_h"
    "R_EXCHG_ecoli2_h2o_p_to_h2o"
    "R_EXCHG_ecoli2_6pgc_c_to_6pgc"
    "R_EXCHG_ecoli2_pi_p_to_pi"
    "R_EXCHG_ecoli2_o2_p_to_o2"
    "R_EXCHG_ecoli2_dha_c_to_dha"
    "R_ATPS4rpp_ecoli2"
    "R_CYTBO3_4pp_ecoli2"
    "R_F6PA_ecoli2"
    "R_FBA_ecoli2"
    "R_G3PD2_ecoli2"
    "R_G3PD5_ecoli2"
    "R_G6PDH2r_ecoli2"
    "R_H2Otpp_ecoli2"
    "R_NADH16pp_ecoli2"
    "R_NADTRHD_ecoli2"
    "R_O2tpp_ecoli2"
    "R_PFK_ecoli2"
    "R_PGI_ecoli2"
    "R_PGL_ecoli2"
    "R_PIt2rpp_ecoli2"
    "R_TPI_ecoli2"
    "R_ATPS4rpp_ecoli1"
    "R_CYTBO3_4pp_ecoli1"
    "R_F6PA_ecoli1"
    "R_FBA_ecoli1"
    "R_G3PD2_ecoli1"
    "R_G3PD5_ecoli1"
    "R_G6PDH2r_ecoli1"
    "R_H2Otpp_ecoli1"
    "R_NADH16pp_ecoli1"
    "R_NADTRHD_ecoli1"
    "R_O2tpp_ecoli1"
    "R_PFK_ecoli1"
    "R_PGI_ecoli1"
    "R_PGL_ecoli1"
    "R_PIt2rpp_ecoli1"
    "R_TPI_ecoli1"
];
load("./0Astherisc/publication_runs/ecoli_models/models_and_dG0_data/ecolicore2double.mat", 'cnap')
allowedReactionsIndices = [];
for reactionName = allowedReactionsNames'
    index = PSBCNAFindReaction(reactionName, cnap);
    allowedReactionsIndices = [allowedReactionsIndices index];
end
for reactionIndex = 1:cnap.numr
    if ~ismember(reactionIndex, allowedReactionsIndices)
        cnap.reacMin(reactionIndex) = 0;
        cnap.reacMax(reactionIndex) = 0;
    end
end
modelPath = './0Astherisc/publication_runs/ecoli_models/models_and_dG0_data/ecolicore2double_kdo8p_only.mat';
save(modelPath, 'cnap')

astherisc(targetSubstrateReaction, targetMaxSubstrateUptake,...
minimalYieldFactor, minAbsoluteYield,...
minimalMdfAdvantage, minimalMdfWithCommunity,...
numMaxExchanges, temperature, vEpsilon,...
modelPath, reportPath, reportHeader,...
showCommunityModelBottlenecks, showSingleSpeciesBottlenecks, calculateIndirectBottlenecks,...
reactionsWithSetMinimumZero, ignoredMetabolites, fixedRatiosNameMat, exchangeReactionSelection,...
maximalMilpRunTime, verboseSettings, printValFiles, getDetailedOriginalSingleStrainSolution)

delete './0Astherisc/publication_runs/ecoli_models/models_and_dG0_data/ecolicore2double_kdo8p_only.mat'
