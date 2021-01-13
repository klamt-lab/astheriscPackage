generalPublicationRunSettings
% "Original" OptMDFPathway-paper based formulation of *fixed*
% metabolite ratio
% fixedRatios = [
%     PSBCNAFindMetabolite("M_atp_c_ecoli1", cnap) PSBCNAFindMetabolite('M_adp_c_ecoli1', cnap) 10;
%     PSBCNAFindMetabolite('M_adp_c_ecoli1', cnap) PSBCNAFindMetabolite('M_amp_c_ecoli1', cnap) 1;
%     PSBCNAFindMetabolite('M_nad_c_ecoli1', cnap) PSBCNAFindMetabolite('M_nadh_c_ecoli1', cnap) 10;
%     PSBCNAFindMetabolite('M_nadph_c_ecoli1', cnap) PSBCNAFindMetabolite('M_nadp_c_ecoli1', cnap) 10;
%     PSBCNAFindMetabolite('M_atp_c_ecoli2', cnap) PSBCNAFindMetabolite('M_adp_c_ecoli2', cnap) 10;
%     PSBCNAFindMetabolite('M_adp_c_ecoli2', cnap) PSBCNAFindMetabolite('M_amp_c_ecoli2', cnap) 1;
%     PSBCNAFindMetabolite('M_nad_c_ecoli2', cnap) PSBCNAFindMetabolite('M_nadh_c_ecoli2', cnap) 10;
%     PSBCNAFindMetabolite('M_nadph_c_ecoli2', cnap) PSBCNAFindMetabolite('M_nadp_c_ecoli2', cnap) 10;
% ];
% Derived formulation of metabolite ratio *ranges*
fixedRatiosNameMat = [
       "M_atp_c_ecoli1"   "M_adp_c_ecoli1" 3 10;
       "M_adp_c_ecoli1"   "M_amp_c_ecoli1" .5 2;
       "M_nad_c_ecoli1"   "M_nadh_c_ecoli1" 3 10;
       "M_nadph_c_ecoli1" "M_nadp_c_ecoli1" 3 10;
       "M_atp_c_ecoli2"   "M_adp_c_ecoli2" 3 10;
       "M_adp_c_ecoli2"   "M_amp_c_ecoli2" .5 2;
       "M_nad_c_ecoli2"   "M_nadh_c_ecoli2" 3 10;
       "M_nadph_c_ecoli2" "M_nadp_c_ecoli2" 3 10;
       "M_atp_c_ecoli3"   "M_adp_c_ecoli3" 3 10;
       "M_adp_c_ecoli3"   "M_amp_c_ecoli3" .5 2;
       "M_nad_c_ecoli3"   "M_nadh_c_ecoli3" 3 10;
       "M_nadph_c_ecoli3" "M_nadp_c_ecoli3" 3 10;
];
exchangeReactionSelection = [];

reactionsWithSetMinimumZero = ["R_ATPM_ecoli1", "R_ATPM_ecoli2", "R_ATPM_ecoli3"];
ignoredMetabolites = [
    "M_h2o_c_ecoli1", "M_h2o_c_ecoli2", "M_h2o_c_ecoli3",...
    "M_h2o_p_ecoli1", "M_h2o_p_ecoli2", "M_h2o_p_ecoli3",...
    "M_h2o_e_ecoli1", "M_h2o_e_ecoli2", "M_h2o_e_ecoli3", "M_h2o_exchg",...
    "M_h_c_ecoli1", "M_h_c_ecoli2", "M_h_c_ecoli3",...
    "M_h_p_ecoli1", "M_h_p_ecoli2", "M_h_p_ecoli3",...
    "M_h_e_ecoli1", "M_h_e_ecoli2", "M_h_e_ecoli3", "M_h_exchg"
];
