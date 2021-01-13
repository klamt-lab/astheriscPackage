function [mdf, v, conc, dfs]= CNAcomputeOptMDFpathway_ratio_range(maximalMilpRunTime, cnap, RT, G0, Cmin, Cmax, D, d, fixed_ratios)
%
% This function calculates a pathway (flux vector v) with the maximal
% max-min driving force (MDF),i.e., a pathway with associated metabolite
% concentrations where the minimum driving force of all participating
% reactions is maximal.
% (for MDF see Noor et al. (2014), PLOS Comp Biol, 10:e1003483)
%
% Usage: [mdf, v, conc, dfs]= CNAcomputeOptMDFpathway(cnap, RT, G0, Cmin, Cmax, D, d, fixed_ratios, stoichMat)
%
% The calculated results fulfill
%
%       cnap.stoichMat * v = 0
%       cnap.reacMin <= v <= cnap.reacMax
%       D * v <= d
%       Cmin <= exp(conc) <= Cmax
%
% and in addition the driving forces (dfs) of all reactions that
% participate in v is >= mdf (with dfs equal to mdf for at least one
% participating reaction).
% The thermodynamic parameters for the calculation of the driving forces
% are RT (product of the gas constant R and the temperature) and G0 (Gibbs
% energy of the reactions). Cmin and Cmax are lower/upper bounds for the
% metabolite concentrations in the network.
%
% In order to run this function it is necessary that both the MATLAB CPLEX and Java CPLEX
% interfaces work. If CPLEX is installed under /cluster/apps/cplex-124 the commands for
% this are (see also startup.m file in CNA's home directory):
%      addpath('/cluster/apps/cplex-124/cplex/matlab/');
%      javaaddpath('/cluster/apps/cplex-124/cplex/lib/cplex.jar');
% Additionally, the MATLAB JVM needs to have the CPLEX shared library on its library path
% which must be set up before (!!) starting MATLAB. For MATLAB versions up to 2010b
% this can be achieved by adding
%       /cluster/apps/cplex-124/cplex/bin/x86-64_sles10_4.1
% to Matlab's librarypath.txt (or javalibrarypath.txt from version 2012 onwards) configuration file
%(see also manual).
%
% Inputs (the first 5 arguments are mandatory):
% ---------------------------------------------
%
%   cnap: (mandatory) is a CellNetAnalyzer (mass-flow) project variable representing
%         a (metabolic) reaction network. You may easily generate such a structure
%         by the CNAgenerateMFNetwork function.
%
%     The function accesses the following fields of cnap (see manual):
%
%     cnap.stoichMat: the stoichiometric matrix of the network
%     cnap.specInternal: vector with the indices of the internal species
%     cnap.reacMin: lower boundaries of reaction rates; -Inf will be replaced by -1000
%     cnap.reacMax: upper boundaries of reaction rates; Inf will be replaced by 1000
%     cnap.epsilon: can be used to specify the tolerance used for solving
%       the MILP (if set to a smaller value than 1E-9 the tolerance will
%       still be 1E-9 since this is the minimal tolerance for the solver)
%
%   RT: product of the gas constant R and the temperature
%
%   G0: column vector of length cnap.numr containing the Gibbs energy of
%       each reaction; if NaN is given as value for a reaction then this
%       reaction may participate in the pathway v but no driving force
%       constraint constraint for this reactions needs to be fulfilled;
%       if G0 is a cnap.numr x 2 matrix then the second column specifies
%       the uncetainty of the Gibbs energy which means that the Gibbs
%       energy G0(i, 1) of reaction i may vary in the interval
%       [G0(i, 1) - G0(i, 2) ... G0(i, 1) + G0(i, 2)]
%
%   Cmin: lower bounds for the metabolite concentrations (in [M] = mol/L)
%
%   Cmax: upper bounds for the metabolite concentrations (in [M] = mol/L)
%
%   D: the matrix specifying (with vector d) the desired flux vectors as given above.
%      D has Dimension numDesiredConst x cnap.numr. D and d can be empty
%
%   d: the vector specifying (with matrix D) the desired flux vectors as given above.
%      d has Dimension numDesiredConst x 1. D and d can be empty
%
%  fixed_ratios (optional): matrix with three columns; can be used to fix
%      the ratios of metabolites to fixed values; each such constraint is a
%      row in fixed_ratios with the i-th row setting the constraint
%      metabolite(fixed_ratios(i, 1))/metabolite(fixed_ratios(i, 2)) = fixed_ratios(i, 3)
%
%  stoichMat (optional): if given is used instead of cnap.stoichMat to set
%      up the driving force constraints
%
%
% Results:
% --------
%
%  mdf: max-min driving force of the calculated pathway (= min(dfs)).
%
%  v: flux vector describing the pathway
%
%  conc: natural logarithm of the metabolite concentrations
%
%  dfs: driving forces of the reactions
%

% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2019 by Steffen Klamt and Axel von Kamp,
% Max Planck Institute for Dynamics of Complex Technical Systems, Magdeburg, Germany.
%
% Contributors are listed in CONTRIBUTORS.txt.
%
% This software can be used under the terms of our CellNetAnalyzer License.
% A copy of the license agreement is provided in the file named "LICENSE.txt"
% included with this software distribution. The license is also available online at
% http://www2.mpi-magdeburg.mpg.de/projects/cna/license.html
%
% For questions please contact: cellnetanalyzer@mpi-magdeburg.mpg.de

if nargin < 9+1
    stoichMat= cnap.stoichMat;
    if nargin < 8+1
        fixed_ratios= [];
        if nargin < 6+1
            D= [];
            d= [];
        end
    end
end

% Old precision: max(1e-9, cnap.epsilon)
[mdf, v, conc, dum1, dum2, dfs]= max_min_driving_force_pathway_ratio_range(maximalMilpRunTime, stoichMat,...
    cnap.stoichMat(cnap.specInternal, :), cnap.reacMin, cnap.reacMax, D, d,...
    1e-9, [], 1000, RT, G0, Cmin, Cmax, fixed_ratios,...
    isnan(G0(:, 1)), false, cnap.reacMin, cnap.reacMax);
