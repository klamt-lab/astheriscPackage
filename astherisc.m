function [finalReport, mdfsWithCommunity, mdfsWithoutCommunity] = astherisc(targetSubstrateReaction, targetMaxSubstrateUptake,...
    minimalYieldFactor, minAbsoluteYield,...
    minimalMdfAdvantage, minimalMdfWithCommunity,...
    numMaxExchanges, temperature, vEpsilon,...
    modelPath, reportPath, reportHeader,...
    showCommunityModelBottlenecks, showSingleSpeciesBottlenecks, calculateIndirectBottlenecks,...
    reactionsWithSetMinimumZero, ignoredMetabolites, fixedRatiosNameMat, exchangeReactionSelection,...
    maximalMilpRunTime, verboseSettings, printValFiles)
    %% Algorithmic Search of THERmodynamic advantages In stoichiometric Single-species Community models
    % This is the ASTHERISC package's main function astherisc. For more about ASTHERISC itself, read the "README.md" in the same folder as this script.
    %
    % ==ARGUMENTS==
    % In order to use the astherisc function, you have to provide the following arguments in the given order:
    % targetSubstrateReaction : String ~ The exchange reaction with which the target (carbon) substrate can be taken up, e.g. "R_EX_C_glc__D_exchg".
    % targetMaxSubstrateUptake : Float ~ The maximal flux for the substrate uptake reaction given in targetSubstrateReaction for a single species. The negative value
    %     will be used as the actual lower limit.
    % minimalYieldFactor : Float ~ The percentage of the maximal yield reached without community that must be reached in a community solution.
    % minAbsoluteYield : Float ~ The necessary minimal absolute yield value.
    % minimalMdfAdvantage : Float ~ The necessary minimal absolute optimal MDF advantage of the community solution over a single-species solution.
    % minimalMdfWithCommunity : Float ~ The absolute minimal MDF value that must be reached by the solution with community.
    % numMaxExchanges : Int ~ The maximal number of allowed extra exchanges, can be a fixed positive integer or inf (in the latter case, there is no
    %     restriction of the maximal number of allowed extra exchanges).
    % temperature : Float ~ The assumed temperature in Kelvin. For reference: 25°C (77°F) equal 298.15 K,
    % vEpsilon : Float ~ The minimal absolute flux value that a reaction must have in order to be considered active.
    % modelPath : String ~ The (relative to CellNetAnalyzer's main folder or absolute) path to the .mat file containing the community metabolic model.
    % reportPath : String ~ The (relative to CellNetAnalyzer's main folder or absolute) path to the text report of the ASTHERISC run.
    %     If it is "", no report file will be written (although the finalReport variable will be still generateed,.
    % showCommunityModelBottlenecks : bool ~ Whether or not thermodynamic bottleneck reactions in the whole community model with its community solution shall be
    %     calculated and shown in the text report.
    % showSingleSpeciesBottlenecks : bool ~ Whether or not thermodynamic bottleneck reactions in the single species solutions with its community solution shall be
    %     calculated and shown in the text report.
    % calculateIndirectBottlenecks : bool ~ CURRENTLY DYSFUNCTIONAL, just set it to "false"
    % reactionsWithSetMinimumZero : List[String] ~ List of reaction IDs for which the minimal alloweed flux of the associated reactions shall
    %     be set to zero.
    % ignoredMetabolites : List[String] ~ List of metabolite IDs for which their associated metabolites a) have a fixed concentration range of [1;1] M i.e.,
    %     their concentrations are fixed to 1 M and b) are not unused, i.e., they are meant to have this concentration, and instead of
    %     "[unused]", "[1;1]" will be printed in the final text report for their metabolite concentrations.
    % fixedRatiosNameMat : Matrix ~ A matrix describing enforced metabolite concentration ratio ranges of the following structure:
    % [
    %    "$METABOLITE_ID_1" "METABOLITE_ID_2" $minimal_concentration_ratio $maximal_concentration_ratio
    % ]
    % exchangeReactionSelection : List[String] or [Empty] ~ If empty, all target metabolites as defined with commodepy in the model generation
    %     are being used. If not empty, on the the given target metabolites are used in ASTHERISC runs.
    % verboseSettings: bool ~ If true, all dG'0 of the reactions as well as all metabolite concentration ranges are printed in the settings. If false,
    %                         they will not be printed.
    % printValFiles: bool ~ If true, CellNetAnalyzer .val file strings (which can be e.g. directly copied into a separate text file) will be added to the text reports for each target metabolite with a community advantage,
    %                       including .val file strings for the single-species solution using the active community reactions as well as the community solution itself,
    %                       in the latter case separated for each single species of the community.
    %
    % ==OUTPUTS==
    % finalReport : String ~ The final text report. It will be also exported as text file to reportPath.
    % mdfsWithCommunity: Vector[float] ~ The list of calculated mdfs with
    %    the community in cases where a community solution provides
    %    an MDF advantage.
    % mdfsWithoutCommunity: Vector[float] ~ The list of calculated mdfs
    %    without the community in cases where a community solution provides
    %    an MDF advantage.

    fullRunTime = tic; % Start measuring the time

    %% Check if report path is for an existing file
    if exist(reportPath, "file") == 2
        disp("File "+reportPath+" already exists and will be overwritten!");
        userChoice = input('Do you want to continue? Y/N [default: N]: ', 's');
        if isempty(userChoice)
            userChoice = 'N';
        end
        if userChoice ~= 'Y'
            return
        end
    end

    %% Load model
    % [cnap] = TGLoadCommunityModel();
    % save(modelPath, 'cnOptimal MDF with communityap');
    load(modelPath, 'cnap');
    cnapReactionIdsCellstring = cellstr(cnap.reacID); % Used in PSBCNAFindReactionInCellstring


    %% Set up text report
    finalReport = reportHeader;
    finalReport = finalReport + "\n=SETTINGS=\n";
    finalReport = finalReport + "Substrate reaction: " + targetSubstrateReaction + "\n";
    finalReport = finalReport + "Maximal substrate uptake for single species (in units of the model): " + targetMaxSubstrateUptake + "\n";
    finalReport = finalReport + "Minimal yield factor (i.e., minimal necessary percentage of maximal possible yield): " + minimalYieldFactor + "\n";
    finalReport = finalReport + "Minimal absolute yield: " + minAbsoluteYield + "\n";
    finalReport = finalReport + "Minimal absolute MDF community advantage: " + minimalMdfAdvantage + "\n";
    finalReport = finalReport + "Temperature [K] (for comparison: 25°C=77°F=298.15 K): "+ temperature + "\n";
    finalReport = finalReport + "Minimal absolute mdf with community: "+ minimalMdfWithCommunity+ "\n";
    finalReport = finalReport + "Minimal absolute flux (in units of the model) for a reaction in order to be considered 'active' ('v epsilon'): " + num2str(vEpsilon) + "\n";
    finalReport = finalReport + "Model path: " + modelPath + "\n";
    finalReport = finalReport + "Report path: " + reportPath + "\n";
    finalReport = finalReport + "Calculate and show community solution bottlenecks? " + showCommunityModelBottlenecks + "\n";
    finalReport = finalReport + "Calculate and show single-species bottlenecks with active reactions of the community? " + showSingleSpeciesBottlenecks + "\n";
    finalReport = finalReport + "Calculate indirect bottlenecks for selected bottleneck types? " + calculateIndirectBottlenecks + "\n";
    finalReport = finalReport + "Number of maximal allowed exchanges: " + numMaxExchanges + "\n";
    finalReport = finalReport + "Reactions with set minimal flux of 0:\n";
    if isempty(reactionsWithSetMinimumZero)
        finalReport = finalReport + " None\n";
    else
        for reactionId = reactionsWithSetMinimumZero
            finalReport = finalReport + " " + reactionId + "\n";
        end
    end
    finalReport = finalReport + "'Ignored' metabolites (i.e., metabolites which are deemed to play no role thermodynamically as their concentration range is fixed to 1, so that instead of the default output '[unused]' in this case, '[1;1]' is printed instead:\n";
    if isempty(ignoredMetabolites)
        finalReport = finalReport + " None\n";
    else
        for metaboliteId = ignoredMetabolites
            finalReport = finalReport + " " + metaboliteId + "\n";
        end
    end
    finalReport = finalReport + "Fixed ratio ranges (ID of metabolite A, ID of metabolite B, minimal concentration ratio of A:B, maximal concentration ratio of A:B:\n";
    if isempty(fixedRatiosNameMat)
        finalReport = finalReport + " None\n";
    else
        for currentLine = 1:size(fixedRatiosNameMat, 1)
            ratiosNameMatLine = fixedRatiosNameMat(currentLine, :);
            finalReport = finalReport + " " + ratiosNameMatLine(1) + ", " + ratiosNameMatLine(2) + ", ";
            finalReport = finalReport + ratiosNameMatLine(3) + ", " + ratiosNameMatLine(4) + "\n";
        end
    end
    finalReport = finalReport + "Target exchange reaction selection:\n";
    if isempty(exchangeReactionSelection)
        finalReport = finalReport + " All possible ones\n";
    else
        for reactionId = exchangeReactionSelection
            finalReport = finalReport + " " + reactionId + "\n";
        end
    end

    %% Get species IDs
    [speciesIds] = getSpeciesIdsFromCommunityModel(cnap);

    %% Get species-specific substrate uptake reaction IDs
    targetSubstrate = strrep(targetSubstrateReaction, "R_EX_C_", "");
    targetSubstrate = strrep(targetSubstrate, "_exchg", "");
    speciesSubstrateReactions = [];
    for speciesId = speciesIds
        try
            reactionId = "R_EXCHG_" + speciesId + "_" + targetSubstrate + "_p_to_" + targetSubstrate;
            PSBCNAFindReactionInCellstring(reactionId, cnapReactionIdsCellstring); % Throws error if ID is not correct
        catch
            try
                reactionId = "R_EXCHG_" + speciesId + "_" + targetSubstrate + "_c_to_" + targetSubstrate;
                PSBCNAFindReactionInCellstring(reactionId, cnapReactionIdsCellstring); % Throws error if ID is not correct
            catch
                reactionId = "R_EXCHG_" + speciesId + "_" + targetSubstrate + "_to_" + targetSubstrate;
                PSBCNAFindReactionInCellstring(reactionId, cnapReactionIdsCellstring); % Throws error if ID is not correct
            end
        end
        speciesSubstrateReactions = [speciesSubstrateReactions reactionId];
    end

    %% Get all default community input and output metabolites
    defaultMetabolites = [];
    for i = 1:numel(cellstr(cnap.reacID))
        reactionId = getIdAsString(cnap.reacID(i,:));
        if ~startsWith(reactionId, "R_EX_C_")
            continue
        end
        if (cnap.reacMin(i) ~= 0) || (cnap.reacMax(i) ~= 0)
            reactionId = strrep(reactionId, "R_EX_C_", "");
            reactionId = strrep(reactionId, "_exchg", "");

            defaultMetabolites = [defaultMetabolites reactionId];
        end
    end

    %% Get non-default exchange reaction IDs
    [eligibleExchanges, allNondefaultExchangeReactionIndices] = getNondefaultExchangeReactionIdsFromCommunityModel(cnap, speciesIds);

    %% Get all non-default community exchange reactions
    [allCommunityExchangeReactions] = getAllCommunityExchangeReactionsFromCommunityModel(eligibleExchanges, speciesIds);

    %% Get all potential products
    if isempty(exchangeReactionSelection)
        potentialProductReactions = [];
        for i = 1:numel(cnap.reacMax)
            reactionId = getIdAsString(cnap.reacID(i,:));
            if ~startsWith(reactionId, "R_EX_C_")
                continue
            end
            potentialProductReactions = [potentialProductReactions reactionId];
        end
    else
        potentialProductReactions = exchangeReactionSelection;
    end

    %% Get general MDF MILP arguments for the given model (used for all versions of the MILP)
    [RT, dG0sAndUncertainties, minConcentrationsMat, maxConcentrationsMat, fixedConcentrationRatioRanges] = getGeneralMdfMilpParameters(cnap, temperature, fixedRatiosNameMat);

    if verboseSettings
        finalReport = finalReport + "Set dG0 values (reaction ID, dG0):\n";
        for currentReactionIndex = 1:length(dG0sAndUncertainties)
            reactionId = strtrim(convertCharsToStrings(cnap.reacID(currentReactionIndex, :)));
            dG0 = dG0sAndUncertainties(currentReactionIndex);
            finalReport = finalReport + " " + reactionId + ", " + pCatchNanInStr(dG0) + "\n";
        end

        finalReport = finalReport + "Set metabolite concentrations (metabolite ID, minimal concentration, maximal concentration):\n";
        for currentMetaboliteIndex = 1:length(minConcentrationsMat)
            metaboliteId = strtrim(convertCharsToStrings(cnap.specID(currentMetaboliteIndex, :)));
            minConcentration = minConcentrationsMat(currentMetaboliteIndex);
            maxConcentration = maxConcentrationsMat(currentMetaboliteIndex);
            finalReport = finalReport + " " + metaboliteId + ", " + minConcentration + ", " + maxConcentration + "\n";
        end
    end
    finalReport = finalReport + "\n=RUN RESULTS=";
    
    %% Set all given reaction's flux minimum to 0 (especially useful for ATPM reactions)
    for reactionWithSetMinimumZero = reactionsWithSetMinimumZero
        cnap.reacMin(PSBCNAFindReactionInCellstring(reactionWithSetMinimumZero, cnapReactionIdsCellstring)) = 0;
        cnap.reacMin(PSBCNAFindReactionInCellstring(reactionWithSetMinimumZero, cnapReactionIdsCellstring)) = 0;
    end

    %% Deactivate dG0=NaN reactions except of special transporters (as identified by BiGG IDs)
    cnapWithNanReactions = cnap;
    counter = 1;
    while counter < length(dG0sAndUncertainties)
        dG0 = dG0sAndUncertainties(counter);
        reactionId = convertCharsToStrings(cnap.reacID(counter, :));

        if isnan(dG0)
            if contains(reactionId, "tpp_") || contains(reactionId, "t1pp_") || contains(reactionId, "tex_") || endsWith(reactionId, "tex")
                % tpp: Facilitated transport or (proton) symport
                % t1pp: Facilitated transport or (proton) symport
                % tex: "Via diffusion"
                counter = counter + 1;
                continue
            end
            if ~startsWith(reactionId, "R_EX_C_") && ~startsWith(reactionId, "R_EXCHG_")
                disp(reactionId)
                cnap.reacMin(counter) = 0;
                cnap.reacMax(counter) = 0;
            end
        end
        counter = counter + 1;
    end

    %% Perform actual MDF community benefit analysis
    % -> Set up major result variables
    uniqueProductsWithCommunityBenefit = [];
    mdfsWithCommunity = [];
    mdfsWithoutCommunity = [];
    numbersOfExtraExchanges = [];
    extraExchangesList = [];
    forbiddenReactions = [];
    mdfAdvantages = [];
    maxyieldDifferencesSolutionsOnly = [];
    maxyieldDifferencesAll = [];
    activeExchangeCounter = zeros(1, length(cnap.reacMin));
    worksWithoutNan = 0;
    worksWithNan = 0;
    % -> Go through each potential product
    currentRound = 0; % Counter for current round
    targetProductReactionIndex = 1; % find(potentialProductReactions == "R_EX_C_clpn160_exchg")+1;  %%%
    while targetProductReactionIndex <= length(potentialProductReactions)
        targetProductReaction = potentialProductReactions(targetProductReactionIndex);

        % -> Print current round and target product reaction
        currentRound = currentRound + 1;
        disp("Round "+num2str(currentRound)+"/"+num2str(length(potentialProductReactions))+": "+targetProductReaction);

        % -> Get current target exchanges
        targetExchanges = [];
        targetExchangesIndices = [];
        metaboliteId_part = strrep(targetProductReaction, "R_EX_C_", "");
        metaboliteId = strrep(metaboliteId_part, "_exchg", "");
        for j = 1:numel(speciesIds)
            speciesId = speciesIds(j);
            targetExchange = "R_EXCHG_"+speciesId+"_"+metaboliteId+"_p_to_"+metaboliteId;

            try
                PSBCNAFindReactionInCellstring(targetExchange, cnapReactionIdsCellstring);
            catch
                try
                    targetExchange = "R_EXCHG_"+speciesId+"_"+metaboliteId+"_c_to_"+metaboliteId;
                    PSBCNAFindReactionInCellstring(targetExchange, cnapReactionIdsCellstring);
                catch
                    targetExchange = "R_EXCHG_"+speciesId+"_"+metaboliteId+"_to_"+metaboliteId;
                end
            end

            targetExchanges = [targetExchanges targetExchange];

            targetExchangesIndices = [targetExchangesIndices PSBCNAFindReactionInCellstring(targetExchange, cnapReactionIdsCellstring)];
        end

        % Starting ASTHERISC algorithm steps, with steps named as in
        % ASTHERISC's publication Note that some substeps (denoted with "A"
        % or "B" after the step number) may occur after some numbered steps
        % (e.g. step 3B after 2A) since the numbers reflect the order logic
        % as given in the publication which is different here due to the
        % described usage of MDF and yield approximations which are used in
        % order to reduce the computational complexity.
        % -> STEP 0A: Get maximal yield for the current product with
        %             NaN reactions without community capability
        % --> Delete any previous dYield as we use another "dYield" variable for the next MILP :D (without clear chaos would arise D:)
        clear dYield
        % --> Set yield function (c*r/d*r)
        cYield = zeros([1, cnap.numr]);
        cYield(1, PSBCNAFindReactionInCellstring(targetProductReaction, cnapReactionIdsCellstring)) = 1;

        dYield = zeros([1, cnap.numr]);
        dYield(1, PSBCNAFindReactionInCellstring(targetSubstrateReaction, cnapReactionIdsCellstring)) = -1; % Substrate uptake has a negative flux :O

        % --> Run \o/
        cnapWithNanReactionsWithoutCommunity = getCommunityModelWithSetExchanges(cnapWithNanReactions, targetSubstrateReaction, targetMaxSubstrateUptake, targetProductReaction, targetExchanges, speciesIds);
        cnapWithNanReactionsWithoutCommunity = pDeactivateAllExchangesOutsideSpecies1(cnapWithNanReactionsWithoutCommunity, speciesIds);
        [maxyieldWithNaNReactions, ~, success, ~] = CNAoptimizeYield(cnapWithNanReactionsWithoutCommunity,...
            cYield,...
            dYield,...
            [],...
            cnap.macroDefault,...
            2);
        % --> Catch potential errors
        if success ~= 1
            finalReport = finalReport + "\nWARNING: In NaN reactions model: No maximal yield could be calculated for "+convertCharsToStrings(targetProductReaction)+" and the model without allowed community\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end

        % -> STEP 0B: Get maximal yield for the current product without NaN reactions
        %    and without community capability
        % --> Get model without community
        cnapWithoutCommunity = getCommunityModelWithSetExchanges(cnap, targetSubstrateReaction, targetMaxSubstrateUptake, targetProductReaction, targetExchanges, speciesIds);
        cnapWithoutCommunity = pDeactivateAllExchangesOutsideSpecies1(cnapWithoutCommunity, speciesIds);
        % --> Run \o/
        [maxyield, ~, success, ~] = CNAoptimizeYield(cnapWithoutCommunity,...
            cYield,...
            dYield,...
            [],...
            cnap.macroDefault,...
            2);
        % --> Catch potential errors
        if success ~= 1
            finalReport = finalReport + "\nWARNING: No maximal yield could be calculated for "+convertCharsToStrings(targetProductReaction)+" and the model without allowed community (i.e., an error occurred in CNAoptimizeYield)\n";
            finalReport = finalReport + "Maxyield with NaN reactions: " + maxyieldWithNaNReactions + "\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end
        maxyieldDifferencesAll = [maxyieldDifferencesAll maxyieldWithNaNReactions-maxyield];
        % --> Collect statistics
        if (maxyieldWithNaNReactions > 0.0)
            if (maxyield <= 0.0)
                worksWithNan = worksWithNan + 1;
            else
                worksWithoutNan = worksWithoutNan + 1;
                worksWithNan = worksWithNan + 1;
            end
        end
        % --> Catch too low yield
        if maxyield < minAbsoluteYield
            finalReport = finalReport + "\nINFO: Maximal yield of " + maxyield + " is smaller than minimal yield for "+convertCharsToStrings(targetProductReaction)+" and the model without allowed community\n";
            finalReport = finalReport + "Maxyield with NaN reactions: " + maxyieldWithNaNReactions + "\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end

        % -> STEP 0C: Check if the maximal production flux is high
        %    enough in the model without NaN reactions
        cnapWithoutCommunityFBACheck = cnapWithoutCommunity;
        for currentObjectiveReactionIndex = 1:length(cnap.objFunc)
            cnapWithoutCommunityFBACheck.objFunc(currentObjectiveReactionIndex) = 0;
        end
        currentTargetIndex = PSBCNAFindReactionInCellstring(targetProductReaction, cnapReactionIdsCellstring);
        cnapWithoutCommunityFBACheck.objFunc(currentTargetIndex) = -1;
        [~, success, ~, optval] = CNAoptimizeFlux(cnapWithoutCommunityFBACheck, [], cnap.macroDefault, 2);
        % --> Catch no success
        if success ~= 1
            finalReport = finalReport + "\nWARNING: No maximal product flux could be calculated for "+convertCharsToStrings(targetProductReaction)+" and the model without allowed community\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end
        % --> Catch too low flux
        if abs(optval) < vEpsilon
            finalReport = finalReport + "\nINFO: Maximal target product production reaction rate of " + optval + " is smaller than vEpsilon (set to "+vEpsilon+") for "+convertCharsToStrings(targetProductReaction)+" and the model without allowed community\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end

        % -> Set new yield minimum for subsequent MDF MILPs with the
        % current product
        targetMinYield = maxyield * minimalYieldFactor; % A bit lower in order to prevent numeric problems

        % -->Set desired phenotype for OptMDFPathway-based MILPs with
        %    the minimal yield
        numDesiredConst = 1;
        D = zeros([numDesiredConst, cnap.numr]);
        % Product/Substrate >= minYield <=> minYield*Substrate - Product <= 0
        D(1, PSBCNAFindReactionInCellstring(targetSubstrateReaction, cnapReactionIdsCellstring)) = -targetMinYield; % Substrate uptake has a negative flux :O
        D(1, PSBCNAFindReactionInCellstring(targetProductReaction, cnapReactionIdsCellstring)) = -1;
        d(1) = 0;
        d = d';


        % -> STEP 1: Get optimal MDF without community using
        % OptMDFPathway, "without community" in the definition that all
        % exchange reactions are blocked.
        [optmdfWithoutCommunity, vWithoutCommunity, ~, ~] = CNAcomputeOptMDFpathway_ratio_range(maximalMilpRunTime, cnapWithoutCommunity,...
            RT,...
            dG0sAndUncertainties,...
            minConcentrationsMat,...
            maxConcentrationsMat,...
            D,...
            d,...
            fixedConcentrationRatioRanges);

        % --> Catch potential NaN result
        if isnan(optmdfWithoutCommunity)
            finalReport = finalReport + "\nINFO: No result (probably infeasible or timeout) with OptMDFPathway and no allowed community for "+convertCharsToStrings(targetProductReaction)+"\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end

        reachedYieldWithoutCommunity = vWithoutCommunity(PSBCNAFindReactionInCellstring(targetProductReaction, cnapReactionIdsCellstring)) / -vWithoutCommunity(PSBCNAFindReactionInCellstring(targetSubstrateReaction, cnapReactionIdsCellstring));


        % -> STEP 2A:
        % Check if there is a higher MDF (at least higher or equal to
        % minimalMdf) with community than
        % the optimal one without community using a "cut" OptMDFPathway.
        % --> Get community exchange reactions
        if isempty(forbiddenReactions)
            usedCommunityExchangeReactions = allCommunityExchangeReactions;
        else
            usedCommunityExchangeReactions = allCommunityExchangeReactions(~ismember(allCommunityExchangeReactions, forbiddenReactions));
        end
        usedCommunityExchangeReactionsIndices = [];
        for reactionCounter = 1:length(usedCommunityExchangeReactions)
            currentReactionId = usedCommunityExchangeReactions(reactionCounter);
            currentReactionIndex = PSBCNAFindReactionInCellstring(currentReactionId, cnapReactionIdsCellstring);
            usedCommunityExchangeReactionsIndices = [usedCommunityExchangeReactionsIndices currentReactionIndex];
        end
        % --> Get model with community capability
        cnapWithCommunityAndAllExchanges = getCommunityModelWithSetExchanges(cnap, targetSubstrateReaction, targetMaxSubstrateUptake, targetProductReaction, usedCommunityExchangeReactions, speciesIds);

        % --> Run \o/
        minimalMdf = optmdfWithoutCommunity + minimalMdfAdvantage;
        if minimalMdf < minimalMdfWithCommunity
            minimalMdf = minimalMdfWithCommunity;
        end
        measureMilp2Time = tic;
        [mdfCommunityAndHigherMdf, vCommunityAndHigherMdf, ~, dfsCommunity] = CNAcomputeOptMDFpathway_higher_mdf_max_exchanges(maximalMilpRunTime, numMaxExchanges, allNondefaultExchangeReactionIndices, minimalMdf, false,...
            cnapWithCommunityAndAllExchanges,...
            RT,...
            dG0sAndUncertainties,...
            minConcentrationsMat,...
            maxConcentrationsMat,...
            D,...
            d,...
            fixedConcentrationRatioRanges);
        timeUsed = toc(measureMilp2Time);
        % --> Catch potential NaN result
        if isnan(mdfCommunityAndHigherMdf)
            if timeUsed >= 999.9999
                finalReport = finalReport + "\nWARNING: MILP time excess with OptMDFPathway-based MILP2 and allowed community for "+convertCharsToStrings(targetProductReaction)+"\n";
            else
                finalReport = finalReport + "\nINFO: MILP infeasible with OptMDFPathway-based MILP2 and allowed community for "+convertCharsToStrings(targetProductReaction)+"\n";
            end
            finalReport = finalReport + "Reached MDF without community is " + optmdfWithoutCommunity + "\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end
        % --> Recheck if the MDF with community is actually higher than the one without community
        %    If this is not the case, proceed with the next target metabolite :O
        if (mdfCommunityAndHigherMdf) <= optmdfWithoutCommunity - abs(optmdfWithoutCommunity)*.0001
            finalReport = finalReport + "\nINFO: No higher MDF found with community compared to no community for "+convertCharsToStrings(targetProductReaction)+"\n";
            finalReport = finalReport + "Reached MDF without community is " + optmdfWithoutCommunity + "\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end


        % Test for numeric problem: In the community solution, check if *both*
        % species are actually used.
        numberActiveReactionsPerSpecies = zeros(1, length(speciesIds));
        for currentIndex = 1:length(cnap.reacMin)
           if abs(vCommunityAndHigherMdf(currentIndex)) > vEpsilon
               for currentSpecies = 1:length(speciesIds)
                   reactionIdAsStr = getIdAsString(cnap.reacID(currentIndex,:));
                   if endsWith(reactionIdAsStr, "_"+speciesIds(currentSpecies))
                       numberActiveReactionsPerSpecies(currentSpecies) = numberActiveReactionsPerSpecies(currentSpecies) + 1;
                   end
               end
           end
        end
        if sum(numberActiveReactionsPerSpecies == 0) == (length(numberActiveReactionsPerSpecies) - 1)
            finalReport = finalReport + "\nWARNING: Numeric problem! Higher MDF found with community compared to no community even though only one species is used for "+convertCharsToStrings(targetProductReaction)+"\n";
            finalReport = finalReport + "Reached MDF without community is " + optmdfWithoutCommunity + "\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            pWriteReport(reportPath, finalReport);
            continue
        end



        % -> STEP 2B: Perform optmdf approximation with the community model
        %    using an approximation loop
        stepSizes = [1, 0.5, 0.1];
        optmdfWithCommunity = NaN;
        for stepSize = stepSizes
            if isnan(optmdfWithCommunity)
                currentMinimalMdf = optmdfWithoutCommunity;
            else
                currentMinimalMdf = optmdfWithCommunity;
            end
            while 1
                currentMinimalMdf = currentMinimalMdf + stepSize;
                [newOptmdfWithCommunity, newVOptMdf, ~, dfsCommunity] = CNAcomputeOptMDFpathway_higher_mdf_max_exchanges(maximalMilpRunTime, numMaxExchanges, allNondefaultExchangeReactionIndices, currentMinimalMdf, false,...
                    cnapWithCommunityAndAllExchanges,...
                    RT,...
                    dG0sAndUncertainties,...
                    minConcentrationsMat,...
                    maxConcentrationsMat,...
                    D,...
                    d,...
                    fixedConcentrationRatioRanges);
                if isnan(newOptmdfWithCommunity)
                    break
                else
                    optmdfWithCommunity = newOptmdfWithCommunity;
                    mdfConstraint = optmdfWithCommunity;
                    vOptMdf = newVOptMdf;
                end
            end
        end
        % --> Catch numeric problem
        if (optmdfWithCommunity < minimalMdf - 0.0001)
            finalReport = finalReport + "\nWARNING: Numeric problem (community OptMDF is lower than single species OptMDF, still MILP2 works) for "+convertCharsToStrings(targetProductReaction)+"\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            continue
        end
        % --> Catch NaN
        if isnan(optmdfWithCommunity)
            errorOptMdfCommunity = true;
            vOptMdf = vCommunityAndHigherMdf;
            mdfConstraint = mdfCommunityAndHigherMdf - 0.00001;
            optmdfWithCommunity = mdfCommunityAndHigherMdf - 0.00001;
        else
            errorOptMdfCommunity = false;
        end

        % -> Check for further numeric problems in solution (sometimes,
        % reactions with a *very* low flux occur, even though they are
        % shown with a flux of "0" in the resulting OptMDFPathway flux
        % vector D:)
        cnapTest = cnapWithCommunityAndAllExchanges;
        dG0sAndUncertaintiesTest = dG0sAndUncertainties;
        for i = 1:length(cnap.reacMin)
            if vOptMdf(i) == 0
                cnapTest.reacMin(i) = 0;
                cnapTest.reacMax(i) = 0;
            end
        end
        [optmdfWithCommunityTest, ~, ~, ~] = CNAcomputeOptMDFpathway_ratio_range_max_exchanges(maximalMilpRunTime, numMaxExchanges, allNondefaultExchangeReactionIndices, cnapTest,...
                    RT,...
                    dG0sAndUncertaintiesTest,...
                    minConcentrationsMat,...
                    maxConcentrationsMat,...
                    D,...
                    d,...
                    fixedConcentrationRatioRanges);
        if (optmdfWithCommunityTest) - abs(optmdfWithCommunityTest)*.000001 <= optmdfWithoutCommunity
            finalReport = finalReport + "\nINFO: Numeric problem! Higher MDF found with community compared to no community cannot be reproduced for "+convertCharsToStrings(targetProductReaction)+"\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            pWriteReport(reportPath, finalReport);
            continue
        end


        % -> STEP 3A: Perform yield optimization at the maximal
        %  community MDF
        stepSizes = [0.5, 0.1, 0.05, 0.01];
        maxyieldCommunityAtOptmdf = NaN;
        for stepSize = stepSizes
            if isnan(maxyieldCommunityAtOptmdf)
                currentMinimalYield = targetMinYield;
            else
                currentMinimalYield = maxyieldCommunityAtOptmdf;
            end
            while 1
                currentMinimalYield = currentMinimalYield + stepSize;

                % -->Set desired phenotype for OptMDFPathway-based MILPs with
                %    the minimal yield
                D_maxyield = zeros([1, cnap.numr]);
                % Product/Substrate >= minYield <=> minYield*Substrate - Product <= 0
                D_maxyield(1, PSBCNAFindReactionInCellstring(targetSubstrateReaction, cnapReactionIdsCellstring)) = -currentMinimalYield; % Substrate uptake has a negative flux :O
                D_maxyield(1, PSBCNAFindReactionInCellstring(targetProductReaction, cnapReactionIdsCellstring)) = -1;
                d_maxyield(1) = 0;
                d_maxyield = d_maxyield';

                [newOptmdfWithCommunity, newVOptMdf, ~, dfsCommunity] = CNAcomputeOptMDFpathway_higher_mdf_max_exchanges(maximalMilpRunTime, numMaxExchanges, allNondefaultExchangeReactionIndices, mdfConstraint, false,...
                    cnapWithCommunityAndAllExchanges,...
                    RT,...
                    dG0sAndUncertainties,...
                    minConcentrationsMat,...
                    maxConcentrationsMat,...
                    D_maxyield,...
                    d_maxyield,...
                    fixedConcentrationRatioRanges);
                if isnan(newOptmdfWithCommunity)
                    break
                else
                    optmdfWithCommunity = newOptmdfWithCommunity;
                    mdfConstraint = newOptmdfWithCommunity;
                    vOptMdf = newVOptMdf;
                    maxyieldCommunityAtOptmdf = currentMinimalYield;
                end
            end
        end

        % -->Set desired phenotype for OptMDFPathway-based MILPs with
        %    the maximal yield at the maximal MDF
        if isnan(maxyieldCommunityAtOptmdf)
            errorApproximatedYield = true;
            maxyieldCommunityAtOptmdf = targetMinYield;
        else
            errorApproximatedYield = false;
        end

        D = zeros([1, cnap.numr]);
        % Product/Substrate >= minYield <=> minYield*Substrate - Product <= 0
        D(1, PSBCNAFindReactionInCellstring(targetSubstrateReaction, cnapReactionIdsCellstring)) = -maxyieldCommunityAtOptmdf; % Substrate uptake has a negative flux :O
        D(1, PSBCNAFindReactionInCellstring(targetProductReaction, cnapReactionIdsCellstring)) = -1;
        d(1) = 0;
        d = d';

        % -> STEP 2C: Calculate exact optimal MDF at the approximated optimal yield
        cnapWithMaxyieldMaxmdfReactionsOnly = pGetCnapWithBlockedReactionsWhereVIsZero(cnapWithCommunityAndAllExchanges, vOptMdf, vEpsilon);
        [mdfWithMaxyieldMaxmdf, vWithMaxyieldMaxmdf, ~, dfsCommunity] = CNAcomputeOptMDFpathway_ratio_range_max_exchanges(maximalMilpRunTime, numMaxExchanges, allNondefaultExchangeReactionIndices,...
            cnapWithMaxyieldMaxmdfReactionsOnly,...
            RT,...
            dG0sAndUncertainties,...
            minConcentrationsMat,...
            maxConcentrationsMat,...
            D,...
            d,...
            fixedConcentrationRatioRanges);
        if ~isnan(mdfWithMaxyieldMaxmdf)
            optmdfWithCommunity = mdfWithMaxyieldMaxmdf - 0.0001;
            mdfConstraint = mdfWithMaxyieldMaxmdf - 0.0001;
            vOptMdf = vWithMaxyieldMaxmdf;
            exactCommunityMdfCalculationError = false;
            mdfAdvantage = mdfConstraint - optmdfWithoutCommunity;
            mdfAdvantages = [mdfAdvantages mdfAdvantage];
        else
            exactCommunityMdfCalculationError = true;
        end


        % -> STEP 4: Get mdf solution with minimal flux at maximal
        % community MDF & yield at this MDF
        cnapWithCommunityAndMinimalFluxReactionsOnly = pGetCnapWithBlockedReactionsWhereVIsZero(cnapWithCommunityAndAllExchanges, vOptMdf, vEpsilon);
        if ~errorOptMdfCommunity
            [mdfMinimalFlux, vMinimalFlux, ~, dfsCommunity] = CNAcomputeOptMDFpathway_higher_mdf_max_exchanges(maximalMilpRunTime, numMaxExchanges, allNondefaultExchangeReactionIndices, mdfConstraint - abs(mdfConstraint)*.000001, true,...
                cnapWithCommunityAndMinimalFluxReactionsOnly,...
                RT,...
                dG0sAndUncertainties,...
                minConcentrationsMat,...
                maxConcentrationsMat,...
                D,...
                d,...
                fixedConcentrationRatioRanges);
            % Catch NaN
            if isnan(mdfMinimalFlux)
                errorOptMdfMinimalFlux = true;
            else
                errorOptMdfMinimalFlux = false;
                vOptMdf = vMinimalFlux;
            end
        else
            errorOptMdfMinimalFlux = true;
        end

        % --> STEP 3B: Calculate exact maximal yield for text report
        cnapWithCommunityAndMinimalFluxReactions = pGetCnapWithBlockedReactionsWhereVIsZero(cnapWithCommunityAndMinimalFluxReactionsOnly, vOptMdf, vEpsilon);
        % --> Set yield function (c*r/d*r) and get optimal yield with minimal active reactions
        [reachedYield, ~, success, ~] = CNAoptimizeYield(cnapWithCommunityAndMinimalFluxReactions,...
            cYield,...
            dYield,...
            [],...
            cnap.macroDefault,...
            2);
        % --> Catch potential errors
        if success ~= 1
            errorMaximalCommunityYield = true;
            reachedYield = -vOptMdf(PSBCNAFindReactionInCellstring(targetProductReaction, cnapReactionIdsCellstring)) / vOptMdf(PSBCNAFindReactionInCellstring(targetSubstrateReaction, cnapReactionIdsCellstring));
        else
            errorMaximalCommunityYield = false;
        end

        % --> Get number of reactions and total flux for each species
        %     with the minimal total flux solution
        substrateUptakeRate = vOptMdf(PSBCNAFindReactionInCellstring(targetSubstrateReaction, cnapReactionIdsCellstring));
        fluxScaleVector = abs(targetMaxSubstrateUptake / substrateUptakeRate);
        numberSpecies = length(speciesIds);
        numberReactionsPerSpecies = zeros(numberSpecies, 1);
        totalFluxPerSpecies = zeros(numberSpecies, 1);
        currentReactionIndex = 1;
        while currentReactionIndex < length(vOptMdf)
            vSelection = vOptMdf(currentReactionIndex);
            if abs(vSelection) < vEpsilon
                currentReactionIndex = currentReactionIndex + 1;
                continue
            end

            current_reactionId = getIdAsString(cnap.reacID(currentReactionIndex, :));
            currentSpeciesIndex = 1;
            for speciesId = speciesIds
                if ~endsWith(current_reactionId, speciesId)
                    currentSpeciesIndex = currentSpeciesIndex + 1;
                    continue
                end

                numberReactionsPerSpecies(currentSpeciesIndex) = numberReactionsPerSpecies(currentSpeciesIndex) + 1;
                totalFluxPerSpecies(currentSpeciesIndex) = totalFluxPerSpecies(currentSpeciesIndex) + abs(vSelection);
                currentSpeciesIndex = currentSpeciesIndex + 1;
            end
            currentReactionIndex = currentReactionIndex + 1;
        end
        speciesSpecificReactionReport = "";
        currentSpeciesIndex = 1;
        for speciesId = speciesIds
            numberReactions = numberReactionsPerSpecies(currentSpeciesIndex);
            total_flux = totalFluxPerSpecies(currentSpeciesIndex) * fluxScaleVector;
            speciesSpecificReactionReport = speciesSpecificReactionReport + speciesId + ": " + num2str(numberReactions) + " | " + num2str(total_flux) + "\n";
            currentSpeciesIndex = currentSpeciesIndex + 1;
        end



        % -> STEP 5: Get active metabolite concentration ranges
        %   Based on (heavily modified) code provided by Axel von Kamp :-)
        % --> Get active metabolites (i.e., metabolites which occur in
        % reactions that are active in the OptMDFPathway solution)
        activeMetaboliteIndices = find(any(cnapWithCommunityAndMinimalFluxReactions.stoichMat(:, (abs(vOptMdf)>vEpsilon)), 2));

        % --> Get concentration range vectors
        [cMins, cMaxs, activeMetaboliteIndices, errorConcentrationRanges] = getConcentrationRanges(numMaxExchanges, allNondefaultExchangeReactionIndices, mdfConstraint,...
            cnapWithCommunityAndMinimalFluxReactions,...
            D, d,...
            dG0sAndUncertainties,...
            RT,...
            minConcentrationsMat,...
            maxConcentrationsMat,...
            fixedConcentrationRatioRanges,...
            activeMetaboliteIndices,...
            maximalMilpRunTime);

        % --> Find non-overlapping metabolite concentration ranges
        occuringMetabolitesStr = "";
        nonOccuringMetabolitesStr = "";
        analyzedMetaboliteIndices = [];
        for metaboliteIndex = 1:length(cMins)
            % Check if there is a metabolite range with NaN
            currentCMin = cMins(metaboliteIndex);
            currentCMax = cMaxs(metaboliteIndex);
            % Ignore if cMin and cMax are equal to 0
            if isnan(currentCMin) || isnan(currentCMax)
                continue
            end

            % Ignore if metabolite is already analyzed
            if ismember(metaboliteIndex, analyzedMetaboliteIndices)
                continue
            end

            % Get the metabolite's ID
            try
                baseMetaboliteId = cnap.specID(metaboliteIndex,:);
            catch
                disp("A")
            end
            baseMetaboliteId = getIdAsString(baseMetaboliteId);
            baseMetaboliteSpecies = pGetSpeciesFromId(baseMetaboliteId);
            % Catch exchanges
            if ~ismember(baseMetaboliteSpecies, speciesIds)
                continue
            end

            % Get metabolite index in all species
            metaboliteIndexInAllSpecies = [];
            metaboliteIdInAllSpecies = [];
            for speciesId = speciesIds
                speciesMetaboliteId = strrep(baseMetaboliteId, "_"+baseMetaboliteSpecies, "_"+speciesId);
                species_metaboliteIndex = PSBCNAFindMetabolite(speciesMetaboliteId, cnap);
                metaboliteIndexInAllSpecies = [metaboliteIndexInAllSpecies species_metaboliteIndex];
                metaboliteIdInAllSpecies = [metaboliteIdInAllSpecies speciesMetaboliteId];
            end

            % Check if the metabolite ranges are overlapping
            firstIndex = metaboliteIndexInAllSpecies(1);
            selectionCMin = cMins(firstIndex);
            selectionCMax = cMaxs(firstIndex);
            foundNan = false;
            for unused_counter_xD = 1:length(metaboliteIndexInAllSpecies)
                for currentSpecies = 1:length(metaboliteIndexInAllSpecies)
                    currentMetaboliteIndexInSpecies = metaboliteIndexInAllSpecies(currentSpecies);
                    currentCMin = cMins(currentMetaboliteIndexInSpecies);
                    currentCMax = cMaxs(currentMetaboliteIndexInSpecies);
                    if isnan(currentCMin) || isnan(currentCMax)
                        foundNan = true;
                    end
                    % Extension below
                    if (currentCMin < selectionCMin) && (currentCMax >= selectionCMin)
                        selectionCMin = currentCMin;
                    end
                    % Extension above
                    if (currentCMax > selectionCMax) && (currentCMin <= selectionCMax)
                        selectionCMax = currentCMax;
                    end
                end
            end
            speciesWideMetaboliteCMin = min(cMins(metaboliteIndexInAllSpecies));
            speciesWideMetaboliteCMax = max(cMaxs(metaboliteIndexInAllSpecies));
            % Catch metabolites which are not active in all species
            if foundNan
                continue
            end

            % Add metabolite with its concentration range if it is
            % non-overlapping to text report
            if (selectionCMin ~= speciesWideMetaboliteCMin) || (selectionCMax ~= speciesWideMetaboliteCMax)
                speciesReport = "";
                foundNaN = false;
                for currentSpecies = 1:length(metaboliteIndexInAllSpecies)
                    metaboliteIndex = metaboliteIndexInAllSpecies(currentSpecies);
                    metaboliteId = metaboliteIdInAllSpecies(currentSpecies);
                    currentCMin = exp(cMins(metaboliteIndex));
                    currentCMax = exp(cMaxs(metaboliteIndex));
                    currentCMinAsStr = num2str(currentCMin);
                    currentCMaxAsStr = num2str(currentCMax);
                    if (currentCMin == currentCMax  == 1)
                        foundNan = true;
                    end
                    speciesReport = speciesReport + metaboliteId + " " + pAddConcentrationString(currentCMinAsStr, currentCMaxAsStr, metaboliteId, ignoredMetabolites) + " & ";
                end
                if foundNan
                    nonOccuringMetabolitesStr = nonOccuringMetabolitesStr + speciesReport + "X";
                else
                    occuringMetabolitesStr = occuringMetabolitesStr + speciesReport + "X";
                end
            end

            analyzedMetaboliteIndices = [analyzedMetaboliteIndices metaboliteIndexInAllSpecies];
        end
        numOverlappingMetabolitesStr = occuringMetabolitesStr + nonOccuringMetabolitesStr;
        % Post-process output string of non-overlapping metabolite
        % concentration ranges
        numOverlappingMetabolitesStr = numOverlappingMetabolitesStr;
        numOverlappingMetabolitesStr = strrep(numOverlappingMetabolitesStr, " & X", "\n");
        if numOverlappingMetabolitesStr == "X"
            numOverlappingMetabolitesStr = "";
        end

        %%%%
        % --> Get all non-default active exchanges from the solution with
        %     community and minimal extra exchanges (5th main step) as a list of strings containing the
        %     exchange IDs
        activeExchanges = cnap.reacID(abs(vOptMdf)>vEpsilon, :);
        activeExchangesAsStrList = [];
        activeExchangesIndices = [];
        for j = 1:length(cellstr(activeExchanges))
            reactionId = getIdAsString(activeExchanges(j,:));
            if ~startsWith(reactionId, "R_EXCHG_")
                continue
            end

            reactionIdSplit = strsplit(reactionId, "_to_");
            metaboliteId = reactionIdSplit(2);
            if ismember(metaboliteId, defaultMetabolites)
                continue
            end

            index = PSBCNAFindReactionInCellstring(reactionId, cnapReactionIdsCellstring);
            activeExchangeCounter(index) = activeExchangeCounter(index) + 1;
            activeExchangesIndices = [activeExchangesIndices index];
            activeExchangesAsStrList = [activeExchangesAsStrList reactionId];
        end
        numberOfExtraExchanges = length(activeExchangesAsStrList);

        % --> Get extra exchanges and default exchanges output string
        extraExchangesAsStr = "";
        extraExchangesAsVector = [];
        defaultExchangesAsStr = "";
        for j = 1:length(cellstr(activeExchanges))
            reactionId = getIdAsString(activeExchanges(j,:));
            if ~startsWith(reactionId, "R_EXCHG_")
                continue
            end

            reactionIdSplit = strsplit(reactionId, "_to_");
            metaboliteId = reactionIdSplit(2);

            reactionIndex = PSBCNAFindReactionInCellstring(reactionId, cnapReactionIdsCellstring);
            scaledFlux = fluxScaleVector * vOptMdf(reactionIndex);
            if abs(scaledFlux) < vEpsilon
                continue
            end
            if ismember(reactionId, speciesSubstrateReactions)
                continue
            end
            exchangeReportAsStr = reactionId + " | SCALED FLUX: " + num2str(scaledFlux);

            % Get exchange report after identifying the exchange's
            % metabolites
            reaction_row = cnap.stoichMat(:, reactionIndex);
            metabolite_indices = find(reaction_row ~= 0);
            for metaboliteIndex = metabolite_indices'
                metaboliteId = getIdAsString(cnap.specID(metaboliteIndex,:));
                metaboliteCRangeAsStr = pAddConcentrationString(exp(cMins(metaboliteIndex)), exp(cMaxs(metaboliteIndex)), metaboliteId, ignoredMetabolites);
                exchangeReportAsStr = exchangeReportAsStr + " | " + metaboliteId + " " + metaboliteCRangeAsStr;
            end

            exchangeReportAsStr = exchangeReportAsStr + "\n";
            rawMetaboliteId = strrep("X"+metaboliteId, "XM_", "");
            rawMetaboliteId = strrep(rawMetaboliteId+"X", "_exchgX", "");
            if ismember(rawMetaboliteId, defaultMetabolites)
                disp("DEFAULT")
                disp(metaboliteId)
                defaultExchangesAsStr = defaultExchangesAsStr + exchangeReportAsStr;
            else
                disp("NON-DEFAULT")
                disp(metaboliteId)
                extraExchangesAsStr = extraExchangesAsStr + exchangeReportAsStr;
                extraExchangesAsVector = [extraExchangesAsVector reactionId];
            end
        end



        % -> STEP 6A: Get active species and calculate bottlenecks
        % with the single-species solutions for active species in the
        % community-free solution with the minimal absolute flux sum
        % Get bottlenecks of each active species
        singleSpeciesBottlenecksStr = "";
        hasNumericError = false;
        if showSingleSpeciesBottlenecks
            for currentSpecies = speciesIds(1) % Should be changed for multi-species communities
                cnapWithOneSpeciesOnly = pGetCnapWithActiveReactionsInBothSpecies(cnapWithCommunityAndAllExchanges, vOptMdf, speciesIds, vEpsilon);
                cnapWithOneSpeciesOnly = getCommunityModelWithSetExchanges(cnapWithOneSpeciesOnly, targetSubstrateReaction, targetMaxSubstrateUptake, targetProductReaction, targetExchanges, speciesIds);
                cnapWithOneSpeciesOnly = pDeactivateAllExchangesOutsideSpecies1(cnapWithOneSpeciesOnly, speciesIds);
                for currentReactionIndex = 1:length(cnapWithOneSpeciesOnly.reacMin)
                    current_reactionId = getIdAsString(cnapWithOneSpeciesOnly.reacID(currentReactionIndex,:));
                    if ~startsWith(current_reactionId, "R_EXCHG_")
                        continue
                    end
                    if ~startsWith(current_reactionId, "R_EXCHG_"+currentSpecies+"_")
                        cnapWithOneSpeciesOnly.reacMin(currentReactionIndex) = 0;
                        cnapWithOneSpeciesOnly.reacMax(currentReactionIndex) = 0;
                    end
                end

                [mdfWithSingleSpecies, vWithSingleSpecies, ~, dfsWithSingleSpecies] = CNAcomputeOptMDFpathway_ratio_range_max_exchanges(maximalMilpRunTime, numMaxExchanges, allNondefaultExchangeReactionIndices, cnapWithOneSpeciesOnly,...
                    RT,...
                    dG0sAndUncertainties,...
                    minConcentrationsMat,...
                    maxConcentrationsMat,...
                    D,...
                    d,...
                    fixedConcentrationRatioRanges);
                if isnan(mdfWithSingleSpecies)
                    vWithSingleSpecies = [];
                    continue
                end
                if (mdfWithSingleSpecies-0.001) > optmdfWithoutCommunity
                    finalReport = finalReport + "WARNING: Inconsistent single-species solution for " + targetProductReaction + "! This is probably caused by a numeric problem\n";
                end

                cnapWithOneSpeciesOnlyAndActiveReactionsOnly = pGetCnapWithBlockedReactionsWhereVIsZero(cnapWithOneSpeciesOnly, vWithSingleSpecies, vEpsilon);

                [~, directBottlenecksStr, indirectBottlenecksStr, connectionMetaboliteIndices] = getDirectAndIndirectBottlenecks(numMaxExchanges, allNondefaultExchangeReactionIndices, cnapWithOneSpeciesOnlyAndActiveReactionsOnly,...
                    dfsWithSingleSpecies, RT, dG0sAndUncertainties,...
                    minConcentrationsMat, maxConcentrationsMat,...
                    D, d, fixedConcentrationRatioRanges,...
                    mdfWithSingleSpecies, speciesIds, vWithSingleSpecies, vEpsilon, calculateIndirectBottlenecks,...
                    maximalMilpRunTime);

                [cMinsSingleSpecies, cMaxsSingleSpecies, ~, ~] = getConcentrationRanges(numMaxExchanges, allNondefaultExchangeReactionIndices, mdfWithSingleSpecies,...
                    cnapWithOneSpeciesOnlyAndActiveReactionsOnly,...
                    D, d,...
                    dG0sAndUncertainties,...
                    RT,...
                    minConcentrationsMat, maxConcentrationsMat,...
                    fixedConcentrationRatioRanges,...
                    connectionMetaboliteIndices,...
                    maximalMilpRunTime);

                connectionMetaboliteStr = "";
                for connectionMetaboliteIndex = connectionMetaboliteIndices
                    connectionMetaboliteId = getIdAsString(cnap.specID(connectionMetaboliteIndex, :));

                    singleSpeciesMinConcentration = exp(cMinsSingleSpecies(connectionMetaboliteIndex));
                    singleSpeciesMinConcentration = pCatchNanInStr(singleSpeciesMinConcentration);

                    singleSpeciesMaxConcentration = exp(cMaxsSingleSpecies(connectionMetaboliteIndex));
                    singleSpeciesMaxConcentration = pCatchNanInStr(singleSpeciesMaxConcentration);

                    connectionMetaboliteStr = connectionMetaboliteStr + "~" + connectionMetaboliteId + "~" + "\nIn single species:\n";
                    connectionMetaboliteStr = connectionMetaboliteStr + pAddConcentrationString(singleSpeciesMinConcentration, singleSpeciesMaxConcentration, connectionMetaboliteId, ignoredMetabolites);
                    if (~isstring(singleSpeciesMinConcentration)) && (~isstring(singleSpeciesMaxConcentration))
                        singleSpeciesConcentrationRange = singleSpeciesMaxConcentration - singleSpeciesMinConcentration;
                        connectionMetaboliteStr = connectionMetaboliteStr + " (range " + singleSpeciesConcentrationRange + ")";
                    end
                    connectionMetaboliteStr = connectionMetaboliteStr + "\nIn community:\n";

                    metaboliteSpeciesId = pGetSpeciesFromId(connectionMetaboliteId);
                    metaboliteSpeciesId = metaboliteSpeciesId + "XXX";
                    for singleSpeciesId = speciesIds
                        tempMetaboliteSpeciesId = strrep(connectionMetaboliteId+"XXX", "_"+metaboliteSpeciesId, "_"+singleSpeciesId);
                        tempMetaboliteIndex = PSBCNAFindMetabolite(tempMetaboliteSpeciesId, cnap);

                        communityMinConcentration = exp(cMins(tempMetaboliteIndex));
                        communityMinConcentration = pCatchNanInStr(communityMinConcentration);

                        communityMaxConcentration = exp(cMaxs(tempMetaboliteIndex));
                        communityMaxConcentration = pCatchNanInStr(communityMaxConcentration);

                        connectionMetaboliteStr = connectionMetaboliteStr + tempMetaboliteSpeciesId +  ": " + pAddConcentrationString(communityMinConcentration, communityMaxConcentration, tempMetaboliteSpeciesId, ignoredMetabolites);
                        if (~isstring(communityMaxConcentration)) && (~isstring(communityMinConcentration))
                            communityConcentrationRange = communityMaxConcentration - communityMinConcentration;
                            connectionMetaboliteStr = connectionMetaboliteStr + " (range " + communityConcentrationRange + ")";
                        end
                        if singleSpeciesId ~= speciesIds(end)
                            connectionMetaboliteStr = connectionMetaboliteStr + "\n";
                        end
                    end
                    connectionMetaboliteStr = connectionMetaboliteStr + "\n";
                end

                difference = optmdfWithCommunity-mdfWithSingleSpecies;
                if difference < 0
                    hasNumericError = true;
                end
                singleSpeciesBottlenecksStr = singleSpeciesBottlenecksStr + "Bottleneck reactions (driving force=OptMDF) of single-species solution (with " + currentSpecies + " only):\n";
                singleSpeciesBottlenecksStr = singleSpeciesBottlenecksStr + "MDF in single species with reactions occuring in community solution: " + mdfWithSingleSpecies + " (OptMDF in community was "+optmdfWithCommunity+", difference is "+difference+")\n";
                singleSpeciesBottlenecksStr = singleSpeciesBottlenecksStr + ">Bottleneck reactions (driving force=OptMDF):\n" + directBottlenecksStr;
                singleSpeciesBottlenecksStr = singleSpeciesBottlenecksStr + ">Connecting metabolite ranges:\n" + connectionMetaboliteStr;
                if calculateIndirectBottlenecks
                    singleSpeciesBottlenecksStr = singleSpeciesBottlenecksStr + ">Indirect bottleneck reactions:\n" + indirectBottlenecksStr;
                end
            end
        else
            vWithSingleSpecies = [];
        end

        if hasNumericError
            finalReport = finalReport + "\nINFO: Numeric problem! Higher MDF found with community compared to single species cannot be reproduced for "+convertCharsToStrings(targetProductReaction)+"\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            pWriteReport(reportPath, finalReport);
            continue
        end

        if optmdfWithoutCommunity > optmdfWithCommunity
            finalReport = finalReport + "\nINFO: Numeric problem! Higher MDF wihtout community than with community for "+convertCharsToStrings(targetProductReaction)+"\n";
            targetProductReactionIndex = targetProductReactionIndex + 1;
            pWriteReport(reportPath, finalReport);
            continue
        end

        % -> STEP 6B: Community model bottleneck calculation
        if showCommunityModelBottlenecks
            [allBottleneckReactionsStr, directBottlenecksStr, indirectBottlenecksStr] = getDirectAndIndirectBottlenecks(numMaxExchanges, allNondefaultExchangeReactionIndices, cnapWithCommunityAndMinimalFluxReactions,...
                                                                    dfsCommunity, RT, dG0sAndUncertainties,...
                                                                    minConcentrationsMat, maxConcentrationsMat,...
                                                                    D, d, fixedConcentrationRatioRanges,...
                                                                    optmdfWithCommunity, speciesIds, vOptMdf, vEpsilon, calculateIndirectBottlenecks,...
                                                                    maximalMilpRunTime);
        end



        % -> (STEP 7:) Assemble text report for target metabolite
        % --> Get .val file string for single species in community solution
        % and single-species solution
        if printValFiles
            valFileString = "CellNetAnalyzer .val files for different solutions (formatted for e.g. iJO1366 model map):";
            % Community solution
            valFileString = valFileString + "\n~Community solution .val files~\n";
            valFileString = valFileString + getValFileString(cnap, vOptMdf, vEpsilon, speciesIds);
            % Single-species solution
            valFileString = valFileString + "\n~Single-species solution .val file~\n";
            firstSpecies = speciesIds(1);
            valFileString = valFileString + getValFileString(cnap, vWithSingleSpecies, vEpsilon, firstSpecies);
        else
            valFileString = "";
        end

        % --> Get solution reaction data string
        detailedSolutionString = "\nDetailed reaction-wise solution reports:";
        % Community solution
        detailedSolutionString = detailedSolutionString + "\n~Community solution detailed solutions~\n";
        detailedSolutionString = detailedSolutionString + getDetailedSolutionString(cnap, vOptMdf, vEpsilon, speciesIds, dG0sAndUncertainties);
        % Single-species solution
        detailedSolutionString = detailedSolutionString + "\n~Single-species solution detailed solution~\n";
        firstSpecies = speciesIds(1);
        if ~isempty(vWithSingleSpecies)
            detailedSolutionString = detailedSolutionString + getDetailedSolutionString(cnap, vWithSingleSpecies, vEpsilon, firstSpecies, dG0sAndUncertainties);
        end
        % --> Get species-specific substrate uptake string
        speciesUptakeAsStr = "";
        for speciesSubstrateReaction = speciesSubstrateReactions
            scaledFlux = fluxScaleVector * vOptMdf(PSBCNAFindReactionInCellstring(speciesSubstrateReaction, cnapReactionIdsCellstring));
            speciesUptakeAsStr = speciesUptakeAsStr + speciesSubstrateReaction + " | SCALED FLUX " + num2str(scaledFlux) + "\n";
        end

        % --> Set output variables
        uniqueProductsWithCommunityBenefit = [uniqueProductsWithCommunityBenefit targetProductReaction];
        mdfsWithCommunity = [mdfsWithCommunity optmdfWithCommunity];
        mdfsWithoutCommunity = [mdfsWithoutCommunity optmdfWithoutCommunity];
        numbersOfExtraExchanges = [numbersOfExtraExchanges numberOfExtraExchanges];
        extraExchangesList = [extraExchangesList extraExchangesAsStr];

        % --> Expand final report text for the current target product and
        %     write the expanded report to the text file
        finalReport = finalReport + "\n==TARGET WITH COMMUNITY ADVANTAGE: " + convertCharsToStrings(targetProductReaction) + "==\n";
        if errorOptMdfCommunity
            finalReport = finalReport + "WARNING: The optimal MDF with the calculated minimal exchanges could not be calculated due to status unknown (e.g. timeout), showing solution of MILP2 instead\n";
        end
        if errorOptMdfMinimalFlux
            finalReport = finalReport + "WARNING: The minimal absolute flux sum with the optimal MDF could not be calculated due to status unknown (e.g. timeout), showing solution of OptMDFPathway with (minimal, if no error occured) exchanges instead\n";
        end
        if errorMaximalCommunityYield
            finalReport = finalReport + "WARNING: No maximal yield could be calculated for "+convertCharsToStrings(targetProductReaction)+" and the model with allowed community, showing yield with minimal flux solution instead\n";
        end
        if errorApproximatedYield
            finalReport = finalReport + "WARNING: Approximated optimal yield at approximated optimal MDF for the community could not be calculated, this is possible e.g. due to timeouts.\n";
        end
        if exactCommunityMdfCalculationError
            finalReport = finalReport + "WARNING: Exact maximal MDF calculation with community at maximal yield failed (e.g. due to timeout)\n";
        end
        if errorConcentrationRanges
            finalReport = finalReport + "WARNING: Metabolite concentration range calculation failed (e.g. due to numeric problem or previous warning-inducing errors)\n";
        end
        finalReport = finalReport + "Optimal yield without community and with deactivated dG0=NaN reactions without any MDF constraint: " + num2str(maxyieldWithNaNReactions) + "\n";
        finalReport = finalReport + "Optimal yield without community without any MDF constraint: " + num2str(maxyield) + "\n";
        finalReport = finalReport + "Optimal MDF without community with minimal necessary yield as minimal yield constraint: " + num2str(optmdfWithoutCommunity) + "\n";
        finalReport = finalReport + "Reached yield without community at optimal MDF solution: " + reachedYieldWithoutCommunity + "\n";
        finalReport = finalReport + "Optimal MDF with community with minimal necessary yield as minimal yield constraint (if no warning given): " + num2str(optmdfWithCommunity) + "\n";
        finalReport = finalReport + "Approximated optimal yield with community at maximal community MDF (if no warning given): " + num2str(reachedYield) + "\n";
        finalReport = finalReport + "Number of active metabolites: " + num2str(length(activeMetaboliteIndices)) + "\n";

        finalReport = finalReport + "\nNumber of active reactions | Absolute total flux per species; Both with minimal absolute flux sum solution and substrate uptake scaled to maximum:\n" + speciesSpecificReactionReport + "\n";
        finalReport = finalReport + "Extra exchanges and scaled flux with minimal absolute flux solution (negative flux means uptake to species, positive flux secretion from species):\n" + extraExchangesAsStr + "\n";
        finalReport = finalReport + "Used default exchanges and scaled flux with minimal absolute flux solution (negative flux means uptake to species, positive flux secretion from species):\n" + defaultExchangesAsStr + "\n";
        finalReport = finalReport + "Species-specific uptake:\n" + speciesUptakeAsStr + "\n";
        if showCommunityModelBottlenecks
            finalReport = finalReport + "All bottleneck reactions:\n" + allBottleneckReactionsStr + "\n";
            finalReport = finalReport + "Direct bottleneck reactions:\n" + directBottlenecksStr + "\n";
            finalReport = finalReport + "Indirect bottleneck reactions:\n" + indirectBottlenecksStr + "\n";
        end
        finalReport = finalReport + "Metabolites with non-overlapping concentration ranges between species:\n" + numOverlappingMetabolitesStr + "\n";
        finalReport = finalReport + singleSpeciesBottlenecksStr + "\n";
        finalReport = finalReport + valFileString;
        finalReport = finalReport + detailedSolutionString;
        pWriteReport(reportPath, finalReport);

        targetProductReactionIndex = targetProductReactionIndex + 1;
    end

    %% --> Expand final report with statistics from all found community benefit metabolite
    maxyieldDifferencesSolutionsOnly = [maxyieldDifferencesSolutionsOnly maxyieldWithNaNReactions-maxyield];
    finalReport = finalReport + "\n==FINAL STATISTICS==\n";
    finalReport = finalReport + "Total number analyzed target metabolites:\n" + length(potentialProductReactions) + "\n";
    finalReport = finalReport + "Number (fraction) producible target metabolites with dG0=NaN reactions:\n";
    finalReport = finalReport + + worksWithNan + " (" + num2str(100 * (worksWithNan / length(potentialProductReactions))) + " pct.)\n";
    finalReport = finalReport + "Number (fraction) producible target metabolites without dG0=NaN reactions:\n";
    finalReport = finalReport + worksWithoutNan + " (" + num2str(100 * (worksWithoutNan / length(potentialProductReactions))) + " pct.)\n";
    if ~isempty(numbersOfExtraExchanges)
        finalReport = finalReport + "Maximal extra exchanges:\n" + max(numbersOfExtraExchanges) + "\n";
    end
    if ~isempty(numbersOfExtraExchanges)
        finalReport = finalReport + "Minimal extra exchanges:\n" + min(numbersOfExtraExchanges) + "\n";
    end
    finalReport = finalReport + "Mean extra exchanges:\n";
    if isempty(numbersOfExtraExchanges)
        finalReport = finalReport + "N/A\n";
    else
        finalReport = finalReport + mean(numbersOfExtraExchanges) + "\n";
    end
    finalReport = finalReport + "Mean maximal possible yield change due to NaN reactions in MDF community advantage solutions:\n";
    if isempty(maxyieldDifferencesSolutionsOnly)
        finalReport = finalReport + "N/A\n";
    else
        finalReport = finalReport + mean(maxyieldDifferencesSolutionsOnly) + "\n";
    end
    finalReport = finalReport + "Mean maximal possible yield change due to NaN reactions in all solutions: \n";
    if isempty(maxyieldDifferencesAll)
        finalReport = finalReport + "N/A\n";
    else
        finalReport = finalReport + mean(maxyieldDifferencesAll) + "\n";
    end
    finalReport = finalReport + "Number metabolites with found community MDF benefit:\n" +  length(uniqueProductsWithCommunityBenefit) + "\n";
    finalReport = finalReport + "Number metabolites without found community MDF benefit:\n" +  num2str(length(potentialProductReactions) - length(uniqueProductsWithCommunityBenefit)) + "\n";
    finalReport = finalReport + "Mean MDF advantage (found MDF with community minus MDF without community in all cases where the community showed an advantage and where OptMDFPathway worked with the community):\n";
    if isempty(mdfAdvantages)
        finalReport = finalReport + "N/A\n";
    else
        finalReport = finalReport + mean(mdfAdvantages) + "\n";
    end
    finalReport = finalReport + "List of metabolites with found community MDF benefit:\n";
    for benefitProduct = uniqueProductsWithCommunityBenefit
        finalReport = finalReport + benefitProduct + "\n";
    end
    finalReport = finalReport + "List of used exchange reactions, ordered by number of community advantage solutions where they occur:\n";
    currentValue = max(activeExchangeCounter);
    while currentValue >= 1
        reactionsWithCurrentValue = find(activeExchangeCounter == currentValue);
        for reactionIndex = reactionsWithCurrentValue
            reactionId = cnap.reacID(reactionIndex, :);
            finalReport = finalReport + getIdAsString(reactionId) + ": ";
            finalReport = finalReport + currentValue + "\n";
        end
        currentValue = currentValue - 1;
    end
    pWriteReport(reportPath, finalReport);
    fullRunTimeInSeconds = toc(fullRunTime);
    finalReport = finalReport + "\n\nFull run time (min): " + fullRunTimeInSeconds/60 + "\n";
    finalReport = finalReport + "Full run time (h): " + fullRunTimeInSeconds/3600 + "\n";
    finalReport = finalReport + "Full run time (days): " + fullRunTimeInSeconds/86400 + "\n";
    pWriteReport(reportPath, finalReport);
end


%% Small subfunctions used in astherisc()
function concentrationString = pAddConcentrationString(minConcentration, maxConcentration, metaboliteId, ignoredMetabolites)
    minConcentrationConversion = str2double(minConcentration);
    if ~isnan(minConcentrationConversion)
        minConcentration = minConcentrationConversion;
    end
    maxConcentrationConversion = str2double(maxConcentration);
    if ~isnan(maxConcentrationConversion)
        maxConcentration = maxConcentrationConversion;
    end
    if minConcentration == maxConcentration == 1
        if ~isempty(ignoredMetabolites)
            if ismember(metaboliteId, ignoredMetabolites)
                concentrationString = "[1;1]";
            else
                concentrationString = "[unused]";
            end
        else
            concentrationString = "[unused]";
        end
    else
        concentrationString = "[" + minConcentration + ";" + maxConcentration + "]";
    end
end


function str = pCatchNanInStr(str)
    if isnan(str)
        str = "NaN";
    end
end


function cnap = pDeactivateAllExchangesOutsideSpecies1(cnap, speciesIds)
    species1Id = speciesIds(1);
    species1ExchangeStart = "R_EXCHG_"+species1Id+"_";
    for reactionCounter = 1:length(cnap.reacMin)
        currentReactionId = getIdAsString(cnap.reacID(reactionCounter, :));
        if startsWith(currentReactionId, "R_EXCHG_")
            if ~startsWith(currentReactionId, species1ExchangeStart)
                cnap.reacMin(reactionCounter) = 0;
                cnap.reacMax(reactionCounter) = 0;
            end
        end
    end
end


function [cnap, deactivatedReactions] = pGetCnapWithActiveReactionsInBothSpecies(cnap, v, speciesIds, vEpsilon)
    allAnalyzedIndices = [];
    deactivatedReactions = [];
    cnapReactionIdsCellstring = cellstr(cnap.reacID);
	for currentReactionIndex = 1:length(v)
        if abs(v(currentReactionIndex)) < vEpsilon
            if ismember(currentReactionIndex, allAnalyzedIndices)
                continue
            end

            reactionId = getIdAsString(cnap.reacID(currentReactionIndex, :));
            currentSpecies = pGetSpeciesFromId(reactionId);
            speciesReactionIndices = [];
            if ismember(currentSpecies, speciesIds)
                for speciesId = speciesIds

                    speciesReactionId = strrep(reactionId+"X", currentSpecies+"X", speciesId);
                    speciesReactionIndex =  PSBCNAFindReactionInCellstring(speciesReactionId, cnapReactionIdsCellstring);
                    speciesReactionIndices = [speciesReactionIndices speciesReactionIndex];
                end
            end

            allAnalyzedIndices = [allAnalyzedIndices speciesReactionIndices];
            if sum(abs(v(speciesReactionIndices))) < (vEpsilon * length(speciesReactionIndices))
                if ~isempty(speciesReactionIndices)
                    deactivatedReactions = [deactivatedReactions speciesReactionId];
                    cnap.reacMin(speciesReactionIndices) = 0;
                    cnap.reacMax(speciesReactionIndices) = 0;
                end
            end
        end
	end
end


function [cnap, deactivatedReactions] = pGetCnapWithBlockedReactionsWhereVIsZero(cnap, v, vEpsilon)
    deactivatedReactions = [];
    for currentReactionIndex = 1:length(v)
        if abs(v(currentReactionIndex)) < vEpsilon
            cnap.reacMin(currentReactionIndex) = 0;
            cnap.reacMax(currentReactionIndex) = 0;
            deactivatedReactions = [deactivatedReactions getIdAsString(cnap.reacID(currentReactionIndex, :))];
        end
    end
end


function species = pGetSpeciesFromId(id)
    idSplit = strsplit(id, "_");
    species = idSplit(end);
end


function pWriteReport(reportPath, finalReport)
    if reportPath ~= ""
        file = fopen(reportPath, "wt");
        fprintf(file, finalReport);
        fclose(file);
    end
end
