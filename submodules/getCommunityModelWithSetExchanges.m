function cnap = getCommunityModelWithSetExchanges(cnap, targetSubstrateReaction, targetMaxSubstrateUptake, targetProductReaction, openExchangeReactions, speciesIds)
    % This function contains all community-advantage-search-related changes
    % which are applied to the iJO1366 community model.
    reactionIdsCellstring = cellstr(cnap.reacID);

    num_species = numel(speciesIds);

    %% Set maximal uptake for whole community and single species
    targetSubstrateReactionIndex = PSBCNAFindReactionInCellstring(targetSubstrateReaction, reactionIdsCellstring);
    substrateReactionMetabolites = cnap.stoichMat(:, targetSubstrateReactionIndex);
    substrateReactionMetaboliteIndex = find(substrateReactionMetabolites<0);
    substrateMetaboliteReactions = cnap.stoichMat(substrateReactionMetaboliteIndex, :);
    substrateMetaboliteReactionIndices = find(substrateMetaboliteReactions~=0);
    for substrateMetaboliteReactionIndex = substrateMetaboliteReactionIndices
        if substrateMetaboliteReactionIndex ~= targetSubstrateReactionIndex
            cnap.reacMin(substrateMetaboliteReactionIndices) = -inf;
            cnap.reacMax(substrateMetaboliteReactionIndices) = 0;
        end
    end
    cnap.reacMin(targetSubstrateReactionIndex) = -targetMaxSubstrateUptake;
    cnap.reacMax(targetSubstrateReactionIndex) = -targetMaxSubstrateUptake;

    %% Open product reaction
    if targetProductReaction ~= ""
        cnap.reacMax(PSBCNAFindReactionInCellstring(targetProductReaction, reactionIdsCellstring)) = inf;
    end

    %% Set exchange reaction upper and lower bounds
    for exchangeReactionId = openExchangeReactions
        exchangeReactionIndex = PSBCNAFindReactionInCellstring(exchangeReactionId, reactionIdsCellstring);
        cnap.reacMin(exchangeReactionIndex) = -inf;
        cnap.reacMax(exchangeReactionIndex) = inf;
    end
end
