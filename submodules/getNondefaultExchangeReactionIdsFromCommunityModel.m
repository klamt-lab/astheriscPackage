function [eligibleExchanges, eligibleExchangesIndices] = getNondefaultExchangeReactionIdsFromCommunityModel(cnap, speciesIds)
    eligibleExchanges = [];
    eligibleExchangesIndices = [];
    for i = 1:numel(cellstr(cnap.reacID))
        reactionId = convertCharsToStrings(cnap.reacID(i,:));

        if ~startsWith(reactionId, "R_EXCHG_")
            continue
        end
        if ~((cnap.reacMin(i) == 0) && (cnap.reacMax(i) == 0))
            continue
        end

        reactionIdPart = reactionId;
        for j = 1:numel(speciesIds)
            reactionIdPart = strrep(reactionIdPart, "R_EXCHG_"+speciesIds(j)+"_", "");
        end
        reactionIdPart_split = strsplit(reactionIdPart, " ");
        usedReactionIdPart = reactionIdPart_split(1);

        eligibleExchanges = [eligibleExchanges usedReactionIdPart];
        eligibleExchangesIndices = [eligibleExchangesIndices i];
    end
end

