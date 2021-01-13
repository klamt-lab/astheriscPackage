function [allCommunityExchangeReactions] = getAllCommunityExchangeReactionsFromCommunityModel(eligibleExchanges, speciesIds)
    allCommunityExchangeReactions = [];
    for i = 1:numel(eligibleExchanges)
        eligibleExchangeId = eligibleExchanges(i);

        for j = 1:numel(speciesIds)
            fullReactionId = "R_EXCHG_"+speciesIds(j)+"_"+eligibleExchangeId;
            allCommunityExchangeReactions = [allCommunityExchangeReactions fullReactionId];
        end
    end
end
