function [reactionSpeciesIndexMapping] = getReactionSpeciesIndexMapping(cnap, speciesIds, exchangeCompartmentId)
    numReactions = length(cnap.reacID);
    reactionSpeciesIndexMapping = zeros(1, numReactions);
    for currentReaction = 1:numReactions
        idAsString = getIdAsString(cnap.reacID(currentReaction,:));
        idSplit = strsplit(idAsString, "_");
        species = idSplit(end);
        if species ~= exchangeCompartmentId
            try
                reactionSpeciesIndexMapping(currentReaction) = find(species == speciesIds);
            catch % Nothing found
                continue
            end
        end
    end
end
