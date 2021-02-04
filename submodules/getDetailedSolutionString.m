function [fileContent] = getDetailedSolutionString(cnap, vCommunity, vEpsilon, dfCommunity, speciesIds, dG0sAndUncertainties)
    % Creates CellNetAnalyzer val files for the given pyredcom community
    % model and the given associated flux distribution vector. These val
    % files can be used to visualize flux pathways on CellNetAnalyzer maps.
    %
    % INPUT
    % cnap: CNAModel ~ The pyredcom community model
    % vCommunity: List[float] ~ The flux distribution vector, ordered after the order
    %  of reactions in cnap
    % vEpsilon: Float ~ The minimal neccessary flux value in order to be regarded
    %  in the .val file
    % valFilePathBase: String ~ Base path (i.e., folder and base file name
    %  for the created .val files. For each organism, a .val file is
    %  created, the single organism names are added to the base path as
    %  well as the ending ".val".
    %
    % OUTPUT
    % CellNetAnalyzer .val files for every species with the given file path
    if isempty(vCommunity)
        fileContent = "";
        return
    end

    fileContent = "";
    exchangeString = "Active reactions of exchange compartment<->environment exchange reactions:\n";
    for singleSpecies = speciesIds
        fileContent = fileContent + "Active reactions associated with species " + singleSpecies + " (Reaction ID; Flux; Direction; dG0 (adjusted to direction); dG; Reaction string (adjusted to direction):\n";
        activeReactionIDs = cnap.reacID(abs(vCommunity) > vEpsilon, :);
        for activeReactionIDIndex = 1:length(cellstr(activeReactionIDs))
            activeReactionIDString = getIdAsString(activeReactionIDs(activeReactionIDIndex,:));
            activeReactionIndex = PSBCNAFindReaction(activeReactionIDString, cnap);
            activeReactionIDString = "___" + activeReactionIDString;
            activeReactionIDString = strrep(activeReactionIDString, "___R_", "");
            activeReactionIDString = strrep(activeReactionIDString, "___", "");

            activeReactionIDStringSplit = strsplit(activeReactionIDString, "_");
            isEXCHGReactionOfSpecies = (activeReactionIDStringSplit(1) == "EXCHG") && (activeReactionIDStringSplit(2) == singleSpecies);
            isEXCReaction = (activeReactionIDStringSplit(1) == "EX") && (activeReactionIDStringSplit(2) == "C");
            isFirstSpecies = (singleSpecies == speciesIds(1));
            if ~(activeReactionIDStringSplit(end) == singleSpecies) && ~isEXCHGReactionOfSpecies && ~(isEXCReaction && isFirstSpecies)
                continue
            end

            fileLine = "";
            dG0 = dG0sAndUncertainties(activeReactionIndex);
            if isnan(dG0)
                dG0 = "NaN";
            end
            flux = vCommunity(activeReactionIndex);
            if flux < 0
                direction = "reverse";
            else
                direction = "forward";
            end
            df = -dfCommunity(activeReactionIndex);
            if isnan(df)
                df = "NaN";
            end
            fileLine = fileLine + activeReactionIDString + "; " + flux + "; " + direction + "; " + dG0 + "; " + df + "; ";

            % Create reaction string
            metaboliteMatrixColumn = cnap.stoichMat(:, activeReactionIndex);
            if vCommunity(activeReactionIndex) > 0
                eductMetaboliteIndices = find(metaboliteMatrixColumn < 0)';
                productMetaboliteIndices = find(metaboliteMatrixColumn > 0)';
            else
                eductMetaboliteIndices = find(metaboliteMatrixColumn > 0)';
                productMetaboliteIndices = find(metaboliteMatrixColumn < 0)';
                metaboliteMatrixColumn = (-1) * metaboliteMatrixColumn;
            end
            reactionString = "";
            for eductMetaboliteIndex = eductMetaboliteIndices
                reactionString = reactionString + abs(metaboliteMatrixColumn(eductMetaboliteIndex)) + " " + getIdAsString(cnap.specID(eductMetaboliteIndex,:)) + " ";
                if eductMetaboliteIndex ~= eductMetaboliteIndices(end)
                    reactionString = reactionString + "+ ";
                end
            end
            if cnap.reacMin(activeReactionIndex) < 0
                reactionString = reactionString + "(<)-> ";
            else
                reactionString = reactionString + "-> ";
            end
            for productMetaboliteIndex = productMetaboliteIndices
                reactionString = reactionString + abs(metaboliteMatrixColumn(productMetaboliteIndex)) + " " + getIdAsString(cnap.specID(productMetaboliteIndex,:)) + " ";
                if productMetaboliteIndex ~= productMetaboliteIndices(end)
                    reactionString = reactionString + "+ ";
                end
            end

            fileLine = fileLine + reactionString + "\n";

            isInternalSpeciesReaction = (sum(speciesIds == activeReactionIDStringSplit(end)) == 0);
            isSpeciesExchangeReaction = (activeReactionIDStringSplit(1) == "EXCHG") && (activeReactionIDStringSplit(2) == singleSpecies);
            if isInternalSpeciesReaction && ~isSpeciesExchangeReaction
                exchangeString = exchangeString + fileLine;
            else
                fileContent = fileContent + fileLine;
            end
        end
    end
    fileContent = fileContent + exchangeString;
end
