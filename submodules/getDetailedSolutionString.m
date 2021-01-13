function [fileContent] = getDetailedSolutionString(cnap, vCommunity, vEpsilon, speciesIds, dG0sAndUncertainties)
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
    for singleSpecies = speciesIds
        fileContent = fileContent + "Detailed solution for species " + singleSpecies + " (Reaction ID; Flux; dG0; Reaction string):\n";
        activeReactionIDs = cnap.reacID(abs(vCommunity) > vEpsilon, :);
        for activeReactionIDIndex = 1:length(cellstr(activeReactionIDs))
            activeReactionIDString = getIdAsString(activeReactionIDs(activeReactionIDIndex,:));
            activeReactionIndex = PSBCNAFindReaction(activeReactionIDString, cnap);
            activeReactionIDString = "___" + activeReactionIDString;
            activeReactionIDString = strrep(activeReactionIDString, "___R_", "");
            activeReactionIDString = strrep(activeReactionIDString, "___", "");
            activeReactionIDStringSplit = strsplit(activeReactionIDString, "_");

            if (activeReactionIDStringSplit(end) == singleSpecies) || ((singleSpecies == speciesIds(1)) && (sum(speciesIds == activeReactionIDStringSplit(end)) == 0))
                activeReactionIDString = strrep(activeReactionIDString, "_" + singleSpecies, "");
                if isnan(dG0sAndUncertainties(activeReactionIndex))
                    dG0string = "NaN";
                else
                    dG0string = dG0sAndUncertainties(activeReactionIndex);
                end
                fileContent = fileContent + activeReactionIDString + "; " + vCommunity(activeReactionIndex) + "; " + dG0string + "; ";
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
                fileContent = fileContent + reactionString + "\n";
            end
        end
        % valFilePath = valFilePathBase + singleSpecies + ".val";
        % file = fopen(valFilePath, "wt");
        % fprintf(file, fileContent);
        % fclose(file);
    end
end
