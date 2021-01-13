function [fileContent] = getValFileString(cnap, vCommunity, vEpsilon, speciesIds)
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
    
    fileContent = "";
    for singleSpecies = speciesIds
        fileContent = fileContent + ".val file for species " + singleSpecies + ":\n";
        activeReactionIDs = cnap.reacID(abs(vCommunity) > vEpsilon, :);
        for activeReactionIDIndex = 1:length(cellstr(activeReactionIDs))
            activeReactionIDString = getIdAsString(activeReactionIDs(activeReactionIDIndex,:));
            activeReactionIndex = PSBCNAFindReaction(activeReactionIDString, cnap);
            activeReactionIDString = "___" + activeReactionIDString;
            activeReactionIDString = strrep(activeReactionIDString, "___R_", "");
            activeReactionIDString = strrep(activeReactionIDString, "___", "");
            activeReactionIDStringSplit = strsplit(activeReactionIDString, "_");
            if activeReactionIDStringSplit(end) == singleSpecies
                activeReactionIDString = strrep(activeReactionIDString, "_"+singleSpecies, "");
                fileContent = fileContent + activeReactionIDString + " " + vCommunity(activeReactionIndex) + "\n";
            end
        end
        fileContent = fileContent + "\n## Dont delete this line\n";
        
        % valFilePath = valFilePathBase + singleSpecies + ".val";
        % file = fopen(valFilePath, "wt");
        % fprintf(file, fileContent);
        % fclose(file);
    end
end
