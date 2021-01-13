function [allBottleneckReactionsStr, directBottlenecksStr, indirectBottlenecksStr, connectionMetaboliteIndices] = getDirectAndIndirectBottlenecks(numMaxExchanges, exchangeIndices, cnap, dfs,...
                                                                RT, dG0sAndUncertainties,...
                                                                minConcentrationsMat, maxConcentrationsMat,...
                                                                D, d, fixedConcentrationRatioRanges,...
                                                                unchangedMdf, speciesIds, v, vEpsilon,...
                                                                calculateIndirectBottlenecks,...
                                                                maximalMilpRunTime)
    for vIndex = 1:length(v)
        if abs(v(vIndex)) < vEpsilon
            dfs(vIndex) = NaN;
        end
    end

    cnapReactionIdCellstring = cellstr(cnap.reacID);

    % --> Get reactions with df at mdf
    reactionWithDfAtMdfIndices = find((dfs >= unchangedMdf - .01) & (dfs <= unchangedMdf + .01)); % Some room in order to prevent numeric problems
    % --> Create text report of all reactions with df at mdf
    allBottleneckReactionsStr = "";
    for currentBottleNeckIndexI = 1:length(reactionWithDfAtMdfIndices)
        currentBottleneckIndex = reactionWithDfAtMdfIndices(currentBottleNeckIndexI);
        bottleneckReactionId = getIdAsString(cnap.reacID(currentBottleneckIndex, :));
        allBottleneckReactionsStr = allBottleneckReactionsStr + bottleneckReactionId + "\n";
    end

    % --> Get all *direct* bottlenecks
    allMetaboliteColumns = [];
    directBottleneckIndices = [];
    directBottlenecksStr = ""; % Text report variable
    for changedReactionIndexI = 1:length(reactionWithDfAtMdfIndices)
        changedReactionIndex = reactionWithDfAtMdfIndices(changedReactionIndexI);

        % Check mdf change with changes dG0s of the selected single reaction,
        % if the mdf is higher with this change we found direct bottleneck \o/
        is_direct_bottleneck = pIsMdfHigherWithChangedDeltaG0s(numMaxExchanges, exchangeIndices, cnap,...
                                                               RT, dG0sAndUncertainties,...
                                                               minConcentrationsMat, maxConcentrationsMat,...
                                                               D, d, fixedConcentrationRatioRanges,...
                                                               changedReactionIndex, unchangedMdf,...
                                                               maximalMilpRunTime);
        if is_direct_bottleneck
            directBottleneckIndices = [directBottleneckIndices changedReactionIndex];
            changed_reaction_id = getIdAsString(cnap.reacID(changedReactionIndex, :));

            % Add reaction ID to output string
            directBottlenecksStr = directBottlenecksStr + changed_reaction_id + " | ";

            % Get all affected metabolites and add to output string
            metaboliteMatrixColumn = cnap.stoichMat(:, changedReactionIndex);
            if v(changedReactionIndex) > 0
                eductMetaboliteIndices = find(metaboliteMatrixColumn < 0)';
                productMetaboliteIndices = find(metaboliteMatrixColumn > 0)';
            else
                eductMetaboliteIndices = find(metaboliteMatrixColumn > 0)';
                productMetaboliteIndices = find(metaboliteMatrixColumn < 0)';
                metaboliteMatrixColumn = (-1) * metaboliteMatrixColumn;
            end
            allMetaboliteColumns = [allMetaboliteColumns metaboliteMatrixColumn];
            for eductMetaboliteIndex = eductMetaboliteIndices
                directBottlenecksStr = directBottlenecksStr + abs(metaboliteMatrixColumn(eductMetaboliteIndex)) + " " + getIdAsString(cnap.specID(eductMetaboliteIndex,:)) + " ";
                if eductMetaboliteIndex ~= eductMetaboliteIndices(end)
                    directBottlenecksStr = directBottlenecksStr + "+ ";
                end
            end
            if cnap.reacMin(changedReactionIndex) < 0
                directBottlenecksStr = directBottlenecksStr + "(<)-> ";
            else
                directBottlenecksStr = directBottlenecksStr + "-> ";
            end
            for productMetaboliteIndex = productMetaboliteIndices
                directBottlenecksStr = directBottlenecksStr + abs(metaboliteMatrixColumn(productMetaboliteIndex)) + " " + getIdAsString(cnap.specID(productMetaboliteIndex,:)) + " ";
                if productMetaboliteIndex ~= productMetaboliteIndices(end)
                    directBottlenecksStr = directBottlenecksStr + "+ ";
                end
            end

            % Add newline to output string
            directBottlenecksStr = directBottlenecksStr + "\n";
        end
    end


    % Get connection status
    connectionMetaboliteIndices = [];
    directBottlenecksStr = directBottlenecksStr + ">Connecting metabolites:\n";
	if length(directBottleneckIndices) >= 2
        for metaboliteCounter = 1:length(allMetaboliteColumns)
            reactionLine = allMetaboliteColumns(metaboliteCounter, :);
            numProducts = sum(reactionLine < 0);
            numEducts = sum(reactionLine > 0);
            if (numProducts >= 1) && (numEducts >= 1)
                eductReactionIndexInMatrix = find(reactionLine < 0);
                productReactionIndexInMatrix = find(reactionLine > 0);

                eductReactionIndex = directBottleneckIndices(eductReactionIndexInMatrix);
                productReactionIndex = directBottleneckIndices(productReactionIndexInMatrix);

                eductReactionId = "";
                for currentEductReactionIndex = eductReactionIndex
                    eductReactionId = eductReactionId + getIdAsString(cnap.reacID(currentEductReactionIndex, :));
                    if currentEductReactionIndex ~= eductReactionIndex(end)
                        eductReactionId = eductReactionId + "/";
                    end
                end
                productReactionId = "";
                for currentProductReactionIndex = productReactionIndex
                    productReactionId = productReactionId + getIdAsString(cnap.reacID(currentProductReactionIndex, :));
                    if currentProductReactionIndex ~= productReactionIndex(end)
                        productReactionId = productReactionId + "/";
                    end
                end

                directBottlenecksStr = directBottlenecksStr + productReactionId;
                directBottlenecksStr = directBottlenecksStr + "-> ";
                directBottlenecksStr = directBottlenecksStr + eductReactionId;
                directBottlenecksStr = directBottlenecksStr + "; ";
                directBottlenecksStr = directBottlenecksStr + getIdAsString(cnap.specID(metaboliteCounter, :));
                directBottlenecksStr = directBottlenecksStr + " \n";
                connectionMetaboliteIndices = [connectionMetaboliteIndices metaboliteCounter];
            end
        end
    end

    % --> Get all species-wide *indirect* bottlenecks
    indirectBottleneckIndices = [];
    indirectBottlenecksStr = ""; % Text report variable
    if calculateIndirectBottlenecks
        for currentReactionIndexI = 1:length(reactionWithDfAtMdfIndices)
            currentReactionIndex = reactionWithDfAtMdfIndices(currentReactionIndexI);

            % Exclude direct bottlenecks
            if ismember(currentReactionIndex, directBottleneckIndices)
                continue
            end
            % Exclude previously found indirect bottlenecks
            if ismember(currentReactionIndex, indirectBottleneckIndices)
                continue
            end

            % Get the reaction's indices and ids in all species
            currentReactionId_full = getIdAsString(cnap.reacID(currentReactionIndex, :));
            currentReactionId = currentReactionId_full;
            if startsWith(currentReactionId, "R_EXCHG_")
                currentReactionId = strrep(currentReactionId_full, "R_EXCHG_", "");
                for currentSpeciesId = speciesIds
                    currentReactionId = strrep(currentReactionId, currentSpeciesId + "_", "R_EXCHG_");
                end
            else
                for currentSpeciesId = speciesIds
                    currentReactionId = strrep(currentReactionId, "_"+currentSpeciesId, "");
                end
            end

            reactionIndexInAllSpecies = [];
            reactionIdInAllSpecies = [];
            for currentSpeciesId = speciesIds
                if startsWith(currentReactionId, "R_EXCHG_")
                    reactionIdInSpecies = strrep(currentReactionId, "R_EXCHG_", "R_EXCHG_"+currentSpeciesId+"_");
                else
                    reactionIdInSpecies = currentReactionId+"_"+currentSpeciesId;
                end
                reactionIdInAllSpecies = [reactionIdInAllSpecies reactionIdInSpecies];
                reaction_index_of_species = PSBCNAFindReactionInCellstring(reactionIdInSpecies, cnapReactionIdCellstring);
                reactionIndexInAllSpecies = [reactionIndexInAllSpecies reaction_index_of_species];
            end
            if sum(ismember(reactionIndexInAllSpecies, directBottleneckIndices)) > 0
                continue
            end

            % Check mdf change with changed dG0s of the selected reaction in all species,
            % if the mdf is higher with this change we found an indirect bottleneck \o/
            isIndirectBottleneck = pIsMdfHigherWithChangedDeltaG0s(numMaxExchanges, exchangeIndices, cnap,...
                                                                   RT, dG0sAndUncertainties,...
                                                                   minConcentrationsMat, maxConcentrationsMat,...
                                                                   D, d, fixedConcentrationRatioRanges,...
                                                                   reactionIndexInAllSpecies, unchangedMdf,...
                                                                   maximalMilpRunTime);
            if isIndirectBottleneck
                indirectBottleneckIndices = [indirectBottleneckIndices reactionIndexInAllSpecies];
                for single_species_reaction_id = reactionIdInAllSpecies
                    indirectBottlenecksStr = indirectBottlenecksStr + single_species_reaction_id + " & ";
                end
                indirectBottlenecksStr = indirectBottlenecksStr + "X\n";
            end
        end
        indirectBottlenecksStr = strrep(indirectBottlenecksStr, " & X", "");
    end
end


function isMdfHigher = pIsMdfHigherWithChangedDeltaG0s(numMaxExchanges, exchangeIndices, cnap, RT, dG0sAndUncertainties, minConcentrationsMat, maxConcentrationsMat, D, d, fixedConcentrationRatioRanges, changedReactionIndices, oldMdf, maximalMilpRunTime)
    % Check with higher MDF
    isMdfHigher = false;
    higherMdfWithHigherDeltaG0 = pGetHigherMdfWithChangedDeltaG0s(numMaxExchanges, exchangeIndices, cnap, RT, dG0sAndUncertainties,...
                                                                  minConcentrationsMat, maxConcentrationsMat,...
                                                                  D, d, fixedConcentrationRatioRanges,...
                                                                  changedReactionIndices, oldMdf,...
                                                                  +25, maximalMilpRunTime);
    if ~isnan(higherMdfWithHigherDeltaG0)
        isMdfHigher = true;
        return
    end

    % Check with lower MDF
    higherMdfWithLowerDeltaG0 = pGetHigherMdfWithChangedDeltaG0s(numMaxExchanges, exchangeIndices, cnap, RT, dG0sAndUncertainties,...
                                                                 minConcentrationsMat, maxConcentrationsMat,...
                                                                 D, d, fixedConcentrationRatioRanges,...
                                                                 changedReactionIndices, oldMdf,...
                                                                 -25, maximalMilpRunTime);
    if ~isnan(higherMdfWithLowerDeltaG0)
        isMdfHigher = true;
        return
    end
end

function higherMdfWithChanges = pGetHigherMdfWithChangedDeltaG0s(numMaxExchanges, exchangeIndices, cnap, RT, dG0sAndUncertainties, minConcentrationsMat, maxConcentrationsMat, D, d, fixedConcentrationRatioRanges, changedReactionIndices, oldMdf, changeValue, maximalMilpRunTime)
    dG0sAndUncertaintiesCopy = dG0sAndUncertainties;
    for current_mdf_index = changedReactionIndices
        dG0sAndUncertaintiesCopy(current_mdf_index,1) = dG0sAndUncertainties(current_mdf_index,1) + changeValue;
    end
    higherMdfWithChanges = CNAcomputeOptMDFpathway_higher_mdf_max_exchanges(maximalMilpRunTime, numMaxExchanges, exchangeIndices, oldMdf*1.01, false,...
        cnap,...
        RT,...
        dG0sAndUncertaintiesCopy,...
        minConcentrationsMat,...
        maxConcentrationsMat,...
        D,...
        d,...
        fixedConcentrationRatioRanges);
end
