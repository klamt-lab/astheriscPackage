function [cMins, cMaxs, activeMetaboliteIndices, errorConcentrationRanges] = getConcentrationRanges(numMaxExchanges, exchangeIndices, minMdf, cnap, D, d, dG0sAndUncertainties, RT, minConcentrationsMat, maxConcentrationsMat, fixedConcentrationRatioRanges, activeMetaboliteIndices, maximalMilpRunTime)
    % This OptMDFpathway version (where minimal and maximal MDF constraints are calculated under
    % the given minimal MDF constraint) is called OptMDF5 in ASTHERISC's publication.

    % --> Get CPLEX obj
    [~, ~, ~, ~, ~, ~, ~, ~, ~, obj, ~, ~, cVars] = max_min_driving_force_pathway_higher_mdf_max_exchanges(maximalMilpRunTime, numMaxExchanges, exchangeIndices, minMdf - abs(minMdf)*.0001, false, cnap.stoichMat,...
        cnap.stoichMat(cnap.specInternal, :), cnap.reacMin, cnap.reacMax, D, d,...
        max(1e-9, cnap.epsilon), [], 1000, RT, dG0sAndUncertainties, minConcentrationsMat, maxConcentrationsMat, fixedConcentrationRatioRanges,...
        isnan(dG0sAndUncertainties(:, 1)), false, cnap.reacMin, cnap.reacMax);

    % --> Get metabolite concentration ranges
    metabolitesSize = size(cnap.specID);
    numMetabolites = metabolitesSize(1);
    cMins = zeros(1, numMetabolites);
    cMaxs = zeros(1, numMetabolites);
    obj.cpx.setOut([]);
    errorConcentrationRanges = false;
    for currentActiveMetaboliteIndex = 1:length(activeMetaboliteIndices)
        activeMetaboliteIndex = activeMetaboliteIndices(currentActiveMetaboliteIndex);
        obj.cpx.getObjective().setExpr(cVars(activeMetaboliteIndex));
        obj.cpx.solve();
        try
            cMin = obj.cpx.getObjValue();
        catch
            errorConcentrationRanges = true;
            break;
        end
        obj.cpx.getObjective().setExpr(obj.cpx.prod(-1, cVars(activeMetaboliteIndex)));
        obj.cpx.solve();
        try
            cMax = -obj.cpx.getObjValue();
        catch
            errorConcentrationRanges = true;
            break;
        end

        cMins(activeMetaboliteIndex) = cMin;
        cMaxs(activeMetaboliteIndex) = cMax;
    end
end
