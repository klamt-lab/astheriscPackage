function [RT, dG0sAndUncertainties, minConcentrationsMat, maxConcentrationsMat, fixedRatios] = getGeneralMdfMilpParameters(cnap, temperature, fixedRatiosNameMat)
    %% Set RT
    R = 8.314e-3; % kJ⋅K⁻1⋅mol⁻1 (standard value is in J⋅K⁻1⋅mol⁻1)
    RT = R * temperature;
    
    %% Set dG0s and uncertainties
    [dG0s, errval] = CNAgetGenericReactionData_as_array(cnap, 'dG0');
    if errval
        error("ERROR: dG0s could not be retrieved");
    end
    [uncertainties, errval] = CNAgetGenericReactionData_as_array(cnap, 'uncertainty');
    if errval
        error("ERROR: dG0 uncertainties could not be retrieved");
    end

    dG0sAndUncertainties = NaN(length(dG0s), 1);
    for i = 1:length(dG0s)
        if uncertainties{i} < 100 %%%
            dG0sAndUncertainties(i,1) = dG0s{i};
        end
        % dG0sAndUncertainties(i,2) = 0; %uncertainties{i};
    end

    %% Set metabolite concentration ranges
    [maxConcentrations, errval] = CNAgetGenericSpeciesData_as_array(cnap, 'maxConcentration');
    if errval
        error("ERROR: Maximal metabolite concentrations could not be retrieved");
    end
    [minConcentrations, errval] = CNAgetGenericSpeciesData_as_array(cnap, 'minConcentration');
    if errval
        error("ERROR: Minimal metabolite concentrations could not be retrieved");
    end
    minConcentrationsMat = [];
    maxConcentrationsMat = [];
    for i = 1:length(minConcentrations)
        minConcentrationsMat = [minConcentrationsMat minConcentrations{i}];
        maxConcentrationsMat = [maxConcentrationsMat maxConcentrations{i}];
    end

    %% Add fixed metabolite ratios    
    fixedRatiosSize = size(fixedRatiosNameMat);
    numFixedRatios = fixedRatiosSize(1);
    fixedRatios = [];
    for currentFixedRatio = 1:numFixedRatios
        singleFixedRatio = [];
        singleFixedRatio = [singleFixedRatio PSBCNAFindMetabolite(fixedRatiosNameMat(currentFixedRatio, 1), cnap)];
        singleFixedRatio = [singleFixedRatio PSBCNAFindMetabolite(fixedRatiosNameMat(currentFixedRatio, 2), cnap)];
        singleFixedRatio = [singleFixedRatio str2double(fixedRatiosNameMat(currentFixedRatio, 3))];
        singleFixedRatio = [singleFixedRatio str2double(fixedRatiosNameMat(currentFixedRatio, 4))];
        fixedRatios = [fixedRatios; singleFixedRatio];
    end
end
