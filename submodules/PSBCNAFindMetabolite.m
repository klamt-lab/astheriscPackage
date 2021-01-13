function index = PSBCNAFindMetabolite(metaboliteName, model)
    % Arguments:
    % 1: metaboliteName: string
    % 2: model: CNAModel
    % Causes an error and prints an explanation if a reaction is not found
    % Otherwise, it returns the position of the reactionName
    % in the model's metabolite matrix (i.e., cnap.specID).
    
    % If a string is given, convert is to char as this is the CNA model's
    % format for names
    metaboliteName = convertStringsToChars(metaboliteName);
    
    % Get a vector with 1 for every name's occurrence, otherwise 0
    metaboliteIndexVector = strcmp(metaboliteName, cellstr(model.specID));

    if sum(metaboliteIndexVector) == 0
        disp("No index found for metabolite ")
        disp(metaboliteName)
        disp(" D:\n")
        error("ERROR: See message above");
    elseif sum(metaboliteIndexVector) > 1
        disp("Multiple indices found for metabolite ")
        disp(metaboliteName)
        disp(" D:\n")
        error("ERROR: See message above");
    else
        index = find(metaboliteIndexVector);
    end
end
