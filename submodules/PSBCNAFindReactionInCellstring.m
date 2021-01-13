function index = PSBCNAFindReaction(reactionName, reactionIdsCellstring)
    % Arguments:
    % 1: reactionName: string
    % 2: model: CNAModel
    % Causes an error and prints an explanation if a reaction is not found
    % Otherwise, it returns the position of the reactionName
    % in the model's reaction matrix (i.e., in cnap.reacID).
    
    % If a string is given, convert is to char as this is the CNA model's
    % format for names
    reactionName = convertStringsToChars(reactionName);
    
    % Get a vector with 1 for every name's occurrence, otherwise 0
    reactionIndexVector = strcmp(reactionName, reactionIdsCellstring);

    if sum(reactionIndexVector) == 0
        disp("No index found for reaction ")
        disp(reactionName)
        disp(" D:\n")
        error("ERROR: See message above");
    elseif sum(reactionIndexVector) > 1
        disp("Multiple indices found for reaction ")
        disp(reactionName)
        disp(" D:\n")
        error("ERROR: See message above");
    else
        index = find(reactionIndexVector);
    end
end
