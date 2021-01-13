function [cnap] = loadAndSaveCommModelPyCommunityModel(sbmlFilePath, dG0FilePath, savePath, defaultMinConcentration, defaultMaxConcentration, nondefaultConcentrationsMat, exchangeMinConcentration, exchangeMaxConcentration)

ext_comparts = {};

% Load metabolic model SBML file
[cnap, ~] = CNAsbmlModel2MFNetwork(sbmlFilePath, ext_comparts);

% Load dG0 JSON file
text = fileread(dG0FilePath);
jsondata = jsondecode(text);
jsonReactionKeys = fieldnames(jsondata);

% Set up result matrices
numReactions = numel(cnap.reacID(:,1));
dG0s = cell(1, numReactions);
uncertainties = cell(1, numReactions);

for currentReaction = 1:numReactions
   reactionName = cnap.reacID(currentReaction,:);
   % reverse = contains(string(reactionName), "_TG_reverse");

   reactionNameTemp = "\b"+ string(reactionName);
   reactionNameAsKey = strrep(reactionNameTemp, '\bR_', '');
   reactionNameAsKey = strtrim(reactionNameAsKey);

   if numel(reactionNameAsKey{1}) >= 1
      isXToAdd = isstrprop(reactionNameAsKey, "digit");
      if isXToAdd(1)
          reactionNameAsKey = "x"+string(reactionNameAsKey);
      end
   end

   isJsonKey = 0;
   for j = 1:numel(jsonReactionKeys)
       jsonReactionKey = jsonReactionKeys{j};
       if reactionNameAsKey == jsonReactionKey
           isJsonKey = 1;
           break
       end
   end

   if isJsonKey
       thermodynamicReactionData =  eval("jsondata."+reactionNameAsKey);
       dG0 = thermodynamicReactionData.dG0;
       uncertainty = thermodynamicReactionData.uncertainty;

       dG0s{currentReaction} = dG0;
       uncertainties{currentReaction} = uncertainty;
   else
       dG0s{currentReaction} = NaN;
       uncertainties{currentReaction} = NaN;
   end
end


% METABOLITES
numMetabolites = numel(cnap.specID(:,1));
minConcentrations = cell(1, numMetabolites);
maxConcentrations = cell(1, numMetabolites);
% Set default concentration ranges
for currentMetabolite = 1:length(minConcentrations)
    metaboliteId = strtrim(convertCharsToStrings(cnap.specID(currentMetabolite, :)));
    metaboliteIdSplit = strsplit(metaboliteId, "_");
    metaboliteSpecies = metaboliteIdSplit(end);
    
    if metaboliteSpecies == "exchg"
        minConcentrations{currentMetabolite} = exchangeMinConcentration;
        maxConcentrations{currentMetabolite} = exchangeMaxConcentration;
    else
        minConcentrations{currentMetabolite} = defaultMinConcentration;
        maxConcentrations{currentMetabolite} = defaultMaxConcentration;
    end
end
% Set non-default concentration ranges
nondefaultMatSize = size(nondefaultConcentrationsMat);
numNondefaultConcentrations = nondefaultMatSize(1);
for currentNondefault = 1:numNondefaultConcentrations
    metaboliteIndex = PSBCNAFindMetabolite(nondefaultConcentrationsMat(currentNondefault, 1), cnap);
    nondefaultMinConcentration = str2double(nondefaultConcentrationsMat(currentNondefault, 2));
    nondefaultMaxConcentration = str2double(nondefaultConcentrationsMat(currentNondefault, 3));
    minConcentrations{metaboliteIndex} = nondefaultMinConcentration;
    maxConcentrations{metaboliteIndex} = nondefaultMaxConcentration;
end


% ADD GENERIC DATA
[cnap, ~] = CNAsetGenericReactionData_with_array(cnap, 'dG0', dG0s);
[cnap, ~] = CNAsetGenericReactionData_with_array(cnap, 'uncertainty', uncertainties);
[cnap, ~] = CNAsetGenericSpeciesData_with_array(cnap, 'minConcentration', minConcentrations);
[cnap, ~] = CNAsetGenericSpeciesData_with_array(cnap, 'maxConcentration', maxConcentrations);

save(savePath, 'cnap')
end
