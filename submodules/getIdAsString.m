function stringifiedId = getIdAsString(id)
    stringifiedId = convertCharsToStrings(id);
    stringifiedId = strrep(stringifiedId, " ", "");
end
