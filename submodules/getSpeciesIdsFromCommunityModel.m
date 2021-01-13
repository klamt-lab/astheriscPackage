function [species_ids] = getSpeciesIdsFromCommunityModel(cnap)
    species_ids = [];
    for species_index = 1:numel(cnap.specID(:,1))
        species_id = convertCharsToStrings(cnap.specID(species_index,:));
        if ~startsWith(species_id, "M_community_biomass_")
            continue
        end
        split_species_id = strsplit(species_id, "_");
        species_id_with_spaces = split_species_id(numel(split_species_id));
        species_id_with_spaces_split = strsplit(species_id_with_spaces, " ");
        final_species_id = species_id_with_spaces_split(1);
        species_ids = [species_ids final_species_id];
    end
end
