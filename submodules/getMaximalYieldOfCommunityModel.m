function [maxyield, success] = getMaximalYieldOfCommunityModel(cnap, target_substrate_reaction, target_max_substrate_uptake,...
                                             target_product_reaction, community_exchange_reactions, species_ids)
    %% Set yield function (c*r/d*r)
    c = zeros([1, cnap.numr]);
    c(1, PSBCNAFindReaction(target_product_reaction, cnap)) = 1;
    d = zeros([1, cnap.numr]);
    d(1, PSBCNAFindReaction(target_substrate_reaction, cnap)) = -1; % Substrate uptake has a negative flux :O

    %% Set maximal uptake
    cnap.reacMin(PSBCNAFindReaction(target_substrate_reaction, cnap)) = -target_max_substrate_uptake * numel(species_ids);
    cnap.reacMax(PSBCNAFindReaction(target_substrate_reaction, cnap)) = 0;

    %% Open product reaction
    cnap.reacMin(PSBCNAFindReaction(target_product_reaction, cnap)) = 0;
    cnap.reacMax(PSBCNAFindReaction(target_product_reaction, cnap)) = inf;

    %% Set exchange reactions
    for i = 1:numel(community_exchange_reactions)
        exchange_reaction_id = community_exchange_reactions(i);
        cnap.reacMin(PSBCNAFindReaction(exchange_reaction_id, cnap)) = -inf;
        cnap.reacMax(PSBCNAFindReaction(exchange_reaction_id, cnap)) = inf;
    end

    %% Perform yield optimization
    fixedFluxes = [];
    c_macro = cnap.macroDefault;
    solver = 2;
    [maxyield, ~, success, ~] = CNAoptimizeYield(cnap, c, d, fixedFluxes, c_macro, solver);
end
