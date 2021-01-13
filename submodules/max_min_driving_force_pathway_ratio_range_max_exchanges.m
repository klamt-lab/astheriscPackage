function [mmdf, fv, x, min_df, max_df, dfs, fva_lb, fva_ub, reac_map, obj, z_vars, flux_vars, c_vars, first_z, first_fv, deltaGf]= max_min_driving_force_pathway_ratio_range_max_exchanges(maximalMilpRunTime, numMaxExchanges, exchange_reaction_indices, stoichMat, st, reacMin,...
  reacMax, des, db, lp_tol, min_active_flux, default_flux_limit, RT, G0, c_min, c_max, fixed_ratios,...
  ignore, thermo_fva, fva_lb, fva_ub, max_reac_with_mdf, shortest_path_mdf, max_reac_dfs, setup_only, obj)

% stoichMat: includes external metabolites, may contain modifications for
% special metabolites, e.g. only zeros for H2O
% st: only internal metabolites, no modifications of the coefficients | st * fv = 0
% reacMin, reacMax: lower/upper bounds for the fluxes, not used when fva_lb, fva_ub are provided
% des, db: inhomogeneous constraints | des * fv <= db
% lp_tol: value for optimality and feasibility tolerance
% min_active_flux: if not empty used as threshold in indicator constrains; if empty bigM constraints are used
% default_flux_limit: if not empty used as limit for unbounded fluxes
% RT: R * T
% G0: free Gibbs energy of reactions or metabolites (latter if cell array)
% c_min, c_max: concentration limits for the metabolites; order as in stoichMat
% fixed ratios:
% ignore: indices of reactions or logical vector; for the marked reactions no thermodynamic constraints are added
% thermo_fva: if true run thermodynamic FVA with MDF >= 0 for all reactions;
% when cell array: thermo_fva{1} is the MDF lower limit, thermo_fva{2} a row vector of reaction indices for which to run the FVA
% fva_lb, fva_ub: used directly as lower/upper bounds for the fluxes; no
% structural FVA will be performed; when fva_lb = fva_ub then only the MDF for the path
% described by this flux distribution is calculated
% shortest_path_mdf: find the shortest path with MDF >= shortest_path_mdf
% max_reac_with_mdf: maximize the flux through reaction max_reac_with_mdf(1) with MDF >= max_reac_with_mdf(2)
% max_reac_dfs: iteratively maximize the driving force of each reaction with MDF >= max_reac_dfs(1);
% if max_reac_dfs is a vector of length 2 then max_reac_dfs(2) then the driving force will not be
% maximized all the way if it was already found to be >= max_reac_dfs(2)
% setup_only: if true set up the MILP only but do not solve
% obj: optional CplexJava object to use for setting up the MILP
%
% x: vector of logarithmized concentrations when G0 are reaction deltaG or vector of deltaGf when G0 are metabolite deltaG

% This file is part of CellNetAnalyzer. Please visit
% http://www.mpi-magdeburg.mpg.de/projects/cna/cna.html
% for more information and the latest version of CellNetAnalyzer.
%
% Copyright (C) 2000-2019 by Steffen Klamt and Axel von Kamp,
% Max Planck Institute for Dynamics of Complex Technical Systems, Magdeburg, Germany.
%
% Contributors are listed in CONTRIBUTORS.txt.
%
% This software can be used under the terms of our CellNetAnalyzer License.
% A copy of the license agreement is provided in the file named "LICENSE.txt"
% included with this software distribution. The license is also available online at
% http://www2.mpi-magdeburg.mpg.de/projects/cna/license.html
%
% For questions please contact: cellnetanalyzer@mpi-magdeburg.mpg.de

if nargin < 25+1
  obj= [];
  if nargin < 24+1
    setup_only= false;
    if nargin < 23+1
      max_reac_dfs= NaN;
      if nargin < 22+1
        shortest_path_mdf= []; % find shortest path with at least this MDF
        if nargin < 21+1
          max_reac_with_mdf= []; % maximize reaction (1) with at at least this (2) MDF
          if nargin < 18+1
            thermo_fva= false;
          end
        end
      end
    end
  end
end

use_reaction_uncertainty= ~iscell(G0) && size(G0, 2) == 2;

[num_internal_met, num_reac]= size(st);
mmdf= NaN;
fv= NaN(num_reac, 1);
x= NaN(1, num_internal_met);
min_df= fv;
max_df= fv;
dfs= [];
reac_map= [];

if isempty(min_active_flux) && isempty(default_flux_limit)
  error('You have to set the ''min_active_flux'' or ''default_flux_limit'' parameter.');
end

if isempty(des)
  des= zeros(0, num_reac);
end

c_min= log(c_min(:));
c_max= log(c_max(:));

tmp= false(num_reac, 1);
tmp(ignore)= true;
ignore= tmp | ~any(stoichMat, 1)'; % also ignore empty reaction columns

cplex_inner= setup_cplex_inner_class_access();
fvaPars= cplex_inner.ParameterSet.constructor.newInstance([]);
fvaPars.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
fvaPars.setParam(cplex_inner.DoubleParam.EpOpt, lp_tol);
fvaPars.setParam(cplex_inner.DoubleParam.EpRHS, lp_tol);
fvaPars.setParam(cplex_inner.IntParam.RootAlg, 2);
num_met_zv= zeros(num_internal_met, 1);
if iscell(db)
  lhs= [num_met_zv; db{1}];
  rhs= [num_met_zv; db{2}];

else
  lhs= [num_met_zv; -Inf(length(db), 1)];
  rhs= [num_met_zv; db];
end
if nargin < 18 || isempty(fva_lb) || isempty(fva_ub)
  disp('Running FVA without thermodynamic constraints.');
  [ind_i, ind_j, val]= find([st; des]);
  fvaj= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, reacMin, reacMax, fvaPars);
  clear ind_i ind_j val;
  if fvaj(3)
    fva_lb= [];
    fva_ub= [];
    disp('FVA infeasible; no flux possible in this network.');
    return;
  else
    fva_lb= fvaj(1);
    fva_ub= fvaj(2);
    fva_lb(abs(fva_lb) < lp_tol)= 0;
    fva_ub(abs(fva_ub) < lp_tol)= 0;
  end
end
flux_lb= fva_lb;
flux_ub= fva_ub;
if ~isempty(default_flux_limit)
  flux_lb(isinf(flux_lb))= -default_flux_limit;
  flux_ub(isinf(flux_ub))= default_flux_limit;
end

if isequal(flux_lb, flux_ub) && ~iscell(thermo_fva) && ~thermo_fva && isnan(max_reac_dfs)
  mmdf_active_reactions_only= true;
  ignore= ignore | flux_lb == 0;
else
  mmdf_active_reactions_only= false;
end

% split reversible reactions so that all flux_lb >= 0 (could only split effectively reversible reactions)
to_split= find(flux_lb < -lp_tol);
if nargout > 8
  reac_map= [(1:num_reac)'; -to_split];
end
to_split= int32(to_split);
% even ignored reactions are currently split because of association with z_vars
st= [st, -st(:, to_split)];
stoichMat= [stoichMat, -stoichMat(:, to_split)];
des= [des, -des(:, to_split)];
n= size(st, 2);
flux_lb(end+1:n)= max(0, -flux_ub(to_split));
flux_ub(end+1:n)= -flux_lb(to_split);
flux_lb(to_split)= 0;
split_off= flux_ub < flux_lb;
flux_lb(split_off)= 0;
flux_ub(split_off)= 0;
if ~iscell(G0)
  if use_reaction_uncertainty
    G0(end+1:n, :)= [-G0(to_split, 1), G0(to_split, 2)];
  else
    G0(end+1:n, :)= -G0(to_split, :);
  end
end
ignore(end+1:n)= ignore(to_split);

blocked= (flux_lb >= -lp_tol & flux_ub <= lp_tol)';
essential= (flux_lb > lp_tol)';
variable= ~blocked & ~essential;
essential= find(essential);
if ~use_reaction_uncertainty
  % would mess with split u_vars
  ignore(blocked)= true;
end

% for the purpose of finding min_df/max_df to calculate bigM the ratio
% constraints can be ignored; this might miss some essential reactions with
% driving force < 0 which is unproblematic here
[min_df, max_df]= driving_forces(stoichMat, RT, G0, c_min, c_max, ignore, [], lp_tol, to_split);
min_min_df= min(min_df);
if isnan(min_min_df)
  warning('No thermodynamic constraints are active in this system.');
  return;
end
max_max_df= max(max_df);
fprintf('Minimal driving force: %.2f, maximal driving force: %.2f.\n', min_min_df, max_max_df);
ind= essential(max_df(essential) < -lp_tol);
if ~isempty(ind)
  fprintf('Essential reactions with maximal driving force < 0 found:\n');
  fprintf('%d ', ind);
  fprintf('\nNo solution with max min driving force > 0 possible.\n');
end

if isempty(obj)
  obj= CplexJava();
end
[obj, c_vars, driving_force, u_vars, deltaGf]= setup_driving_force_constraints_ratio_ranges(stoichMat, RT, G0, c_min, c_max,...
  ignore, fixed_ratios, to_split, use_reaction_uncertainty, obj);

if mmdf_active_reactions_only
  z_vars= [];
  flux_vars= [];
  first_z= NaN;
  first_fv= NaN;
else
  sel= flux_lb ~= flux_ub;
  flux_lb(sel)= flux_lb(sel).*(1 - sign(flux_lb(sel))*2*lp_tol); % give some wiggle room
  flux_ub(sel)= flux_ub(sel).*(1 + sign(flux_ub(sel))*2*lp_tol);
  flux_vars= obj.cpx.numVarArray(n, flux_lb, flux_ub);
  first_fv= obj.mipmat.addCols(flux_vars);
  obj.mipmat.addRows(lhs, rhs, [], []);
  [i, j, val]= find([st; des]);
  obj.mipmat.setNZs(int32(i - 1), first_fv + int32(j - 1), val); %A Java indices
  clear i j val

  z_vars= obj.cpx.boolVarArray(n); % z_vars of essential/blocked reactions can be fixed to 1/0
  for i= find(blocked)
    z_vars(i).setUB(0);
  end
  for i= essential
    z_vars(i).setLB(1);
  end

  first_z= int32(obj.mipmat.getNcols());
  m= int32(obj.mipmat.getNrows());
  obj.mipmat.addCols(z_vars);
  if isempty(min_active_flux)
    % this variant allows zero flux even if the corresponding z_var is true
    % which is OK in conjunction with the maximization of mmdf_var
    ind= find(variable);
    if ~isempty(ind)
      ind= int32(ind);
      r_ind= m:length(ind)+m-1;
      obj.mipmat.addRows(zeros(length(ind), 1), Inf(length(ind), 1), [], []);
      obj.mipmat.setNZs(r_ind, first_fv + ind - 1, -ones(1, length(ind)));
      obj.mipmat.setNZs(r_ind, first_z + ind - 1, flux_ub(ind));
    end
  end

  % only allow one direction of a reversible reaction to be active
  m= int32(obj.mipmat.getNrows());
  if length(to_split) > 0 %#ok for older Matlab versions isempty does not work on int32
    r_ind= m:length(to_split)+m-1;
    obj.mipmat.addRows(zeros(length(to_split), 1), ones(length(to_split), 1), [], []);
    obj.mipmat.setNZs(r_ind, first_z + to_split - 1, ones(length(to_split), 1));
    obj.mipmat.setNZs(r_ind, first_z+int32(num_reac):first_z+int32(n-1), ones(length(to_split), 1));
  end
end

if iscell(thermo_fva) || thermo_fva
  if iscell(thermo_fva)
    mmdf= thermo_fva{1};
    if length(thermo_fva) > 1
      [a, b]= ismember(thermo_fva{2}, to_split); % thermo_fva{2} must be a row vector
      fva_ind= [thermo_fva{2}, num_reac + b(a)];
    else
      fva_ind= 1:n;
    end
    thermo_fva= true;
  else
    mmdf= 0;
    fva_ind= 1:n;
  end
  fva_ind= setdiff(fva_ind, find(blocked));
  mmdf_var= obj.cpx.constant(mmdf);
  obj.cpx.addMinimize(flux_vars(1));
  minimize= obj.cpx.getObjective().getSense();
  obj.cpx.getObjective().clearExpr();
  maximize= obj.cpx.maximize(flux_vars(1)).getSense();
elseif ~isnan(max_reac_dfs)
  mmdf_var= obj.cpx.constant(max_reac_dfs(1));
  obj.cpx.addMaximize(flux_vars(1)); % only to initiate the 'objective' variable
  objective= obj.cpx.getObjective();
else
  if isempty(max_reac_with_mdf) && isempty(shortest_path_mdf)
    mmdf_var= obj.cpx.numVar(min_min_df - lp_tol, max_max_df + lp_tol, 'MDF');
    obj.cpx.addMaximize(mmdf_var);
  elseif ~isempty(max_reac_with_mdf)
    mmdf_var= obj.cpx.numVar(max_reac_with_mdf(2), max_max_df + lp_tol, 'MDF');
    obj.cpx.addMaximize(flux_vars(max_reac_with_mdf(1))); % must be irreversible forward
  elseif ~isempty(shortest_path_mdf)
    mmdf_var= obj.cpx.numVar(shortest_path_mdf, max_max_df + lp_tol, 'MDF');
    obj.cpx.addMinimize(obj.cpx.sum(z_vars));
  else
    error('Invalid parameter configuration');
  end
  obj.mipmat.addColumn(mmdf_var);
end

obj.cpx.setParam(cplex_inner.DoubleParam.EpInt, 1e-10); % AvK/PSB change: from lp_tol to 1e-10
obj.cpx.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
obj.cpx.setParam(cplex_inner.DoubleParam.EpOpt, lp_tol);
obj.cpx.setParam(cplex_inner.DoubleParam.EpRHS, lp_tol);


%%% PSB CHANGE START
[~, exchange_reaction_indices]= ismember([exchange_reaction_indices, -exchange_reaction_indices], reac_map);
exchange_reaction_indices = exchange_reaction_indices(exchange_reaction_indices~=0);
if numMaxExchanges ~= inf
    obj.cpx.addLe((obj.cpx.sum(z_vars(exchange_reaction_indices))), numMaxExchanges);
end
%%% PSB CHANGE END

if ~mmdf_active_reactions_only
  bigM= max_max_df - min_df; % reaction-specific bigM
  for i= 1:n
    if ~isempty(min_active_flux) && variable(i)
      obj.cpx.add(obj.cpx.ifThen(obj.cpx.le(z_vars(i), int32(0)), obj.cpx.le(flux_vars(i), 0)));
      obj.cpx.add(obj.cpx.ifThen(obj.cpx.ge(z_vars(i), int32(1)), obj.cpx.ge(flux_vars(i), min_active_flux)));
    end

    if ~ignore(i)
      obj.cpx.addGe(obj.cpx.sum(obj.cpx.prod(bigM(i), obj.cpx.diff(1, z_vars(i))), driving_force(i)), mmdf_var); % ranges(count)=
    end
  end
else
  for i= 1:n
    if ~ignore(i)
      if flux_lb(i) ~= 0 && ~isempty(driving_force(i))
        obj.cpx.addGe(driving_force(i), mmdf_var);
      elseif use_reaction_uncertainty
        obj.cpx.addGe(driving_force(i), -Inf); % kludge to include associated u_var
      end
    end
  end
end

if setup_only
  return;
end

if thermo_fva
  fprintf('Running FVA for %d reactions with thermodynamic constraints.\n', length(fva_ind));
  err_msg= 'No flux feasible under given thermodynamic constraints.';
  flux_lb(:)= NaN; %x
  flux_ub(:)= NaN; %v
  flux_lb(blocked)= 0;
  flux_ub(blocked)= 0;

% run tFVA externally
%   fname= tempname();
%   mipfile= [fname, '.sav'];
%   prmfile= [fname, '.prm'];
%   obj.cpx.exportModel(mipfile);
%   obj.cpx.writeParam(prmfile);
%   obj.cpx.clearModel();
% %   fvaj= CplexFVA.sav_fva(mipfile, prmfile, first_fv:(first_fv+length(flux_vars)-1));
%   fvaj= CplexFVA.sav_fva(mipfile, prmfile, first_fv + fva_ind - 1);
%   if fvaj(3)
%     error(err_msg);
%   end
%   tmp= fvaj(1);
%   flux_lb(fva_ind)= tmp(first_fv + fva_ind);
%   tmp= fvaj(2);
%   flux_ub(fva_ind)= tmp(first_fv + fva_ind);
% %   cgp= Cplex();
% %   cgp.readModel(mipfile);
% %   flux_lb= cgp.Model.lb;
% %   flux_ub= cgp.Model.ub;
%   delete(mipfile);
%   delete(prmfile);

  obj.cpx.setOut([]);
  obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, 10); % better make this an explicit thermo_fva parameter
  objective= obj.cpx.getObjective();
  fprintf('1');
  for i= fva_ind
    objective.setExpr(flux_vars(i));
    if flux_vars(i).getLB() == 0
      flux_lb(i)= 0;
    else
      objective.setSense(minimize); % !! leaves status Optimal, only clears solution !!
      while true
        obj.cpx.solve();
        if ~obj.cpx.getStatus().equals(cplex_inner.Status.Unknown)
          break;
        end
      end
      if obj.cpx.getStatus().equals(cplex_inner.Status.Optimal)
        flux_lb(i)= obj.cpx.getObjValue();
      elseif obj.cpx.getStatus().equals(cplex_inner.Status.Feasible)
        fprintf('\nFeasible only for LB %d.\n', i);
        flux_lb(i)= obj.cpx.getBestObjValue();
      else
        error(err_msg);
      end
    end
    objective.setSense(maximize); % !! leaves status Optimal, only clears solution !!
    while true
      obj.cpx.solve();
      if ~obj.cpx.getStatus().equals(cplex_inner.Status.Unknown)
        break;
      end
    end
    if obj.cpx.getStatus().equals(cplex_inner.Status.Optimal)
      flux_ub(i)= obj.cpx.getObjValue();
    elseif obj.cpx.getStatus().equals(cplex_inner.Status.Feasible)
      fprintf('\nFeasible only for UB %d.\n', i);
      flux_ub(i)= obj.cpx.getBestObjValue();
    else
      error(err_msg);
    end
    if mod(i, 10) == 0
      fprintf(' %d', i);
      save tFVAtmp.mat flux_lb flux_ub reac_map
      if mod(i, 100) == 0
        fprintf('\n');
      end
    end
  end
  fprintf('\n');
  save tFVAtmp.mat flux_lb flux_ub reac_map

  if any(flux_lb(to_split) > lp_tol & flux_lb(num_reac+1:n) > lp_tol)
    warning('Inconsistent bounds for split reversible reactions.');
  end
  flux_lb(to_split)= -flux_ub(num_reac+1:n);
  ind= find(flux_ub(to_split) < lp_tol);
  flux_ub(to_split(ind))= -flux_lb(num_reac+ind);
  fv= flux_lb(1:num_reac);
  x= flux_ub(1:num_reac);
elseif isnan(max_reac_dfs(1))
  if isempty(shortest_path_mdf)
    obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, maximalMilpRunTime); % PSB; prevent hang-ups, better make this an explicit parameter
  else
    obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, maximalMilpRunTime); % PSB; prevent hang-ups, better make this an explicit parameter
  end
  obj.cpx.setParam(cplex_inner.DoubleParam.SolnPoolGap, lp_tol);
  obj.solve();
  status= obj.getStatus();
  if status.equals(cplex_inner.Status.Optimal) || status.equals(cplex_inner.Status.Feasible)
    num_sol= 1; %obj.cpx.getSolnPoolNsolns();
    fv= zeros(n, num_sol);
    for i= 1:num_sol
    if mmdf_active_reactions_only
      fv= flux_lb;
      zv= logical(fv);
    else
      fv(:, i)= obj.getValues(flux_vars);
      zv= logical(round(obj.getValues(z_vars)));
      if ~isempty(min_active_flux) && (~all(fv(zv~=0)>0) || ~all(fv(zv==0)==0))
        warning('CAUTION: Could not properly distinguish between active and inactive fluxes;\nwiden the gap between ''lp_tol'' and ''min_active_flux''.\n');
        return;
      end
    end
    mmdf= obj.getValue(mmdf_var);
    [dfs, x]= calc_dfs(obj, num_reac, n, c_vars, u_vars, G0, RT, stoichMat, ignore, to_split, deltaGf);
    if mmdf >= max_max_df - lp_tol && all(isnan(dfs(zv)));
      warning('The solution does not involve any reactions with thermodynamic constraints.');
      mmdf= NaN;
    end
    rev_act= find(zv(num_reac+1:n) & ~ignore(num_reac+1:n));
    dfs(to_split(rev_act))= dfs(num_reac+rev_act);
    dfs= dfs(1:num_reac);
    fv(to_split, i)= fv(to_split, i) - fv(num_reac+1:n, i); % merge split reversible reactions
    end
    fv= fv(1:num_reac, :);
  end
  if ~status.equals(cplex_inner.Status.Optimal)
    warning('MILP terminated with status %s.\n', char(status.toString()));
    mmdf= NaN; % can still retrieve current best solution from dfs
  end
else
  dfs= NaN(n, 1);
  obj.cpx.setOut([]);

%   for i= 1:n
%     if ~ignore(i) && (length(max_reac_dfs) == 1 || ~(dfs(i) >= max_reac_dfs(2))) % latter is true also when isnan(dfs(i))
%       if length(max_reac_dfs) > 1
%         ti_lim= 0.1;
%         obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, ti_lim);
%       end
%       objective.setExpr(driving_force(i));
%       while (1)
%         obj.cpx.solve();
%         status= obj.cpx.getStatus();
%         if status.equals(cplex_inner.Status.Optimal) ||...
%             (length(max_reac_dfs) > 1 && status.equals(cplex_inner.Status.Feasible) && obj.cpx.getObjValue() >= max_reac_dfs(2))
%           current_dfs= calc_dfs(obj, lp_tol, num_reac, n, obj.cpx.getValues(c_vars)',...
%             u_vars, G0, RT, stoichMat, ignore, to_split, 0);
%           sel= ~ignore & ~(current_dfs <= dfs);
%           dfs(sel)= current_dfs(sel);
%           break; % while
%         elseif status.equals(cplex_inner.Status.Infeasible)
%           error('Unexpected solver status %s.\n', char(status.toString()));
%         else
% %           disp(char(status.toString()));
%           ti_lim= 2*ti_lim;
%           obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, ti_lim);
%         end
%       end
%     end
%   end

  % here the second solve() has to finish with status 'Optimal'
  for i= 1:n
    if ~ignore(i) && (length(max_reac_dfs) == 1 || ~(dfs(i) >= max_reac_dfs(2))) % latter is true also when isnan(dfs(i))
      objective.setExpr(driving_force(i));
      if length(max_reac_dfs) > 1
        obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, 2);
        obj.cpx.solve();
        status= obj.cpx.getStatus();
        if status.equals(cplex_inner.Status.Optimal) ||...
            (status.equals(cplex_inner.Status.Feasible) && obj.cpx.getObjValue() >= max_reac_dfs(2))
          current_dfs= calc_dfs(obj, num_reac, n, c_vars,...
            u_vars, G0, RT, stoichMat, ignore, to_split, deltaGf);
          sel= ~ignore & ~(current_dfs <= dfs);
          dfs(sel)= current_dfs(sel);
          continue;
        end
        obj.cpx.setParam(cplex_inner.DoubleParam.TiLim, 1e75);
      end
      obj.cpx.solve();
      status= obj.cpx.getStatus();
      if ~status.equals(cplex_inner.Status.Optimal)
        error('Unexpected solver status %s.\n', char(status.toString()));
      end
      dfs(i)= obj.cpx.getObjValue();
      if length(max_reac_dfs) > 1
        current_dfs= calc_dfs(obj, num_reac, n, c_vars,...
          u_vars, G0, RT, stoichMat, ignore, to_split, deltaGf);
        sel= ~ignore & ~(current_dfs <= dfs);
        dfs(sel)= current_dfs(sel);
      end
    end
  end

end

end % function

function [dfs, x]= calc_dfs(obj, num_reac, n, c_vars, u_vars, G0, RT, stoichMat, ignore, to_split, deltaGf)
  if iscell(G0)
    x= obj.cpx.getValues(deltaGf);
    dfs= -(stoichMat'*x);
  else
    x= obj.cpx.getValues(c_vars)';
    rev_fwd_ind= to_split(~ignore(to_split));
    rev_bck_ind= find(~ignore(num_reac+1:n))' + num_reac;
    if isempty(u_vars)
      dfs= -G0 - RT*(x*stoichMat)';
    else
      tmp= obj.getValues(u_vars);
      uv= zeros(n, 1);
      uv(~ignore(1:num_reac))= tmp;
      uv(rev_bck_ind)= -uv(rev_fwd_ind);
      dfs= -G0(:, 1) + uv.*G0(:, 2) - RT*(x*stoichMat)';
    end
  end
end
