function [obj, c_vars, driving_force, u_vars1, deltaGf]= setup_driving_force_constraints_ratio_ranges(stoichMat, RT,...
  G0, log_cMin, log_cMax, ignore, fixed_ratios, to_split, use_reaction_uncertainty, obj)
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

if nargin < 10
  obj= CplexJava();
end

[m, n]= size(stoichMat);
numsm= nums(m);
c_vars= obj.cpx.numVarArray(m, log_cMin, log_cMax, strcat('C', numsm));
obj.mipmat.addCols(c_vars);
for i= 1:size(fixed_ratios, 1)
    %%% <PSB: metabolite ratio *ranges* :D>
    obj.cpx.add(obj.cpx.ge(obj.cpx.sum(c_vars(fixed_ratios(i, 1)), obj.cpx.prod(-1, c_vars(fixed_ratios(i, 2)))),...
                           log(fixed_ratios(i, 3))));
    obj.cpx.add(obj.cpx.le(obj.cpx.sum(c_vars(fixed_ratios(i, 1)), obj.cpx.prod(-1, c_vars(fixed_ratios(i, 2)))),...
                           log(fixed_ratios(i, 4))));
    %%% </PSB>
end

if use_reaction_uncertainty
  num_reac= n - length(to_split);
  u_vars_ind= find(~ignore(1:num_reac));
  u_vars1= obj.cpx.numVarArray(length(u_vars_ind), -1, 1);
  u_vars= javaArray('ilog.concert.IloNumExpr', n);
  u_vars(u_vars_ind)= u_vars1;
  if ~isempty(to_split)
    rev_ind= find(~ignore(num_reac+1:n))' + num_reac;
    if ~isempty(rev_ind)
      u_vars(rev_ind)= u_vars(to_split(~ignore(to_split)));
      for i= rev_ind
        u_vars(i)= obj.cpx.negative(u_vars(i));
      end
    end
  end
else
  u_vars1= [];
end

driving_force= javaArray('ilog.concert.IloNumExpr', n);
if iscell(G0)
  Gf0= G0{1};
  if length(G0) > 1
    default_uncertainty= G0{2};
  else
    default_uncertainty= 10000;
  end
  if size(Gf0, 2) == 1
    Gf0(end, 2)= 0; % add second column
  end
  % metabolites with NaN for deltaGf0 are set up with a deltaGf0 of 0 and default_uncertainty
  sel= isnan(Gf0(:, 1));
  Gf0(sel, 1)= 0;
  Gf0(sel, 2)= default_uncertainty;
  lb= Gf0(:, 1) - Gf0(:, 2) + RT*log_cMin;
  ub= Gf0(:, 1) + Gf0(:, 2) + RT*log_cMax;
  Gf0(sel, 1)= NaN; % explicit constraints not needed for these, already implicit in lb and ub

  deltaGf= obj.cpx.numVarArray(m, lb, ub, strcat('Gf', numsm));
  for i= 1:m
    if isempty(fixed_ratios) || isnan(Gf0(i, 1))
      obj.cpx.addGe(deltaGf(i), deltaGf(i).getLB()); % kludge to make all deltaGf accessible
    else
      % force this branch to explicitly include the c_vars
      if Gf0(i, 2) == 0
        obj.cpx.addEq(deltaGf(i), obj.cpx.sum(Gf0(i, 1), obj.cpx.prod(RT, c_vars(i))));
      else % with uncertainties
        obj.cpx.addGe(deltaGf(i), obj.cpx.sum(obj.cpx.constant(Gf0(i, 1)),...
          obj.cpx.prod(RT, c_vars(i)), obj.cpx.constant(-Gf0(i, 2))));
        obj.cpx.addLe(deltaGf(i), obj.cpx.sum(obj.cpx.constant(Gf0(i, 1)),...
          obj.cpx.prod(RT, c_vars(i)), obj.cpx.constant(Gf0(i, 2))));
      end
    end
  end
  for i= 1:n
    if ~ignore(i)
      ind= find(stoichMat(:, i));
%       if any(isnan(Gf0(ind, 1))) % skip driving force if reaction deltaG is not fully determined
%         continue;
%       end
      if length(ind) == 1
        driving_force(i)= obj.cpx.prod(-stoichMat(ind, i), deltaGf(ind));
      else
        driving_force(i)= obj.cpx.scalProd(-stoichMat(ind, i), deltaGf(ind));
      end
    end
  end
else
  for i= 1:n
    if ~ignore(i)
      ind= find(stoichMat(:, i));
      if use_reaction_uncertainty
        if length(ind) == 1
          driving_force(i)= obj.cpx.diff(obj.cpx.diff(obj.cpx.prod(G0(i, 2), u_vars(i)), G0(i, 1)),...
            obj.cpx.prod(RT*stoichMat(ind, i), c_vars(ind)));
        elseif length(ind) > 1
          driving_force(i)= obj.cpx.diff(obj.cpx.diff(obj.cpx.prod(G0(i, 2), u_vars(i)), G0(i, 1)),...
            obj.cpx.scalProd(RT*stoichMat(ind, i), c_vars(ind)));
        end
      else
        if length(ind) == 1
          driving_force(i)= obj.cpx.diff(-G0(i), obj.cpx.prod(RT*stoichMat(ind, i), c_vars(ind)));
        elseif length(ind) > 1
          driving_force(i)= obj.cpx.diff(-G0(i), obj.cpx.scalProd(RT*stoichMat(ind, i), c_vars(ind)));
        end
      end
    end
  end
  deltaGf= [];
end

end

function res= nums(n)
res= strtrim(cellstr(num2str((1:n)')));
end
