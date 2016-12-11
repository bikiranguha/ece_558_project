function [success, alphamn, alphamx] = gen_dispatch(Ltotal, flow1, flow2, ...
                                                  CapAva1, CapAva2, MWmax)
% [success, alphamn, alphamx] = gen_dispatch (Ltotal, flow1, flow2, ...
%   CapAva1, CapAva2, MWmax)
%
% inputs:
%   Ltotal: total load (MW)
%   flow1: vector of branch flows, based on G1 supplying Ltotal (MW)
%   flow2: vector of branch flows, based on G2 supplying Ltotal (MW)
%   CapAva1: G1 capacity available (MW)
%   CapAva2: G2 capacity available (MW)
%   MWmax: vector of branch flow limits (MW)
%
% outputs:
%   success: set to 0 if no dispatch exists within branch flow and generation
%       capacity constraints; set to 1 if dispatch exists
%   alphamn: minimum alpha value that satisfies constraints
%   alphamx: maximum alpha value that satisfies constraints

% Copyright 2012: Alexander J. Flueck, Illinois Institute of Technology

errf = 0;
[nr1,nc1] = size(flow1);
[nr2,nc2] = size(flow2);
[nrx,ncx] = size(MWmax);
if nr1 ~= nr2 | nr1 ~= nrx,
  fprintf('ERROR: # of rows must be the same for flow1, flow2 and MWmax\n');
  errf = 1;
end
if nc1 ~= 1 | nc2 ~= 1 | ncx ~= 1,
  fprintf('ERROR: # of columns for flow1, flow2 and MWmax must be one\n');
  errf = 1;
end
if CapAva1 < 0 | CapAva2 < 0 | Ltotal < 0,
  fprintf('ERROR: CapAva1, CapAva2 and Ltotal must be non-negative\n');
  errf = 1;
end
if errf == 1,
  error('Exiting "gen_dispatch"');
end

success = 1;
fprintf('\ngen_dispatch: Ltotal=%.1f, CapAva1=%.1f, CapAva2=%.1f\n',...
        Ltotal, CapAva1, CapAva2);

flowdiff = flow1-flow2;
zndx = find(flowdiff == 0);
if Ltotal > CapAva1+CapAva2 || sum(abs(flow1(zndx)) > MWmax(zndx)),
  success = 0;
  alphamn = inf;
  alphamx = -inf;
  fprintf(['gen_dispatch: FAILED - could not find dispatch satisfying ' ...
           'constraints\n']);
  return
end

nzndx = find(flowdiff ~= 0);
Llim = -(MWmax+flow2);
Rlim = MWmax-flow2;
Llim(nzndx) = Llim(nzndx)./flowdiff(nzndx);
Rlim(nzndx) = Rlim(nzndx)./flowdiff(nzndx);
Llim(zndx) = -inf;
Rlim(zndx) = inf;
ndx = find(flowdiff < 0);
if length(ndx) > 0,
  tmp = Llim(ndx);
  Llim(ndx) = Rlim(ndx);
  Rlim(ndx) = tmp;
end
Llim = [Llim; (Ltotal-CapAva2)/Ltotal];
Rlim = [Rlim; CapAva1/Ltotal];

alphamn = max(Llim);
alphamx = min(Rlim);

if alphamn > alphamx,
  success = 0;
end

if success == 0,
  fprintf(['gen_dispatch: FAILED - could not find dispatch satisfying ' ...
           'constraints\n']);
else
  fprintf('gen_dispatch: SUCCESS - found dispatch satisfying constraints\n');
end
