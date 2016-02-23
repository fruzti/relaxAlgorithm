% function x = mls(K,reps)
%
% Generates a K'th order mls signal up to an order of 20.
%
% The recursion indices are taken from
% W. Stahnke. Primitive binary polynomials.
% Mathematics of Computation, 27:977-980, 1973. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K: MLS order ( 0 < K < 21)
%    size: 1x1
% reps: number of repetitions of the MLS sequence
%    size: 1x1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x: Data column vector
%    size: Nx1 where N = 2^K-1
%
% Created:        2013-10-08  by Jesper Kjær Nielsen (jkn@es.aau.dk)
% Last modified:  2013-11-07  by Jesper Kjær Nielsen (jkn@es.aau.dk)

function x = mls(K,reps)
  % validate the MLS order
  if K < 1 || K > 20
    error('K is not in the range 1 < K < 21');
  end
  % the length of a single MLS sequence of order K
  N = 2^K-1;
  x = zeros(N,1);
  a = zeros(K,1);
  % generate a random initial point
  while sum(a) == 0
    a = round(rand(K,1));
  end

  % compute the MLS sequence
  idx = recursion_idx(K);
  for ii=1:N
    x(ii) = a(1);
    a(1:K-1) = a(2:K);
    a(K) = mod(sum([x(ii);a(idx)]),2);
  end
  x = -2*x+1;
  % repeat the sequence reps times
  x = repmat(x,reps,1);
end

function idx = recursion_idx(k)
  if k == 1
    idx = [];
  elseif sum((k-[2,3,4,6,7,15]) == 0)
    idx = 1;
  elseif sum((k-[5,11]) == 0)
    idx = 2;
  elseif sum((k-[8,19]) == 0)
    idx = [1,5,6];
  elseif k == 9
    idx = 4;
  elseif k == 12
    idx = [3,4,7];
  elseif k == 13
    idx = [1,3,4];
  elseif k == 14
    idx = [1,11,12];
  elseif k == 16
    idx = [2,3,5];
  elseif sum((k-[10,17,20]) == 0)
    idx = 3;
  elseif k == 18
    idx = 7;
  elseif k == 19
    idx = [1,5,6];
  else
    idx = nan;
  end
end
