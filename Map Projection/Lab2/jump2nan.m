function [y,b] = jump2nan(x,thr,a)

% JUMP2NAN inserts NaN's when differences within a vector
% exceed a certain threshold (in absolute sense)
%
% HOW: 
%    y = jump2nan(x,thr)
%    [y,b] = jump2nan(x,thr,a)
%
% IN:
%    x ...... input vector (where threshold condition is checked)
%    thr .... threshold
%    a ...... matrix, in which NaN's will be inserted as well 
%
% OUT:
%    y ...... copy of vector X, including NaN's
%    b ...... copy of matrix A, including NaN's
%
% REMARKS:
%    The NaN's in matrix A are inserted following the jumps in vector X.

% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Nico SNEEUW (NS), IAPG, TU-Munich
%    <bundle@gis.uni-stuttgart.de>
% -------------------------------------------------------------------------
% revision history:
%    1996-11-13: NS, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

% diagnostics and preliminaries
if min(size(x)) > 1,   error('X must be vector.'), end
if max(size(thr)) > 1, error('THReshold must be scalar'), end
isup   = diff(size(x)) < 0;			% is X a column vector?
x      = x(:);					% put X upright
nx     = length(x);

% determine index vector IDX
jumps  = (abs(diff(x)) > thr); 			% mask of jump occurences
idx    = cumsum([0;jumps]) + (1:nx)';		% index into output vector

% now the output Y vector can be created.
y      = NaN*ones(nx+sum(jumps),1);		% initialize output vector
y(idx) = x;					% do the mapping
if ~isup, y = y'; end				% lay down Y, if necessary

% handle the matrix A, if supplied
if nargin == 3
   [n,m] = size(a);
   if n == nx
      b = NaN*ones(nx+sum(jumps),m);
      b(idx,:) = a;
   elseif m == nx
      b = NaN*ones(n,nx+sum(jumps));
      b(:,idx) = a;
   else
      error('Dimensions of A don''t correspond to X')
   end
end
