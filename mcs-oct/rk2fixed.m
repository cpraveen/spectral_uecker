function [tout,xout] = rk2fixed(FUN,tspan,x0,Nsteps,ode_fcn_format,trace,count)

% Copyright (C) 2001, 2000 Marc Compere
% This file is intended for use with Octave.
% rk2fixed.m is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% rk2fixed.m is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details at www.gnu.org/copyleft/gpl.html.
%
% --------------------------------------------------------------------
%
% rk2fixed (v1.11) integrates a system of ordinary differential equations using a
% 2nd order Runge-Kutta formula called Ralston's method.
% This choice of 2nd order coefficients provides a minimum bound on truncation error.
% For more, see Ralston & Rabinowitz (1978) or 
% Numerical Methods for Engineers, 2nd Ed., Chappra & Cannle, McGraw-Hill, 1985
%
% rk2fixed() requires 2 function evaluations per integration step.
%
% Usage:
%         [tout, xout] = rk2fixed(FUN, tspan, x0, Nsteps, ode_fcn_format, trace, count)
%
% INPUT:
% FUN    - String containing name of user-supplied problem derivatives.
%          Call: xprime = fun(t,x) where FUN = 'fun'.
%          t      - Time (scalar).
%          x      - Solution column-vector.
%          xprime - Returned derivative COLUMN-vector; xprime(i) = dx(i)/dt.
% tspan  - [ tstart, tfinal ]
% x0     - Initial value COLUMN-vector.
% Nsteps - number of steps used to span [ tstart, tfinal ]
% ode_fcn_format - this specifies if the user-defined ode function is in
%          the form:     xprime = fun(t,x)   (ode_fcn_format=0, default)
%          or:           xprime = fun(x,t)   (ode_fcn_format=1)
%          Matlab's solvers comply with ode_fcn_format=0 while
%          Octave's lsode() and sdirk4() solvers comply with ode_fcn_format=1.
% trace  - If nonzero, each step is printed. (optional, default: trace = 0).
% count  - if nonzero, variable 'rhs_counter' is initalized, made global
%          and counts the number of state-dot function evaluations
%          'rhs_counter' is incremented in here, not in the state-dot file
%          simply make 'rhs_counter' global in the file that calls rk2fixed
%
% OUTPUT:
% tout  - Returned integration time points (row-vector).
% xout  - Returned solution, one solution column-vector per tout-value.
%
% The result can be displayed by: plot(tout, xout).
%
% Marc Compere
% CompereM@asme.org
% created : 06 October 1999
% modified: 17 January 2001

if nargin < 7, count = 0; end
if nargin < 6, trace = 0; end
if nargin < 5, Nsteps = 400/(tspan(2)-tspan(1)); end % <-- 400 is a guess for a default,
                                                %  try verifying the solution with rk4fixed
if nargin < 4, ode_fcn_format = 0; end

if count==1,
 global rhs_counter
 if ~exist('rhs_counter'),rhs_counter=0;,end
end % if count

% Initialization
t = tspan(1);
h = (tspan(2)-tspan(1))/Nsteps;
xout(1,:) = x0';           
tout(1)   = t;
x = x0(:);

if trace
 clc, t, h, x
end

% The main loop
h = (tspan(2)-tspan(1))/Nsteps;

for i=1:Nsteps,
     if (ode_fcn_format==0),
      k1 = feval(FUN,t,x);
      k2 = feval(FUN,t+3/4*h,x+3/4*h*k1);
     else,
      k1 = feval(FUN,x,t);
      k2 = feval(FUN,x+3/4*h*k1,t+3/4*h);
     end % if (ode_fcn_format==0)

     % increment rhs_counter
     if count==1,
      rhs_counter = rhs_counter + 2;
     end % if

     t = t + h;
     x = (x+h*(1/3*k1+2/3*k2));
     tout = [tout; t];
     xout = [xout; x.'];

     if trace,
      home, t, h, x
     end

end
