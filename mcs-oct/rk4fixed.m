function [tout,xout] = rk4fixed(FUN,tspan,x0,Nsteps,ode_fcn_format,trace,count)

% Copyright (C) 2001, 2000 Marc Compere
% This file is intended for use with Octave.
% rk4fixed.m is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% rk4fixed.m is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details at www.gnu.org/copyleft/gpl.html.
%
% --------------------------------------------------------------------
%
% rk4fixed (v1.11) is a 4th order Runge-Kutta numerical integration routine.
% It requires 4 function evaluations per step.
%
% Usage:
%         [tout, xout] = rk4fixed(FUN, tspan, x0, Nsteps, ode_fcn_format, trace, count)
%
% INPUT:
% FUN    - String containing name of user-supplied problem derivatives.
%          Call: xprime = fun(t,x) where FUN = 'fun'.
%          t      - Time or independent variable (scalar).
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
%          simply make 'rhs_counter' global in the file that calls rk4fixed
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
if nargin < 5, Nsteps = 100/(tspan(2)-tspan(1)); end % <-- 100 is a guess for a default,
                                                %  try verifying the solution with rk8fixed
if nargin < 4, ode_fcn_format = 0; end

if count==1,
 global rhs_counter
 if ~exist('rhs_counter'),rhs_counter=0;,end
end % if count

% Initialization
t = tspan(1);
h = (tspan(2)-tspan(1))/Nsteps;
xout(1,:) = x0';
tout(1) = t;
x = x0(:);
halfh = 0.5*h;

if trace
 clc, t, h, x
end

for i=1:Nsteps,
     if (ode_fcn_format==0),
      RK1 = feval(FUN,t,x);
      thalf = t+halfh;
      xtemp = x+halfh*RK1;
      RK2 = feval(FUN,thalf,xtemp);
      xtemp = x+halfh*RK2;
      RK3 = feval(FUN,thalf,xtemp);
      tfull = t+h;
      xtemp = x+h*RK3;
      RK4 = feval(FUN,tfull,xtemp);
     else,
      RK1 = feval(FUN,x,t);
      thalf = t+halfh;
      xtemp = x+halfh*RK1;
      RK2 = feval(FUN,xtemp,thalf);
      xtemp = x+halfh*RK2;
      RK3 = feval(FUN,xtemp,thalf);
      tfull = t+h;
      xtemp = x+h*RK3;
      RK4 = feval(FUN,xtemp,tfull);
     end % if (ode_fcn_format==0)

     % increment rhs_counter
     if count==1,
      rhs_counter = rhs_counter + 4;
     end % if

     t = t + h;
     x = (x+h/6*(RK1+2.0*(RK2+RK3)+RK4));
     tout = [tout; t];
     xout = [xout; x.'];

     if trace,
      home, t, h, x
     end

end
