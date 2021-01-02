function [tout,xout] = rk8fixed(FUN,tspan,x0,Nsteps,ode_fcn_format,trace,count)

% Copyright (C) 2001, 2000 Marc Compere
% This file is intended for use with Octave.
% rk8fixed.m is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% rk8fixed.m is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details at www.gnu.org/copyleft/gpl.html.
%
% --------------------------------------------------------------------
%
% rk8fixed (v1.11) is an 8th order Runge-Kutta numerical integration routine.
% It requires 13 function evaluations per step.  This is not the most
% efficient 8th order implementation.  It was just the easiest to put
% together as a variant from ode78.m.
%
% Usage:
%         [tout, xout] = rk8fixed(FUN, tspan, x0, Nsteps, ode_fcn_format, trace, count)
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
if nargin < 5, Nsteps = 50/(tspan(2)-tspan(1)); end % <-- 50 is a guess for a default,
                                                %  try verifying the solution with ode78
if nargin < 4, ode_fcn_format = 0; end

if count==1,
 global rhs_counter
 if ~exist('rhs_counter'),rhs_counter=0;,end
end % if count

alpha_ = [ 2./27., 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1 ]';
beta_  = [ 2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
          1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
          1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
          5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0, 0 ;
          0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0, 0 ;
          -25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0, 0 ;
          31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0, 0 ;
          2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0, 0 ;
          -91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0, 0 ;
          2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0, 0 ;
          3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0, 0 ;
          -1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1, 0 ]';
chi_   = [ 0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840]';

% Initialization
t = tspan(1);
h = (tspan(2)-tspan(1))/Nsteps;
xout(1,:) = x0';
tout(1) = t;
x = x0(:);
f = x*zeros(1,13);

if trace
 clc, t, h, x
end

for i=1:Nsteps,

     % Compute the slopes
     if (ode_fcn_format==0),
      f(:,1) = feval(FUN,t,x);
      for j = 1:12
         f(:,j+1) = feval(FUN, t+alpha_(j)*h, x+h*f*beta_(:,j));
      end
     else,
      f(:,1) = feval(FUN,x,t);
      for j = 1:12
         f(:,j+1) = feval(FUN, x+h*f*beta_(:,j), t+alpha_(j)*h);
      end
     end % if (ode_fcn_format==0)

     % increment rhs_counter
     if count==1,
      rhs_counter = rhs_counter + 13;
     end % if

     t = t + h;
     x = x + h*f*chi_;
     tout = [tout; t];
     xout = [xout; x.'];

     if trace,
      home, t, h, x
     end

end
