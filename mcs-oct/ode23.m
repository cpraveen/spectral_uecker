function [tout,xout] = ode23(FUN,tspan,x0,ode_fcn_format,tol,trace,count)

% Copyright (C) 2001, 2000 Marc Compere
% This file is intended for use with Octave.
% ode23.m is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% ode23.m is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details at www.gnu.org/copyleft/gpl.html.
%
% --------------------------------------------------------------------
%
% ode23 (v1.11) Integrates a system of ordinary differential equations using
% 2nd & 3rd order Runge-Kutta formulas.  This particular 3rd-order method reduces
% to Simpson's 1/3 rule and uses the 3rd order estimate for xout.
%
% 3rd-order RK methods have local and global errors of O(h^4) and O(h^3),
% respectively and yield exact results when the solution is a cubic.
%
% The order of the RK method is the order of the local *truncation* error, d.
% The local truncation error is defined as the principle error term in the
% portion of the Taylor series expansion that gets dropped.  This portion of
% the Taylor series exapansion is within the group of terms that gets multipled
% by h in the solution definition of the general RK method.  Therefore, the
% order-p solution created by the RK method will be roughly accurate to
% within O(h^(p+1)).  The difference between two different-order solutions
% is the definition of the "local error," l.  This makes the local error as
% large as the error in the lower order method, which is the truncation error
% times h, resulting in O(h^(p+1)).
% Summary:   For an RK method of order-p,
%            - the local truncation error is O(h^p)
%            - the local error (used for stepsize adjustment) is O(h^(p+1))
%
% This requires 3 function evaluations per integration step.
%
% Relevant discussion on step size choice can be found on pp.90,91 in
% U.M. Ascher, L.R. Petzold, Computer Methods for  Ordinary Differential Equations
% and Differential-Agebraic Equations, Society for Industrial and Applied Mathematics
% (SIAM), Philadelphia, 1998
%
% The error estimate formula and slopes are from
% Numerical Methods for Engineers, 2nd Ed., Chappra & Cannle, McGraw-Hill, 1985
%
% Usage:
%         [tout, xout] = ode23(FUN, tspan, x0, ode_fcn_format, tol, trace, count)
%
% INPUT:
% FUN   - String containing name of user-supplied problem description.
%         Call: xprime = fun(t,x) where FUN = 'fun'.
%         t      - Time (scalar).
%         x      - Solution column-vector.
%         xprime - Returned derivative COLUMN-vector; xprime(i) = dx(i)/dt.
% tspan - [ tstart, tfinal ]
% x0    - Initial value COLUMN-vector.
% ode_fcn_format - this specifies if the user-defined ode function is in
%         the form:     xprime = fun(t,x)   (ode_fcn_format=0, default)
%         or:           xprime = fun(x,t)   (ode_fcn_format=1)
%         Matlab's solvers comply with ode_fcn_format=0 while
%         Octave's lsode() and sdirk4() solvers comply with ode_fcn_format=1.
% tol   - The desired accuracy. (optional, default: tol = 1.e-6).
% trace - If nonzero, each step is printed. (optional, default: trace = 0).
% count - if nonzero, variable 'rhs_counter' is initalized, made global
%         and counts the number of state-dot function evaluations
%         'rhs_counter' is incremented in here, not in the state-dot file
%         simply make 'rhs_counter' global in the file that calls ode23
%
% OUTPUT:
% tout  - Returned integration time points (column-vector).
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
if nargin < 5, tol = 1.e-3; end
if nargin < 4, ode_fcn_format = 0; end

pow = 1/4; % see p.91 in the Ascher & Petzold reference for more infomation.

% The 2(3) coefficients:
 a(1,1)=0;
 a(2,1)=1/2;
 a(3,1)=-1; a(3,2)=2;
 % 2nd order b-coefficients
 b2(1)=0; b2(2)=1;
 % 5th order b-coefficients
 b3(1)=1/6; b3(2)=2/3; b3(3)=1/6;
 for i=1:3
  c(i)=sum(a(i,:));
 end

% Initialization
t0 = tspan(1);
tfinal = tspan(2);
t = t0;
hmax = (tfinal - t)/2.5;
hmin = (tfinal - t)/1e12;
h = (tfinal - t)/200; % initial guess at a step size
x = x0(:);            % this always creates a column vector, x
tout = t;             % first output time
xout = x.';           % first output solution

if count==1,
 global rhs_counter
 if ~exist('rhs_counter'),rhs_counter=0; end
end % if count

if trace
 clc, t, h, x
end

% The main loop
   while (t < tfinal) & (h >= hmin)
      if t + h > tfinal, h = tfinal - t; end

      % compute the slopes
      if (ode_fcn_format==0), % (default)
         k(:,1)=feval(FUN,t,x);
         k(:,2)=feval(FUN,t+c(2)*h,x+h*(a(2,1)*k(:,1)));
         k(:,3)=feval(FUN,t+c(3)*h,x+h*(a(3,1)*k(:,1)+a(3,2)*k(:,2)));
      else, % ode_fcn_format==1
         k(:,1)=feval(FUN,x,t);
         k(:,2)=feval(FUN,x+h*(a(2,1)*k(:,1)),t+c(2)*h);
         k(:,3)=feval(FUN,x+h*(a(3,1)*k(:,1)+a(3,2)*k(:,2)),t+c(3)*h);
      end % if (ode_fcn_format==1)

      % increment rhs_counter
      if count==1, rhs_counter = rhs_counter + 3; end

      % compute the 2nd order estimate
      x2=x + h*b2(2)*k(:,2);
      % compute the 3rd order estimate
      x3=x + h*(b3(1)*k(:,1) + b3(2)*k(:,2) + b3(3)*k(:,3));

      % estimate the local truncation error
      gamma1 = x3 - x2;

      % Estimate the error and the acceptable error
      delta = norm(gamma1,'inf');
      tau = tol*max(norm(x,'inf'),1.0);

      % Update the solution only if the error is acceptable
      if delta <= tau
         t = t + h;
         x = x3;    % <-- using the higher order estimate is called 'local extrapolation'
         tout = [tout; t];
         xout = [xout; x.'];
      end
      if trace
         home, t, h, x
      end

      % Update the step size
      if delta == 0.0
       delta = 1e-16;
      end
      h = min(hmax, 0.8*h*(tau/delta)^pow);

   end;

   if (t < tfinal)
      disp('Step size grew too small.')
      t, h, x
   end
