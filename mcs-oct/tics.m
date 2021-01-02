## tics(axis,[pos1,pos2,...],['lab1';'lab2';...])
##
## Explicitly set the tic positions and labels for the given axis.
##
## If no positions or labels are given, then restore the default.
## If positions are given but no labels, use those positions with the
## normal labels.  If positions and labels are given, each position
## labeled with the corresponding row from the label matrix.
##
## Axis is 'x', 'y' or 'z'.

## This program is in the public domain
## Author: Paul Kienzle <pkienzle@users.sf.net>

## Modified to use new gnuplot interface in octave > 2.9.0
## Dmitri A. Sergatskov <dasergatskov@gmail.com>
## April 18, 2005

## Modifications which makes the y, z axis tics work. It was set to
## always do x axis before.
## 2007-08-12 Russel Valentine and Peter Gustafson

function tics (axis,pos,lab)

  if (nargin == 0)
    usage ("tics(axis,[pos1,pos2,...],['lab1';'lab2';...])");
  endif

  t = lower (axis);
  if (t ~= "x" && t ~= "y" && t ~= "z")
    error ("First input argument must be one of 'x', 'y' or 'z'");
  endif

  if (nargin == 1)
    set (gca(), [t, "tick"], []);
    set (gca(), [t, "tickmode"], "auto");
    set (gca(), [t, "ticklabel"], "");
    set (gca(), [t, "ticklabelmode"], "auto");
  elseif (nargin == 2)
    set (gca(), [t, "tick"], pos);
    set (gca(), [t, "ticklabel"], "");
    set (gca(), [t, "ticklabelmode"], "auto");
  elseif (nargin == 3)
    set (gca(), [t, "tick"], pos);
    set (gca(), [t, "ticklabel"], lab);
  else
    usage ("tics(axis,[pos1,pos2,...],['lab1';'lab2';...])");
  endif

endfunction
