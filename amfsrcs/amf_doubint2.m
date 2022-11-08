function [ xx ] = amf_doubint2( nx, xx, der )
% modified from software Madagascar (University of Texas at Austin)

%%
%   Copyright (C) 2016 Delft University of Technology -- Delphi consortium - Shan Qu
%   
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%%
% integrate forward
t = 0.0;
for i = 0 : nx-1
    t = t + xx(i+1);
    xx(i+1) = t;
end


if(der)
    return;
end

% integrate backward
t = 0.0;
for i = nx-1 : -1 : 0
    t = t + xx(i+1);
    xx(i+1) = t;
end



end

