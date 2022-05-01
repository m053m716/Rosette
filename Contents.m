%ROSETTE       This object creates a Rosette centered at some coordinate
%
% Public Properties:
%   Parent - Parent axes handle object.
%   Vector - A table that tracks the associated graphics objects.
%
% Protected Properties:
%   el    - Eventlistener that deletes Rosette if parent is deleted.
%   basis - Array of handles to line objects that form basis of Rosette
%   label - Array of handles to text objects that label the Rosette.
%   ring  - Line object that indicates values along basis.
%
% Methods:
%   Rosette - Constructor
%   add_vector - (Public) method to add vector(s) and labels to basis
%   get_scaled_value - (Public) Returns the value for indexed basis, including center offset and basis scaling.
%   modify_label - (Public) method that modifies the text label associated to a given indexed basis vector.
%   modify_value - (Public) method to modify value(s) of indexed basis vectors
%   modify_vector - (Public) method to modify indexed existing basis vectors on the Rosette.
%   move - (Public) method to move the center of the Rosette.
%
% Syntax:
%   r = Rosette(ax, xc, yc); % Create a Rosette with no vectors
%   r = Rosette(ax, xc, yc, V); % Create a Rosette with unlabeled vectors
%   r = Rosette(ax, xc, yc, V, label); % Create a Rosette with labeled vectors.
%   r = Rosette(ax, xc, yc, V, label, value); % Create a Rosette with labeled vectors and values on those vectors.
%   r = Rosette(__, prop_vals); % Specify <'Name', Value> pairs for each vector.
%
% Inputs:
%   ax - Parent axes handle.
%   xc - Center x-coordinate on axes (default: 0).
%   yc - Center y-coordinate on axes (default: 0).
%   V  - nVector x 2 array of <x, y> coordinates specifying the
%        zero-centered vector for each Rosette "spoke" or
%        "axis".
%   labels - nVector x 1 array of strings or cell array of
%            character vectors labeling each axis.
%   value - nVector x 1 array of values indicating scaling
%            along each axis.
%   prop_val - nVector x 1 array of cells, where each cell
%               element is itself either an empty cell array,
%               or a cell array with an even number of
%               elements, where the first element is a `Line`
%               property name and the second is the
%               corresponding value to set for that `Line`
%               object, which corresponds to a "spoke" on the
%               Rosette.
%
% Outputs:
%   r - Rosette object (handle subclass)
    
