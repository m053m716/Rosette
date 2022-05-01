classdef Rosette < handle
    %ROSETTE  This object creates a Rosette centered at some coordinate
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
    %
    % See also: Contents
    
    properties (GetAccess = public, SetAccess = protected)
        Parent      matlab.graphics.axis.Axes
        Position    double = [0, 0]
        Vector      table
    end
    
    properties (Hidden, Access = public)
        label_offset    double = 1.25;  % The scalar for offsetting labels from tip of the basis indicator 
    end
    
    properties (Access = protected)
        basis   matlab.graphics.primitive.Line  % Array of line objects indicating basis in 2D plane
        label   matlab.graphics.primitive.Text  % Primitive text object array labeling each basis vector
        ring    matlab.graphics.primitive.Line  % Primitive line object indicating values along basis
        r       double = [] % Array of radii to simplify scaling along this basis.
        theta   double = [] % Array of thetas to simplify scaling along the basis.
    end
    
    methods
        function self = Rosette(ax, xc, yc, V, label, value, prop_vals)
            %Rosette  Constructor for Rosette graphics object
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
            %   label - nVector x 1 array of strings or cell array of
            %            character vectors labeling each axis.
            %   value - nVector x 1 array of values indicating scaling
            %            along each axis.
            %   prop_vals - nVector x 1 array of cells, where each cell
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
            %
            % See also: Contents
            self.Parent = ax;
            if nargin >= 3
                self.Position = [xc, yc];
            end
            self.Vector = table('Size', [0, 5], ...
                'VariableTypes', {'double', 'double', 'double', 'string', 'double'}, ...
                'VariableNames', {'index', 'x', 'y', 'label', 'value'});
            if nargin < 4
                return;
            end
            if nargin < 5
                label = strings(size(V, 1), 1);
            end
            if nargin < 6
                value = nan(size(V, 1), 1);
            end
            if nargin < 7
                prop_vals = repmat({{}}, size(V, 1), 1);
            end
            self.ring = line(ax, nan, nan, ...
                'LineWidth', 1, 'Color', 'k', 'LineStyle', ':', ...
                'Marker', 'o', 'MarkerFaceColor', 'b');
            self.ring.Annotation.LegendInformation.IconDisplayStyle = 'off';
            self.add_vector(V, label, value, prop_vals);
            if isempty(self.Parent.DeleteFcn)
                self.Parent.UserData.chain_callback = {@(~, ~)delete(self)};
                self.Parent.DeleteFcn = @Rosette.callback_chainer;
            else
                self.Parent.UserData.chain_callback{end+1} = @(~, ~)delete(self);
            end
        end
        
        function add_vector(self, V, label, value, prop_vals)
            %ADD_VECTOR  Add new vector(s) to the Rosette
            %
            % Syntax:
            %   r.add_vector(V);
            %   r.add_vector(V, label, value, prop_vals);
            %
            % Inputs:
            %   V - nVector2add x 2 double array, first element is x center
            %       and second element is y center for array.
            %   label - nVector2add x 1 string or character array of labels
            %           corresponding to rows of V.
            %   value - nVector2add x 1 vector of values (optional) to set
            %           the value along the new bases.
            %   prop_vals - nVector2add x 1 cell array of cell arrays of
            %               {'Name',value} property value pairs for MATLAB
            %               `line` builtin.
            %
            % See also: Contents, Rosette, modify_value, modify_vector
            if nargin < 2
                label = strings(size(V, 1), 1);
            end
            if nargin < 3
                value = nan(size(V, 1), 1);
            end
            if nargin < 4
                prop_vals = repmat({{}}, size(V, 1), 1);
            end
            if size(V, 1) > 1
                if numel(label) ~= size(V, 1)
                    error("Dimension mismatch: number of rows in V must match number of elements of label.");
                end
                
                for ii = 1:size(V, 1)
                    self.add_vector(V(ii, :), label(ii), value(ii), prop_vals{ii});
                end
                return;
            end
            index = size(self.Vector, 1) + 1;
            x = V(1);
            y = V(2);
            X = [self.Position(1), self.Position(1) + x];
            Y = [self.Position(2), self.Position(2) + y];
            
            T = table(index, x, y, label, value);
            self.Vector = vertcat(self.Vector, T);
            h = line(self.Parent, X, Y, ...
                'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '-', ...
                prop_vals{:});
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            self.basis(index) = h;
            self.theta(index) = atan2(y, x);
            self.r(index) = sqrt((X(2) - X(1))^2 + (Y(2) - Y(1))^2);
            self.label(index) = text(self.Parent, ...
                self.label_offset .* self.r(index) .* cos(self.theta(index)) + self.Position(1), ...
                self.label_offset .* self.r(index) .* sin(self.theta(index)) + self.Position(2), label, ...
                'FontName', 'Tahoma', 'Color', 'k');
            xv = value .* self.r(index) .* cos(self.theta(index)) + self.Position(1);
            yv = value .* self.r(index) .* sin(self.theta(index)) + self.Position(2);
            self.ring.XData(index) = xv;
            self.ring.YData(index) = yv;
            self.ring.XData(index+1) = self.ring.XData(1);
            self.ring.YData(index+1) = self.ring.YData(1); % So it "loops back" on itself.
        end
        
        function [xv, yv] = get_scaled_value(self, index)
            %GET_SCALED_VALUE  Return scaled, translated value of Rosette along indexed basis vector.
            %
            % Syntax:
            %   [xv, yv] = r.get_scaled_value(index);
            %
            % Inputs:
            %   index - Scalar or vector indexing back into self.Vector
            %           table matching the index variable there to
            %           reference the corresponding row (and elements from
            %           self.r and self.theta).
            %
            % Output:
            %   [xv, yv] - Coordinates in 2D plane including the offset of
            %               the center of the Rosette and the scaling along
            %               the indexed basis axis.
            %
            % See also: Contents            
            if numel(self) > 1
                xv = cell(size(self));
                yv = cell(size(self));
                for ii = 1:numel(self)
                    [xv{ii}, yv{ii}] = self(ii).get_scaled_value(index{ii});
                end
                return;
            end
            
            xv = nan(size(index));
            yv = nan(size(index));
            for ii = 1:numel(index)
                k = self.Vector.index == index(ii);
                xv = self.Vector.value(k) .* self.r(k) .* cos(self.theta(k)) + self.Position(1);
                yv = self.Vector.value(k) .* self.r(k) .* sin(self.theta(k)) + self.Position(2);
            end
        end
        
        function modify_label(self, index, label, varargin)
            %MODIFY_LABEL  Change values of 'Name', value pairs for label
            %
            % Syntax:
            %   r.modify_label(index, label, 'Name', value, ...);
            %       % Apply same formatting to all label elements
            %   r.modify_label(index, label, {'Name', value, ...},{'Name'...});
            %       % Match each specific cell subset list to indexed
            %       %   label elements.
            %
            % Inputs:
            %   index - The index from self.Vector to use when matching
            %           indexing for which basis vectors to change labels
            %           on.
            %   label - String or string array or cell array of character
            %           vectors that is the new label for any indexed basis
            %           vectors of the Rosette. Can leave this empty to
            %           skip it and just use varargin 'Name', value syntax.
            %   varargin - {'Name', value} pairs for the indexed label,
            %               can use any modifiable {'Name',value} pair from
            %               MATLAB text object builtin.
            %
            % See also: Contents, modify_ring, text
            if islogical(index)
                vec = 1:numel(index);
                index = vec(index);
            end
            if isempty(label)
                label = self.Vector.label(index);
            end
            if numel(label) ~= numel(index)
                error("Number of elements in value and index should match.");
            end
            if numel(varargin) ~= numel(index)
                varargin = repmat({varargin}, size(index)); 
            end
            for ii = 1:numel(index)
                v = varargin{ii};
                k = self.Vector.index == index(ii);
                set(self.label(k), 'String', label, v{:});
            end
        end
        
        function modify_ring(self, varargin)
            %MODIFY_RING  Change values of 'Name', value pairs for "value" ring
            %
            % Syntax:
            %   r.modify_ring('Name', value, ...);
            %
            % Accepts any 'Name', value argument pairs that MATLAB builtin
            % `line` would accept.
            %
            % See also: Contents, line
            set(self.ring, varargin{:});
        end
        
        function modify_value(self, index, value)
            %MODIFY_VALUE  Change the value of some indexed point on basis
            %
            % Syntax:
            %   r.modify_value(index, value);
            %
            % Inputs:
            %   index - Index to the basis value to be modified, see the
            %           `Vector` property. If this is a vector then number
            %           of elements in value should match number of
            %           elements in index.
            %   value - The new value(s) to update to.
            %
            % See also: Contents, Rosette, modify_vector
            if islogical(index)
                vec = 1:numel(index);
                index = vec(index);
            end
            if numel(value) ~= numel(index)
                error("Number of elements in value and index should match.");
            end
            for ii = 1:numel(index)
                k = find(self.Vector.index == index(ii), 1, 'first');
                xv = value(k) .* self.r(k) .* cos(self.theta(k)) + self.Position(1);
                yv = value(k) .* self.r(k) .* sin(self.theta(k)) + self.Position(2);
                self.ring.XData(k) = xv;
                self.ring.YData(k) = yv;
                self.Vector.value(k) = value(ii);
                if k == 1
                    self.ring.XData(end) = xv;
                    self.ring.YData(end) = yv;
                end
            end            
        end
        
        function modify_vector(self, index, V, varargin)
            %MODIFY_VECTOR  Change properties or labels on basis
            %
            % Syntax:
            %   r.modify_vector(index, V, 'Name', value, ...);
            %   % Modifies each vector with same 'Name', value, ... options
            %
            %   r.modify_vector(index, V, {'Name', value, ...}, {'Name', value, ...})
            %   % Modify each basis vector with matched-specific 'Name',
            %   % value pairs.
            %   
            %   r.modify_vector(index, [], ...);
            %   % Does not change bases.
            %
            %   r.modify_vector([], [], ...);
            %   % Apply changes to all bases.
            %
            % Inputs:
            %   index - Index to the basis vector(s) to be modified, see
            %           `Vector` property. If this is a vector then number
            %           of rows of V, number of elements in label, and cell
            %           elements of prop_vals, should all match number of
            %           elements in index.
            %   V - nIndex x 2 new zero-centered coordinates to point
            %           indexed basis vector(s) to. Can also leave this
            %           empty to keep it the same and just modify varargin
            %           instead.
            %   varargin - nIndex x 1 cell array of {'Name', value} cell
            %               arrays, used to update indexed primitive line
            %               objects according to MATLAB builtin properties
            %               for line primitive. OR standard 'Name', value
            %               syntax and it will be generalized to all the
            %               indexed lines.
            %
            % See also: Contents, Rosette, modify_value, line
            if isempty(index)
                index = 1:size(self.Vector, 1); 
            end
            if islogical(index)
                vec = 1:numel(index);
                index = vec(index);
            end
            if isempty(V)
                V = [self.Vector.x(index), self.Vector.y(index)];
            end
            if size(V, 1) ~= numel(index)
                error("Number of rows in V and index should match.");
            end
            if numel(varargin) ~= size(V, 1)
                varargin = repmat({varargin}, size(V, 1), 1);
            end
            for ii = 1:numel(index)
                k = self.Vector.index == index(ii);
                v = varargin{ii};
                x = V(ii, 1);
                y = V(ii, 2);
                self.Vector.x(k) = x;
                self.Vector.y(k) = y;
                X = [self.Position(1), self.Position(1) + x];
                Y = [self.Position(2), self.Position(2) + y];
                self.theta(k) = atan2(y, x);
                self.r(k) = sqrt((X(2) - X(1))^2 + (Y(2) - Y(1))^2);
                set(self.basis(k), 'XData', X, 'YData', Y, v{:});
                self.label(k).Position(1:2) = [X(1) + self.label_offset * self.r(k) .* cos(self.theta(k)), ...
                                               Y(1) + self.label_offset * self.r(k) .* sin(self.theta(k))];
                [self.ring.XData(k), self.ring.YData(k)] = self.get_scaled_value(index(ii));
            end 
        end
        
        function move(self, xc, yc)
            %MOVE  Move center coordinate of Rosette.
            %
            % Inputs:
            %   xc - X-center position
            %
            % See also: Contents, Rosette
            
            if numel(self) > 1
                if numel(xc) == 1
                    xc = ones(size(self)) .* xc;
                end
                if numel(yc) == 1
                    yc = ones(size(self)) .* yc;
                end
                for ii = 1:numel(self)
                    self(ii).move(xc(ii), yc(ii));
                end
                return;
            end
            xd = xc - self.Position(1);
            yd = yc - self.Position(2);
            for ii = 1:numel(self.basis)
                set(self.basis(ii), 'XData', self.basis(ii).XData + xd, ...
                    'YData', self.basis(ii).YData + yd);
                self.label(ii).Position(1) = self.label(ii).Position(1) + xd;
                self.label(ii).Position(2) = self.label(ii).Position(2) + yd;
            end
            set(self.ring, ...
                'XData', self.ring.XData + xd, ...
                'YData', self.ring.YData + yd);
            self.Position = [xc, yc];
        end
    end
    
    methods (Static)
        function callback_chainer(src, ~)
            %CALLBACK_CHAINER  Hack around issue with not being able to string together callbacks easily.
            for ii = 1:numel(src.UserData.chain_callback)
                feval(src.UserData.chain_callback{ii}); 
            end
        end
    end
end

