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
    %   modify_cloud - (Public) method that changes properties of the "point cloud" along one or more rosette axes.
    %   modify_label - (Public) method that modifies the text label associated to a given indexed basis vector.
    %   modify_label_offset - (Public) Change value of protected `label_offset` property
    %   modify_value - (Public) method to modify value(s) of indexed basis vectors
    %   modify_vector - (Public) method to modify indexed existing basis vectors on the Rosette.
    %   move - (Public) method to move the center of the Rosette.
    %
    % Syntax:
    %   r = Rosette(ax, xc, yc); % Create a Rosette with no vectors
    %   r = Rosette(ax, xc, yc, V); % Create a Rosette with unlabeled vectors
    %   r = Rosette(ax, xc, yc, V, label); % Create a Rosette with labeled vectors.
    %   r = Rosette(ax, xc, yc, V, label, value); % Create a Rosette with labeled vectors and values on those vectors.
    %   r = Rosette(__, {'Name', value}, ...); % Specify <'Name', Value> pairs for each vector.
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
        
        CData           uint8           % Colormap used to key the labels and vectors.
        CDataMode       string = "auto";% Mode for assigning colors. Default is "auto" and distributes colorspace along CData. "flat" just increases color index by 1 for each additional rosette vector.
        Jitter          double = 0.020; % The jitter on 2D Gaussian applied to points in data cloud along each axes
        LabelOffset     double = 1.25;  % The scalar for offsetting labels from tip of the basis indicator 
    end
    
    properties (Access = protected)
        basis   matlab.graphics.primitive.Line           % Array of line objects indicating basis in 2D plane
        label   matlab.graphics.primitive.Text           % Primitive text object array labeling each basis vector
        cloud   matlab.graphics.chart.primitive.Scatter  % Primitive scatter object indicating values along basis
        c       uint8  = []                              % The actual colors that are used (depends on CDataMode and CData)
        r       double = []                              % Array of radii to simplify scaling along this basis.
        theta   double = []                              % Array of thetas to simplify scaling along the basis.
    end
    
    methods
        function self = Rosette(ax, xc, yc, V, label, value, varargin)
            %Rosette  Constructor for Rosette graphics object
            %
            % Syntax:
            %   r = Rosette(ax, xc, yc); % Create a Rosette with no vectors
            %   r = Rosette(ax, xc, yc, V); % Create a Rosette with unlabeled vectors
            %   r = Rosette(ax, xc, yc, V, label); % Create a Rosette with labeled vectors.
            %   r = Rosette(ax, xc, yc, V, label, value); % Create a Rosette with labeled vectors and values on those vectors.
            %   r = Rosette(__, {'Name', value, ...}, {'Name', value, ...}, 'Name', value, ...); 
            %       
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
            %   varargin - Specify {'Name', Value, ...} pairs for each
            %                   nVector x 1 array element in `label` and
            %                   `value`. If left empty they will use
            %                   default values.
            %
            %              Specify <'Name', value> pairs (outside of cell 
            %               arrays) to set global Rosette properties:
            %                   * CData: Colormap used to key the labels 
            %                   * CDataMode: "auto" or "fixed"
            %                   * Jitter: Jitter on 2D Gaussian applied to 
            %                             points in data cloud along each 
            %                             axes.
            %                   * LabelOffset: The scalar for offsetting 
            %                                  labels
            %
            % Outputs:
            %   r - Rosette object (handle subclass)
            %
            % See also: Contents
            self.CData = cm.map('rosette'); 
            if nargin < 1
                ax = gca;
            end
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

            k = 1;
            while k <= numel(varargin)
                if iscell(varargin{k})
                    k = k + 1;
                else
                    try
                        self.(varargin{k}) = varargin{k+1};
                    catch
                        if (k + 1) > numel(varargin)
                            error("Bad number of 'Name', value pairs, check varargin input.");
                        end
                        warning("Error assigning varargin{%d} (varargin{%d}): ", k, k+1);
                        disp(varargin{k});
                        disp(varargin{k+1});
                    end
                    varargin([k, k+1]) = [];
                end
            end

            if numel(varargin) == 1
                varargin = repmat(varargin, size(V, 1), 1); 
            elseif numel(varargin) == 0
                varargin = repmat({{}}, size(V, 1), 1);
            end
            self.add_vector(V, label, value, varargin);
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
                value = num2cell(zeros(size(V, 1), 1));
            end
            if nargin < 4
                prop_vals = repmat({{}}, size(V, 1), 1);
            else
                if numel(prop_vals) == 1
                    prop_vals = repmat(prop_vals, size(V, 1), 1);
                end
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
            self.modify_colormap(self.CData, self.CDataMode);
            h = line(self.Parent, X, Y, ...
                'Color', self.c(index, :), ...
                'LineWidth', 0.5, 'LineStyle', ':', ...
                'Tag', sprintf('%s.basis', label), ...
                'MarkerIndices', 1, 'Marker', 'o', ...
                'MarkerSize', 5, 'MarkerFaceColor', 'k', ...
                prop_vals{:});
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            self.basis(index) = h;
            self.theta(index) = atan2(y, x);
            self.r(index) = sqrt((X(2) - X(1))^2 + (Y(2) - Y(1))^2);
            rot = rad2deg(self.theta(index));
            if (rot > 90) && (rot <= 270)
                rot = rot - 180;
                ha = 'right';
            elseif (rot < -90) && (rot >= -270)
                rot = rot + 180;
                ha = 'right';
            else
                ha = 'left';
            end
            self.label(index) = text(self.Parent, ...
                self.LabelOffset .* self.r(index) .* cos(self.theta(index)) + self.Position(1), ...
                self.LabelOffset .* self.r(index) .* sin(self.theta(index)) + self.Position(2), label, ...
                'FontName', 'Tahoma', 'FontWeight', 'bold', ...
                'Color', self.c(index, :), ...
                'HorizontalAlignment', ha, ...
                'VerticalAlignment', 'middle', ...
                'Rotation', rot, 'Tag', sprintf('%s.label', label));
            z = [cos(self.theta(index)) -sin(self.theta(index)); sin(self.theta(index)) cos(self.theta(index))] * randn(2, numel(value{1})) .* self.Jitter .* self.r(index);
            xv = value{1} .* self.r(index) .* cos(self.theta(index)) + self.Position(1) + z(1, :)';
            yv = value{1} .* self.r(index) .* sin(self.theta(index)) + self.Position(2) + z(2, :)';
            self.cloud(index) = scatter(self.Parent, xv, yv, ...
                'MarkerFaceAlpha', 0.20, 'MarkerEdgeColor', 'none', 'SizeData', 4, ...
                'Marker', 'o', 'MarkerFaceColor', self.c(index, :), ...
                'Tag', sprintf('%s.cloud', label));
            self.cloud(index).Annotation.LegendInformation.IconDisplayStyle = 'off';
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
            if numel(index) > 1
                xv = cell(size(index));
                yv = cell(size(index));
                for ii = 1:numel(index)
                    k = self.Vector.index == index(ii);
                    z = [cos(self.theta(k)) -sin(self.theta(k)); sin(self.theta(k)) cos(self.theta(k))] * randn(2, numel(self.Vector.value{k})) .* self.Jitter .* self.r(k);
                    xv{ii} = self.Vector.value{k} .* self.r(k) .* cos(self.theta(k)) + self.Position(1) + z(1, :)';
                    yv{ii} = self.Vector.value{k} .* self.r(k) .* sin(self.theta(k)) + self.Position(2) + z(2, :)';
                end
            else
                k = self.Vector.index == index;
                z = [cos(self.theta(k)) -sin(self.theta(k)); sin(self.theta(k)) cos(self.theta(k))] * randn(2, numel(self.Vector.value{k})) .* self.Jitter .* self.r(k);
                xv = self.Vector.value{k} .* self.r(k) .* cos(self.theta(k)) + self.Position(1) + z(1, :)';
                yv = self.Vector.value{k} .* self.r(k) .* sin(self.theta(k)) + self.Position(2) + z(2, :)';
            end
        end
        
        function modify_colormap(self, CData, CDataMode)
            %MODIFY_COLORMAP  Update the colormap and way it is used.
            %
            % Syntax:
            %   r.modify_colormap(CData); % Only changes CData
            %   r.modify_colormap(CData, CDataMode); % Change both
            %   r.modify_colormap([], CDataMode); % Only changes CDataMode
            %
            % Inputs:
            %   CData - New nColors x 3 array to use for label color keys.
            %               Specify as [] to skip this input argument and
            %               continue to use the colormap assigned to the
            %               Rosette object.
            %   CDataMode - (Optional) "auto" or "fixed". If not specified,
            %                   the initial default setting is "auto" which
            %                   sets CData indices to be distributed
            %                   equally along the elements of the CData
            %                   rows. "fixed" sets it so that each new
            %                   vector is directly indexed to a new row of
            %                   CData.
            %
            % See also: Contents, Rosette, add_vector, modify_vector
            if isempty(CData)
                CData = self.CData;
            end
            self.CData = CData;
            if nargin > 1
                if ~ismember(CDataMode, ["auto", "fixed"])
                    error("CDataMode must be 'auto' or 'fixed' (value given was <strong>%s</strong>)\n", ...
                        CDataMode);
                end
                self.CDataMode = CDataMode; 
            end
            switch self.CDataMode
                case "auto"
                    self.c = self.CData(1:size(self.Vector, 1), :);
                case "fixed"
                    iC = round(linspace(1, size(self.CData, 1), size(self.Vector, 1)));
                    self.c = self.CData(iC, :);
                otherwise
                    error("Invalid value of CDataMode property: %s", self.CDataMode);
            end
            for ii = 1:size(self.Vector, 1)
                if (numel(self.basis) >= ii) && ~isa(self.basis(ii), 'matlab.graphics.GraphicsPlaceholder')
                    self.basis(ii).Color = self.c(ii, :);
                    self.label(ii).Color = self.c(ii, :);
                    self.cloud(ii).MarkerFaceColor = self.c(ii, :);
                end
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
            % See also: Contents, modify_cloud, text
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
        
        function modify_label_offset(self, offset)
            %MODIFY_LABEL_OFFSET  Change value of protected `LabelOffset` property
            %
            % Syntax:
            %   r.modify_label_offset(offset);
            %       % Applies same offset to all labels.
            %
            % Inputs:
            %   offset - The new value of LabelOffset property.
            %
            % See also: Contents, modify_cloud, text
            self.LabelOffset = offset;
            self.modify_vector();
        end
        
        function modify_cloud(self, varargin)
            %MODIFY_CLOUD  Change values of 'Name', value pairs for "value" ring
            %
            % Syntax:
            %   r.modify_cloud('Name', value, ...);
            %
            % Accepts any 'Name', value argument pairs that MATLAB builtin
            % `line` would accept.
            %
            % See also: Contents, line
            set(self.cloud, varargin{:});
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
                z = [cos(self.theta(k)) -sin(self.theta(k)); sin(self.theta(k)) cos(self.theta(k))] * randn(2, numel(value{ii})) .* self.Jitter .* self.r(k);
                xv = value{k} .* self.r(k) .* cos(self.theta(k)) + self.Position(1) + z(1,:)';
                yv = value{k} .* self.r(k) .* sin(self.theta(k)) + self.Position(2) + z(2,:)';
                self.cloud(k).XData = xv;
                self.cloud(k).YData = yv;
                self.Vector.value{k} = value{ii};
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
            if nargin < 2
                index = [];
                V = [];
            end
            
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
                self.label(k).Position(1:2) = [X(1) + self.LabelOffset * self.r(k) .* cos(self.theta(k)), ...
                                               Y(1) + self.LabelOffset * self.r(k) .* sin(self.theta(k))];
                [self.cloud(k).XData, self.cloud(k).YData] = self.get_scaled_value(index(ii));
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
            set(self.cloud, ...
                'XData', self.cloud.XData + xd, ...
                'YData', self.cloud.YData + yd);
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

