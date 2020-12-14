classdef (Sealed) MatrixRestriction
    % Matrix restriction with multiple representations.

    properties (GetAccess = private, SetAccess = private)
        data_;
    end

    methods
        function self = MatrixRestriction(Null)
            % TODO: Help.

            if isa(Null, 'MatrixRestriction')
                self.data_ = Null.data_;
            else
                self.data_ = containers.Map();
                if issparse(Null)
                    assert(all(nonzeros(Null) == 1), ...
                        'Expected all nonzeros to be one.');
                    assert(all(sum(Null, 1) == 1), ...
                        'Expected all columns to have exactly one entry.');
                    assert(all(sum(Null, 2) <= 1), ...
                        'Expected all rows to have at most one entry.');
                    [rows, ~] = find(Null);
                    assert(issorted(rows), ...
                        'Expected restriction matrix to preserve order.');
                    self.data_('size') = size(Null);
                    self.data_('Null') = Null;
                elseif islogical(Null)
                    assert(iscolumn(Null), ...
                        'Expected logical vector to be a column vector.');
                    self.data_('size') = [length(Null), sum(Null)];
                    self.data_('mask') = Null;
                else
                    error('Unsupported parameter type.');
                end
            end
        end

        function varargout = size(val, dim)
            % TODO: Help.

            sz = val.data_('size');

            if nargin == 1
                if nargout <= 1
                    varargout = {sz};
                else
                    varargout(1:2) = {sz(1), sz(2)};
                    [varargout{3:nargout}] = deal(1);
                end
            else
                if dim <= 2
                    varargout = {sz(dim)};
                else
                    varargout = {1};
                end
            end
        end

        function Null = getNull(val)
            % TODO: Help.

            if ~val.data_.isKey('Null')
                sz = val.data_('size');
                mask = val.data_('mask');
                new_Null = sparse(find(mask), 1:sz(2), 1, sz(1), sz(2));
                val.data_('Null') = new_Null;
            end

            Null = val.data_('Null');
        end

        function mask = getMask(val)
            % TODO: Help.

            if ~val.data_.isKey('mask')
                sz = val.data_('size');
                Null = val.data_('Null');
                [rows, ~] = find(Null);
                new_mask = false(sz(1), 1);
                new_mask(rows) = true;
                val.data_('mask') = new_mask;
            end

            mask = val.data_('mask');
        end
    end
end
