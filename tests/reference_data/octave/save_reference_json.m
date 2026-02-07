function save_reference_json(filepath, data)
% SAVE_REFERENCE_JSON Save reference data to a JSON file with complex number support.
%   Complex numbers are stored as [real, imag] pairs.
%   Matrices are stored as arrays of row arrays.
%
%   save_reference_json(filepath, data)
%
%   data is a struct with fields:
%     .metadata  - struct with string/numeric metadata
%     .matrices  - struct of named complex matrices
%     .vectors   - struct of named complex vectors
%     .scalars   - struct of named scalar values

    fid = fopen(filepath, 'w');
    if fid == -1
        error('Cannot open file: %s', filepath);
    end

    fprintf(fid, '{\n');

    % Write metadata
    fprintf(fid, '  "metadata": ');
    write_value(fid, data.metadata, 2);
    fprintf(fid, ',\n');

    sections = {};
    if isfield(data, 'matrices')
        sections{end+1} = 'matrices';
    end
    if isfield(data, 'vectors')
        sections{end+1} = 'vectors';
    end
    if isfield(data, 'scalars')
        sections{end+1} = 'scalars';
    end

    for si = 1:length(sections)
        sec = sections{si};
        fprintf(fid, '  "%s": {\n', sec);
        names = fieldnames(data.(sec));
        for i = 1:length(names)
            name = names{i};
            val = data.(sec).(name);
            fprintf(fid, '    "%s": ', name);

            if strcmp(sec, 'matrices')
                write_complex_matrix(fid, val);
            elseif strcmp(sec, 'vectors')
                write_complex_vector(fid, val);
            else
                write_value(fid, val, 4);
            end

            if i < length(names)
                fprintf(fid, ',');
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '  }');
        if si < length(sections)
            fprintf(fid, ',');
        end
        fprintf(fid, '\n');
    end

    fprintf(fid, '}\n');
    fclose(fid);
end

function write_complex_matrix(fid, M)
    [rows, cols] = size(M);
    fprintf(fid, '[\n');
    for r = 1:rows
        fprintf(fid, '      [');
        for c = 1:cols
            write_complex_number(fid, M(r, c));
            if c < cols
                fprintf(fid, ', ');
            end
        end
        fprintf(fid, ']');
        if r < rows
            fprintf(fid, ',');
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '    ]');
end

function write_complex_vector(fid, v)
    n = length(v);
    fprintf(fid, '[');
    for i = 1:n
        write_complex_number(fid, v(i));
        if i < n
            fprintf(fid, ', ');
        end
    end
    fprintf(fid, ']');
end

function write_complex_number(fid, z)
    fprintf(fid, '[%.17g, %.17g]', real(z), imag(z));
end

function write_value(fid, val, indent)
    prefix = repmat(' ', 1, indent);
    if isstruct(val)
        fprintf(fid, '{\n');
        names = fieldnames(val);
        for i = 1:length(names)
            fprintf(fid, '%s  "%s": ', prefix, names{i});
            write_value(fid, val.(names{i}), indent + 2);
            if i < length(names)
                fprintf(fid, ',');
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '%s}', prefix);
    elseif ischar(val) || isstring(val)
        fprintf(fid, '"%s"', val);
    elseif isnumeric(val) && isscalar(val)
        if val == floor(val) && abs(val) < 1e15
            fprintf(fid, '%d', val);
        else
            fprintf(fid, '%.17g', val);
        end
    elseif isnumeric(val) && isvector(val)
        fprintf(fid, '[');
        for i = 1:length(val)
            if val(i) == floor(val(i)) && abs(val(i)) < 1e15
                fprintf(fid, '%d', val(i));
            else
                fprintf(fid, '%.17g', val(i));
            end
            if i < length(val)
                fprintf(fid, ', ');
            end
        end
        fprintf(fid, ']');
    elseif islogical(val)
        if val
            fprintf(fid, 'true');
        else
            fprintf(fid, 'false');
        end
    else
        fprintf(fid, 'null');
    end
end
