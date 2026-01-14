
function result = HHH(n)
    result = [];
    for i = 1:2^(n*n)
        array = dec2bin(i-1, n*n) - '0';
        matrix = reshape(array, n, n);
                result(:, :, end+1) = matrix; %#ok<AGROW>
        doStuff(matrix);
    end
end


function doStuff(matrix)
    matrix;
end