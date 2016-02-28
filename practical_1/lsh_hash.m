function [ hash ] = lsh_hash( x, b, T )
% x = [n x p] matrix of p-dimensional points.
% b = number of bits used for hashing.
% T = number of hash tables to build.
%
% hash = cell array of length T. Each entry contains a struct with the 
% following fields
%   - keys = a vector where each entry represents the decimal representation
%   of the bitstring key. 
%   - buckets = a cell array of length length(key), where each entry
%   contains all points that have the corresponding key. 

hash = cell(T,1);
for i = 1 : T
    P = randn(size(x,2), b);
    H = (sign(x*P) + 1)/2;
    d = bin2dec(num2str(H));

    hash{i}.keys = unique(d);
    hash{i}.buckets = cell(length(hash.keys),1);
    for j = 1 : length(hash.keys)
        hash{i}.buckets{j} = x(d == j,:);
    end
end





end

