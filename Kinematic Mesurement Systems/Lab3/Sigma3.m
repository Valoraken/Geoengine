function [p, b] = Sigma3(a)
% [p, a] = Sigma3(a)
% SIGMA3 aims to delete gross errors of a according to 3 sigma principle
% Input:
% a:       array                      (column vector)      [length(a)*1]
% Output:
% p:       position of gross errors   (column vector)      [length(p)*1]
% b:       updated array              (column vector)      [(length(a)-length(p))*1]
% Name: YI HONG
% Matriculation number:3294211
m = mean(a);
s = std(a);
p = find(abs(a - m)>3*s);
a(p) = [];
b = a;
end