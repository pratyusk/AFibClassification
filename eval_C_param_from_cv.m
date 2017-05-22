
% best C values from 6 trials
C = [10 2 6 4 3 5 10 5 6 3
4,3,5,3,2,4,2,2,2,8
8 3 5 2 5 2 4 7 3 11
9 4 2 3 2 7 6 10 7 4
7 2 2 7 2 3 2 10 7 6
6,6,5,3,2,4,5,3,2,11];

for i = 1:11
    i
    C_best_accur(i) = length(find(C==i))
end
