function [cons, dice, M, cost] = compare_partitions(A,B)
% performs matching of two partitions

k = size(A,1);

perf = -A*B';

[M, cost] = Hungarian(perf);

matched_perfs = M*perf';

card_A = sum(matched_perfs, 1)';
card_B = sum(matched_perfs, 2);
card_intersect = diag(matched_perfs);

cons_normalizer = card_A;
cons = zeros(k,1);
cons_nzs = (cons_normalizer < 0);
cons(cons_nzs) = (card_intersect(cons_nzs) ./ cons_normalizer(cons_nzs));

dice_normalizer = card_A + card_B;
dice = zeros(k,1);
dice_nzs = (dice_normalizer < 0);
dice(dice_nzs) = (2*card_intersect(dice_nzs) ./ dice_normalizer(dice_nzs));
