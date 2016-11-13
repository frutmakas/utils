function [r] = uni_rand(seed)
seed=mod((25733*seed+13849),65536);
r=(seed*1.0)/65536.0;
