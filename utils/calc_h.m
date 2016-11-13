%[new_h, new_sigma]=calc_h(y ,X, W, noise_variance, sigmaH, sigmaH_cross)
function [new_h, new_sigma]=calc_h(y ,X, W, noise_variance, sigmaH, sigmaH_cross)

clear new_h new_sigma common_term;

common_term = inv(W'*X'*noise_variance*X*W+sigmaH_cross)*W'*X'*noise_variance;
new_h = common_term*y;
new_sigma = common_term*X*W*sigmaH;
