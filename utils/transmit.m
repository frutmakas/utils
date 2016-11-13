% [y, h_est, sigma_h_chapeau] = transmit(X, h, W, noise, sigma_h, sigma_h_cross)
function [y, h_est, sigma_h_chapeau] = transmit(X, h, W, noise, sigma_h, sigma_h_cross)

y = X*W*h+noise;
h_est = inv(W'*X'*X*W+sigma_h_cross)*W'*X'*y;
sigma_h_chapeau = sigma_h - inv(W'*X'*X*W+sigma_h_cross)*W'*X'*X*W*sigma_h;

