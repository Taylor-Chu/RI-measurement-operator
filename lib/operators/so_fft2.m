function X = so_fft2(x, No, scale)
% Computes and oversampled and scaled FFT2 with the scale factors
% precomputed from the ``nufft_init`` function.
%
% Parameters
% ----------
% x : double[:, :]
%     Input image.
% No : int[2]
%     Overscale size.
% scale : double[:, :]
%     Scale parameters precomputed by ``nufft_init``.
%
% Returns
% -------
% X : complex[:]
%     2D FFT coefficients.
%

% if x is 1D, reshape to square image 
if size(x, 2) == 1
    N = sqrt(numel(x));
    x = reshape(x, [N, N]);
end

%%
% apply scaling factors
x = x .* scale;

% oversampled FFT
X = fft2(x, No(1), No(2));
X = X(:);

end
