% test_all_util.m

if 1 % test jf_protected_names methods
    pn = jf_protected_names;
    pn.prctile('test');
end

% this list requires mex files so test it second
list2 = { ...
'bspline_1d_coef test', ...
'bspline_1d_interp test', ...
'bspline_1d_synth test'
};

% test these basic things first
list1 = { ...
'arg_pair test', ...
'cpu test', ...
'dft_sym_check test', ...
'downsample2 test', ...
'downsample3 test', ...
'embed test', ...
'fld_write test', ...
'fractional_delay test', ...
'fwhm1 test', ...
'fwhm2 test', ...
'fwhm_match test', ...
'hist_equal test', ...
'ifft_sym test', ...
'interp1_jump test', ...
'interp1_lagrange test', ...
'ir_apply_tridiag_inv test', ...
'ir_conv test', ...
'ir_idct2 test', ...
'ir_dwt_filters test', ...
'ir_im2col test', ...
'ir_interpft test', ...
'ir_odwt1 test', ...
'ir_odwt2 test', ...
'ir_pad_into_center test', ...
'ir_patch_avg test', ...
'ir_poly2_fun test', ...
'ir_project_k_sparse test', ...
'kde_pmf1 test', ...
'kde_pmf2 test', ...
'kde_pmf_width test', ...
'lloyd_max_hist test', ...
'masker test', ...
'min_cos_quad test', ...
'padn test', ...
'poisson test', ...
'reale test', ...
'reshaper test', ...
'sinc_periodic test', ...
'strum test', ...
'unpadn test', ...
'upsample_rep test', ...
'vararg_pair test'
};

% run_mfile_local(list1, 'pause', false)
% run_mfile_local(list2)
run_mfile_local({list1{:}, list2{:}}, 'pause', 0, 'abort', 1);

% octave todo:
% lloyd_max_hist test
% min_cos_quad test
% bspline* ... mex (list1)
