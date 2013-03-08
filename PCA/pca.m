% Octave script for computing PCA using an input covariance matrix
% The covariance matrix should be an ASCII file.
%
% Usage:
%     octave -q pca.m  cov_mat_file
% The -q is optional, but it supresses annoying output from Octave.
% Note that octave is an open source matlab clone, and this could
% can work in Matlab with little, if any, modification.
%
% Try running:
%    octave -q pca.m --help
% for more info.

% Constants -------------------------------------------------------------------
EPS_PRINT_ARGUMENTS = '-depsc';  % Args to print function for saving .eps images.
LINE_WIDTH           = 5;
TITLE_FONT_SZ        = 22;
AXES_LABEL_FONT_SIZE = 18;



% Parse command line arguments ------------------------------------------------
args = argv();
if nargin != 0 && strcmp('--help', args{1})
    printf('Usage\n');
    printf('\toctave -q pca.m cov_matrix_file_name\n', program_name());
    printf('\t\tcov_matrix_file_name = Path to ASCII covariance matrix.\n');
    printf('Will output files to disk like <cov_matrix_file_name>_XXX.eps\n');\
    exit(0)
end
mat_file_name = '';
if nargin >= 1
    mat_file_name = args{1};
else
    printf('See %s --help\n', program_name());
    exit(1)
end

% Load the covariance matrix --------------------------------------------------
printf('Loading covariance matrix from %s.\n', mat_file_name);
cov_mat = load(mat_file_name);

% Generate heatmap ------------------------------------------------------------
printf('Generating heatmap...\n');
heatmap_out = sprintf('%s_heatmap.eps', mat_file_name);
imagesc(cov_mat);
t = title('Covariance Matrix Heatmap');
set(t,  'FontSize' , TITLE_FONT_SZ);
colorbar;
print(heatmap_out, EPS_PRINT_ARGUMENTS);


% Perform PCA -----------------------------------------------------------------
[U,D,pc] = svd(cov_mat);
sv = diag(D);
h = plot(sv);

% Plot eigenvalue magnitude ---------------------------------------------------
figure; %New figure
printf('Generating eigenvalue magnitude plot...\n');
eig_mag_out = sprintf('%s_eig_mag.eps', mat_file_name);

figure; %New figure

%sqrt of singular values = eigenvalue mag.
h = plot(sqrt(sv));
grid on;
t = title('Plot of eigenvalue magnitude resulting from PCA');
ax = xlabel('Eigenvalue Index');
ay = ylabel('Eigenvalue Magnitude');
set(ax, 'FontSize' , AXES_LABEL_FONT_SIZE);
set(ay, 'FontSize' , AXES_LABEL_FONT_SIZE);
set(h,  'linewidth', LINE_WIDTH);
set(t,  'FontSize' , TITLE_FONT_SZ);

print(eig_mag_out, EPS_PRINT_ARGUMENTS);
printf('Done.\n');
exit(0)
