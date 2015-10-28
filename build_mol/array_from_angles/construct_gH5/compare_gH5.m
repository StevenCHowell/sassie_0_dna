
%loadxy('1zbb_tetra.pdb.dat');
t1zbb.raw = loadxy('1zbb_tetra_combined.pdb.dat');
c11.raw = loadxy('c11_r.pdb.dat');
gh5_initial.raw = loadxy('manual_gH5x4.pdb.dat');
gh5_protein.raw = loadxy('gH5x4_allProteins.pdb.dat');
gh5_complete.raw = loadxy('complete_gH5x4.pdb.dat');

% normalize the data to 1
% t1zbb.n = [t1zbb.raw(:,1), t1zbb.raw(:,2:end)/t1zbb.raw(1,2)];
% c11.n = [c11.raw(:,1), c11.raw(:,2:end)/c11.raw(1,2)];
% gh5_initial.n = [gh5_initial.raw(:,1), gh5_initial.raw(:,2:end)/gh5_initial.raw(1,2)];
% gh5_protein.n = [gh5_protein.raw(:,1), gh5_protein.raw(:,2:end)/gh5_protein.raw(1,2)];
% gh5_complete.n = [gh5_complete.raw(:,1), gh5_complete.raw(:,2:end)/gh5_complete.raw(1,2)];

%data = loadxy('/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/1406data/chess/iqdata/c000_4x167_h5_k010.i0q');
% data.raw = loadxy('/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/1406data/chess/iqdata/c000_4x167_h5_mg1.i0q');
% data.iq = [c11.raw(:,1), interp1(data.raw(:,1), data.raw(:,2), c11.raw(:,1)), interp1(data.raw(:,1), data.raw(:,3), c11.raw(:,1))];
% match_type = 'all';
% t1zbb.iq = scale_offset(t1zbb.n, data.iq, match_type);
% c11.iq = scale_offset(c11.n, data.iq, match_type);
% gh5_initial.iq = scale_offset(gh5_initial.n, data.iq, match_type);
% gh5_protein.iq = scale_offset(gh5_protein.n, data.iq, match_type);
% gh5_complete.iq = scale_offset(gh5_complete.n, data.iq, match_type);

data.raw = loadxy('/home/schowell/Dropbox/gw_phd/paper_tetranucleosome/1406data/iqdata/coarse_grid/c000_4x167_h5_mg1_33.wi0');
data.iq = data.raw;
i0 = data.iq(1,2);

% normalize the data to I(0) from the experimental data
t1zbb.iq        = [t1zbb.raw(:,1),        i0 * t1zbb.raw(:,2:end)/t1zbb.raw(1,2)];
c11.iq          = [c11.raw(:,1),          i0 * c11.raw(:,2:end)/c11.raw(1,2)];
gh5_initial.iq  = [gh5_initial.raw(:,1),  i0 * gh5_initial.raw(:,2:end)/gh5_initial.raw(1,2)];
gh5_protein.iq  = [gh5_protein.raw(:,1),  i0 * gh5_protein.raw(:,2:end)/gh5_protein.raw(1,2)];
gh5_complete.iq = [gh5_complete.raw(:,1), i0 * gh5_complete.raw(:,2:end)/gh5_complete.raw(1,2)];

figure;
hold all
xyerror(data.iq, 's');
xyplot(t1zbb.iq);
xyplot(c11.iq);
xyplot(gh5_initial.iq);
xyplot(gh5_complete.iq);

% legend('1. Exp 4x167 gH5', '2. PDB:1ZBB', '3. All DNA & Proteins added to 2', '4. 4 gH5 NCPs', '5. All Proteins added to 4', '6. All DNA added to 6', 'location', 'southwest');
lh = legend('1. Exp 4x167 gH5 in 1 mM MgCl_2', '2. PDBID:1ZBB', '3. Missing residues added to 2', '4. Positioned gH5 NCPs', '5. Missing residues added to 4', 'location', 'southwest');
legend boxoff

logxy;
axis tight
zoomout(0.1);
xlim([0.0055 0.21])
shift_underscore_legend_entries(lh)
iqlabel
% set(ax,'YTickLabel',[])

saveps(gcf, 'compare_gH5.eps')
savepng(gcf, 'compare_gH5.png')
