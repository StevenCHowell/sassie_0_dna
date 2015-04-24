% analyze the tetramer angles (want to plot them)
ncp1_axes = loadxy('gH5c11_r_n1_axes.txt');
ncp2_axes = loadxy('gH5c11_r_n2_axes.txt');
origin = loadxy('gH5c11_r_orig.txt');
h = loadxy('gH5c11_r_h.txt');

figure
hold all
axis equal
r = 70;
for i=1:3
    xyplot3(cat(1, origin(1,:), r*ncp1_axes(i,:)+origin(1,:)), 'linewidth', 2)
end
hold off
hold all
for i=1:3
    %     xyplot3(cat(1, origin(2,:), r*ncp2_axes(i,:)+origin(2,:)), 'linewidth', 2)
    offset = origin(2,:) - h(1) * ncp1_axes(3,:);
    xyplot3(cat(1, offset, r*ncp2_axes(i,:)+offset), 'linewidth', 2)
end

legend('x1','y1','z1','x2','y2','z2')
% set(gca, 'plotBoxAspectRatio', [1 1 1])