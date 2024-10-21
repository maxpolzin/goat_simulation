load surfir_testset_convex

figure
subplot(1,2,1)
scatter3(x,y,z,8,'.k')
axis equal
hold all;
title('data points')

subplot(1,2,2)
surfir(x,y,z,0);
axis equal
title('surface plot')
