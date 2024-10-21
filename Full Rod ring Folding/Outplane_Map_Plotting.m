%%
%draw the inplane and outplane map
load('Rod_Origin.mat')
load('Rod_Inplane.mat')
load('Rod_Outplane.mat')

X0=Q_sec_M0(1,1);
disp_x=[0:5:X0,X0];% sampling displacement x
N_x=length(disp_x)-1

X_1=[]
Y_2=[]
Z_1=[]
Z_1_flag=[]
for mm=1:N_x+1
    X_1=[X_1;reshape(Q_sec_M2{mm}(1,1,:),[],1)];
    Z_1=[Z_1;reshape(Q_sec_M2{mm}(3,1,:),[],1)];
    Y_2=[Y_2;reshape(Q_sec_M2{mm}(2,end,:),[],1)];
    
end
Z_1_flag= abs(Z_1)<5
figure(1)
scatter3(X_1,Y_2,Z_1,[],100*Z_1_flag)
hold on
scatter3(Y_2,X_1,Z_1,[],100*Z_1_flag)
colorbar
caxis([0 100])

figure(2)
surfir(X_1,Y_2,100*Z_1_flag,0.5)
hold on
surfir(Y_2,X_1,100*Z_1_flag,0.5)
hold on
colorbar
caxis([0 100])
axis equal
view([0 0 1])