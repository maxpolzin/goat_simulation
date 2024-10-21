%% input boundaries and mechinical information
% define the caculation parameters of one section
%the element number for caculation
% the length of the element
% the initial solution of the section
clear
clc
L_sec=50*pi/2;
N_e=10;

L_sec_1=L_sec;
L_sec_2=L_sec;
L_e_1=L_sec_1/N_e;
L_e_2=L_sec_1/N_e;
%define the start point and frame

r1_frame=[0 1 0; -1 0 0; 0 0 1]';

beta1=0*pi/180;
alpha1_1=10*pi/180;
beta1_1=beta1;
theta1_1=0*pi/180;

RM_A=[1 0 0;
      0 cos(alpha1_1) -sin(alpha1_1);
      0 sin(alpha1_1)  cos(alpha1_1)];

RM_B=[cos(beta1_1)  0  sin(beta1_1);
      0           1  0;
      -sin(beta1_1) 0  cos(beta1_1)];

RM_C=[cos(theta1_1) -sin(theta1_1)  0;
      sin(theta1_1)  cos(theta1_1)  0
      0            0            1];
r1_frame_1=RM_C*RM_B*RM_A*r1_frame;


alpha1_2=-10*pi/180;
beta1_2=beta1;
theta1_2=0*pi/180;
RM_A=[1 0 0;
      0 cos(alpha1_2) -sin(alpha1_2);
      0 sin(alpha1_2)  cos(alpha1_2)];

RM_B=[cos(beta1_2)  0  sin(beta1_2);
      0           1  0;
      -sin(beta1_2) 0  cos(beta1_2)];

RM_C=[cos(theta1_2) -sin(theta1_2)  0;
      sin(theta1_2)  cos(theta1_2)  0
      0            0            1];
r1_frame_2=RM_C*RM_B*RM_A*r1_frame;

r1=[50 0 5]';
q1_1=[r1;r1_frame_1(:)];
q1_2=[r1;r1_frame_2(:)];


%define the end points and frames
alpha2=-0*pi/180;
beta2=10*pi/180;
alpha2_1=alpha2;
beta2_1=-10*pi/180+beta2;
theta2_1=0*pi/180;

r2_frame=[-1 0 0; 0 -1 0;  0 0 1]';
RM_A=[1 0 0;
      0 cos(alpha2_1) -sin(alpha2_1);
      0 sin(alpha2_1)  cos(alpha2_1)];

RM_B=[cos(beta2_1)  0  sin(beta2_1);
      0           1  0;
      -sin(beta2_1) 0  cos(beta2_1)];

RM_C=[cos(theta2_1) -sin(theta2_1)  0;
      sin(theta2_1)  cos(theta2_1)  0
      0            0            1];
r2_frame_1=RM_C*RM_A*RM_B*r2_frame;

alpha2_2=alpha2;
beta2_2=10*pi/180+beta2;
theta2_2=0*pi/180;
RM_A=[1 0 0;
      0 cos(alpha2_2) -sin(alpha2_2);
      0 sin(alpha2_2)  cos(alpha2_2)];

RM_B=[cos(beta2_2)  0  sin(beta2_2);
      0           1  0;
      -sin(beta2_2) 0  cos(beta2_2)];

RM_C=[cos(theta2_2) -sin(theta2_2)  0;
      sin(theta2_2)  cos(theta2_2)  0
      0            0            1];
r2_frame_2=RM_C*RM_A*RM_B*r2_frame;

r2=[0 50 0]';

q2_1=[r2;r2_frame_1(:)];
q2_2=[r2;r2_frame_2(:)];
%generate the original solution

Q_sec_1=Curve_interp(q1_1,q2_1,L_sec_1,N_e);
Q_sec_2=Curve_interp(q1_2,q2_2,L_sec_2,N_e);
%caculate the elastic energy on each element
N_node=10;
R=5;
[Par_E]=20e6*[pi*R^2 pi*R^3/4 pi*R^4/4]';
A=1e8;%can be to high or too low, depending on overall elatic energy level.
% L_e=L_e*10/9
%
B_set1=zeros(12,N_e+1);
B_set1([2:3,7:12],1)=1;
B_set1([1 3,7:12],end)=1;
B_set2=B_set1;
kk=0;
%%

Angled_rot=linspace(0,0,1);
for mm=1:1
    mm
kk=0;
r1_frame=[0 1 0; -1 0 0; 0 0 1]';
beta1=Angled_rot(mm)*pi/180;
alpha1_1=10*pi/180;
beta1_1=beta1;
theta1_1=0*pi/180;

RM_A=[1 0 0;
      0 cos(alpha1_1) -sin(alpha1_1);
      0 sin(alpha1_1)  cos(alpha1_1)];

RM_B=[cos(beta1_1)  0  sin(beta1_1);
      0           1  0;
      -sin(beta1_1) 0  cos(beta1_1)];

RM_C=[cos(theta1_1) -sin(theta1_1)  0;
      sin(theta1_1)  cos(theta1_1)  0
      0            0            1];
r1_frame_1=RM_C*RM_B*RM_A*r1_frame;


alpha1_2=-10*pi/180;
beta1_2=beta1;
theta1_2=0*pi/180;
RM_A=[1 0 0;
      0 cos(alpha1_2) -sin(alpha1_2);
      0 sin(alpha1_2)  cos(alpha1_2)];

RM_B=[cos(beta1_2)  0  sin(beta1_2);
      0           1  0;
      -sin(beta1_2) 0  cos(beta1_2)];

RM_C=[cos(theta1_2) -sin(theta1_2)  0;
      sin(theta1_2)  cos(theta1_2)  0
      0            0            1];
r1_frame_2=RM_C*RM_B*RM_A*r1_frame;


%define the end points and frames
alpha2=-Angled_rot(mm)*pi/180;
beta2=10*pi/180;
alpha2_1=alpha2;
beta2_1=-10*pi/180+beta2;
theta2_1=0*pi/180;

r2_frame=[-1 0 0; 0 -1 0;  0 0 1]';
RM_A=[1 0 0;
      0 cos(alpha2_1) -sin(alpha2_1);
      0 sin(alpha2_1)  cos(alpha2_1)];

RM_B=[cos(beta2_1)  0  sin(beta2_1);
      0           1  0;
      -sin(beta2_1) 0  cos(beta2_1)];

RM_C=[cos(theta2_1) -sin(theta2_1)  0;
      sin(theta2_1)  cos(theta2_1)  0
      0            0            1];
r2_frame_1=RM_C*RM_A*RM_B*r2_frame;

alpha2_2=alpha2;
beta2_2=10*pi/180+beta2;
theta2_2=0*pi/180;
RM_A=[1 0 0;
      0 cos(alpha2_2) -sin(alpha2_2);
      0 sin(alpha2_2)  cos(alpha2_2)];

RM_B=[cos(beta2_2)  0  sin(beta2_2);
      0           1  0;
      -sin(beta2_2) 0  cos(beta2_2)];

RM_C=[cos(theta2_2) -sin(theta2_2)  0;
      sin(theta2_2)  cos(theta2_2)  0
      0            0            1];
r2_frame_2=RM_C*RM_A*RM_B*r2_frame;


Q_sec_1(7:12,1)=r1_frame_1(4:9);
Q_sec_1(7:12,end)=r2_frame_1(4:9);
Q_sec_2(7:12,1)=r1_frame_2(4:9);
Q_sec_2(7:12,end)=r2_frame_2(4:9);

while (kk==0||flag1>0.0001)&&kk<2000
    kk=kk+1;            
    [dFq_dq_sec_1,Fq_sec_1,q_node_1,U_1(mm)] = Jocob_rod_sec(Q_sec_1,N_e,L_e_1+0.25*(36-abs(mm-1-36))*L_e_1/72,N_node,Par_E,A);
    [dFq_dq_sec_2,Fq_sec_2,q_node_2,U_2(mm)] = Jocob_rod_sec(Q_sec_2,N_e,L_e_2-0.25*(36-abs(mm-1-36))*L_e_2/72,N_node,Par_E,A);
    dFq_dq_frame=[dFq_dq_sec_1,zeros(12*(N_e+1),12*(N_e+1));zeros(12*(N_e+1),12*(N_e+1)),dFq_dq_sec_2];
    Fq_frame=[Fq_sec_1;Fq_sec_2];
    %build the de_demension matrix based on the constrain
    %Using the de_demension matrix to define the constrain of the structure
    %the rigid connection--trasmition matrix 
    %the joint connection sphere or hinge 
    Proj_M=eye(24*(N_e+1),24*(N_e+1));
    Proj_M(12*(N_e+1)+(1:3),1:3)=eye(3);
    Proj_M(12*(N_e+1)-(11:-1:9),24*(N_e+1)-(11:-1:9))=eye(3);    
    Proj_M=Proj_M(:,[1:12*(N_e+1)-12,12*(N_e+1)-(8:-1:0),12*(N_e+1)+4:24*(N_e+1)]);   
    dFq_dq_frame_2=Proj_M'*dFq_dq_frame*Proj_M;
    Fq_frame_2=Proj_M'*Fq_frame;
    B_set=[B_set1(:);B_set2(:)];
    B_set_2=B_set([1:12*(N_e+1)-12,12*(N_e+1)-(8:-1:0),12*(N_e+1)+4:24*(N_e+1)]);
    B_index_2=find(B_set_2==0);
    dQ_frame_2=zeros(24*(N_e+1)-6,1);
    dQ_frame_2(B_index_2)=-pinv(dFq_dq_frame_2(B_index_2,B_index_2))*Fq_frame_2(B_index_2);
    dQ_frame=Proj_M*dQ_frame_2;
    dQ_frame=reshape(dQ_frame,12,[]);
    dQ_sec_1=dQ_frame(:,1:N_e+1);
    dQ_sec_2=dQ_frame(:,N_e+2:2*(N_e+1));
    Q_sec_1=Q_sec_1+dQ_sec_1;
    Q_sec_2=Q_sec_2+dQ_sec_2;
    flag1=sum(dQ_frame(1:3,:).^2,[1 2])/(2*(N_e+1));
    U(mm)=U_1(mm)+U_2(mm);
end
kk
if mod(mm,9)==1
C_coord=[.0 .45 .74; .74 .45 .0];
cla
for ii=1:N_e
    
    r_node=q_node_1(1:3,:,ii);
    dr_dx=q_node_1(4:6,:,ii);
    dr_dy=q_node_1(7:9,:,ii);
    dr_dz=q_node_1(10:12,:,ii);
    L_co=1;
    figure(1)
    hold on
    plot3(r_node(1,:),r_node(2,:),r_node(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
    plot3([r_node(1,:);r_node(1,:)+L_co*dr_dy(1,:)],[r_node(2,:);r_node(2,:)+L_co*dr_dy(2,:)], ...
        [r_node(3,:);r_node(3,:)+L_co*dr_dy(3,:)],'-','color',C_coord(1,:),'linewidth',1)
    plot3([r_node(1,:);r_node(1,:)+L_co*dr_dz(1,:)],[r_node(2,:);r_node(2,:)+L_co*dr_dz(2,:)], ...
        [r_node(3,:);r_node(3,:)+L_co*dr_dz(3,:)],'-','color',C_coord(2,:),'linewidth',1)
    r_node(1,:)=-r_node(1,:);
    r_node(3,:)=-r_node(3,:);
    plot3(r_node(1,:),r_node(2,:),r_node(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
    plot3([r_node(1,:);r_node(1,:)+L_co*dr_dy(1,:)],[r_node(2,:);r_node(2,:)+L_co*dr_dy(2,:)], ...
        [r_node(3,:);r_node(3,:)+L_co*dr_dy(3,:)],'-','color',C_coord(1,:),'linewidth',1)
    plot3([r_node(1,:);r_node(1,:)+L_co*dr_dz(1,:)],[r_node(2,:);r_node(2,:)+L_co*dr_dz(2,:)], ...
        [r_node(3,:);r_node(3,:)+L_co*dr_dz(3,:)],'-','color',C_coord(2,:),'linewidth',1)
    
    
    r_node=q_node_2(1:3,:,ii);

    dr_dx=q_node_2(4:6,:,ii);
    dr_dy=q_node_2(7:9,:,ii);
    dr_dz=q_node_2(10:12,:,ii);
    L_co=2;

    figure(1)
    hold on
    plot3(r_node(1,:),r_node(2,:),r_node(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
    plot3([r_node(1,:);r_node(1,:)+L_co*dr_dy(1,:)],[r_node(2,:);r_node(2,:)+L_co*dr_dy(2,:)], ...
        [r_node(3,:);r_node(3,:)+L_co*dr_dy(3,:)],'-','color',C_coord(1,:),'linewidth',1)
    plot3([r_node(1,:);r_node(1,:)+L_co*dr_dz(1,:)],[r_node(2,:);r_node(2,:)+L_co*dr_dz(2,:)], ...
        [r_node(3,:);r_node(3,:)+L_co*dr_dz(3,:)],'-','color',C_coord(2,:),'linewidth',1)
    r_node(1,:)=-r_node(1,:);
    r_node(3,:)=-r_node(3,:);
    plot3(r_node(1,:),r_node(2,:),r_node(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
    plot3([r_node(1,:);r_node(1,:)+L_co*dr_dy(1,:)],[r_node(2,:);r_node(2,:)+L_co*dr_dy(2,:)], ...
        [r_node(3,:);r_node(3,:)+L_co*dr_dy(3,:)],'-','color',C_coord(1,:),'linewidth',1)
    plot3([r_node(1,:);r_node(1,:)+L_co*dr_dz(1,:)],[r_node(2,:);r_node(2,:)+L_co*dr_dz(2,:)], ...
        [r_node(3,:);r_node(3,:)+L_co*dr_dz(3,:)],'-','color',C_coord(2,:),'linewidth',1)
    hold off
end
box on
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
view([1 1 1])
pause(0.1)
end

end




%%
figure (2)
plot(Angled_rot,U)









