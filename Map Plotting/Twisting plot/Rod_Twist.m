%% input boundaries and mechinical information
% define the caculation parameters of one section
%the element number for caculation
% the length of the element
% the initial solution of the section
clear
clc
L_sec=50*pi/2;
N_e=8;

L_sec_1=L_sec;
L_sec_2=L_sec;
L_e=L_sec_1/N_e;
%define the start point and frame

r1_frame=[0 1 0; -1 0 0; 0 0 1]';


alpha1=0*pi/180;
beta1=0*pi/180;
theta1=0*pi/180;

RM_A=[1 0 0;
      0 cos(alpha1) -sin(alpha1);
      0 sin(alpha1)  cos(alpha1)];

RM_B=[cos(beta1)  0  sin(beta1);
      0           1  0;
      -sin(beta1) 0  cos(beta1)];

RM_C=[cos(theta1) -sin(theta1)  0;
      sin(theta1)  cos(theta1)  0
      0            0            1];
r1_frame=RM_C*RM_B*RM_A*r1_frame;

r1=[50 0 0]';
q1=[r1;r1_frame(:)];



%define the end points and frames
alpha2=0*pi/180;
beta2=0*pi/180;
theta2=0*pi/180;

r2_frame=[-1 0 0; 0 -1 0;  0 0 1]';
RM_A=[1 0 0;
      0 cos(alpha2) -sin(alpha2);
      0 sin(alpha2)  cos(alpha2)];

RM_B=[cos(beta2)  0  sin(beta2);
      0           1  0;
      -sin(beta2) 0  cos(beta2)];

RM_C=[cos(theta2) -sin(theta2)  0;
      sin(theta2)  cos(theta2)  0
      0            0            1];
r2_frame=RM_C*RM_A*RM_B*r2_frame;

r2=[0 50 0]';

q2=[r2;r2_frame(:)];

%generate the original solution
Q_sec=Curve_interp(q1,q2,L_sec_1,N_e);

%caculate the elastic energy on each element
N_node=5;
R=5;
E=20e6
G=E/2*(1+0.4)
[Par_E]=[E*pi*R^2, E*pi*R^4/64 0.2*G*pi*R^4/32]';
A=0.1e8;%can be to high or too low, depending on overall elatic energy level.
% L_e=L_e*10/9
%
B_set=zeros(12,N_e+1);
B_set([2:3,7:12],1)=1;
B_set([1 3,7:12],end)=1;

kk=0;
%%
%boundary conditions
%BC1 symetrical compressing
% r1
% position:
% control: X
% Fixed: Y Z
% Rostation
% Fixed: X Y Z

% r2
% position:
% free: Y
% Fixed: X Z
% Rostation
% Fixed: X Y Z
%build the de_demension matrix based on the constrain
%Using the de_demension matrix to define the constrain of the structure
Proj_M=eye(12*(N_e+1),12*(N_e+1));%dedemesionalize the variables

%BC2 Antisymmetric twisting
% r1 start point
% position:
% control: X
% Fixed: Y Z
% Rotation
% Free: X axis (how to constrain this)
% r1_frame(2,1)-r1_frame(1,2)=0
% r1_frame(2,3)-r1_frame(3,2)=0
% considering the orthogonality, the above equals to:
% r1_frame(2,1)=0
% r1_frame(1,2)=0
% r1_frame(2,3)=0
% r1_frame(3,2)=0
D_fix=0*Q_sec;%fixed dimensions/varibales
D_fix(1,1)=1;%control the X disp
D_fix(2:3,1)=1;%fix the Y Z disp

D_fix(4,1)=1;%fix the rotation
D_fix(8:9,1)=1;%fix the rotation
D_fix(10,1)=1;%fix the rotation
%
% r2
% position:
% Free: Y
% Fixed: X Z
% Rotation
% Free: Y axis (how to constrain this)


% D_fix([1 3],end)=1;%fix the Y Z disp
% D_fix([7  9],end)=1;%fix the rotation
% D_fix([5 ],end)=1;%fix the rotation
% D_fix([11],end)=1;%fix the rotation

% Fixed: X Y Z axis
D_fix([1 3],end)=1;%fix the Y Z disp
D_fix([7 8 9],end)=1;%fix the rotation
D_fix([4 5 6],end)=1;%fix the rotation
D_fix([10 11 12],end)=1;%fix the rotation

B_index=find(D_fix==0);
%
Proj_M=Proj_M(:,B_index);

% displacement x
N_disp_x=50;
disp_x=linspace(0,50,N_disp_x+1);
% Angle_y=linspace(0,90,19)*pi/180;
%define the start points and frames
U=0*disp_x;
for mm=1:N_disp_x+1
    mm
    r1=[50-disp_x(mm) 0 0]';
    Q_sec(1:3,1)=r1(1:3);
    Q_sec(4:12,end)=r2_frame(:);
    %define the end points and frames
    kk=0;
    flag1=1
    while (kk==0||flag1>0.00001||flag2>0.10)&&kk<2000
        kk=kk+1;
        [dFq_dq_sec,Fq_sec(:,mm),q_node,U(mm)] = Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A);
        dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
        dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,mm);
        dQ_sec=Proj_M*dQ_sec_2;
        dQ_sec=reshape(dQ_sec,12,[]);
        Q_sec=Q_sec+dQ_sec;
        flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e+1);
        flag2=mean(abs(Proj_M'*Fq_sec(:,mm)));
    end

    Angle_T(mm)=atan2d(Q_sec(6,end),-Q_sec(4,end));
    Q_sec_M(:,:,mm)=Q_sec;
    figure(1)
    C_coord=[.0 .45 .74; .74 .45 .0];
    % cla
    if mod(mm,5)==1
        for ii=1:N_e

            r_node_1=q_node(1:3,:,ii);
            dr_dx_1=q_node(4:6,:,ii);
            dr_dy_1=q_node(7:9,:,ii);
            dr_dz_1=q_node(10:12,:,ii);
            L_co=2;
            figure(1)
            hold on
            plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

            %     r_node_2=[-1 1 -1]'.*r_node_1;
            %     dr_dx_2=[-1 1 -1]'.*dr_dx_1;
            %     dr_dy_2=[-1 1 -1]'.*dr_dy_1;
            %     dr_dz_2=[-1 1 -1]'.*dr_dz_1;
            %     plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            %     plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
            %         [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            %     plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
            %         [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

            r_node_3=[1 -1 -1]'.*r_node_1;
            dr_dx_3=[1 -1 -1]'.*dr_dx_1;
            dr_dy_3=[1 -1 -1]'.*dr_dy_1;
            dr_dz_3=[1 -1 -1]'.*dr_dz_1;
            plot3(r_node_3(1,:),r_node_3(2,:),r_node_3(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dy_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dy_3(2,:)], ...
                [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dy_3(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dz_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dz_3(2,:)], ...
                [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dz_3(3,:)],'-','color',C_coord(2,:),'linewidth',1)

            %     r_node_4=[-1 -1 1]'.*r_node_1;
            %     dr_dx_4=[-1 -1 1]'.*dr_dx_1;
            %     dr_dy_4=[-1 -1 1]'.*dr_dy_1;
            %     dr_dz_4=[-1 -1 1]'.*dr_dz_1;
            %     plot3(r_node_4(1,:),r_node_4(2,:),r_node_4(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
            %     plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dy_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dy_4(2,:)], ...
            %         [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dy_4(3,:)],'-','color',C_coord(1,:),'linewidth',1)
            %     plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dz_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dz_4(2,:)], ...
            %         [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dz_4(3,:)],'-','color',C_coord(2,:),'linewidth',1)

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

figure (2)
hold on
% plot3(disp_x(1:7)*0+Angle_T(1:7),disp_x(1:7),U(1:7),'-o')
plot3(disp_x*0+Angle_T,disp_x,U,'-o')
grid on
box on
xlabel('TAngle')
ylabel('disp_x')
zlabel('U')

figure (3)
hold on
% plot3(disp_x(1:7)*0+Angle_T(1:7),disp_x(1:7),-Fq_sec(1,1:7),'-o')
plot3(disp_x*0+Angle_T,disp_x,-Fq_sec(1,:),'-o')
grid on
box on
xlabel('TAngle')
ylabel('disp_x')
zlabel('F_x')



%%
load("Rod_Twist_2.mat")
%BC3 Antisymmetric twisting
% r1 start point
% position:
% control: X
% Fixed: Y Z
% Rotation
% Free: X axis in Y-Z
D_fix=0*Q_sec;%fixed dimensions/varibales
D_fix(1,1)=1;%control the X disp
D_fix(2:3,1)=1;%fix the Y Z disp

% free x local axis
D_fix(4,1)=1;%fixed dimensions/varibales
D_fix(8:9,1)=1;%fix the rotation
D_fix(10,1)=1;%fix the rotation
%
% r2
% position:
% Free: Y
% Fixed: X Z
% Rotation
% % Free: Y axis  
% D_fix([1 3],end)=1;%fix the Y Z disp
% D_fix([7  9],end)=1;%fix the rotation
% D_fix([5 ],end)=1;%fix the rotation
% D_fix([11],end)=1;%fix the rotation

% Fixed: X Y Z axis 
D_fix([1 3],end)=1;%fix the Y Z disp
D_fix([7 8 9],end)=1;%fix the rotation
D_fix([4 5 6],end)=1;%fix the rotation
D_fix([10 11 12],end)=1;%fix the rotation

Proj_M=eye(12*(N_e+1),12*(N_e+1));%dedemesionalize the variables
B_index=find(D_fix==0);%
Proj_M=Proj_M(:,B_index);

for nn=28:size(disp_x,2)
    nn
    if (18+nn*2)>90
        Beta_y=[0:1:90];
    else
        Beta_y=linspace(0,18+nn*2,91);
    end

    U=[];
    for mm=1:size(Beta_y,2)
        mm
        %define the end frame
        alpha2=0*pi/180;
        beta2=Beta_y(mm)*pi/180;
        theta2=0*pi/180;

        r2_frame=[-1 0 0; 0 -1 0;  0 0 1]';
        RM_A=[1 0 0;
            0 cos(alpha2) -sin(alpha2);
            0 sin(alpha2)  cos(alpha2)];

        RM_B=[cos(beta2)  0  sin(beta2);
            0           1  0;
            -sin(beta2) 0  cos(beta2)];

        RM_C=[cos(theta2) -sin(theta2)  0;
            sin(theta2)  cos(theta2)  0
            0            0            1];
        r2_frame=RM_C*RM_A*RM_B*r2_frame;

        Q_sec=Q_sec_M(:,:,nn);
        Q_sec(4:12,end)=r2_frame(:);
        kk=0;
        flag1=1;
        while (kk==0||flag1>0.00001||flag2>0.0010)&&kk<2000
            kk=kk+1;
            [dFq_dq_sec,Fq_sec(:,mm),q_node,U(mm)] = Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A);
            dFq_dq_sec_2=Proj_M'*dFq_dq_sec*Proj_M;
            dQ_sec_2=-pinv(dFq_dq_sec_2)*Proj_M'*Fq_sec(:,mm);
            dQ_sec=Proj_M*dQ_sec_2;
            dQ_sec=reshape(dQ_sec,12,[]);
            Q_sec=Q_sec+dQ_sec;
            flag1=sum(dQ_sec(1:3,:).^2,[1 2 3])/(N_e+1);
            flag2=mean(abs(Proj_M'*Fq_sec(:,mm)));
        end

        if mm==size(Beta_y,2)&&mod(nn,10)==1
            figure(1)
            C_coord=[.0 .45 .74; .74 .45 .0];
            % cla
            for ii=1:N_e

                r_node_1=q_node(1:3,:,ii);
                dr_dx_1=q_node(4:6,:,ii);
                dr_dy_1=q_node(7:9,:,ii);
                dr_dz_1=q_node(10:12,:,ii);
                L_co=2;
                figure(1)
                hold on
                plot3(r_node_1(1,:),r_node_1(2,:),r_node_1(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dy_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dy_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dy_1(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_1(1,:);r_node_1(1,:)+L_co*dr_dz_1(1,:)],[r_node_1(2,:);r_node_1(2,:)+L_co*dr_dz_1(2,:)], ...
                    [r_node_1(3,:);r_node_1(3,:)+L_co*dr_dz_1(3,:)],'-','color',C_coord(2,:),'linewidth',1)

                %             r_node_2=[-1 1 -1]'.*r_node_1;
                %             dr_dx_2=[-1 1 -1]'.*dr_dx_1;
                %             dr_dy_2=[-1 1 -1]'.*dr_dy_1;
                %             dr_dz_2=[-1 1 -1]'.*dr_dz_1;
                %             plot3(r_node_2(1,:),r_node_2(2,:),r_node_2(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                %             plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dy_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dy_2(2,:)], ...
                %                 [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dy_2(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                %             plot3([r_node_2(1,:);r_node_2(1,:)+L_co*dr_dz_2(1,:)],[r_node_2(2,:);r_node_2(2,:)+L_co*dr_dz_2(2,:)], ...
                %                 [r_node_2(3,:);r_node_2(3,:)+L_co*dr_dz_2(3,:)],'-','color',C_coord(2,:),'linewidth',1)

                r_node_3=[1 -1 -1]'.*r_node_1;
                dr_dx_3=[1 -1 -1]'.*dr_dx_1;
                dr_dy_3=[1 -1 -1]'.*dr_dy_1;
                dr_dz_3=[1 -1 -1]'.*dr_dz_1;
                plot3(r_node_3(1,:),r_node_3(2,:),r_node_3(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dy_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dy_3(2,:)], ...
                    [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dy_3(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                plot3([r_node_3(1,:);r_node_3(1,:)+L_co*dr_dz_3(1,:)],[r_node_3(2,:);r_node_3(2,:)+L_co*dr_dz_3(2,:)], ...
                    [r_node_3(3,:);r_node_3(3,:)+L_co*dr_dz_3(3,:)],'-','color',C_coord(2,:),'linewidth',1)

                %             r_node_4=[-1 -1 1]'.*r_node_1;
                %             dr_dx_4=[-1 -1 1]'.*dr_dx_1;
                %             dr_dy_4=[-1 -1 1]'.*dr_dy_1;
                %             dr_dz_4=[-1 -1 1]'.*dr_dz_1;
                %             plot3(r_node_4(1,:),r_node_4(2,:),r_node_4(3,:),'-','color',[0.5 0.5 0.5],'linewidth',1)
                %             plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dy_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dy_4(2,:)], ...
                %                 [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dy_4(3,:)],'-','color',C_coord(1,:),'linewidth',1)
                %             plot3([r_node_4(1,:);r_node_4(1,:)+L_co*dr_dz_4(1,:)],[r_node_4(2,:);r_node_4(2,:)+L_co*dr_dz_4(2,:)], ...
                %                 [r_node_4(3,:);r_node_4(3,:)+L_co*dr_dz_4(3,:)],'-','color',C_coord(2,:),'linewidth',1)
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
        Q_sec_MM{nn}(:,:,mm)=Q_sec;
    end
    % write down the data
    U_M(:,nn)=U;
    TAngle_M(:,nn)=Beta_y;
    disp_x_M(:,nn)=disp_x(nn)+0*Beta_y;

    % caculate the minimal energy state
    [U_min(nn), b]=min(U_M(:,nn));
    TAngle_min(nn)=TAngle_M(b,nn);
    disp_x_min(nn)=disp_x(nn);


    save('Rod_twist_2','disp_x','disp_x_M','TAngle_M',"U_M",'disp_x_min',"TAngle_min",'U_min','Q_sec_MM')

    figure (2)
    hold on
    scatter3(Beta_y,disp_x_M(:,nn),U_M(:,nn), [],U_M(:,nn),'o');
    colormap;
    grid on
    box on
    xlabel('TAngle')
    ylabel('disp_x')
    zlabel('U')    

end






%%

figure(4)
% for ii=1:51
%     U_M_2(:,nn)=U_M{ii};
%     TAngle_M_2(:,nn)=TAngle_M{ii};
%     disp_x_M_2(:,nn)=disp_x_M{ii};
% end
hold on
surf(TAngle_M,disp_x_M,U_M,'FaceAlpha',0.7)
shading interp
% contour3(TAngle_M,disp_x_M,U_M,15,'LineWidth',1,'LineStyle','--','LineColor',[0.5 0.5 0.5])
plot3(TAngle_M(:,1:5:51),disp_x_M(:,1:5:51),U_M(:,1:5:51),'LineWidth',1,'LineStyle','--','Color',[0.5 0.5 0.5])
colormap turbo
grid on
box on
figure (4)
hold on
plot3(TAngle_min,disp_x_min,U_min,'o-','LineWidth',2,'Color',[0 0 0],'MarkerFaceColor','w')
grid on
xlabel('Twisting Angle(°)')
ylabel('Displacement(cm)')
zlabel('Elastic Energy()')

%%









