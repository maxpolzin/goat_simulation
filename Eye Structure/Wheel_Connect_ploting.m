function Wheel_Connect_ploting(P_wheel,RM_wheel,R_wheel,N_spoke,Duplication_M,Options)
%%
%Ploting the spokes of the wheel
% eg.：
% P_wheel=[50,50,0]'
% RM_wheel=eye(3)
% R_wheel=120;
% N_spoke=12;
% Duplication_M=[1 1 1;-1 1 1;-1 -1 1;1 -1 1;];
TM_wheel=[RM_wheel,P_wheel(:,3); 0 0 0 1];
TM_wheel=TM_wheel*[eye(3),[0 20 0]'; 0 0 0 1];
Spoke_angle=linspace(0,2*pi,N_spoke+1);
Spoke_XYZ(1,:)=R_wheel*[0*Spoke_angle,cos(Spoke_angle)];
Spoke_XYZ(2,:)=R_wheel*[0*Spoke_angle,0*Spoke_angle];
Spoke_XYZ(3,:)=R_wheel*[0*Spoke_angle,sin(Spoke_angle)];
Spoke_XYZ(4,:)=1+R_wheel*[0*Spoke_angle,0*Spoke_angle];
Spoke_XYZ=TM_wheel*Spoke_XYZ;
Spoke_X=reshape(Spoke_XYZ(1,:),[],2);
Spoke_Y=reshape(Spoke_XYZ(2,:),[],2);
Spoke_Z=reshape(Spoke_XYZ(3,:),[],2);
Circle_X=Spoke_X(:,2)
Circle_Y=Spoke_Y(:,2)
Circle_Z=Spoke_Z(:,2)
Arrow_wheel=[P_wheel(:,3),P_wheel(:,3)+R_wheel*RM_wheel(:,2)]

for ii=1:size(Duplication_M,1)
    if nargin<6||Options==0||Options==1
        plot3(Duplication_M(ii,1)*reshape(P_wheel(1,[1,2,4,5,7,8]),2,[]),...
            Duplication_M(ii,2)*reshape(P_wheel(2,[1,2,4,5,7,8]),2,[]),...
            Duplication_M(ii,3)*reshape(P_wheel(3,[1,2,4,5,7,8]),2,[]), ...
            'k-',LineWidth=1)
    end
    if nargin<6||Options==0||Options==2
    plot3(Duplication_M(ii,1)*Spoke_X', ...
        Duplication_M(ii,2)*Spoke_Y', ...
        Duplication_M(ii,3)*Spoke_Z','-','color',[0.5 0.5 0.5],LineWidth=1)
    plot3(Duplication_M(ii,1)*Circle_X', ...
        Duplication_M(ii,2)*Circle_Y', ...
        Duplication_M(ii,3)*Circle_Z','--','color',[0.5 0.5 0.5],LineWidth=1)
    plot3(Duplication_M(ii,1)*Arrow_wheel(1,:), ...
        Duplication_M(ii,2)*Arrow_wheel(2,:), ...
        Duplication_M(ii,3)*Arrow_wheel(3,:),'r-',LineWidth=1,color=[1 0.4 0.16])
    end
end

% box on
% grid on
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% axis equal
% axis([-600 600 -700 700 -100 100])
% view([1 1 1])
%%
end