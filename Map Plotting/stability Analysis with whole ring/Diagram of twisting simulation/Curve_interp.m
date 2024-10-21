function [Q_sec] = Curve_interp(q1_sec,q2_sec,L_sec,N_e)
% interplot the curve with boundaries

r1=q1_sec(1:3);
r1_x=q1_sec(4:6);
r1_y=q1_sec(7:9);
r1_z=q1_sec(10:12);


r2=q2_sec(1:3);
r2_x=q2_sec(4:6);
r2_y=q2_sec(7:9);
r2_z=q2_sec(10:12);

q_sec=[q1_sec;q2_sec];
q_sec=reshape(q_sec,3,[]);

X=linspace(0,L_sec,N_e+1);
Y=[0];
Z=[0];
 
Sc(1,:)=1-3*X.^2/(L_sec^2)+2*X.^3/(L_sec^3);
Sc(2,:)=X-2*X.^2/(L_sec^1)+X.^3/(L_sec^2);
Sc(3,:)=Y-Y*X/L_sec;
Sc(4,:)=Z-Z*X/L_sec;
Sc(5,:)=3*X.^2/(L_sec^2)-2*X.^3/(L_sec^3);
Sc(6,:)=-X.^2/(L_sec^1)+X.^3/(L_sec^2);
Sc(7,:)=Y*X/L_sec;
Sc(8,:)=Z*X/L_sec;

dSc_dx(1,:)=-6*X/(L_sec^2)+6*X.^2/(L_sec^3);
dSc_dx(2,:)=1-4*X/L_sec+3*X.^2/(L_sec^2);
dSc_dx(3,:)=-Y/L_sec;
dSc_dx(4,:)=-Z/L_sec;
dSc_dx(5,:)=6*X/(L_sec^2)-6*X.^2/(L_sec^3);
dSc_dx(6,:)=-2*X/L_sec+3*X.^2/(L_sec^2);
dSc_dx(7,:)=Y/L_sec;
dSc_dx(8,:)=Y/L_sec;

dSc_dy=zeros(8,N_e+1);
dSc_dy(3,:)=1-X/L_sec;
dSc_dy(7,:)=X/L_sec;

dSc_dz=zeros(8,N_e+1);
dSc_dz(4,:)=1-X/L_sec;
dSc_dz(8,:)=X/L_sec;

%Y=0 Z=0
r_x_c=q_sec*Sc;
dr_dx_c=q_sec*dSc_dx;
dr_dy_c=q_sec*dSc_dy;
dr_dz_c=q_sec*dSc_dz;

Q_sec=[r_x_c;dr_dx_c;dr_dy_c;dr_dz_c];



end