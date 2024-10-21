function [U_cal,Delta_U_cal,q_node,U_v,Delta_U_E,Dens_E_node] =...
    Element_energy(qe,L_e,N_node,Par_E,A)

% interplot the curve with boundaries
qe=reshape(qe,3,[]);

r1=qe(:,1);
r1_x=qe(:,2);
r1_y=qe(:,3);
r1_z=qe(:,4);

r2=qe(:,5);
r2_x=qe(:,6);
r2_y=qe(:,7);
r2_z=qe(:,8);



EA=Par_E(1);
EI=Par_E(2);
GI=Par_E(3);

EA=Par_E(1);
EI=Par_E(2);
GI=Par_E(3);

X=linspace(0,L_e,N_node+1);
Y=[0];
Z=[0];

Sc(1,:)=1-3*X.^2/(L_e^2)+2*X.^3/(L_e^3);
Sc(2,:)=X-2*X.^2/(L_e^1)+X.^3/(L_e^2);
Sc(3,:)=Y-Y*X/L_e;
Sc(4,:)=Z-Z*X/L_e;
Sc(5,:)=3*X.^2/(L_e^2)-2*X.^3/(L_e^3);
Sc(6,:)=-X.^2/(L_e^1)+X.^3/(L_e^2);
Sc(7,:)=Y*X/L_e;
Sc(8,:)=Z*X/L_e;



dSc_dx(1,:)=-6*X/(L_e^2)+6*X.^2/(L_e^3);
dSc_dx(2,:)=1-4*X/L_e+3*X.^2/(L_e^2);
dSc_dx(3,:)=-Y/L_e;
dSc_dx(4,:)=-Z/L_e;
dSc_dx(5,:)=6*X/(L_e^2)-6*X.^2/(L_e^3);
dSc_dx(6,:)=-2*X/L_e+3*X.^2/(L_e^2);
dSc_dx(7,:)=Y/L_e;
dSc_dx(8,:)=Y/L_e;

dSc_dy=zeros(8,N_node+1);
dSc_dy(3,:)=1-X/L_e;
dSc_dy(7,:)=X/L_e;

dSc_dz=zeros(8,N_node+1);
dSc_dz(4,:)=1-X/L_e;
dSc_dz(8,:)=X/L_e;

dSc_2dx=zeros(8,N_node+1);
dSc_2dx(1,:)=-6/(L_e^2)+12*X/(L_e^3);
dSc_2dx(2,:)=-4/L_e+6*X/(L_e^2);
dSc_2dx(5,:)=6/(L_e^2)-12*X/(L_e^3);
dSc_2dx(6,:)=-2/L_e+6*X/(L_e^2);

dSc_dxdy=zeros(8,N_node+1);
dSc_dxdy(3,:)=-1/L_e;
dSc_dxdy(7,:)=1/L_e;

dSc_dxdz=zeros(8,N_node+1);
dSc_dxdz(4,:)=-1/L_e;
dSc_dxdz(8,:)=1/L_e;

r_x_c=qe*Sc;
dr_dx_c=qe*dSc_dx;
dr_dy_c=qe*dSc_dy;
dr_dz_c=qe*dSc_dz;
dr_2dx_c=qe*dSc_2dx;
dr_dxdy_c=qe*dSc_dxdy;
dr_dxdz_c=qe*dSc_dxdz;

q_node=[r_x_c;dr_dx_c;dr_dy_c;dr_dz_c];

V_b=cross(dr_dx_c,dr_2dx_c)./sum(dr_dx_c.^2);% bending vector
k_b=sqrt(sum(V_b.^2,1));
k_t=(dot(dr_dz_c,dr_dxdy_c)-dot(dr_dy_c,dr_dxdz_c))/2;%torsion rate
Strian_n=0.5*(sum(dr_dx_c.^2)-1);
% Strian_n=sqrt(sum(dr_dx_c.^2))-1;


% constrains to r_y  r_z 
F_ec=[r1_y'*r1_y-1, r1_z'*r1_z-1, r1_y'*r1_z, (sqrt(sum(r1_x.^2)).*cross(r1_y,r1_z)-r1_x)',...
      r2_y'*r2_y-1,r2_z'*r2_z-1,r2_y'*r2_z,(sqrt(sum(r2_x.^2)).*cross(r2_y,r2_z)-r2_x)']';



dUe_b_dx=0.5*EI*sum(V_b.^2,1);
dUe_t_dx=0.5*GI*k_t.^2;
dUe_n_dx=0.5*EA*Strian_n.^2;
dUe_dx=0.5*(EI*sum(V_b.^2,1)+GI*k_t.^2+EA*Strian_n.^2+A*F_ec'*F_ec);% A is the the penalty factor

Dens_E_node=dUe_b_dx+dUe_t_dx+dUe_n_dx;

Ue_b=L_e/N_node*sum(dUe_b_dx(1:end-1)+dUe_b_dx(2:end))/2;
Ue_t=L_e/N_node*sum(dUe_t_dx(1:end-1)+dUe_t_dx(2:end))/2;
Ue_n=L_e/N_node*sum(dUe_n_dx(1:end-1)+dUe_n_dx(2:end))/2;
Ue=Ue_b+Ue_t+Ue_n;
U_penalty=0.5*A*F_ec'*F_ec*L_e;
U_cal=Ue+U_penalty;


Delta_qe=eye(3*8,3*8);
Delta_qe=reshape(Delta_qe,[3 8 24]);
Delta_dr_dx_c=pagemtimes(Delta_qe,dSc_dx);%Dqe
Delta_dr_dy_c=pagemtimes(Delta_qe,dSc_dy);%Dqe
Delta_dr_dz_c=pagemtimes(Delta_qe,dSc_dz);%Dqe

Delta_dr_2dx_c=pagemtimes(Delta_qe,dSc_2dx);
Delta_dr_dxdy_c=pagemtimes(Delta_qe,dSc_dxdy);
Delta_dr_dxdz_c=pagemtimes(Delta_qe,dSc_dxdz);

Delta_V_b=cross(Delta_dr_dx_c,repmat(dr_2dx_c,1,1,24))./sum(dr_dx_c.^2)...
-cross(dr_dx_c,dr_2dx_c).*(2*dot(Delta_dr_dx_c,repmat(dr_dx_c,1,1,24),1))./(sum(dr_dx_c.^2).^2)...
+ cross(repmat(dr_dx_c,1,1,24),Delta_dr_2dx_c)./sum(dr_dx_c.^2);          

Delta_dUe_b_dx=EI*dot(repmat(V_b,1,1,24),Delta_V_b);
Delta_dUe_b_dx=reshape(Delta_dUe_b_dx,N_node+1,[],1)';

Delta_k_t=0.5*dot(Delta_dr_dz_c,repmat(dr_dxdy_c,1,1,24)) ...    
+0.5*dot(repmat(dr_dz_c,1,1,24),Delta_dr_dxdy_c)...
-0.5*dot(Delta_dr_dy_c,repmat(dr_dxdz_c,1,1,24))...
-0.5*dot(repmat(dr_dy_c,1,1,24),Delta_dr_dxdz_c);%torsion rate
Delta_k_t=reshape(Delta_k_t,N_node+1,[],1)';
Delta_dUe_t_dx=GI*k_t.*Delta_k_t;


Delta_Strian_n=dot(repmat(dr_dx_c,1,1,24),Delta_dr_dx_c);
% Delta_Strian_n=dot(repmat(dr_dx_c,1,1,24),Delta_dr_dx_c)./sqrt(sum(dr_dx_c.^2));
Delta_Strian_n=reshape(Delta_Strian_n,N_node+1,[],1)';
Delta_dUe_n_dx=EA*repmat(Strian_n,24,1).*Delta_Strian_n;

Delta_Ue_b=L_e/N_node*sum(Delta_dUe_b_dx(:,1:end-1)+Delta_dUe_b_dx(:,2:end),2)/2;
Delta_Ue_t=L_e/N_node*sum(Delta_dUe_t_dx(:,1:end-1)+Delta_dUe_t_dx(:,2:end),2)/2;
Delta_Ue_n=L_e/N_node*sum(Delta_dUe_n_dx(:,1:end-1)+Delta_dUe_n_dx(:,2:end),2)/2;

Delta_F_ec=zeros(12,24);
Delta_F_ec(1,7:9)=2*r1_y';
Delta_F_ec(2,10:12)=2*r1_z';
Delta_F_ec(3,7:12)=[r1_z',r1_y'];
Delta_F_ec(4:6,4)=r1_x(1)./sqrt(sum(r1_x.^2)).*cross(r1_y,r1_z)-[1 0 0]';
Delta_F_ec(4:6,5)=r1_x(2)./sqrt(sum(r1_x.^2)).*cross(r1_y,r1_z)-[0 1 0]';
Delta_F_ec(4:6,6)=r1_x(3)./sqrt(sum(r1_x.^2)).*cross(r1_y,r1_z)-[0 0 1]';
Delta_F_ec(4:6,7:9)=[0,r1_z(3),-r1_z(2);-r1_z(3),0,r1_z(1); r1_z(2),-r1_z(1),0]*sqrt(sum(r1_x.^2));
Delta_F_ec(4:6,10:12)=[0,-r1_y(3),r1_y(2);r1_y(3),0,-r1_y(1); -r1_y(2),r1_y(1),0]*sqrt(sum(r1_x.^2));
Delta_F_ec(7,19:21)=2*r2_y';
Delta_F_ec(8,22:24)=2*r2_z';
Delta_F_ec(9,19:24)=[r2_z',r2_y'];
Delta_F_ec(10:12,16)=r2_x(1)./sqrt(sum(r2_x.^2)).*cross(r2_y,r2_z)-[1 0 0]';
Delta_F_ec(10:12,17)=r2_x(2)./sqrt(sum(r2_x.^2)).*cross(r2_y,r2_z)-[0 1 0]';
Delta_F_ec(10:12,18)=r2_x(3)./sqrt(sum(r2_x.^2)).*cross(r2_y,r2_z)-[0 0 1]';
Delta_F_ec(10:12,19:21)=[0,r2_z(3),-r2_z(2);-r2_z(3),0,r2_z(1); r2_z(2),-r2_z(1),0]*sqrt(sum(r2_x.^2));
Delta_F_ec(10:12,22:24)=[0,-r2_y(3),r2_y(2);r2_y(3),0,-r2_y(1); -r2_y(2),r2_y(1),0]*sqrt(sum(r2_x.^2));

Delta_U_penalty=L_e*A*dot(repmat(F_ec,1,24),Delta_F_ec)';

Delta_Ue=Delta_Ue_b+Delta_Ue_t+Delta_Ue_n;

Delta_U_cal=Delta_Ue+Delta_U_penalty;

U_v=[Ue,Ue_b,Ue_t,Ue_n,U_penalty];
Delta_U_E=[Delta_Ue,Delta_Ue_b,Delta_Ue_t,Delta_Ue_n,Delta_U_penalty];

end