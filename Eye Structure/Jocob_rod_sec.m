function [dFq_dq_sec,Fq_sec,q_node,U_cal_sec,U_v_sec,U_E_ele,Dens_E_node] = Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A)
%caculated the Jcob matrix of a rod section

Fq_sec=zeros(12*(N_e+1),1);
dFq_dq_sec=zeros(12*(N_e+1),12*(N_e+1));
for jj=1:N_e
    qe1=Q_sec(:,jj);
    qe2=Q_sec(:,jj+1);
    qe=[qe1;qe2];
    [U_cal_ele(jj), F_cal_qe,q_node(:,:,jj),U_v_ele(:,jj),~,Dens_E_node(jj,:)]=Element_energy(qe,L_e,N_node,Par_E,A);
    for ii=1:24
        qe_2=qe;
        
        if ismember(ii,1:3)
            var=L_e/2000;
        else
            var=0.01;
        end
        qe_2(ii)=qe(ii)+var;
        [U_2, F_qe_2,q_node_2]=Element_energy(qe_2,L_e,N_node,Par_E,A);
        dFqe_dqe(:,ii)=(F_qe_2-F_cal_qe)/var;
    end
    Fq_sec((jj-1)*12+1:jj*12+12)=Fq_sec((jj-1)*12+1:jj*12+12)+F_cal_qe;
    dFq_dq_sec((jj-1)*12+1:jj*12+12,(jj-1)*12+1:jj*12+12)=...
        dFq_dq_sec((jj-1)*12+1:jj*12+12,(jj-1)*12+1:jj*12+12)+dFqe_dqe;
end
U_cal_sec=sum(U_cal_ele,"all");
U_v_sec=sum(U_v_ele,2);
U_E_ele=U_v_ele(1,:);
end