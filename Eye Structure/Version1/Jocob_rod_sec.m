function [dFq_dq_sec,Fq_sec,q_node,U_sec,E_elastic_sec] = Jocob_rod_sec(Q_sec,N_e,L_e,N_node,Par_E,A)
%caculated the Jcob matrix of a rod section

Fq_sec=zeros(12*(N_e+1),1);
dFq_dq_sec=zeros(12*(N_e+1),12*(N_e+1));
for jj=1:N_e
    qe1=Q_sec(:,jj);
    qe2=Q_sec(:,jj+1);
    qe=[qe1;qe2];
    [U(jj), Fqe,q_node(:,:,jj),E_elastic(:,jj)]=Element_energy(qe,L_e,N_node,Par_E,A);
    for ii=1:24
        qe_2=qe;
        
        if ismember(ii,1:3)
            var=L_e/2000;
        else
            var=0.01;
        end
        qe_2(ii)=qe(ii)+var;
        [U_2, Fq_2,q_node_2]=Element_energy(qe_2,L_e,N_node,Par_E,A);
        dFq_dqe(:,ii)=(Fq_2-Fqe)/var;
    end
    Fq_sec((jj-1)*12+1:jj*12+12)=Fq_sec((jj-1)*12+1:jj*12+12)+Fqe;
    dFq_dq_sec((jj-1)*12+1:jj*12+12,(jj-1)*12+1:jj*12+12)=...
        dFq_dq_sec((jj-1)*12+1:jj*12+12,(jj-1)*12+1:jj*12+12)+dFq_dqe;
end
U_sec=sum(U,"all");
E_elastic_sec=sum(E_elastic,2);
end