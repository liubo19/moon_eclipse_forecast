%%
%����
global iephem ephname km au
au = 149597870.699626200;% de421 value for astronomical unit (kilometers)
emu = 398600.436233;% earth gravitational constant (km**3/sec**2)
ephname = 'de421.bsp';% initialize DE421 algorithm
iephem = 1;% load binary ephemeris data file
km = 1;% state is kilometers and kilometers/second
LightVelocity=299792.458;%���٣���λkm/s
Re=6371*1.01;%����뾶����λkm����һ��0.01������
Rm=1737.4;%����뾶����λkm
Rs=695700;%̫���뾶����λkm
%%
t_start=datetime('2022-01-01 00:00:00');
t_end=datetime('2023-01-01 00:00:00');
JD_start=juliandate(t_start);%������ʼʱ��������գ�������2022-01-01 00:00:00
JD_end=juliandate(t_end);%�������ʱ��������գ�������2022-12-31 0:0:0
iteration=1/(86400);%��������������
% iteration=1;
U1=[];%��ʳ����������ʱ��
JD=JD_start;
LastFlag=0;%��һʱ���Ƿ�����ʳ�ı�־��LastFlag=1��ʾ��һʱ������ʳ
%ѭ������
while JD<JD_end
    DT_last=0;%����ʱ������ֵ
    %��(�棩��ʸ��
    Vec_S2E_real=jpleph_mice (JD-DT_last,3,11);
    R_S2E_real=Vec_S2E_real(1:3);%ǰ��άΪλ��ʸ��
    l_S2E=norm(R_S2E_real);%�յع��Ӵ�������
    DT=l_S2E/LightVelocity/86400;
    while DT-DT_last>=1e-3/86400%�˳������ı�־
        DT_last=DT;
        Vec_S2E_real=jpleph_mice (JD-DT_last,3,11);
        R_S2E_real=Vec_S2E_real(1:3);
        l_S2E=norm(R_S2E_real);
        DT=l_S2E/LightVelocity/86400;
    end
    Vec_S2E_apparent=jpleph_mice (JD-DT,3,11);%DTΪ��ȷ�Ĺ���ʱ
    R_S2E_apparent=Vec_S2E_apparent(1:3);
    %׶�㵽�����λ��ʸ��
    R_O2E=[R_S2E_apparent(1)*(Re/(Re-Rs));...
        R_S2E_apparent(2)*(Re/(Re-Rs));...
        R_S2E_apparent(3)*(Re/(Re-Rs))];
    
    %�µ�ʸ��
    Vec_M2E=jpleph_mice (JD,3,10);
    R_M2E=Vec_M2E(1:3);%ǰ��άΪλ��ʸ��
    %׶�㵽�����λ��ʸ��
    R_O2M=[R_O2E(1)-R_M2E(1);...
        R_O2E(2)-R_M2E(2);...
        R_O2E(3)-R_M2E(3)];
    
    %��Ƕ��оݣ���-׶-�µļн�
    theta_EOM=acos(dot(R_O2E,R_O2M)/(norm(R_O2M)*norm(R_O2E)));
    theta_E=atan(Re/norm(R_O2E));
    theta_M=atan(Rm/norm(R_O2M));
    %��ǶȲ�
    DeltaAngle=(theta_E+theta_M)-theta_EOM;
    %������о�1��׶-�ؾ����׶-�¾����
    LengthFlag1=norm(R_O2E)-norm(R_O2M);
    %�����¾���
    Vec_S2M_real=jpleph_mice (JD,10,11);
    R_S2M_real=Vec_S2M_real(1:3);%ǰ��άΪλ��ʸ��
    %������о�2����-�ؾ�������¾����
    LengthFlag2=norm(R_S2M_real)-norm(R_S2E_real);
    if DeltaAngle>=0 && LengthFlag1>=0 && LengthFlag2>=0 && LastFlag==0
        LastFlag=1;
        U1=[U1,JD-(32.184+37)/86400];%TDBת��ΪUTC
    elseif DeltaAngle>=0 && LengthFlag1>=0 && LengthFlag2>=0
        LastFlag=1;
    else
        LastFlag=0;
    end
    JD=JD+iteration;
end
%%
%���ݴ���
for i=length(U1)
[iy, im, id, fd]=iauJd2cal(U1(i), 0);%������ת������,����SOFA��ĺ���
%���ʱ����
Hour=floor(24*fd);
Minute=floor((24*fd-Hour)*60);
Second=(60*(24*fd-Hour)-Minute)*60;
fprintf('%i��ĵ�%i����ʳ������ʱ��Ϊ��%i��%i��%iʱ%i��%f��',iy,i,im,id,Hour,Minute,Second)
end

    
    
    
