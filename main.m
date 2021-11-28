%%
%常数
global iephem ephname km au
au = 149597870.699626200;% de421 value for astronomical unit (kilometers)
emu = 398600.436233;% earth gravitational constant (km**3/sec**2)
ephname = 'de421.bsp';% initialize DE421 algorithm
iephem = 1;% load binary ephemeris data file
km = 1;% state is kilometers and kilometers/second
LightVelocity=299792.458;%光速，单位km/s
Re=6371*1.01;%地球半径，单位km；有一个0.01的修正
Rm=1737.4;%月球半径，单位km
Rs=695700;%太阳半径，单位km
%%
t_start=datetime('2022-01-01 00:00:00');
t_end=datetime('2023-01-01 00:00:00');
JD_start=juliandate(t_start);%计算起始时间的儒略日，即公历2022-01-01 00:00:00
JD_end=juliandate(t_end);%计算结束时间的儒略日，即公历2022-12-31 0:0:0
iteration=1/(86400);%计算间隔的儒略日
% iteration=1;
U1=[];%月食初亏的所有时刻
JD=JD_start;
LastFlag=0;%上一时刻是否在月食的标志，LastFlag=1表示上一时刻在月食
%循环计算
while JD<JD_end
    DT_last=0;%光行时迭代初值
    %日(真）地矢量
    Vec_S2E_real=jpleph_mice (JD-DT_last,3,11);
    R_S2E_real=Vec_S2E_real(1:3);%前三维为位置矢量
    l_S2E=norm(R_S2E_real);%日地光子传播距离
    DT=l_S2E/LightVelocity/86400;
    while DT-DT_last>=1e-3/86400%退出迭代的标志
        DT_last=DT;
        Vec_S2E_real=jpleph_mice (JD-DT_last,3,11);
        R_S2E_real=Vec_S2E_real(1:3);
        l_S2E=norm(R_S2E_real);
        DT=l_S2E/LightVelocity/86400;
    end
    Vec_S2E_apparent=jpleph_mice (JD-DT,3,11);%DT为精确的光行时
    R_S2E_apparent=Vec_S2E_apparent(1:3);
    %锥点到地球的位置矢量
    R_O2E=[R_S2E_apparent(1)*(Re/(Re-Rs));...
        R_S2E_apparent(2)*(Re/(Re-Rs));...
        R_S2E_apparent(3)*(Re/(Re-Rs))];
    
    %月地矢量
    Vec_M2E=jpleph_mice (JD,3,10);
    R_M2E=Vec_M2E(1:3);%前三维为位置矢量
    %锥点到月球的位置矢量
    R_O2M=[R_O2E(1)-R_M2E(1);...
        R_O2E(2)-R_M2E(2);...
        R_O2E(3)-R_M2E(3)];
    
    %求角度判据：地-锥-月的夹角
    theta_EOM=acos(dot(R_O2E,R_O2M)/(norm(R_O2M)*norm(R_O2E)));
    theta_E=atan(Re/norm(R_O2E));
    theta_M=atan(Rm/norm(R_O2M));
    %求角度差
    DeltaAngle=(theta_E+theta_M)-theta_EOM;
    %求距离判据1：锥-地距离和锥-月距离差
    LengthFlag1=norm(R_O2E)-norm(R_O2M);
    %求日月距离
    Vec_S2M_real=jpleph_mice (JD,10,11);
    R_S2M_real=Vec_S2M_real(1:3);%前三维为位置矢量
    %求距离判据2：日-地距离和日月距离差
    LengthFlag2=norm(R_S2M_real)-norm(R_S2E_real);
    if DeltaAngle>=0 && LengthFlag1>=0 && LengthFlag2>=0 && LastFlag==0
        LastFlag=1;
        U1=[U1,JD-(32.184+37)/86400];%TDB转化为UTC
    elseif DeltaAngle>=0 && LengthFlag1>=0 && LengthFlag2>=0
        LastFlag=1;
    else
        LastFlag=0;
    end
    JD=JD+iteration;
end
%%
%数据处理
for i=length(U1)
[iy, im, id, fd]=iauJd2cal(U1(i), 0);%儒略日转年月日,调用SOFA库的函数
%算出时分秒
Hour=floor(24*fd);
Minute=floor((24*fd-Hour)*60);
Second=(60*(24*fd-Hour)-Minute)*60;
fprintf('%i年的第%i次月食初亏的时间为：%i月%i日%i时%i分%f秒',iy,i,im,id,Hour,Minute,Second)
end

    
    
    
