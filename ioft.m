function [Tilt_Correction,Stepwise_Unfold] = ioft(filename)
%ioft Export the stepwise unfold results to a csv file, plot 1) the data 
%both before and after tilt correction in equal area projection; 2) 
%precision parameter vesus unfold steps; 3) unfolding of both limbs.
%   e.g. Data = ioft_t_new("test")
%
%   Input format: 
%       string format filename of a spreadsheet, no extenion
%       needed. e.g. type in "test" works for "test.xlsx" file
%
%   Input file: 
%       This code can only recognize spreadsheet file with ".xlsx"
%       extension.
%       Headers should be Group, Gg, Ig, Strick_R, and Dip, 
%           e.g.
%
%             | Group | Dg | Ig | Strike_R | Dip |
%             ------------------------------------
%             |   G1  | 10 | 10 |    38    |  72 |
%             ------------------------------------
%             |   G2  | 10 | 10 |    38    |  72 |
%   
%   Output file:
%       1) csv file containing tilt correction data
%       2) csv file containing stepwise unfolding data
%       2) figure 1: equal area plot both before and after tilt correction
%       3) figure 2: precision parameter vesus unfold steps
%       4) figure 3: unfolding of both limbs


T = readtable(strcat(filename,'.xlsx')); %   import xls data
T = t_cor_step(T,1,1);
T_cor = T;
T.Group = string(T.Group);
G = groupsummary(T,"Group");
Group = G.Group;
T1 = T(T.Group==Group(1),:);
T2 = T(T.Group==Group(2),:);
[Dg_1,Ig_1,Kg_1,a95g_1] = a95calc(T1.Dg,T1.Ig);
[Dg_2,Ig_2,Kg_2,a95g_2] = a95calc(T2.Dg,T2.Ig);
[Dg_m,Ig_m,Kg_m,a95g_m] = a95calc(T.Dg,T.Ig);
[Ds_1,Is_1,Ks_1,a95s_1] = a95calc(T1.Ds,T1.Is);
[Ds_2,Is_2,Ks_2,a95s_2] = a95calc(T2.Ds,T2.Is);
[Ds_m,Is_m,Ks_m,a95s_m] = a95calc(T.Ds,T.Is);
Dg = [Dg_1,Dg_2,Dg_m]'; Ig = [Ig_1,Ig_2,Ig_m]'; a95g = [a95g_1,a95g_2,a95g_m]'; 
Ds = [Ds_1,Ds_2,Ds_m]'; Is = [Is_1,Is_2,Is_m]'; a95s = [a95s_1,a95s_2,a95s_m]'; 
Kg = [Kg_1,Kg_2,Kg_m]'; Ks = [Ks_1,Ks_2,Ks_m]';


for ele = 1:41
    step = (ele-11)/20;
    T_cor = t_cor_step(T_cor,step,1);    
    
    inc = T_cor.Is;
    th = 90-inc;    %   co-inclination
    n = length(inc);    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Result 1: Number of samples
    dr = pi/180;
    [mean_rt,mean_rk] = ArithmeticMean(th);
    Mean = 90-mean_rt;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Result 2: Arithmetic mean
    mean_rk;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Result 3: Inverse variance
    [aralevOutput_ainc,aralevOutput_ak,aralevOutput_t63,aralevOutput_a95,...
        aralevOutput_ierr] = ARALEV(inc);
    c2(ele) = aralevOutput_ainc;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Result 4: Mean inclination
    c3(ele) = aralevOutput_ak;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Result 5: Precision parameter
    aralevOutput_t63;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Result 6: Angular standard deviation 
    c4(ele) = aralevOutput_a95;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Result 7: 95% confidence limits 
    aralevOutput_ierr;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Result 8: Error flag
    c1(ele) = strcat(string(step*100),"%");  
end
varNames = {'Stepwise_unfold','I','K','alpha_95'};
Stepwise_Unfold = table(c1',c2',c3',c4','VariableNames',varNames);
writetable(Stepwise_Unfold,strcat(filename,'_Stepwise_Unfold.csv'),'Delimiter',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1 = [119 137 183]/255; % color 1
c2 = [235 105 105]/255; % color 2

figure   % equal area plots
set(gcf,"position",[10,10,1280,720])

%%%%%%%%%%%%%%%%%% equal area plot before tilt corection %%%%%%%%%%%%%%%%%%
subplot(1,2,1) 
r = 90;
rr = sqrt(2) * r; 
t = linspace(0, 2 * pi);
plot(rr * sin(t), rr * cos(t), 'k');
axis off;
hold on;
for t=0:10:350  
    theta=t/180*pi; 
    polar([theta,theta],[123,127.28],'k'); 
end
A = [-153, 153, 0, 0];
B = [0, 0, -153, 153];
plot(A, B, '.w', 'markersize', 1)
I = [0, 10, 20, 30, 40, 50, 60, 70, 80];
R = 2 * r * sind(I/2);
ax = zeros(1,18);
plot([-R, R], ax, '+k', 'markersize', 8);
plot(ax, [-R, R], '+k', 'markersize', 8);
plot([0, 0], 0, '+k', 'markersize',12,'LineWidth',1.5);
text(0, 145, 'N', 'fontsize', 14);
text(138, 0, 'E', 'fontsize', 14);
text(0, -145, 'S', 'fontsize', 14);
text(-152, 0, 'W', 'fontsize', 14);
axis equal
axis tight

% data points
num = length(T.Ig);
eInc = abs(T.Ig);
ea = 0.5.*(90-eInc);
er = 2 * r .* sind(ea);
ex = er .* sind(T.Dg);
ey = er .* cosd(T.Dg);
for i = 1:num
    if T.Ig(i) >= 0
        if T.Group(i)==Group(1)
            plot(ex(i), ey(i), 'o', 'color', c1, 'MarkerSize', 5,...
                'MarkerFaceColor',c1)
        else
            plot(ex(i), ey(i), 'o', 'color', c2, 'MarkerSize', 5,...
                'MarkerFaceColor',c2)
        end
    else
        if T.Group(i)==Group(1)
            plot(ex(i), ey(i), 'o', 'color',c1, 'MarkerSize', 5)
        else
            plot(ex(i), ey(i), 'o', 'color',c2, 'MarkerSize', 5)
        end
    end
end
text(max(ex(T.Group==Group(1)))+5,min(ey(T.Group==Group(1)))-5,...
    Group(1),"FontSize",14,"Color",c1);
text(max(ex(T.Group==Group(2)))+5,min(ey(T.Group==Group(2)))-5,...
    Group(2),"FontSize",14,"Color",c2);

% mean results with a95 circle
ea = 0.5.*(90-abs(Ig));
er = 2 * r .* sind(ea);
ex = er .* sind(Dg);
ey = er .* cosd(Dg);
if Ig(1) >= 0
    plot(ex(1),ey(1),'p','color',c1,'MarkerSize',12,'MarkerFaceColor',c1)
else
    plot(ex(1),ey(1),'p','color',c1,'MarkerSize',12)
end
if Ig(2) >= 0
    plot(ex(2),ey(2),'p','color',c2,'MarkerSize',12,'MarkerFaceColor',c2)
else
    plot(ex(2),ey(2),'p','color',c2,'MarkerSize',12)
end
if Ig(3) >= 0
    plot(ex(3),ey(3),'pr','MarkerSize',12,'MarkerFaceColor','r')
else
    plot(ex(3),ey(3),'pr','MarkerSize',12)
end

% 95% confidence circle
kz = sind(Ig);
kx = cosd(Ig).*cosd(Dg);
ky = cosd(Ig).*sind(Dg);

vx0 = cosd(Ig+a95g).*cosd(Dg);
vy0 = cosd(Ig+a95g).*sind(Dg);
vz0 = sind(Ig+a95g);


theta = linspace(0,360); count = 0;
for i = 1:3
    k = [kx(i),ky(i),kz(i)];
    v = [vx0(i),vy0(i),vz0(i)];
    Vr = zeros(length(theta),3);
    for j = 1:length(theta)
        Vr(j,:)= v.*cosd(theta(j))+cross(k,v)*sind(theta(j))+dot(k,v).*k.*(1-cosd(theta(j)));
    end
    Vr_Dec = atan2d(Vr(:,2),Vr(:,1));
    Vr_Inc = atan2d(Vr(:,3),sqrt(Vr(:,1).^2+Vr(:,2).^2));
    VrT = table(Vr_Dec,Vr_Inc);
    for i = 1:length(theta)-1
        VrT = [VrT;VrT(1,:)];
        VrT(1,:) = [];
        if Vr_Inc(i)*Vr_Inc(i+1)<0
            break
        end
        if i == 99
            VrT = [VrT;VrT(1,:)];
            VrT(1,:) = [];  
        end
    end
    er = 2 * r .* sind(0.5*(90-abs(VrT.Vr_Inc)));
    ex = er .* sind(VrT.Vr_Dec);
    ey = er .* cosd(VrT.Vr_Dec);
    conf_data = [ex,ey];
    conf_data_p = conf_data(VrT.Vr_Inc>=0,:);
    conf_data_n = conf_data(VrT.Vr_Inc<0,:);
    count = count+1;
    if count < 2
        plot(conf_data_p(:,1), conf_data_p(:,2), '-','Color',c1);
        plot(conf_data_n(:,1), conf_data_n(:,2), '--','Color',c1);
    elseif count <3
        plot(conf_data_p(:,1), conf_data_p(:,2), '-','Color',c2);
        plot(conf_data_n(:,1), conf_data_n(:,2), '--','Color',c2);
    else
        plot(conf_data_p(:,1), conf_data_p(:,2), '-r');
        plot(conf_data_n(:,1), conf_data_n(:,2), '--r');
    end
end

%%%%%%%%%%%%%%%%%%% equal area plot after tilt corection %%%%%%%%%%%%%%%%%%
subplot(1,2,2) 
r = 90;
rr = sqrt(2) * r; 
t = linspace(0, 2 * pi);
plot(rr * sin(t), rr * cos(t), 'k');
axis off;
hold on;
for t=0:10:350  
    theta=t/180*pi; 
    polar([theta,theta],[123,127.28],'k'); 
end
A = [-153, 153, 0, 0];
B = [0, 0, -153, 153];
plot(A, B, '.w', 'markersize', 1)
I = [0, 10, 20, 30, 40, 50, 60, 70, 80];
R = 2 * r * sind(I/2);
ax = zeros(1,18);
plot([-R, R], ax, '+k', 'markersize', 8);
plot(ax, [-R, R], '+k', 'markersize', 8);
plot([0, 0], 0, '+k', 'markersize',12,'LineWidth',1.5);
text(0, 145, 'N', 'fontsize', 14);
text(138, 0, 'E', 'fontsize', 14);
text(0, -145, 'S', 'fontsize', 14);
text(-152, 0, 'W', 'fontsize', 14);
axis equal
axis tight

% data points
num = length(T.Is);
eInc = abs(T.Is);
ea = 0.5.*(90-eInc);
er = 2 * r .* sind(ea);
ex = er .* sind(T.Ds);
ey = er .* cosd(T.Ds);
for i = 1:num
    if T.Is(i) >= 0
        if T.Group(i)==Group(1)
            plot(ex(i), ey(i), 'o', 'color', c1, 'MarkerSize', 5,...
                'MarkerFaceColor',c1)
        else
            plot(ex(i), ey(i), 'o', 'color', c2, 'MarkerSize', 5,...
                'MarkerFaceColor',c2)
        end
    else
        if T.Group(i)==Group(1)
            plot(ex(i), ey(i), 'o', 'color',c1, 'MarkerSize', 5)
        else
            plot(ex(i), ey(i), 'o', 'color',c2, 'MarkerSize', 5)
        end
    end
end
text(max(ex(T.Group==Group(1)))+5,min(ey(T.Group==Group(1)))-5,...
    Group(1),"FontSize",14,"Color",c1);
text(max(ex(T.Group==Group(2)))+5,min(ey(T.Group==Group(2)))-5,...
    Group(2),"FontSize",14,"Color",c2);

% mean results with a95 circle
ea = 0.5.*(90-abs(Is));
er = 2 * r .* sind(ea);
ex = er .* sind(Ds);
ey = er .* cosd(Ds);
if Is(1) >= 0
    plot(ex(1),ey(1),'p','color',c1,'MarkerSize',12,'MarkerFaceColor',c1)
else
    plot(ex(1),ey(1),'p','color',c1,'MarkerSize',12)
end
if Is(2) >= 0
    plot(ex(2),ey(2),'p','color',c2,'MarkerSize',12,'MarkerFaceColor',c2)
else
    plot(ex(2),ey(2),'p','color',c2,'MarkerSize',12)
end
if Is(3) >= 0
    plot(ex(3),ey(3),'pr','MarkerSize',12,'MarkerFaceColor','r')
else
    plot(ex(3),ey(3),'pr','MarkerSize',12)
end

% 95% confidence circle
kz = sind(Is);
kx = cosd(Is).*cosd(Ds);
ky = cosd(Is).*sind(Ds);

vx0 = cosd(Is+a95s).*cosd(Ds);
vy0 = cosd(Is+a95s).*sind(Ds);
vz0 = sind(Is+a95s);


theta = linspace(0,360); count = 0;
for i = 1:3
    k = [kx(i),ky(i),kz(i)];
    v = [vx0(i),vy0(i),vz0(i)];
    Vr = zeros(length(theta),3);
    for j = 1:length(theta)
        Vr(j,:)= v.*cosd(theta(j))+cross(k,v)*sind(theta(j))+dot(k,v).*k.*(1-cosd(theta(j)));
    end
    Vr_Dec = atan2d(Vr(:,2),Vr(:,1));
    Vr_Inc = atan2d(Vr(:,3),sqrt(Vr(:,1).^2+Vr(:,2).^2));
    VrT = table(Vr_Dec,Vr_Inc);
    for i = 1:length(theta)-1
        VrT = [VrT;VrT(1,:)];
        VrT(1,:) = [];
        if Vr_Inc(i)*Vr_Inc(i+1)<0
            break
        end
        if i == 99
            VrT = [VrT;VrT(1,:)];
            VrT(1,:) = [];  
        end
    end
    er = 2 * r .* sind(0.5*(90-abs(VrT.Vr_Inc)));
    ex = er .* sind(VrT.Vr_Dec);
    ey = er .* cosd(VrT.Vr_Dec);
    conf_data = [ex,ey];
    conf_data_p = conf_data(VrT.Vr_Inc>=0,:);
    conf_data_n = conf_data(VrT.Vr_Inc<0,:);
    count = count+1;
    if count < 2
        plot(conf_data_p(:,1), conf_data_p(:,2), '-','Color',c1);
        plot(conf_data_n(:,1), conf_data_n(:,2), '--','Color',c1);
    elseif count <3
        plot(conf_data_p(:,1), conf_data_p(:,2), '-','Color',c2);
        plot(conf_data_n(:,1), conf_data_n(:,2), '--','Color',c2);
    else
        plot(conf_data_p(:,1), conf_data_p(:,2), '-r');
        plot(conf_data_n(:,1), conf_data_n(:,2), '--r');
    end
end

figure   % K plot
set(gcf,"position",[10,10,1280,720])

plot(-0.5:0.05:1.5,Stepwise_Unfold.K,'o-',"MarkerFaceColor","#0072BD",'MarkerSize',4,...
    "LineWidth",1)
set(gca,'xtick',[-0.5, 0, 0.5, 1, 1.5],'xticklabel', ...
    Stepwise_Unfold.Stepwise_unfold(1:10:end),'FontSize',14)
xlabel('Stepwise unfold degree')
ylabel('Precision parameter')
grid on

figure  % 2 group directions
set(gcf,"position",[10,10,1280,720])

for ele = 1:41
    step = (ele-11)/20;
    T1_cor = t_cor_step(T1,step,1);    
    
    inc = T1_cor.Is;
    th = 90-inc;    %   co-inclination
    n = length(inc);
    dr = pi/180;
    [aralevOutput_ainc,~,aralevOutput_a95,~] = ARALEV(inc);
    G1_Inc(ele) = aralevOutput_ainc;  
    G1_a95(ele) = aralevOutput_a95;   
    G1_step(ele) = strcat(string(step*100),"%");  
end

for ele = 1:41
    step = (ele-11)/20;
    T2_cor = t_cor_step(T2,step,1);    
    
    inc = T2_cor.Is;
    th = 90-inc;    %   co-inclination
    n = length(inc);
    dr = pi/180;
    [aralevOutput_ainc,~,aralevOutput_a95,~] = ARALEV(inc);
    G2_Inc(ele) = aralevOutput_ainc;
    G2_a95(ele) = aralevOutput_a95;
    G2_step(ele) = strcat(string(step*100),"%");  
end

x = -0.5:0.05:1.5;
y1 = G1_Inc-G1_a95;
y2 = G1_Inc+G1_a95;
y3 = G2_Inc-G2_a95;
y4 = G2_Inc+G2_a95;


% find intersect point
try
    y_diff1 = y1 - y3;
    idx1 = find(y_diff1(1:end-1) .* y_diff1(2:end) < 0);
    x_intersect1 = zeros(size(idx1));
    y_intersect1 = zeros(size(idx1));
    x_intersect1 = interp1(y_diff1(idx1:idx1+1), x(idx1:idx1+1), 0);
    y_intersect1 = interp1(x, y1, x_intersect1);
    y_diff2 = y1 - y4;
    idx2 = find(y_diff2(1:end-1) .* y_diff2(2:end) < 0);
    x_intersect2 = zeros(size(idx2));
    y_intersect2 = zeros(size(idx2));
    x_intersect2 = interp1(y_diff2(idx2:idx2+1), x(idx2:idx2+1), 0);
    y_intersect2 = interp1(x, y1, x_intersect2);
    y_diff3 = y2 - y3;
    idx3 = find(y_diff3(1:end-1) .* y_diff3(2:end) < 0);
    x_intersect3 = zeros(size(idx3));
    y_intersect3 = zeros(size(idx3));
    x_intersect3 = interp1(y_diff3(idx3:idx3+1), x(idx3:idx3+1), 0);
    y_intersect3 = interp1(x, y2, x_intersect3);
    y_diff4 = y2 - y4;
    idx4 = find(y_diff4(1:end-1) .* y_diff4(2:end) < 0);
    x_intersect4 = zeros(size(idx4));
    y_intersect4 = zeros(size(idx4));
    x_intersect4 = interp1(y_diff4(idx4:idx4+1), x(idx4:idx4+1), 0);
    y_intersect4 = interp1(x, y2, x_intersect4);
    
    x_int = [x_intersect1, x_intersect2, x_intersect3, x_intersect4];
    y_int = [y_intersect1, y_intersect2, y_intersect3, y_intersect4];
    [xmin, min_ind] = min(x_int);
    [xmax, max_ind] = max(x_int);
    label1 = strcat(string(round(xmin*100,1)),"%");
    label2 = strcat(string(round(xmax*100,1)),"%");
end

pic01 = fill([x,fliplr(x)],[y1,fliplr(y2)],c1);
set(pic01,'edgealpha', 0, 'facealpha', 0.4)
hold on; grid on;
plot(x,G1_Inc,'-o',"Color",c1,"MarkerFaceColor",c1,"MarkerSize",4,...
    "LineWidth",1);
pic02 = fill([x,fliplr(x)],[y3,fliplr(y4)],c2);
set(pic02,'edgealpha', 0, 'facealpha', 0.4)
plot(x,G2_Inc,'-or',"color",c2,"MarkerFaceColor",c2,"MarkerSize",4,...
    "LineWidth",1);
yy = ylim;
intv = (yy(2)-yy(1))/25;
try
    plot([xmin,xmin],[yy(1),y_int(min_ind)],'--k','LineWidth',1);
    plot([xmax,xmax],[yy(1),y_int(max_ind)],'--k','LineWidth',1);
    text([xmin+0.03,xmax+0.03],[yy(1)+intv,yy(1)+intv],[label1,label2],'FontSize',14)
end
text(0,G1_Inc(11)-2*intv,Group(1),'Color',c1,'FontSize',14)
text(0,G2_Inc(11)-2*intv,Group(2),'Color',c2,'FontSize',14)
xlabel("Stepwise unfolded degree");
ylabel("Inclination (\circ)")
set(gca,'xtick',[-0.5, 0, 0.5, 1, 1.5],'xticklabel', ...
    Stepwise_Unfold.Stepwise_unfold(1:10:end),'FontSize',14)

name = [Group;"average"];
N = [height(T1);height(T2);height(T)];
Tilt_Correction = table(name,N,Dg,Ig,Kg,a95g,Ds,Is,Ks,a95s);
writetable(Tilt_Correction,strcat(filename,'_Tilt_Correction.csv'),'Delimiter',',');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% t_cor_step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function T = t_cor_step(T,step,fn)
    switch fn
        case 1  %   in-situ to tilt correvtion
            D_g = T.Dg;  %   Declination in-situ in degree
            I_g = T.Ig;  %   Inclination in-situ in degree
            strike = T.Strike_R;   %   strike(R) in degree
            dip = T.Dip.*step;  %   dip in degree
            Mx = sind(90-I_g).*cosd(D_g);
            My = sind(90-I_g).*sind(D_g);
            Mz = cosd(90-I_g);
            MXT = Mx.*((cosd(strike)).^2+(1-(cosd(strike).^2)).*cosd(dip))+...
                My.*(cosd(strike).*cosd(90-strike).*(1-cosd(dip)))+...
                Mz.*(-cosd(90-strike).*sind(dip));
            MYT = Mx.*(cosd(strike).*cosd(90-strike).*(1-cosd(dip)))+...
                My.*(cosd(90-strike).^2+(1-cosd(90-strike).^2).*cosd(dip))+...
                Mz.*(cosd(strike).*sind(dip));
            MZT = Mx.*cosd(90-strike).*sind(dip)+...
                My.*(-cosd(strike).*sind(dip))+...
                Mz.*cosd(dip);
            D_s = atan2(MYT,MXT)*180/pi;
            for i = 1:length(D_s)
                if D_s(i) < 0
                    D_s(i) = D_s(i)+360;
                end
            end
            T.Ds = D_s;
            I_s = atan2(MZT,sqrt(MXT.^2+MYT.^2))*180/pi;
            T.Is = I_s;
        case 2  %   tilt correction to in-situ
            D_s = T.Ds;  %   Declination in-situ in degree
            I_s = T.Is;  %   Inclination in-situ in degree
            strike = T.Strike_R;   %   strike(R) in degree
            dip = T.Dip.*step;  %   dip in degree
            Mx = sind(90-I_s).*cosd(D_s);
            My = sind(90-I_s).*sind(D_s);
            Mz = cosd(90-I_s);
            MXT = Mx.*((cosd(strike-180)).^2+(1-(cosd(strike-180).^2)).*cosd(dip))+...
                My.*(cosd(strike-180).*cosd(270-strike).*(1-cosd(dip)))+...
                Mz.*(-cosd(270-strike).*sind(dip));
            MYT = Mx.*(cosd(strike-180).*cosd(270-strike).*(1-cosd(dip)))+...
                My.*(cosd(270-strike).^2+(1-cosd(270-strike).^2).*cosd(dip))+...
                Mz.*(cosd(strike-180).*sind(dip));
            MZT = Mx.*cosd(270-strike).*sind(dip)+...
                My.*(-cosd(strike-180).*sind(dip))+...
                Mz.*cosd(dip);
            D_g = atan2(MYT,MXT)*180/pi;
            for i = 1:length(D_g)
                if D_g(i) < 0
                    D_g(i) = D_g(i)+360;
                end
            end
            T.Dg = D_g;
            I_g = atan2(MZT,sqrt(MXT.^2+MYT.^2))*180/pi;
            T.Ig = I_g;
    end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% t_cor_step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ArithmeticMean %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [ArithmeticMean_rt,ArithmeticMean_rk] = ArithmeticMean(th)
    s = sum(th);
    s2 = sum(th.^2);
    ArithmeticMean_rt = s/n;
    ArithmeticMean_rk = (n-1)/((s2-s^2/n)*dr^2);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ArithmeticMean %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARALEV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [aralevOutput_ainc,aralevOutput_ak,aralevOutput_t63,...
            aralevOutput_a95,aralevOutput_ierr] = ARALEV(xinc)
    t63max = 105.070062145; %   63% of a sphere.
    a95max = 154.158067237; %   95% of a sphere.
    aralevOutput_ierr = 1;

    % Check for illegal use
    if n == 1
        error('ERROR: Only one or none observed inclination in ARALEV')
    end

     %   Check if inc are out of range
    for i = 1:n
        if xinc(i) < -90 || xinc(i)>90
            error('ERROR: Inclination data out of range [-90, +90] in ARALEV')
        end
    end

    %   Check if all incl are identical
    same = true;
    for i = 1:n
        if xinc(i) ~= xinc(1)
            same = false;
        end
    end
    if same
        error('WARNING: All incl identical in ARALEV')
    end

    %   Calculate arithmetic mean to use as first guess
    c = mean(sum(cosd(th)));
    s = sum(th);
    s2 = sum(th.^2);
    rt = s/n;
    x = (s2-s^2/n)*dr^2;
    rk = 1e10;
    if x/(n-1) > 1e-10
        rk = (n-1)/x;
    end
    rt1 = rt;
    rk1 = rk;

    %   Iterate in the interior to find solution to (theta, kappa) 
    % Start iteration at arithmetic mean (theta, kappa)
    rt = rt1;
    rk = rk1;
    ie1 = 0;

    the1 = rt;
    akap1 = rk;
    conv = false;
    for j = 1:10000
        rt = AL1(th,rt,rk);
        rk = AL2(th,rt,rk);
        dt = abs(rt-the1);
        dk = abs((rk-akap1)/rk);
        the1 = rt;
        akap1 = rk;
        if j > 10 && dt < 1e-6 && dk < 1e-6
            conv = true;
            break;
        end
    end
    if(~conv)
        ie1 = 1;
    end
    
    the1 = rt;
    akap1 = rk;
    xl1 = Xlik(th,rt,rk);
    
    %   Find the maximum on the edge (theta = 0)
    rt = 0;
    rk = rk1;
    x = 1-c;
    if x > 1e-10
        rk = 1/x;
    end
    ie2 = 0;
    akap2 = rk;
    conv = false;
    
    for j = 1:10000
        x = Coth(rk) - c;
        if x > 1e-10
            rk = 1/x;
        else
            rk = 1e10;
        end
        dk = abs((rk-akap2)/rk);
        akap2 = rk;
        if j > 4 && dk < 1e-6
            conv = true;
            break;
        end
    end
    if (~conv)
        ie2 = 1;
    end
    the2 = 0;
    akap2 = rk;
    xl2 = Xlik(th,rt,rk);
    
    %   Find the maximum on the edge (theta = 180)
    rt = 180;
	rk = rk1;
	x = 1+c;
    if x > 1e-10
        rk = 1/x;
    end
    ie3 = 0;
    akap3 = rk;
    conv = false;
    for j = 1:10000
        x = Coth(rk)+c;
        if x > 1e-10
            rk = 1/x;
        else
            rk = 1e10;
        end
        dk = abs((rk-akap3)/rk);
        akap3 = rk;
        if j > 4 && dk < 1e-6
            conv = true;
            break;
        end
    end
    if(~conv)
        ie3 = 1;
    end
    the3 = 180;
    akap3 = rk;
    xl3 = Xlik(th,rt,rk);
    
    %   Find the maximum on the edge (kappa = 0)
    rt = 90;
    rk = 0;
    the4 = rt;
    akap4 = rk;
    xl4 = Xlik(th,rt,rk);
    
    %   Use the best solution of the four
    isol = 1;
    aralevOutput_ierr = ie1;
    if xl2 > xl1
        the1 = the2;
        akap1 = akap2;
        xl1 = xl2;
        isol = 2;
        aralevOutput_ierr = 1;
        if ie2 == 0
            aralevOutput_ierr = 0;
        end
    end
    if xl3 > xl1
        the1 = the3;
        akap1 = akap3;
        xl1 = xl3;
        isol = 3;
        aralevOutput_ierr = 1;
        if ie3 == 0
            aralevOutput_ierr = 0;
        end
    end
    if xl4 > xl1
        the1 = the4;
        akap1 = akap4;
        xl1 = xl4;
        isol = isol4;
        aralevOutput_ierr = 0;
    end
    aralevOutput_ainc = 90-the1;
    aralevOutput_ak = akap1;
    if aralevOutput_ierr ~= 0
        disp('WARNING: Convergence problems in ARALEV');
    end
    
    %   Test robustness of solution theta +/- 0.01 deg and kappa +/- 0.1%
    for x = 1:16
        rt = the1 +0.001*cos(22.5*x*dr);
        rk = akap1*(1+0.001*sin(22.5*x*dr));
        if rt >= 0 && rt <= 180
            xl = Xlik(th,rt,rk);
            if xl > xl1
                aralevOutput_ierr = aralevOutput_ierr+2;
                disp('WARNING: Robustness problem in ARALEV');
            end
        end
    end
    
    %   Estimation of angular standard deviation
    % c	Theta-63 calculated from (kappa), same method as Kono (1980)
    if akap1 >= 20
        co = 1+log(1-0.63)/akap1;
    end
    if akap1 > 0.1 && akap1 < 20
        co = 1+log(1-0.63*(1-exp(akap1)))/akap1;
    end
    if akap1 <= 0.1
        co = -0.26+0.4662*akap1;
    end
    aralevOutput_t63 = 0;
    if co < 0
        aralevOutput_t63 = 180;
    end
    if abs(co) < 1
        aralevOutput_t63 = 90-atan(co/sqrt(1-co^2))/dr;
    end
    if aralevOutput_t63 > t63max
        aralevOutput_t63 = t63max;
    end
    
    %   Estimation of 95% (circular) symmetric confidence limit of the mean
    % Alpha-95 calculated from (N, kappa), same method as Kono (1980)
    co = 1-(n-1)*(20^(1/(n-1))-1)/(n*(akap1-1)+1);
    aralevOutput_a95 = 0;
    if co < 0
        aralevOutput_a95 = 180;
    end
    if abs(co) < 1
        aralevOutput_a95 = 90-atan(co/sqrt(1-co^2))/dr;
    end
    if aralevOutput_a95 > a95max
        aralevOutput_a95 = a95max;
    end
        
    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AL1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function al1 = AL1(th,the,ak)
        s = 0;
        c = 0;
        for i = 1:n
            x = ak*sin(the*dr)*sin(th(i)*dr);
            [bi0e,bi1e,bi1i0] = Bessel(x);
            s = s+sin(th(i)*dr)*bi1i0;
            c = c+cos(th(i)*dr);
        end
        al1 = atan2(s,c)/dr;
        if al1 < 1e-6
            al1 = 1e-6;
        end
        if al1 > 180-1e-6
            al1 = 180-1e-6;
        end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AL1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AL2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function al2 = AL2(th,the,ak)
        s = 0;
        c = 0;
        for i = 1:n
            x=ak*sin(the*dr)*sin(th(i)*dr);
            [bi0e,bi1e,bi1i0] = Bessel(x);
            s = s+sin(th(i)*dr)*bi1i0;
            c = c+cos(th(i)*dr);
        end
        x = n*Coth(ak)-cos(the*dr)*c-sin(the*dr)*s;
        al2 = 1e10;
        if x/n >1e-10
            al2 = n/x;
        end
        if x/n<1e-6
            al2 = 1e-6;
        end 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AL2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bessel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [bi0e,bi1e,bi1i0] = Bessel(x)
        p1 = 1;
        p2 = 3.5156229;
        p3 = 3.0899424;
        p4 = 1.2067492;
        p5 = 0.2659732;
        p6 = 0.360768e-1;
        p7 = 0.45813e-2;

        q1 = 0.39894228;
        q2 = 0.1328592e-1;
        q3 = 0.225319e-2;
        q4 = -0.157565e-2;
        q5 = 0.916281e-2;
        q6 = -0.2057706e-1;
        q7 = 0.2635537e-1;
        q8 = -0.1647633e-1;
        q9 = 0.392377e-2;

        u1 = 0.5;
        u2 = 0.87890594;
        u3 = 0.51498869;
        u4 = 0.15084934;
        u5 = 0.2658733e-1;
        u6 = 0.301532e-2;
        u7 = 0.32411e-3;

        v1 = 0.39894228;
        v2 = -0.3988024e-1;
        v3 = -0.362018e-2;
        v4 = 0.163801e-2;
        v5 = -0.1031555e-1;
        v6 = 0.2282967e-1;
        v7 = -0.2895312e-1;
        v8 = 0.1787654e-1;
        v9 = -0.420059e-2;
        if abs(x) < 3.75
            t = (x/3.75)^2;
            b0 = p1+t*(p2+t*(p3+t*(p4+t*(p5+t*(p6+t*p7)))));
            b1 = x*(u1+t*(u2+t*(u3+t*(u4+t*(u5+t*(u6+t*u7))))));
            bi0e = b0/exp(abs(x));
            bi1e = b1/exp(abs(x));
            bi1i0 = b1/b0;
        else
            t = 3.75/abs(x);
            b0 = q1+t*(q2+t*(q3+t*(q4+t*(q5+t*(q6+t*(q7+t*(q8+t*q9)))))));
            b1 = v1+t*(v2+t*(v3+t*(v4+t*(v5+t*(v6+t*(v7+t*(v8+t*v9)))))));
            if x < 0
                b1 = -b1;
            end
            bi0e = b0/sqrt(abs(x));
            bi1e = b1/sqrt(abs(x));
            bi1i0 = b1/b0;
        end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Bessel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function COTH = Coth(x)
        if x == 0
            COTH = 0;
        end
        t = abs(x);
        if t < 0.001
            COTH = 1/t+t/3-t^3/45+t^5*2/945;
        elseif t <= 15
            ep = exp(t);
            em = exp(-t);
            COTH = (ep+em)/(ep-em);
        else
            COTH = 1;
        end
        if x <0
            COTH = -COTH;
        end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Xlik %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function XLIK = Xlik(th,the,ak)
        %   Illeage use
        if n < 1
            XLIK = -1e10;
        end
        if ak < 0
            XLIK = -1e10;
        end
        a1 = 0;
        if ak >= 0 && ak <0.001
            q = -ak*(1-ak*(2/3-ak*(1/3-ak*(2/15-ak*(8/45)))));
            a1 = n*(-log(2)-log(1+q)-ak);
        elseif ak>= 0.01 && ak <= 15
            a1 = n*(log(ak)-log(1-exp(-2*ak))-ak);
        else
            a1 = n*(log(ak)-ak);
        end 
        a2 = 0;
        for i = 1:n
            x = ak*(sin(the*dr)*sin(th(i)*dr));
            [bi0e,bi1e,bi1i0] = Bessel(x);
            a2 = a2 + ak*cos((th(i)-the)*dr)+log(bi0e);
        end
        a3 = 0;
        for i = 1:n
            x = th(i);
            if x < 1e-6
                x = 1e-6;
            end
            if x >180-1e-6
                x = 180-1e-6;
            end
            a3 = a3+log(sin(x*dr));
        end
        XLIK = a1+a2+a3;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Xlik %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARALEV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% a95calc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [D_mean,I_mean,K,a95] = a95calc(Dec,Inc)
        len = length(Dec);
        x = cosd(Dec).*cosd(Inc);
        y = sind(Dec).*cosd(Inc);
        z = sind(Inc);
        x_mean = mean(x);
        y_mean = mean(y);
        z_mean = mean(z);
        D_mean = atan2d(y_mean,x_mean);
        I_mean = atan2d(z_mean,norm([x_mean,y_mean]));
        R = norm([sum(x),sum(y),sum(z)]);
        K = (len-1)/(len-R);
        a95 = acosd(1-((len-R)/R*((1/0.05).^(1/(len-1))-1)));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% a95calc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end