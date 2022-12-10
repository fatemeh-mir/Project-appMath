clear All
clear clc

%% A:Start

%% The Binary Systems:
%% Benzene-Toluene --> S1
%% Methanol-Water --> S2
%% Heptane-ethylBenzene --> S3
%% Methanol-Toluene --> S4
prompt = 'Select The Binary System Between S1, S2, S3, S4 : ';
sys = input(prompt,'s');

%% input F, q, zF, xD, xW
Feed_rate = 'Enter F:';
F = input(Feed_rate); %kmol/h

feed_thermal_condition = 'Enter feed thermal condition:';
q = input(feed_thermal_condition);

fraction_of_more_volatile_component_at_the_feed = 'Enter fraction of more volatile component at the feed:';
zF = input(fraction_of_more_volatile_component_at_the_feed);

fraction_of_more_volatile_component_at_the_top_product = 'Enter The fraction of more volatile component at the top product:';
xD = input(fraction_of_more_volatile_component_at_the_top_product);

fraction_of_more_volatile_component_at_the_bottom_product = 'Enter The fraction of more volatile component at the bottom product:';
xW = input(fraction_of_more_volatile_component_at_the_bottom_product);

%% Selection of the condenser type(Partial or total condenser)
condenser = 'Enter type of condenser(JUST ENTER t OR p) :';
c = input(condenser, 's');

reflux_ratio = 'Enter reflux ratio:'; 
r = input(reflux_ratio);

%%Top product rate
D = (zF * F - xW * F)/(xD - xW);
fprintf('Top product rate is: %f',D);
fprintf('\n')


%%Bottom product rate
W = F - D;
fprintf('Bottom product rate is: %f',W);
fprintf('\n')


%%S2
if sys == 'S2'
    filename = fullfile('This PC/Downloads/UO_PRJ2_96103852','methanol-water.txt');
    fileID = fopen(filename);
    C = textscan(fileID,'%f %f');
    fclose(fileID);
    x1 = C{1};
    y1 = C{2};
    
    for i = 1 : length(x1) - 1
    
        j = 5 * i - 4;
        xeq(j) = x1(i);
        xeq(j+1) =  xeq(j) + (x1(i+1) - x1(i))/5;
        xeq(j+2) =  xeq(j) + 2 * (x1(i+1) - x1(i))/5;
        xeq(j+3) =  xeq(j) + 3 * (x1(i+1) - x1(i))/5;
        xeq(j+4) =  xeq(j) + 4 * (x1(i+1) - x1(i))/5;
    end
    xeq(1, 251) = 1;
    
    for i = 1 : length(y1) - 1
    
        j = 5 * i - 4;
        yeq(j) = y1(i);
        yeq(j+1) =  yeq(j) + (y1(i+1) - y1(i))/5;
        yeq(j+2) =  yeq(j) + 2 * (y1(i+1) - y1(i))/5;
        yeq(j+3) =  yeq(j) + 3 * (y1(i+1) - y1(i))/5;
        yeq(j+4) =  yeq(j) + 4 * (y1(i+1) - y1(i))/5;
    end
    yeq(1, 251) = 1;
    
    %%draw equil curve.
    %%draw x=y curve.
    origin = [0 0];
    destination = [1 1];
    plot (xeq, yeq, 'b', [origin(1), destination(1)], [origin(2), destination(2)], 'r')
    xlabel('x_M_e_t_h_n_o_l(mole fraction)')
    ylabel('y_M_e_t_h_n_o_l(mole fraction)')  
    
end

%%S1
if sys == 'S1'
    filename = fullfile('This PC/Downloads/UO_PRJ2_96103852','benzene-toluene.txt');
    fileID = fopen(filename);
    C = textscan(fileID,'%f %f');
    fclose(fileID);
    x1 = C{1};
    y1 = C{2};
    
    for i = 1 : length(x1) - 1
    
        j = 5 * i - 4;
        xeq(j) = x1(i);
        xeq(j+1) =  xeq(j) + (x1(i+1) - x1(i))/5;
        xeq(j+2) =  xeq(j) + 2 * (x1(i+1) - x1(i))/5;
        xeq(j+3) =  xeq(j) + 3 * (x1(i+1) - x1(i))/5;
        xeq(j+4) =  xeq(j) + 4 * (x1(i+1) - x1(i))/5;
    end
    xeq(1, 251) = 1;
    
    for i = 1 : length(y1) - 1
    
        j = 5 * i - 4;
        yeq(j) = y1(i);
        yeq(j+1) =  yeq(j) + (y1(i+1) - y1(i))/5;
        yeq(j+2) =  yeq(j) + 2 * (y1(i+1) - y1(i))/5;
        yeq(j+3) =  yeq(j) + 3 * (y1(i+1) - y1(i))/5;
        yeq(j+4) =  yeq(j) + 4 * (y1(i+1) - y1(i))/5;
    end
    yeq(1, 251) = 1;
    
    %%draw equil curve.
    %%draw x=y curve.
    origin = [0 0];
    destination = [1 1];
    plot (xeq, yeq, 'b', [origin(1), destination(1)], [origin(2), destination(2)], 'r')
    xlabel('x_B_e_n_z_e_n_e(mole fraction)')
    ylabel('y_B_e_n_z_e_n_e(mole fraction)')
    
end

%%S3
if sys == 'S3'
    filename = fullfile('This PC/Downloads/UO_PRJ2_96103852','nheptane-ebenzene.txt');
    fileID = fopen(filename);
    C = textscan(fileID,'%f %f');
    fclose(fileID);
    x1 = C{1};
    y1 = C{2};
   
    for i = 1 : length(x1) - 1
    
        j = 5 * i - 4;
        xeq(j) = x1(i);
        xeq(j+1) =  xeq(j) + (x1(i+1) - x1(i))/5;
        xeq(j+2) =  xeq(j) + 2 * (x1(i+1) - x1(i))/5;
        xeq(j+3) =  xeq(j) + 3 * (x1(i+1) - x1(i))/5;
        xeq(j+4) =  xeq(j) + 4 * (x1(i+1) - x1(i))/5;
    end
    xeq(1, 251) = 1;
    
    for i = 1 : length(y1) - 1
    
        j = 5 * i - 4;
        yeq(j) = y1(i);
        yeq(j+1) =  yeq(j) + (y1(i+1) - y1(i))/5;
        yeq(j+2) =  yeq(j) + 2 * (y1(i+1) - y1(i))/5;
        yeq(j+3) =  yeq(j) + 3 * (y1(i+1) - y1(i))/5;
        yeq(j+4) =  yeq(j) + 4 * (y1(i+1) - y1(i))/5;
    end
    yeq(1, 251) = 1;
    
    %%draw equil curve.
    %%draw x=y curve.
    origin = [0 0];
    destination = [1 1];
    plot (xeq, yeq, 'b', [origin(1), destination(1)], [origin(2), destination(2)], 'r')
    xlabel('x_H_e_p_t_a_n_e(mole fraction)')
    ylabel('y_H_e_p_t_a_n_e(mole fraction)')
    
end

%%S4
if sys == 'S4'
    filename = fullfile('This PC/Downloads/UO_PRJ2_96103852','methanol-toluene.txt');
    fileID = fopen(filename);
    C = textscan(fileID,'%f %f');
    fclose(fileID);
    x1 = C{1};
    y1 = C{2};
    
    for i = 1 : length(x1) - 1
    
        j = 5 * i - 4;
        xeq(j) = x1(i);
        xeq(j+1) =  xeq(j) + (x1(i+1) - x1(i))/5;
        xeq(j+2) =  xeq(j) + 2 * (x1(i+1) - x1(i))/5;
        xeq(j+3) =  xeq(j) + 3 * (x1(i+1) - x1(i))/5;
        xeq(j+4) =  xeq(j) + 4 * (x1(i+1) - x1(i))/5;
    end
    xeq(1, 251) = 1;
    
    for i = 1 : length(y1) - 1
    
        j = 5 * i - 4;
        yeq(j) = y1(i);
        yeq(j+1) =  yeq(j) + (y1(i+1) - y1(i))/5;
        yeq(j+2) =  yeq(j) + 2 * (y1(i+1) - y1(i))/5;
        yeq(j+3) =  yeq(j) + 3 * (y1(i+1) - y1(i))/5;
        yeq(j+4) =  yeq(j) + 4 * (y1(i+1) - y1(i))/5;
    end
    yeq(1, 251) = 1;
   
    %%draw equil curve.
    %%draw x=y curve.
    origin = [0 0];
    destination = [1 1];
    plot (xeq, yeq, 'b', [origin(1), destination(1)], [origin(2), destination(2)], 'r')
    xlabel('x_M_e_t_h_n_o_l(mole fraction)')
    ylabel('y_M_e_t_h_n_o_l(mole fraction)')
    
end


hold on
%% xD and xW on plot
plot(xD, xD, 'g.', 'markersize', 15);
plot(xW, xW, 'g.', 'markersize', 15);
plot(zF, zF, 'g.', 'markersize', 15);

plot([xD,xD], [xD,0],'g--')
plot([xW,xW], [xW,0],'g--')
plot([zF,zF], [zF,0],'g--')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% q line.

if q==1
    yf = 0.5 : 0.1 :0.9;
    xf = zF * ones(1,401);
    plot([zF, zF],[zF, 0.9])
    for i = 1 : 251
        eror(i) = abs(xeq(i) - zF);
        if eror(i)<= 0.0025
            int_y_min = yeq(i);
            int_x_min = zF;
            break
        end
    end
    plot(int_x_min,int_y_min,'g.', 'markersize', 15)
    %% find R
    s_min = (xD - int_y_min)/(xD - int_x_min);
    r_min = (s_min)/(1 - s_min);
    R = r_min * r;

    %L1 and V1
    L1 = R*D;
    L = R * D + q * F;
%     fprintf('the rate of cooling water as the coolant: %f',L);
%     fprintf('\n')   

    V1 = (R+1)* D;
    V = V1 + (1-q)*F;
%     fprintf('the rate of LPS  %f',V);
%     fprintf('\n')

    %% op and bottom line
    yf_int = (R/(R+1)) * zF + (xD/(R+1));
    plot(zF, yf_int,'g.', 'markersize', 15)
    plot([zF, xD], [yf_int, xD])
    plot([zF,xW],[yf_int,xW])
    
    xp(1) = xD;
    yp(1) = xD;
    tray_num = 1;
    
    for i = 1:200
        n=1;
        for j = 1 : 251
            if(xp(i) >= zF)
                eror(j) = abs(yeq(j) - ((R/(R+1)) * xp(i) + xD/(R+1)));
                if(eror(j) <= 0.0025)&&(n==1)
                    xp(i+1) = xeq(j);
                    if(xp(i+1) >= zF)
                    yp(i+1) = (R/(R+1)) * xp(i+1) + xD/(R+1);
                    plot([xp(i), xp(i+1)],[yp(i), yp(i)])
                    plot([xp(i+1),xp(i+1)],[yp(i+1),yp(i)])
                    feed_tray = i+1;
                    end
                    n = n+1;
                end
            else
              eror(j) = abs(yeq(j) - (((R * D + q * F)/((R+1)*D)) * xp(i) - (xW * W)/((R+1)*D))); 
              if(eror(j) <= 0.005)
                  yp(i) = (((R * D + q * F)/((R+1)*D)) * xp(i) - W*xW /((R+1)*D));
                  xp(i+1) = xeq(j);
                  yp(i+1) = (((R * D + q * F)/((R+1)*D)) * xp(i+1) - W*xW /((R+1)*D));
                  plot([xp(i), xp(i+1)],[yp(i), yp(i)])
                  plot([xp(i+1),xp(i+1)],[yp(i+1),yp(i)])
                  plot([xp(i), xp(i)],[yp(i-1),yp(i)])
                  plot([xp(i-1), xp(i)],[yp(i-1),yp(i-1)])
                  break
              end
            end
        end
        if(xp(i+1)<=xW)
            break
        end
        tray_num = tray_num + 1;
    end
    
end

if q==0 
    yf = zF * ones(1,401);
    xf = 0.5 : 0.9 :0.1;
    plot([zF, 0.1],[zF, zF])
     for i = 1 : 251
        eror(i) = abs(yeq(i) - zF);
        if eror(i)<= 0.0025
            int_y_min = zF;
            int_x_min = xeq(i);
            break
        end
    end
    plot(int_x_min,int_y_min,'g.', 'markersize', 15)
    s_min = (xD - int_y_min)/(xD - int_x_min);
    r_min = (s_min)/(1 - s_min);
    R = r_min * r;

    %L1 and V1
    L1 = R*D;
    L = R * D + q * F;
%     fprintf('the rate of cooling water as the coolant: %f',L);
%     fprintf('\n')

    V1 = (R+1)* D;
    V = V1 + (1-q)*F;
%     fprintf('the rate of LPS  %f',V);
%     fprintf('\n')
    %% op and bottom line
    xf_int = (zF - (xD/(R+1)))/(R/(R+1));
    plot(xf_int, zF,'g.', 'markersize', 15)
    plot([xf_int, xD], [zF, xD])
    plot([xf_int,xW],[zF,xW])
    
    xp(1) = xD;
    yp(1) = xD;
    tray_num = 1;
    for i = 1:200
        n=1;
        for j = 1 : 251
            if xp(i) >= (zF - (xD/(R+1)))/(R/(R+1))
                eror(j) = abs(yeq(j) - ((R/(R+1)) * xp(i) + xD/(R+1)));
                if(eror(j) <= 0.0025)&& (n==1)
                    xp(i+1) = xeq(j);
                    if(xp(i+1) >= (zF - (xD/(R+1)))/(R/(R+1)))
                    yp(i+1) = (R/(R+1)) * xp(i+1) + (xD/(R+1));
                    plot([xp(i), xp(i+1)],[yp(i), yp(i)])
                    plot([xp(i+1),xp(i+1)],[yp(i+1),yp(i)])
                    feed_tray = i+1;
                    end
                    n = n+1;
                end
            else
              eror(j) = abs(yeq(j) - (((R * D + q * F)/((R+1)*D)) * xp(i) - (xW * W)/((R+1)*D))); 
              if(eror(j) <= 0.005)
                  yp(i) = (((R * D + q * F)/((R+1)*D)) * xp(i) - W*xW /((R+1)*D));
                  xp(i+1) = xeq(j);
                  yp(i+1) = (((R * D + q * F)/((R+1)*D)) * xp(i+1) - W*xW /((R+1)*D));
                  plot([xp(i), xp(i+1)],[yp(i), yp(i)])
                  plot([xp(i+1),xp(i+1)],[yp(i+1),yp(i)])
                  plot([xp(i), xp(i)],[yp(i-1),yp(i)])
                  plot([xp(i-1), xp(i)],[yp(i-1),yp(i-1)])
                  break
              end
            end
        end
        if(xp(i+1)<=xW)
            break
        end
        tray_num = tray_num + 1;
    end

end

if q>0 && q<1
    yf = (-zF)/(q-1);
    xf = 0;
    plot([xf,zF],[yf,zF], 'g')
    
    %L1 and V1
    L1 = R*D;
    L = R * D + q * F;


    V1 = (R+1)* D;
    V = V1 + (1-q)*F;

    
    int_x_min = (((-xD)/(r+1))-(zF/(q-1)))/((r/(r+1))-(q/(q-1)));
    int_y_min = (r/(r+1)) * int_x_min + (xD / (r+1));

    plot(int_x_min,int_y_min,'r.', 'markersize', 15)
    plot([int_x_min, xD], [int_y_min, xD])
    plot([int_x_min, xW], [int_y_min, xW])
    
    xp(1) = xD;
    yp(1) = xD;
    tray_num = 1;
    for i = 1:200
        n=1;
        for j = 1 : 251
            if(xp(i) >= int_x_min)
                eror(j) = abs(yeq(j) - ((R/(R+1)) * xp(i) + xD/(R+1)));
                if(eror(j) <= 0.0025)&& (n==1)
                    xp(i+1) = xeq(j);
                    if(xp(i+1) >= int_x_min)
                    yp(i+1) = (R/(R+1)) * xp(i+1) + (xD/(R+1));
                    plot([xp(i), xp(i+1)],[yp(i), yp(i)])
                    plot([xp(i+1),xp(i+1)],[yp(i+1),yp(i)])
                    feed_tray = i+1;
                    end
                    n = n+1;
                end
            else
              eror(j) = abs(yeq(j) - (((R * D + q * F)/((R+1)*D)) * xp(i) - (xW * W)/((R+1)*D))); 
              if(eror(j) <= 0.005)
                  yp(i) = (((R * D + q * F)/((R+1)*D)) * xp(i) - W*xW /((R+1)*D));
                  xp(i+1) = xeq(j);
                  yp(i+1) = (((R * D + q * F)/((R+1)*D)) * xp(i+1) - W*xW /((R+1)*D));
                  plot([xp(i), xp(i+1)],[yp(i), yp(i)])
                  plot([xp(i+1),xp(i+1)],[yp(i+1),yp(i)])
                  plot([xp(i), xp(i)],[yp(i-1),yp(i)])
                  plot([xp(i-1), xp(i)],[yp(i-1),yp(i-1)])
                  break
              end
            end
        end
        if(xp(i+1)<=xW)
            break
        end
        tray_num = tray_num + 1;
    end
    
end
    



if sys == 'S1'
    h_B = 35570; %kj/kmol
    h_T = 37000; %kj/kmol
    m_reboiler = (h_B * xW + h_T * (1-xW)) * V / 2133.8;
    m_condenser = (h_B * xD + h_T * (1-xD)) * V1 / (30 * 4.2);
    if c == 't' 
        C_duty = (h_B * xD + h_T * (1-xD)) * V1;
        fprintf('Number of theoretical stages: %d',tray_num);
        fprintf('\n')
        fprintf('feed tray: %d',feed_tray);
        fprintf('\n')
    end
    if c == 'p'
        C_duty = (h_B * xD + h_T * (1-xD)) * V; 
        fprintf('Number of theoretical stages: %d',tray_num - 1);
        fprintf('\n')
        fprintf('feed tray: %d',feed_tray);
        fprintf('\n')
    end
    RB_duty = (h_B * xW + h_T * (1-xW)) * L;
    
end
if sys == 'S2'
    h_M = 35210; %kj/kmol
    h_W = 4065; %kj/kmol
    m_reboiler = (h_M * xW + h_W * (1-xW)) * V / 2133.8;
    m_condenser = (h_M * xD + h_W * (1-xD)) * V1 / (30 * 4.2);
    if c == 't' 
        C_duty = (h_M * xD + h_W * (1-xD)) * V1;
        fprintf('Number of theoretical stages: %d',tray_num);
        fprintf('\n')
        fprintf('feed tray: %d',feed_tray);
        fprintf('\n')
    end
    if c == 'p'
        C_duty = (h_M * xD + h_W * (1-xD)) * V; 
        fprintf('Number of theoretical stages: %d',tray_num - 1);
        fprintf('\n')
        fprintf('feed tray: %d',feed_tray);
        fprintf('\n')
    end
    RB_duty = (h_M * xW + h_W * (1-xW)) * L;
end
if sys == 'S3'
    h_H = 31666;
    h_EB = 42468;
    m_reboiler = (h_H * xW + h_EB * (1-xW)) * V / 2133.8;
    m_condenser = (h_H * xD + h_EB * (1-xD)) * V1 / (30 * 4.2); 
    if c == 't' 
        C_duty = (h_H * xD + h_EB * (1-xD)) * V1;
        fprintf('Number of theoretical stages: %d',tray_num);
        fprintf('\n')
        fprintf('feed tray: %d',feed_tray);
        fprintf('\n')
    end
    if c == 'p'
        C_duty = (h_H * xD + h_EB * (1-xD)) * V; 
        fprintf('Number of theoretical stages: %d',tray_num - 1);
        fprintf('\n')
        fprintf('feed tray: %d',feed_tray);
        fprintf('\n')
    end
    RB_duty = (h_H * xW + h_EB * (1-xW)) * L;
end
if sys == 'S4'
    h_M = 35210; %kj/kmol
    h_T = 37000; %kj/kmol
    m_reboiler = (h_M * xW + h_T * (1-xW)) * V / 2133.8;
    m_condenser = (h_M * xD + h_T * (1-xD)) * V1 / (30 * 4.2);
    if c == 't' 
        C_duty = (h_M * xD + h_T * (1-xD)) * V1;
        fprintf('Number of theoretical stages: %d',tray_num);
        fprintf('\n')
        fprintf('feed tray: %d',feed_tray);
        fprintf('\n')
    end
    if c == 'p'
        C_duty = (h_M * xD + h_T * (1-xD)) * V; 
        fprintf('Number of theoretical stages: %d',tray_num - 1);
        fprintf('\n')
        fprintf('feed tray: %d',feed_tray);
        fprintf('\n')
    end
    RB_duty = (h_M * xW + h_T * (1-xW)) * L;
end

fprintf('the rate of LPS:  %f',m_reboiler);
fprintf('\n')

fprintf('the rate of cooling water as the coolant: %f',m_condenser);
fprintf('\n')



fprintf('Condenser Duty: %f',C_duty);
fprintf('\n')



fprintf('Reboiler Dut: %f',RB_duty);
fprintf('\n')




