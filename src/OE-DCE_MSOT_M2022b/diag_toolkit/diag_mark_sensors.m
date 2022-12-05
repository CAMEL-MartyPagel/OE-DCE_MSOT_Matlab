function diag_mark_sensors(tfile,sel_left)
figure;
% 
% gap = 0.0006;                       % gap in mm between 2 halfs
% COR = 0.04;
% span = 270;
% begin = -45;
% 
% gapd = atand(gap/COR);              % gap in deg between elements
% perside = span/2 - gapd/2;          % active element in deg per side
% perelement = perside/elements*2;    % active element in deg
% 
% side1 = begin + perelement/2 : perelement : begin+perside;
% side2 = begin + span - perside + perelement / 2 : perelement : begin + span;
% sensor_angle = [side1 side2];

load('-mat',tfile);

%% plot all sensors
x = -[0 cos(angle_sensor)*r_sensor 0];
y = -[0 sin(angle_sensor)*r_sensor 0];
% figure;
hold off;
plot(x,y,'b*-','DisplayName','All sensors');

%% plot left sensor(s)
x = -[cos(angle_sensor(sel_left))*r_sensor];
y = -[sin(angle_sensor(sel_left))*r_sensor];
hold on;
plot(x,y,'r*','DisplayName','Marked sensors');
text(mean(x)+0.002,mean(y),num2str([sel_left(1) sel_left(end)]),'Color','red');


%% annotate graph
legend('show','location','North');
axis image;
xlabel('x [mm]');
ylabel('y [mm]');


% %% plot 5 spheres 
% plot(-0.009,0.0000,'o','Color','black','MarkerSize',2,'MarkerFaceColor','black','DisplayName','s1');
% plot(-0.006,0.0000,'o','Color','black','MarkerSize',2,'MarkerFaceColor','black','DisplayName','s1');
% plot(-0.0028,0.0000,'o','Color','black','MarkerSize',2,'MarkerFaceColor','black','DisplayName','s2');
% plot(0,0.0000,'o','Color','black','MarkerSize',2,'MarkerFaceColor','black','DisplayName','s3');
% plot(0.0028,0.0000,'o','Color','black','MarkerSize',2,'MarkerFaceColor','black','DisplayName','s4');
% plot(0.006,0.000,'o','Color','black','MarkerSize',2,'MarkerFaceColor','black','DisplayName','s5');
% plot(0.009,0.000,'o','Color','black','MarkerSize',2,'MarkerFaceColor','black','DisplayName','s5');

%%
clear side1 side2 x y elements begin span COR gap perelement perside gapd;