%clear
close all

T_window_s = [200,1000];
XSelect = 240;

measuredMag = readtable("../build/Linux/tests/magLog.txt");
accLog = readtable("../build/Linux/tests/accLog.txt");
measuredQ = readtable("../build/Linux/tests/quatLog.txt");
expectedQ = readtable("../build/Linux/tests/quatExpectedLog.txt");
estimatedQ = readtable("../build/Linux/tests/qEstLog.txt");
measuredv = readtable("../build/Linux/tests/vLog.txt");
PLog = readtable("../build/Linux/tests/PLog.txt");
Pest = readtable("../build/Linux/tests/PestLog.txt");
SLog = readtable("../build/Linux/tests/SLog.txt");
invSLog = readtable("../build/Linux/tests/invSLog.txt");
FLog = readtable("../build/Linux/tests/FLog.txt");
WLog = readtable("../build/Linux/tests/WLog.txt");
QLog = readtable("../build/Linux/tests/QLog.txt");
RLog = readtable("../build/Linux/tests/RLog.txt");
KLog = readtable("../build/Linux/tests/KLog.txt");

figure(1)
subplot(1,2,1)
plot(measuredQ.Variables)
hold on
xline(XSelect);


title("Predicted")
xlim(T_window_s)

subplot(1,2,2)
plot(expectedQ.Variables)
hold on
xline(XSelect);
title("Expected")
xlim(T_window_s);

rows = 4;
col = 3;

figure(2)
subplot(rows,col,1)
plot(measuredv.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("measurement-estimation divergence")

subplot(rows,col,2)
plot(measuredMag.Variables)
hold on
plot(accLog.Variables./10, '--')
xline(XSelect);
xlim(T_window_s)
title("measured mag-acc")
legend

subplot(rows,col,3)
plot(PLog.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("P")

subplot(rows,col,4)
plot(FLog.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("F")

subplot(rows,col,5)
plot(SLog.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("S")

subplot(rows,col,6)
plot(Pest.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("P_est")

subplot(rows,col,7)
plot(estimatedQ.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("q_est")

subplot(rows,col,8)
plot(QLog.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("QLOg")

subplot(rows,col,9)
plot(RLog.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("RLog")

subplot(rows,col,10)
plot(WLog.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("WLog")

subplot(rows,col,11)
plot(invSLog.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("invSLog")

subplot(rows,col,12)
plot(KLog.Variables)
hold on
xline(XSelect);
xlim(T_window_s)
title("KLog")