import numpy as np


# 原始数据
t1i = 80    # 热水入口温度 75-85
t1o = 40    # 热水出口温度 40-50
p1 = 2.5    # 热水运行压力 1.0-3.0
t2i = 20    # 冷却水入口温度 20-25
t2o = 40    # 冷却水出口温度 40-45
p2 = 0.15    # 冷却水运行压力 0.1-0.2
G1 = 18000  # 热水流量 18000

# 物性参数计算
t2avg = (t2i + t2o) / 2
print("8.冷却水的定性温度", format(t2avg, ".2f"))
rho_2 = 995.6
print("9.冷却水的密度", format(rho_2, ".1f"))
cp2 = 4.180
print("10.冷却水的比热", format(cp2, ".3f"))
lambda_2 = 0.615
print("11.冷却水的导热系数", format(lambda_2, ".3f"))
mu_2 = 797.2 * 10 ** (-6)
print("12.冷却水的粘度", format(mu_2, ".7f"))
Pr2 = 5.415
print("13.冷却水的普朗特数", format(Pr2, ".3f"))
t1avg = (t1i + t1o) / 2
print("14.热水的定性温度", format(t1avg, ".2f"))
rho_1 = 983.2
print("15.热水密度", format(rho_1, ".1f"))
cp1 = 4.185
print("16.热水比热", format(cp1, ".3f"))
lambda_1 = 0.654
print("17.热水的导热系数", format(lambda_1, ".3f"))
mu_1 = 466.0 * 10 ** (-6)
print("18.热水的粘度", format(mu_1, ".7f"))
Pr1 = 2.981
print("19.热水的普朗特数", format(Pr1, ".3f"))

# 传热量及水流量
eta = 0.95
Q0 = G1 * cp1 * (t1i - t1o) * eta * 1000 / 3600
print("21.设计传热量", format(Q0, ".1f"))
G2 = 3600 * Q0 / (1000 * cp2 * (t2o - t2i))
print("22.冷却水流量", format(G2, ".1f"))

# 有效平均温差
delta_t_big = t1i - t1o
delta_t_small = t2o - t2i
delta_tN = (delta_t_big - delta_t_small) / np.log(delta_t_big / delta_t_small)
print("23.逆流平均温差", format(delta_tN, ".2f"))
P = (t2o - t2i) / (t1i - t2i)
R = (t1i - t1o) / (t2o - t2i)
print(P, R)
phi = 0.8
print("24.温差校正系数", format(phi, ".1f"))
delta_tm = phi * delta_tN
print("25.有效平均温差", format(delta_tm, ".2f"))

# 管程传热系数
K0 = 700   # 500-1200
print("26.试选传热系数", format(K0, ".1f"))
F0 = Q0 / (K0 * delta_tm)
print("27.初选传热面积", format(F0, ".2f"))
do = 0.025  # 25x2.5 19x2 16x1.5
print("28.管子外径", format(do, ".3f"))
di = do - 2 * 0.0025
print("28.管子内径", format(di, ".2f"))
length = 6   # 2 2.5 3 4.5 6
print("29.管子长度", format(length, ".2f"))
Nt = np.ceil(F0 / (np.pi * do * length))
Nt = Nt + 9
print("30.总管子数", format(Nt, ".0f"))
a2 = (Nt / 2) * (np.pi / 4) * (di ** 2)
print("31.管程流通截面", format(a2, ".3f"))
w2 = G2 / (rho_2 * a2 * 3600)
print("32.管程流速", format(w2, ".2f"))
W2 = rho_2 * w2
print(W2)
Re2 = W2 * di / mu_2
print("33.管程雷诺数", format(Re2, ".2f"))
Nu2 = 0.023 * (Re2 ** 0.8) * (Pr2 ** 0.4)
print(Nu2)
h2 = Nu2 * lambda_2 / di
print("34.管程换热系数", format(h2, ".2f"))

# 初选结构
S = 0.032  # 25-32 19-25 16-22
print("36.管间距", format(S, ".3f"))
Nc = np.ceil(1.1 * (Nt ** 0.5))
print("37.管束中心处一排管数", format(Nc, ".0f"))
e = 2 * do
print("38.管束外沿与壳体间距", format(e, ".3f"))
Ds = S * (Nc - 1) + 4 * do
print("39.壳体内径", format(Ds, ".3f"))
print("40.长径比", format(length / Ds, ".3f"))
h = 0.2 * Ds
print("41.弓形折流板弓高", format(h, ".3f"))
B = Ds / 3
print("42.折流板间距", format(B, ".3f"))
nB = np.ceil((length / B) - 1)
print("43.折流板数", format(nB, ".0f"))

# 壳程换热
a1 = B * Ds * (1 - do / S)
print("44.壳程流通截面", format(a1, ".3f"))
w1 = G1 / 3600 / (rho_1 * a1)
print("45.壳程流速", format(w1, ".2f"))
W1 = rho_1 * w1
print("46.壳程质量流速", format(W1, ".2f"))
de = (Ds ** 2 - Nt * do ** 2) / (Ds + Nt * do)
print("47.壳程当量直径", format(de, ".3f"))
Re1 = W1 * de / mu_1
print("48.壳程雷诺数", format(Re1, ".2f"))
print("49.管间距比值", format((2 / (3 ** 0.5)), ".2f"))
# tw = 50
Prw = 3.551
Nu1 = 0.35 * ((2 / (3 ** 0.5)) ** 0.2) * (Re1 ** 0.6) * (Pr1 ** 0.36) * ((Pr1 / Prw) ** 0.25)
print("50.壳程努塞尔数", format(Nu1, ".2f"))
epsilon_n = 0.99
print("51.管排修正系数", format(epsilon_n, ".2f"))
h1 = epsilon_n * Nu1 * lambda_1 / de
print("52.壳程换热系数", format(h1, ".2f"))

# 传热系数
r2 = 17.2 * 10 ** (-5)
print("53.冷却水侧污垢热阻", format(r2, ".6f"))
r1 = 34.4 * 10 ** (-5)
print("54.热水侧污垢热阻", format(r1, ".6f"))
r_sigma = 1 / h1 + r1 + r2 * (do / di) + (1 / h2) * (do / di)
print("56.总传热热阻", format(r_sigma, ".6f"))
Kj = 1 / r_sigma
print("57.传热系数", format(Kj, ".2f"))
print("58.传热系数比值大小", format(Kj/K0, ".2f"))

# 管壁温度
q1 = Q0 / (Nt * np.pi * do * length)
print("59.管外壁热流密度", format(q1, ".2f"))
tw1 = t1avg - q1 * (1 / h1 + r1)
print("60.管外壁温度", format(tw1, ".2f"))

# 管程压降
tw = 50     # 假定管壁温度
nt = 2      # 管程数
mu_w1 = 546.5 * 10 ** (-6)
mu_w2 = 546.5 * 10 ** (-6)
print("1.壁温下水的粘度", format(mu_w2, ".7f"))
phi_2 = (mu_2 / mu_w2) ** 0.14
print("2.管程粘度修正系数", format(phi_2, ".2f"))
xi_i = 0.035
print("3.管程摩擦系数", format(xi_i, ".3f"))
delta_pi = (W2 ** 2 / (2 * rho_2)) * (length * nt / di) * (xi_i / phi_2)
print("4.管子沿程压降", format(delta_pi, ".2f"))
delta_pr = (W2 ** 2 / (2 * rho_2)) * 4 * nt
print("5.回弯压降", format(delta_pr, ".2f"))
WN2 = rho_2 * (3300 / rho_2) ** 0.5
print("6.进、出口管处质量流速", format(WN2, ".2f"))
delta_pN2 = ((WN2 ** 2) / (2 * rho_2)) * 1.5
print("7.进、出口管处压降", format(delta_pN2, ".2f"))
phi_d2 = 1.4
print("8.管程结垢校正系数", format(phi_d2, ".1f"))
delta_p2 = (delta_pi + delta_pr) * phi_d2 + delta_pN2
print("9.管程压降", format(delta_p2, ".2f"))

# 壳程压降
print("10.当量直径", format(de, ".3f"))
print("11.雷诺数", format(Re1, ".2f"))
xi_0 = 0.42     # 折流板缺圆高度20%
print("12.壳程摩擦系数", format(xi_0, ".2f"))
phi_1 = (mu_1 / mu_w1) ** 0.14
print(phi_1)
delta_p0 = ((W1 ** 2) / (2 * rho_1)) * (Ds * (nB + 1) / de) * (xi_0 / phi_1)
print("13.管束压降", format(delta_p0, ".2f"))
WN1 = rho_1 * (2200 / rho_1) ** 0.5
print("14.管嘴处质量流速", format(WN1, ".2f"))
delta_pN1 = ((WN1 ** 2) / (2 * rho_1)) * 1.5
print("15.进、出口管压降", format(delta_pN1, ".2f"))
epsilon_ip = 7.5    # 5-10
print("16.导流板阻力系数", format(epsilon_ip, ".1f"))
delta_pip = ((WN1 ** 2) / (2 * rho_1)) * epsilon_ip
print("17.导流板压降", format(delta_pip, ".2f"))
phi_do = 1.82
print("18.壳程结垢修正系数", format(phi_do, ".2f"))
delta_p1 = delta_p0 * phi_do + delta_pip + delta_pN1
print("19.壳程压降", format(delta_p1, ".2f"))
delta_p1_allow = 180000
delta_p2_allow = 75000

# 校核传热系数
if 500 <= K0 <= 1200:
    if 4 <= (length / Ds) <= 25:
        if 1.1 <= (Kj / K0) <= 1.2:
            print("传热系数校核通过")

# 校核压降
if delta_p1 < delta_p1_allow:
    if delta_p2 < delta_p2_allow:
        print("压强校核通过")
