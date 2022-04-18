import numpy as np
import math as m
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

GM = 3.986004418 * 10**14
c0 = 2.202095 * 10**10
omega = 7.2921
y = 0.51
a = 6378245
f = 298.3
e = (2*f - 1) / (f - 1)**2


def p(h):  # плотность атмосферы в кг\м^3 в зависимост от высоты в метрах
    return 1.41*m.exp(- h / 6910)


def g(x, v):  # вектор полного ускорениятела в ГСПК
    r = m.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    g0 = -GM*(x-3*c0*x/(r**2)-6*c0*np.array([0, 0, x[2]])/(r**2)+15*c0*x*(x[2]**2)/(r**4))/(r**3)
    gCA = (omega**2)*np.array([x[0], x[1], 0])
    gCOR = 2*omega*np.array([v[1], -v[0], 0])
    v_abs = m.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    s = x[2] / r
    h = r - a/m.sqrt(1 + (e**2)*(s**2))
    gAER = -y*p(h)*v_abs*v/2
    return g0 + gCA + gCOR + gAER


def rungekutta4(f, x0, v0, t, args=()):  # метод Рунге-Кутты 4-го порядка
    n = len(t)
    x = np.zeros((n, len(x0)))
    v = np.zeros((n, len(v0)))
    x[0] = x0
    v[0] = v0
    for i in range(n - 1):
        dt = t[i+1] - t[i]
        k1 = f(x[i], v[i])
        q1 = v[i]
        k2 = f(x[i] + q1 * dt / 2., v[i] + k1 * dt / 2.)
        q2 = q1 + k1 * dt / 2.
        k3 = f(x[i] + q2 * dt / 2., v[i] + k2 * dt / 2.)
        q3 = q1 + k2 * dt / 2.
        k4 = f(x[i] + q3 * dt, v[i] + k3 * dt)
        q4 = q1 + k3 * dt
        v[i+1] = v[i] + (dt / 6.) * (k1 + 2*k2 + 2*k3 + k4)
        x[i+1] = x[i] + (dt / 6.) * (q1 + 2*q2 + 2*q3 + q4)
    return x, v


def x_to_bic(x):
    res = []
    fi = 30
    res.append(m.sqrt(x[0]**2 + x[1]**2) * m.sin(fi*m.pi/180) + x[2] * m.cos(fi*m.pi/180))
    res.append(-m.sqrt(x[0]**2 + x[1]**2) * m.cos(fi*m.pi/180) + x[2] * m.sin(fi*m.pi/180))
    res.append(x[1] / x[0])
    return res


def x_to_grinv(x, v, fi, lam, h):
    fi_r = fi * m.pi / 180
    lam_r = lam * m.pi / 180
    r = a / m.sqrt(1 + (1/((1 - 1/f)**2) - 1)*m.sin(fi)**2)
    fi_ = m.atan(m.tan(fi_r)*(1 - 1/f)**2)
    Amg = np.array([[-m.cos(lam_r)*m.sin(fi_r), -m.sin(lam_r)*m.sin(fi_r), m.cos(fi_r)],
                    [m.cos(lam_r)*m.cos(fi_r), m.sin(lam_r)*m.cos(fi_r), m.sin(fi_r)],
                    [-m.sin(lam_r), m.cos(lam_r), 0]])

    x0 = np.array([r*m.cos(fi_)*m.cos(lam_r) + h*m.cos(fi_r)*m.cos(lam_r),
                   r*m.cos(fi_)*m.sin(lam_r) + h*m.cos(fi_r)*m.sin(lam_r),
                   r*m.sin(fi_) + h*m.sin(fi_r)]).T

    x_res = x0 + Amg.T @ x
    v_res = Amg.T @ v
    return x_res, v_res


def noise(x_bic):
    x_bic[0] += random.uniform(-30, 30)
    x_bic[1] += random.uniform(-30, 30)
    x_bic[2] += random.uniform(-0.3, 0.3)
    return x_bic


def g_h_filter(data, x0, dx, g, h, dt=1.):
    x_est = x0
    results = []
    for z in data:
        x_pred = x_est + dx*dt
        dx = dx

        residual = z - x_pred
        dx = dx + h * (residual) / dt
        x_est = x_pred + g * residual
        results.append(x_est)
    return np.array(results)


def main():
    x = np.array([100, 500, 400])
    v = np.array([3, 2, -1])

    fi = 25.0456  # градусы
    lam = 48.553
    h = 1000  # метры

    x0, v0 = x_to_grinv(x, v, fi, lam, h)
    print(x0, v0)

    t = np.linspace(0, 10, 11)
    x1, v1 = rungekutta4(g, x0, v0, t)
    print(x1)

    xx = []
    yy = []
    zz = []
    for i in range(len(x1)):
        xx.append(x1[i][0]/1000)
        yy.append(x1[i][1]/1000)
        zz.append(x1[i][2]/1000)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(xx, yy, zz)
    ax.set_xlabel('x(км)')
    ax.set_ylabel('у(км)')
    ax.set_zlabel('z(км)')
    plt.savefig('graph.png')

    #x_bic = []
    #for i in range(len(x1)):
    #    x_bic.append(x_to_bic(x1[i]))
    #    x_bic[i] = noise(x_bic[i])

    #result = g_h_filter(data=x_bic, x0=x, dx=, g=.2, h=0.05)


if __name__ == '__main__':
    main()
