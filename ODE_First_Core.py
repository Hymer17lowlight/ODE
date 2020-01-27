from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint

"""
x' = ax + by + alpha
y' = cx + dy + betta
"""
def ode(y, t, a, b, c, d, al, bt): # "transmutation" the ode
    theta, omega = y
    if c != 0 and d != 0 and a / c == b / d:
        x0 = 0
        y0 = 0
    else:
        try:
            x0 = (d*al-b*bt)/(a*d-b*c)
            y0 = (al - a*x0)/b
        except ZeroDivisionError:
            try:
                y0 = (c*al-a*bt)/(b*c-a*d)
                x0 = (al-b*y0)/a
            except ZeroDivisionError:
                x0 = 0
                y0 = 0
    dydt = [c * omega + d * theta + y0, a * omega + b * theta + x0]

    return dydt


def get_type(a, b, c, d):  # getting type of critical point(classification), used in plt.title in DrawPhasePortret
    D = a * a - 2 * a * d + d * d + 4 * b * c

    if c != 0 and d != 0 and a / c == b / d:
        return 'parallel lines'
    if (a == c and a == 0) or (b == d and b == 0):
        return 'parallel lines'
    elif b == c and b == 0 and a * d > 0:
        return 'dicritical node'
    else:
        if D >= 0:
            sd = np.sqrt(D)
            if (a + d + sd) * (a + d - sd) < 0:
                return 'saddle'
            else:
                if a + d + sd >= 0:
                    return 'unstable node'
                else:
                    return 'stable node'
        else:
            if a + d == 0:
                return 'centre'
            else:
                if a + d < 0:
                    return 'stable focus'
                else:
                    return 'unstable focus'


def calcODE(args, y0, dy0, ts=10, nt=0.05):  # getting array of points
    y0 = [y0, dy0]
    t = np.arange(0, ts, nt)
    sol = odeint(ode, y0, t, args)
    return sol

fig = plt.figure()                              # creating a canvas
ax = plt.axes(xlim=(-15, 15), ylim=(-15, 15))
line, = ax.plot([], [], lw=0.5)

def drawPhasePortrait(name, args, deltaX=1, deltaDX=1, startX=-15, stopX=15, startDX=-15, stopDX=15, ts=2, nt=0.05):
    sp = []
    for y0 in range(startX, stopX, deltaX):
        for dy0 in range(startDX, stopDX, deltaDX):
            sol = calcODE(args, y0, dy0, ts, nt)
            plt.plot(sol[:, 1], sol[:, 0], 'black', lw=0.5)
            sp.append(sol)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()
    fname = 'fig{}.png'.format(name)
    plt.savefig(fname)
    plt.title(str(get_type(args[0], args[1], args[2], args[3])))
    plt.show()
    return sp


if __name__ == '__main__':
    spcs = ((1, 0, 2, -3, 0, 0), (1, 1, -1, -1, 0, 0), (0, 2, 1, 0, 0, 0), (1, -1, 0, 1, 0, 0), (0, 1, -1, 0, 0, 0),
            (2, -1, 0, 1, 0, 0), (0, 1, 1, 0, 0, 0), (1, 0, 0, 1, 0, 0), (1, 2, -0.5, 1, 0, 0))
    i = 1
    # start test with different points
    for spc in spcs:
        args = spc
        drawPhasePortrait(i, args)
        i += 1
    # working test
    name = '_temp'
    while True:
        a = float(input('a11='))
        b = float(input('a12='))
        c = float(input('a21='))
        d = float(input('a22='))
        alpha = float(input('alpha='))
        betta = float(input('betta='))
        print('\n')
        drawPhasePortrait(name, (a, b, c, d, alpha, betta))
