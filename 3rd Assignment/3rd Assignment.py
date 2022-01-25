import numpy as np
import scipy.linalg as la
from sympy import symbols, Eq, solve, diff
from pylab import plot, axis, title, subplot, show

def firsteq(x, y, z):
    return 3*(x - (x**3/3) + y + z)

def seceq(x, y):
    return (-1)*((x - 0.7 + 0.8*y)/ 3)

def myFunction(z):
    x, y = symbols('x y')

    eq1 = Eq(firsteq(x, y, z), 0)
    eq2 = Eq(seceq(x, y), 0)

    solve((eq1, eq2), (x, y), dict=True)

    sol_dict = solve((eq1,eq2), (x, y))
    a1 = sol_dict[0][0]
    a2 = sol_dict[0][1]
    return a1, a2

def findequil():
    size = [-2.0, -1.0, -0.5, 0,];
    for i in size:
        print("Equilibrium for " + str(i));
        x, y = myFunction(i)
        print('x = ' + str(x))
        print('y = ' + str(y))

findequil();

def partialderiv():
    x, y, z= symbols('x y z')
    fx = diff(firsteq(x, y, z), x)
    fy = diff(firsteq(x, y, z), y)
    gx = diff(seceq(x, y), x)
    gy = diff(seceq(x, y), y)
    print("Partial Derivative of x on the First Equation: " + str(fx))
    print("Partial Derivative of y on the First Equation: " + str(fy))
    print("Partial Derivative of x on the Second Equation: " + str(gx))
    print("Partial Derivative of y on the Second Equation: " + str(gy))
    stability(fx, fy, gx, gy);

def stability(x1, x2, x3, x4):
    values = [-1.33409398, -0.40886583, 0.8048477, 1.99408035]
    z = [-2.0, -1.0, -0.5, 0,]
    for i in range(0, 4):
        fx1 = 3 - 3 * (values[i]**2)
        A = np.array([[fx1, -0.333],[3, -0.26667]])
        eigvals, eigvecs = la.eig(A)
        print("Eigenvalues when z equals " + str(z[i]) + " : " + str(eigvals))
        maxval = abs(eigvals[0])
        col = 0
        for i in range(0, len(eigvals)):
            if(abs(eigvals[i]) > maxval):
                maxval = eigvals[i];
                col = i;
        print("Dominant Eiganvalue: " + str(maxval))

partialderiv();


Dt = 0.01
def initialize():
    global u, uresult, v, vresult
    u = v = 0.01
    uresult = [u]
    vresult = [v]

def observe():
    global u, uresult, v, vresult
    uresult.append(u)
    vresult.append(v)

def update():
    global u, uresult, v, vresult
    nextx = 3*(u - (u**3/3) + v + z) * Dt
    nexty = (-1)*((u - 0.7 + 0.8*v)/ 3) * Dt
    u, v = nextx, nexty

def plot_phase_space():
    initialize()
    for t in range(10000):
        update()
        observe()
    plot(uresult, vresult)
    axis('image')
    axis([-0.1, 0.1, -0.1, 0.1])
    title('z = ' + str(z))

zs = [-2.0, -1.0, -0.5, 0]
for i in range(len(zs)):
    subplot(1, len(zs), i + 1)
    z = zs[i]
    plot_phase_space()

show()