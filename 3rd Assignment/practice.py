import numpy as np
from scipy.optimize import fsolve

def myFunction(z):
   x = z[0]
   y = z[1]

   F = np.empty((2))
   F[0] = 3*(x - (x**3/3) + y + z)
   F[1] = (-1)*((x - 0.7 + 0.8*y)/ 3)
   return F

zGuess = np.array([1,1])
z = fsolve(myFunction,zGuess)
print(z)


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
    nextx = v + 3*(u - (u**3/3) + v + z) * Dt
    nexty = u + (-1)*((u - 0.7 + 0.8*v)/ 3) * Dt
    u, v = nextx, nexty

def plot_phase_space():
    initialize()
    for t in range(10000):
        update()
        observe()
    plot(uresult, vresult)
    axis('image')
    axis([-3, 3, -3, 3])
    title('z = ' + str(z))

zs = [-2.0, -1.0, -0.5, 0]
for i in range(len(zs)):
    subplot(1, len(zs), i + 1)
    z = zs[i]
    plot_phase_space()

show()