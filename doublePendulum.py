from matplotlib import pyplot as plt
from math import sin, cos, pi, sqrt

# Eqns of motion from
# https://web.mit.edu/jorloff/www/chaosTalk/double-pendulum/double-pendulum-en.html

# Parameters
l1=1
l2=1
m1=1
m2=1
g=9.81
Tstandard = 2*pi*sqrt((l1+l2)/g)
print("Approximate time period = ", Tstandard)

def alpha1(th1, th2, om1, om2):
    return \
    (
      - g*(2*m1 + m2)*sin(th1)
      - m2*g*sin(th1 - 2*th2)
      - 2*sin(th1 - th2)*m2*(om2**2*l2 + om1**2*l1*cos(th1 - th2))
    )/(l1*(2*m1 + m2 - m2*cos(2*th1 - 2*th2)))


def alpha2(th1, th2, om1, om2):
    return 2*sin(th1 - th2)* \
    (
        om1**2*l1*(m1 + m2)
      + g*(m1 + m2)*cos(th1)
      + om2**2*l2*m2*cos(th1 - th2)
    )/(l2*(2*m1 + m2 - m2*cos(2*th1 - 2*th2)))

def energy(th1, th2, om1, om2):
    v1 = om1*l1
    v2r = om2*l2
    Y1 =-l1*cos(th1)
    Y2 = Y1 - l2*cos(th2)
   
    KE = 0.5*m1*v1**2\
       + 0.5*m2*\
       (
          (v2r*cos(th2) + v1*cos(th1))**2
        + (v2r*sin(th2) + v1*sin(th1))**2
       )
    GPE = g*(m1*(Y1+2) + m2*(Y2+2))
    return KE + GPE

# Setup
t   = 0
th1 = 3
th2 = 3
om1 = -2
om2 = -2
ts  = 0.02
nSteps = 500
offc = 0.5
nIters = 2

# Store the energy as a function of time
E = [0.] * (nSteps+1)
E[0] = energy(th1, th2, om1, om2)

plt.gca().set_aspect('equal')
# Time steps
for i in range (nSteps):
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)

    # Update previous time steps
    om1prev = om1
    om2prev = om2
    th1prev = th1
    th2prev = th2
   
    # Iterations each time step
    for it in range(nIters):
        # Update pendulum 1
        al1 = alpha1\
        (
            offc*th1 + (1-offc)*th1prev,
            offc*th2 + (1-offc)*th2prev,
            offc*om1 + (1-offc)*om1prev,
            offc*om2 + (1-offc)*om2prev
        )
        om1 = om1prev + ts*al1
        th1 = (th1prev + ts*(offc*om1 + (1-offc)*om1prev))%(2*pi)
        
        # Update pendulum 2
        al2 = alpha2\
        (
            offc*th1 + (1-offc)*th1prev,
            offc*th2 + (1-offc)*th2prev,
            offc*om1 + (1-offc)*om1prev,
            offc*om2 + (1-offc)*om2prev
        )
        om2 = om2prev + ts*al2
        th2 = (th2prev + ts*(offc*om2 + (1-offc)*om2prev))%(2*pi)
   
    # Plotting
    X1 = l1*sin(th1)
    X2 = X1 + l2*sin(th2)
    Y1 =-l1*cos(th1)
    Y2 = Y1 - l2*cos(th2)
   
    x=[0, X1, X2]
    y=[0, Y1, Y2]
   
    plt.plot(x,y, 'o-', color = "purple")
   
    E[i+1] = energy(th1, th2, om1, om2)
    print(i, ": ", E[i+1])
   
    t = t + ts
    plt.pause(0.01)
    plt.cla()

plt.figure(2)
plt.plot(E)
plt.show()
