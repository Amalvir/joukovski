import numpy as np
import matplotlib.pyplot as plt

def f(z,U0,zc,alpha,R,gamma):
        return U0*((z-zc)*np.exp(-1j*alpha) + R**2/(z - zc)*np.exp(1j*alpha)) - 1j*gamma/(2*np.pi)*np.log(z - zc)

x = np.arange(-5,5,0.01)
y = np.arange(-5,5,0.01)
X, Y = np.meshgrid(x,y)
Z = X + 1j*Y
F = f(Z,1,0,0,1,0)
plt.contour(X,Y,F.imag,np.arange(-5,5,0.2)) # Tracé des lignes de -5 à 5 avec un pas de 0.2
plt.contour(X,Y,F.imag,[0],colors='red',linewidths=3) # Tracé du 0 en rouge
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
