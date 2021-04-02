import numpy as np
from numpy import cos,sin
import matplotlib.pyplot as plt

def f(z,U0,zc,alpha,R,gamma):
	""" Potentiel d'un écoulement autour d'un cylindre """
        return U0*((z-zc)*np.exp(-1j*alpha) + R**2/(z - zc)*np.exp(1j*alpha)) - 1j*gamma/(2*np.pi)*np.log(z - zc)

def Joukovski(z,k):
	""" Transformation de Joukovski """
	return z + k**2/z

# Définition des données
U0 = 1
zc = 0.08 - 0.01j
alpha = np.pi/12
k = 1
R = np.sqrt((k-zc.real)**2 + zc.imag)
beta = np.angle((zc - k)/R)
gamma = -4*np.pi*R*U0*np.sin(alpha - beta)

t = np.arange(0,2*np.pi,0.01)
cercle = R*np.exp(1j*t) + zc

# Création du maillage
x = np.arange(-5,5,0.01)
y = np.arange(-5,5,0.01)
X, Y = np.meshgrid(x,y)

Z = X + 1j*Y
F = f(Z,U0,zc,alpha,R,gamma) # Potentiel de l'écoulement
TJ = Joukovski(F,k) # Transformé de l'écoulement 


# Trace les figures
plt.figure()
plt.contour(X,Y,F.imag,np.arange(-5,5,0.2)) # Tracé des lignes de -5 à 5 avec un pas de 0.2
plt.plot(cercle.real,cercle.imag,color='red',linewidth=3) # Tracé du 0 en rouge
plt.axis('equal')
plt.xlabel('Réel')
plt.ylabel('Imaginaire')
plt.title("Écoulement autour d'un cylindre avant transformation")

plt.figure()
plt.contour(X,Y,TJ.imag,np.arange(-5,5,0.2)) # Tracé des lignes de -5 à 5 avec un pas de 0.2
cerclej = Joukovski(cercle,k)
plt.plot(cerclej.real,cerclej.imag,color='red',linewidth=3)
plt.axis('equal')
plt.xlabel('Réel')
plt.ylabel('Imaginaire')
plt.title("Écoulement autour d'un cylindre après transformation")

plt.figure()
plt.plot(cercle.real, cercle.imag, label='Non transformé')
TJ = Joukovski(cercle,k)
plt.plot(TJ.real, TJ.imag, label='Transformé')
plt.legend()
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title("Transformé d'un cercle qui prouve que la transfo fonctionne")

plt.show()
