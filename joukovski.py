import numpy as np
from numpy import cos,sin
import matplotlib.pyplot as plt

def f(z,U0,zc,alpha,R,gamma):
	""" Potentiel d'un ecoulement autour d'un cylindre """
	return U0*((z-zc)*np.exp(-1j*alpha) + R**2/(z - zc)*np.exp(1j*alpha)) - 1j*gamma/(2*np.pi)*np.log(z - zc)

def Joukovski(z,k):
	""" Transformation de Joukovski """
	return z + k**2/z

# Definition des donnees
U0 = 1.0
xc = -0.08
yc = 0.01
zc = xc + yc*1j
alpha = np.pi/12
k = 1.0
R = np.sqrt((k-zc.real)**2 + zc.imag**2)
#beta = np.angle((zc - k)/R)
beta = np.arctan(yc/(k-xc))
gamma = -4*np.pi*R*U0*np.sin(alpha - beta)

t = np.arange(0,2*np.pi,0.01)
cercle = R*np.exp(1j*t) + zc

print ("xc = ", xc ," yc=", yc ,", R=", R, ", alpha=", alpha, ", beta=", beta)

# Creation du maillage
xx = np.arange(-5,5,0.01)
yy = np.arange(-5,5,0.01)
x,y = np.meshgrid(xx,yy)
#Z = []
#for i in x:
#	for j in y:
#		if abs(i + 1j*j) >= 1.:
#			Z.append(i + 1j*j)
#Z = np.array(Z)

for i in range(xx.shape[0]):
	for j in range(yy.shape[0]):
		if abs((x[i,j] -xc) + 1j*(y[i,j]-yc)) < R:
			x[i,j] = np.nan 
			y[i,j] = np.nan

z = x + 1j*y
F = f(z,U0,zc,alpha,R,gamma) # Potentiel de l'ecoulement (potentiel generateur et transforme identiques)
Z = Joukovski(z,k) #espace transforme

# Trace les figures
plt.figure()
plt.plot(cercle.real,cercle.imag,color='red',linewidth=3) # Trace du 0 en rouge
plt.contour(x,y,F.imag,np.arange(-5,5,0.2)) # Trace des lignes de -5 a 5 avec un pas de 0.2
plt.axis('equal')
plt.xlim([-3, 3])
plt.ylim([-3, 3])
plt.xlabel('x')
plt.ylabel('y')
plt.title("ecoulement autour d'un cylindre, domaine generateur")

plt.figure()
#plt.contour(ZJ.real,ZJ.imag,TJ.imag,np.arange(-5,5,0.2)) # Trace des lignes de -5 a 5 avec un pas de 0.2
#plt.contour(Z.real,Z.imag,TJ.imag,np.arange(-5,5,0.2)) # Trace des lignes de -5 a 5 avec un pas de 0.2
plt.contour(Z.real,Z.imag,F.imag,np.arange(-5,5,0.2)) # Trace des lignes de -5 a 5 avec un pas de 0.2
cerclej = Joukovski(cercle,k)
plt.plot(cerclej.real,cerclej.imag,color='red',linewidth=3)
plt.axis('equal')
plt.xlim([-3, 3])
plt.ylim([-3, 3])
plt.xlabel('X')
plt.ylabel('Y')
plt.title("ecoulement autour d'un profil, domaine transforme")

plt.figure()
plt.plot(cercle.real, cercle.imag, label='Non transforme')
TJ = Joukovski(cercle,k)
plt.plot(TJ.real, TJ.imag, label='Transforme')
plt.legend()
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title("Transforme d'un cercle qui prouve que la transfo fonctionne")

plt.show()
