import numpy as np
import matplotlib.pyplot as plt
import time
from pylab import *
from astropy import constants as cs
from scipy import integrate
UA=1.49597860*10**13

##Pregunta 1

#Guardammos los datos según columna en variables
longitud= np.loadtxt('sun_AM0.dat', usecols = [0])
flujo= np.loadtxt('sun_AM0.dat', usecols = [1])

#Hay que cambiar las unidades de medida de los datos
longitudb=longitud*10
flujob=flujo*100

#Graficamos
semilogx(longitudb,flujob)
xlabel('$Longitud\; de\; onda\; (\lambda) \;[\AA]$')
ylabel('$Flujo\; [ergs\cdot s^{-1} \cdot cm^{-2} \cdot \AA^{-1}]$')
title('Espectro solar')
grid(True)
savefig("EspectroSolar.png")
show()

##Pregunta 2

#Método de integración
def cuadrilatero(x,y):
    s=0
    for i in range(len(x)-1):
        s=s+(x[i+1]-x[i])*y[i]
    return s

#Cálculo de luminosidad
Luminosidad2=cuadrilatero(longitudb,flujob)*4*np.pi*UA**2

print "Luminosidad parte 2=" + str(Luminosidad2) + "Erg/s"


##Pregunta 3

T=5777
h=cs.h.cgs.value
c=cs.c.cgs.value
k=cs.k_B.cgs.value

def fun(x):
    #Función luego de hacer el cambio de variable sugerido
    f=((np.tan(x)**3)*(1+np.tan(x)**2))/(np.exp(np.tan(x))-1)
    return f

def medio(f,a,b):
    #Método sugerido en clases para calcular integrales. Utiliza el valor medio
    m=(a+b)/2
    s=(b-a)*f(m)
    return s

def linspace(a,b,step):
    #Función análoga a la de MatLab
    while a<=b:
            yield a
            a+=step

def simpson(f,a,b):
    #Función que permite calcular integrales siguiendo el método de Simpson
    n=100000
    s=f(a)+f(b)
    h=(b-a)/n
    for i in linspace(1,n-1,2):
        s=s+4*f(a+i*h)
    for j in linspace(2,n-2,2):
        s=s+2*f(a+j*h)
    return s*h/3

constante=((2*pi*h/c**2)*(k*T/h)**4)
t0=time.time()
Energia3=constante*(medio(fun,0.0,0.25)+simpson(fun,0.25,(np.pi/2)-0.25)+medio(fun,(np.pi/2)-0.25,(np.pi/2)))
tf=time.time()
#Se calcula el valor de la integral para todos sus valores, luego de aplicar el c.v.
#El primer y último término equivalen a los extremos donde antes no se podía calcular directamente

Energia3a=constante*((np.pi**4)/15) #Valor analítico usado para comparar


print "Energia parte 3=" + str(Energia3) + "Erg"
print "Se demora" + str(tf-t0) + "segundos en correr"
print "Energia calculada analíticamente="+str(Energia3a)

r=np.sqrt(Luminosidad2/(4*np.pi*Energia3))
print "Radio solar estimado="+str(r)+" metros"


##Pregunta 4

#Luminosidad calculada con el método "trapz"
t0=time.time()
Luminosidad4=np.trapz(flujob,x=longitudb)*4*np.pi*UA**2
tf=time.time()

print "Luminosidad parte 4=" + str(Luminosidad4) + "Erg/s"


#Energía total calculada con el método "quadz"
t0=time.time()
F=lambda x:(x**3)/(np.exp(x)-1) #Función a la que hay que aplicar la integral
integrando4=((2*pi*h/c**2)*(k*T/h)**4)
i4=integrate.quad(F,0,np.inf)
Energia4=integrando4*i4[0]
tf=time.time()

print "Energia parte 4=" + str(Energia4) + "Erg"
print "Se demora" + str(tf-t0) + "segundos en correr"
