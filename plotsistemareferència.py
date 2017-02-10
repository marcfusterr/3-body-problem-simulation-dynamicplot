import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D

#AUTHORS: DANI GONCALVES, MARTÍ BERENGUER, MARC FUSTER
 
###################################################################
##--------------------Declaració de variables--------------------##
###################################################################
 
m1=1.989E30
m2=5.972E24
m3=7.347673E22
m4=1E40
G=6.67E-11

h=50 #pas d'integració temporal
n=10*10**4
divisor=5*10**2
distanciats=1.49597E11    #tot en SI
distanciatl=3.844E8

#sol
xo1=0.0000
yo1=0.0000
zo1=0.0000
vox1=0.000
voy1=0.000
voz1=0.000
#terra
xo2=distanciats
yo2=0.0000
zo2=0.0000
vox2=0.000
voy2=30.3e3
voz2=0.000
#lluna
xo3=distanciats-distanciatl
yo3=0.0000
zo3=0.0000
vox3=0.000
voy3=voy2-1022    
voz3=0.000
#asteroide
xo4=1.001E11
yo4=0.026E11
zo4=0
vox4=3E4
voy4=2.975E4
voz4=0

P1=(xo1,yo1,zo1)
P2=(xo2,yo2,zo2)
P3=(xo3,yo3,zo3)
P4=(xo4,yo4,zo4)
 
P1=np.array(P1)
P2=np.array(P2)
P3=np.array(P3)
P4=np.array(P4)
 
V1=(vox1,voy1,voz1)
V2=(vox2,voy2,voz2)
V3=(vox3,voy3,voz3)
V4=(vox4,voy4,voz4)
 
V1=np.array(V1)
V2=np.array(V2)
V3=np.array(V3)
V4=np.array(V4)

#diccionari de les posicions
posicionsit=dict()
 
 
def mod(x): #Definim el mòdul d'un vector com mod(x)
  return (sum(x*x))**0.5


#A les equacions he definit les i com les masses m_i, les Pi son les posicions de
#les particules i en un temps donat i Vi les velocitats
def Energia(i,Pi,Vi,j,Pj,k,Pk):
    return i*(0.5*mod(Vi)**2-G*(j/mod(Pi-Pj)+k/mod(Pi-Pk)))
 
def Acceleracio(Pi,j,Pj,k,Pk,Pl,l):
    return -G*(j*(Pi-Pj)/(mod(Pi-Pj)**3)+k*(Pi-Pk)/(mod(Pi-Pk)**3))+l*(Pi-Pl)/(mod(Pi-Pl)**3)
     

solucio=np.zeros((n,9))    #ENS crea una matriu on guardarem les solucions de cada iteracio
   
 
###################################################################
##--------------------Programa 3 cossos a ex.--------------------##
###################################################################
 
Pi1=np.zeros(3)
Pi2=np.zeros(3)
Pi3=np.zeros(3)
Pi4=np.zeros(3)
Vi1=np.zeros(3)
Vi2=np.zeros(3)
Vi3=np.zeros(3)
Vi4=np.zeros(3) 
#El subindes i fa referencia a aquells elements que es corresponen al temps n+1 
#a la iteracio n. Es a dir, son les solucions de l'equacio diferencial, els nous
#punts de la trajectoria.
 

for i in range(n):
    
    #Part 1.1:  Noves posicions despres de dt - Mètode de Euler

    Vi1=V1+Acceleracio(P1,m2,P2,m3,P3,m4,P4)*h
    Pi1=P1+Vi1*h
     
    Vi2=V2+Acceleracio(P2,m1,P1,m3,P3,m4,P4)*h
    Pi2=P2+Vi2*h
 
    Vi3=V3+Acceleracio(P3,m1,P1,m2,P2,m4,P4)*h
    Pi3=P3+Vi3*h
    
    Vi4=V4+Acceleracio(P4,m1,P1,m2,P2,m3,P3)*h
    Pi4=P4+Vi4*h
         
    "MOLT IMPORTANT; NOMES GUARDAREM 1 DE CADA 20 PUNTS, JA QUE NO HU HA DIFERÈNCIA"
    "EN ELS GRÀFICS I PERMET MANIPULARLOS MOLT MES FACILMENT"
    
    if i%divisor ==0:
    #per canviar de multiples canviar el divisor
     #Part 1.2.2: guardar en matriu txt els punts de la trajectoria 
     
            vectorposicionsit=np.zeros(13)
            for j in range(3):
                vectorposicionsit[j]=Pi1[j]
                vectorposicionsit[j+3]=Pi2[j]
                vectorposicionsit[j+6]=Pi3[j]
                vectorposicionsit[j+9]=Pi4[j]
                vectorposicionsit[12]=i*h
            posicionsit[i/divisor]=vectorposicionsit
            
            """ DESBLOQUEJAR PER PLOT REAL TIME MENTRES CALCULA
            plt.pause(0.001)
            plt.clf()
            plt.grid()
            plt.xlim([1.3e11,1.5e11])
            plt.ylim([0,1.3e11])
            plt.scatter(Pi1[0],Pi1[1],color='red')
            plt.scatter(Pi2[0],Pi2[1],color='blue')
            plt.scatter(Pi3[0],Pi3[1],color='grey')
            plt.scatter(Pi4[0],Pi4[1],color='green')
            plt.show(1)
            """
      
    #Part 1.3: recuperació dels valors inicials
    P1=Pi1
    P2=Pi2
    P3=Pi3
    P4=Pi4
    V1=Vi1
    V2=Vi2
    V3=Vi3
    V4=Vi4
     
"END DEL FOR"


#aixoensparaelscalculsfinsque apretem una nova tecla si ho desbloquejam
"""
input('press any key to continue')
"""

#CREACIO DEL PLOT UN COP JA GUARDATS TOTS ELS VALORS. 

plt.figure(2)
for i in range(len(posicionsit)):
    k=i+int(round(len(posicionsit))/4)
    for j in range(3):
        
        plt.pause(0.001)
        plt.clf()
        plt.grid()
        plt.xlim([posicionsit[k][3]-0.5E9,posicionsit[k][3]+0.5E9])             #([-2e11,2e11])
        plt.ylim([posicionsit[k][4]-0.5E9,posicionsit[k][4]+0.5E9])             #([-2e11,2e11])
        plt.scatter(posicionsit[k][0],posicionsit[k][1],color='red')
        plt.scatter(posicionsit[k][3],posicionsit[k][4],color='blue',label='Terra')
        plt.scatter(posicionsit[k][6],posicionsit[k][7],color='grey',label='lluna')
        plt.scatter(posicionsit[k][9],posicionsit[k][10],color='grey')
        plt.legend()
        plt.show()
       

     
