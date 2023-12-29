import numpy as np
import matplotlib.pyplot as plt   
import math as ma
#from consti_b3 import*
import pandas as pd
import random
import time
from numpy.linalg import norm
t0 = time.time()


L=1 #Lenght in meter 
n=2001; #odd number of nodes   

arr=np.arange(0,L,L/n)
dl=L/(n-1)/L; 


sig0=250e6
E=2e11
h=0.01                    #depth of beam
K0=(2*sig0)/(E*h)

############################# Use this for auto compute the loading discretization based on the applied residual curvature width and mesh size ##############################
mulp=200 # width of each loading residual curvature = mulp x mesh size 
# f=int(mulp)
# Kl=(L/(n-1))*f


# ndiv=np.round((L/((Kl))),4)

# if ndiv==int(ndiv):
#     print("OK")
#     ndiv=int(ndiv)
   
# else:
#     print("not ok,Change f,I.M of L")

##################################################

ndiv=10 # number of loading discretization 

dl=L/(n-1)/L;


tmax=2.5*np.pi;                               

t=np.linspace(0, tmax, n)

jJ=np.zeros((n,n));   # the jacobian matrix where slope is obtained from curvature : JJ@phi=BB
bB=np.zeros((n,4)); # BB(:,i);i=1 to 4: corresponding to T1,T1+1,T2,T2+1
bB1=np.zeros((n,4)); # BB(:,i);i=1
Kp=np.zeros((n,2))

Kp1=np.zeros((n,2))
ph=np.zeros((n,1))
phy=np.zeros((n,1))
y=np.zeros(n)
u2=np.zeros(n)
u44=np.zeros(n)
R=np.zeros(n)

##################################### Generation of a shape #############################################

###############################################         

# problem-1 n=1001

# for i in range(0,n):
#   ph[i]=ph[i-1]
 
#   # if i<0.3*n:
#   #      ph[i]=0
 
    
#   if i>=.2*n and i<=.4*n:
#         ph[i]=ph[i-1]+0.0005               

#############################################2.Custom #########################################

#Problem -2 n=2001

for i in range(0,n):
  ph[i]=ph[i-1]
 
  # if i<0.3*n:
  #      ph[i]=0
 
    
  if i>=.3*n and i<=.5*n:
        ph[i]=ph[i-1]+0.0003
       

  if i>.5*n and i<=0.7*n:
                  ph[i]=ph[i-1]-0.001




test=np.zeros(n)
############################################# Curvature of the given shape calculated by FDM################

##################################### by using forward and backward difference######################

for i in range(1,int(n-1)):                    # n-1 is for the last node. Because in case of indexing in [], exact position is taken.                                     
    Kp[i,0]=(ph[i+1]-ph[i-1])/(2*dl)/K0   


Kp[n-1,0]=np.round((3*ph[n-1]-4*ph[n-2]+ph[n-3])/(2*dl),3)/K0  # round is used, otherwise a number close to zero is produced 


Kp[0,0]=(-3*ph[0]+4*ph[1]-ph[2])/(2*dl)/K0



###########################################################################################################
maxk=np.zeros((n,ndiv))
mink=np.zeros((n,ndiv))
Kr=np.zeros((n,ndiv))
Krn=np.zeros((n,ndiv))
Krr=np.zeros((n,ndiv))
Kmaxc=np.zeros((n,ndiv))
Mmaxc=np.zeros((n,ndiv))
Kmaxc2=np.zeros((n,ndiv))
Kr1=np.zeros(n)
ipk=np.zeros(ndiv)
urk=np.zeros(ndiv)
sumk=np.zeros(ndiv)
subk=np.zeros(ndiv)
znl=int(n/ndiv)
u=np.zeros(n);

for i in range(0,ndiv):
 
  ipk[0]=0
  ipk[i]=znl*(i)
 
  urk[i]=znl*(i+1)
 
 

###################################### Residual sum of squares #######################################################

#This will obtain the rectagular residual curvature for only on the applied length (d*) 
for j in range(0,ndiv):
  for i in range(int(ipk[j]),int(urk[j])+1):
    if i<=urk[j] and i>=ipk[j] : 
    
      sumk[j]=sumk[j]+(Kp[i,0])
   
for j in range(0,ndiv):
   
 
 
  for i in range(int(ipk[j]),int(urk[j])+1):   
     Kr[i,j]=sumk[j]/(urk[j]-ipk[j]+1)     #sum[{1 to n}(kr1i)^2]=n*1^2=n       

#############################################################################################################


# ################################################## From given curvature to shape generation #######################################################   
nn=1000000  
for j in range(0,ndiv):
  for i in range(int(ipk[j]),int(urk[j]+1)):   # +1 is added to match the M position with Kr
   Kp[i,1]=Kr[i,j]
   Kp1[i,1]=K0*Kp[i,1]
  ####################################################   revision needed
   Kmaxr=1000
   Kmaxr2=1000    
  #print(Kr[i,j]) 
   if Kr[i,j]>0:
    for h in range(0,nn):
       Kmaxr1=Kmaxr-((Kmaxr**(3)-Kr[i,j]*Kmaxr**(2)-1.5*Kmaxr**(2)+0.5)/(3*Kmaxr**(2)-2*Kr[i,j]*Kmaxr-3*Kmaxr))
       if abs(abs(Kmaxr1)-abs(Kmaxr))<=0.001:
            #print('itr=',i)
             Kmaxc[i,j]=Kmaxr1
             Mmaxc[i,j]=0.5*(3-Kmaxc[i,j]**(-2)) 
             
        
             break
       else:
           Kmaxr=Kmaxr1
   if Kr[i,j]<0:
    Krn[i,j]=Kr[i,j]*-1   
    for h in range(0,nn):
       Kmaxr1=Kmaxr2-((Kmaxr2**(3)-Krn[i,j]*Kmaxr2**(2)-1.5*Kmaxr2**(2)+0.5)/(3*Kmaxr2**(2)-2*Krn[i,j]*Kmaxr2-3*Kmaxr2))
     
       if abs(abs(Kmaxr1)-abs(Kmaxr2))<=0.001:
           
             Kmaxc[i,j]=Kmaxr1*-1
             Mmaxc[i,j]=-0.5*(3-Kmaxc[i,j]**(-2))
        
             break
       else:
           Kmaxr2=Kmaxr1       
  
##########################################################################################



################################### using simp. to obtain the phir from Kr ###################
from simp2 import*
u44=simp2(n,Kp1[:,1],dl,ph[0])  
   

       
from simp11 import*
from simp11 import*



[xx2,xy2]=simp11(n,ph,dl/1)

[xx5,xy5]=simp11(n,u44,dl/1)
##################################################################################################


######################################### Error Minimization ##################################################
# Calculating L2 norm =sqrt(e1^2+e2^2+e3^2)
Er=(norm(ph)-norm(u44))**2
print('Error between two shapes=',round(Er,5))



################################################### Plot #############################################


npp=10
xo=np.linspace(0,1,npp)
yo=np.zeros((npp,1))

plt.plot(xx2,xy2,'b--',linewidth=1)

plt.plot(xx5,xy5,'c-.',linewidth=1)
plt.xlabel("x/L",fontsize="12.5")
plt.ylabel("y/L",rotation=0,fontsize="12.5")


plt.legend(["Given shape","Obtained shape"], fontsize="7",)


plt.show()



plt.plot(-Kp[:,0],'k--',linewidth=1)
for i in range(0,ndiv):
  plt.plot(-Kr[:,i],'-',color="red",linewidth=1,)
  plt.plot(-Kmaxc[:,i],'--',color="darkviolet",linewidth=0.7,)
 
plt.ylabel(r'$\bar{K}_{max}$ or $\bar{K}_{rr}$',fontsize="12.5")
plt.xlabel('Nodal position',fontsize="11")
plt.legend([r'$\bar{K}_{rg}$',r'$\bar{K}_{rr}$',r'$\bar{K}_{max}$',], fontsize="10",loc=1)  


plt.show()
plt.plot(ph,'k-',linewidth=1)
plt.ylabel('Slope($\u03A6_{rg}$)',fontsize="11")
plt.xlabel('Nodal position',fontsize="11")

plt.show()
plt.plot(ph,'r--',linewidth=0.9)

plt.plot(u44,'b--',linewidth=1) 
plt.ylabel('Slope($\u03A6$)',fontsize="12.5")
plt.xlabel('Nodal position',fontsize="11")
plt.legend(["$\u03A6_{rg}$","$\u03A6_{rr}$"])
plt.show()
xn=np.zeros(n)
yn=np.zeros(n)


for i in range(0,ndiv):

  plt.plot(-Mmaxc[:,i],'-.',color="blue",linewidth=1)
plt.ylabel(r'$\bar{M}$',fontsize="12.5")
plt.xlabel('Nodal position',fontsize="11")
print('Width of applied curvatures (kl)=',((urk[1]-ipk[1])*dl*L*1000),'mm')
print('Size of mesh interval (dl)=',(dl*L*1000),'mm')
print('Factor (kl/dl)=',((urk[1]-ipk[1])*dl*L*1000)/(dl*L*1000))



t1=time.time()

t=t1-t0
print(t,'seconds')


r'$\bar{x}$'
