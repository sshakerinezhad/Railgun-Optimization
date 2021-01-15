# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 11:14:02 2020

@author: shaya
"""



import subprocess
import scipy as sp

flexcode = """the whole shabang

shabang

shabangerino"""




FlexCode = """TITLE 'RailGun'
COORDINATES cartesian2
VARIABLES
u(threshold=1e-6)	!Displacement in x
v(threshold=1e-6)!Displacement in y

vx(threshold=1e-6)
vy(threshold=1e-6)
rx(threshold=1e-6)
! SELECT
DEFINITIONS 
mag=0.1*globalmax(magnitude(x,y))/globalmax(magnitude(U,V))
!Stress Parameters
E=400e9
nu=0.28

!Railgun
r = 1.2 !rail radius
d = 3.9 !rail seperation
L = 40 !rail length
mu0 = 1.25663706e-6
I = %s
Farm= mu0*(I^2)/(2*Pi)*ABS(ln((d-r)/r))!Force produced as a product of the current (I), distance between rails (d), and rail radius (r)
Fwire = mu0/(2*Pi)*(I^2)/(d)*L

armatureMass = 30
!Kinematics
theta = 0*Pi/180

ax = if rx<L then Farm/armatureMass*cos(theta) else 0
ay = 9.81
vel = sqrt((vx*cos(theta))^2+(vy*sin(theta))^2)

magicShowDist=97.92
flightTime = magicShowDist/vel
dely = 1/2*ay*flightTime^2
delx =ARCTAN (val(v,L,r/2)/L)*180/Pi*magicShowDist
accuracy = 1-(sqrt(dely^2+delx^2)/6)
!Heat
DeltaTemp = 0

alpha=0
alphax = alpha
alphay = alpha

!Stiffness matrix components (isotropic material)
C11 =E*(1-nu)/(1+nu)/(1-2*nu)
C22 = C11
C33 = C11

C12 = E*nu/(1+nu)/(1-2*nu)
C13 = C12
C21 = C12
C23 = C12
C31 = C12
C32 = C12

!! Strain
!Axial Strain
ex=dx(u)
ey=dy(v)

!Engineering Shear Strain
G=E/(2*(1+nu))
gxy=(dx(v)+dy(u))


!mechanical strain
exm=ex - alphax*DeltaTemp
eym=ey - alphay*DeltaTemp

!!Stress via Hooke's law
!Axial Stress
sx = C11*exm+C12*eym
sy = C21*exm+C22*eym

!Shear stress
sxy=G*gxy
! INITIAL VALUES
EQUATIONS     
u:	dx(sx)+dy(sxy)= 0
v:	dx(sxy)+dy(sy)= 0

rx: dt(rx) = vx
vx: dt(vx) = ax
vy: dt(vy) = ay
! CONSTRAINTS 
BOUNDARIES    

  REGION 'Bottom Rail'    
 
    START(0,0)  !BL
    LINE TO (L,0) !BR					
	LINE TO (L,2*r) !TR	
		load(v) = -Fwire/(L*r)				
	LINE TO (0,2*r) !TL
		value(u)=0
		value(v)=0
	LINE TO CLOSE !BL


  REGION 'Top Rail'    
 
    START(0,d+2*r)  !BL
		load(v) = Fwire/(L*r)
    LINE TO (L,d+2*r) !BR					
	LINE TO (L,d+2*r+2*r) !TR					!Quadrilatteral
	LINE TO (0,d+2*r+2*r) !TL
		value(u)=0
		value(v)=0
	LINE TO CLOSE !BL
	
TIME 0 TO 10 halt(rx>L)
MONITORS        
PLOTS           
for t = 0 by endtime/10 to endtime
grid(x+u*mag, y+v*mag)
history (accuracy) at (0,0) PrintOnly Export Format '#t#b#1'
	file='test.txt'
SUMMARY
	report val(v,L,r/2)
	report val(rx,0,0)
	report Farm
	report Fwire
	report I
	report val(vel,0,0)
	report val(dely,0,0)
	report val(delx,0,0)
	report val(accuracy,0,0)

END"""

FlexFileName = "2CM4_DP_Scripting1.pde"

import matplotlib.pyplot as plt
minI = 410000
maxI = 415000
numSteps = 100
iRange = sp.arange(minI,maxI+1,(maxI-minI)/(numSteps-1))
print(iRange)
t = []
for k in range(numSteps):
    t.append(k)
    
acc = []    
for i in iRange:
    with open(FlexFileName, "w") as f:
        print(FlexCode%i, file=f)
    
    #completed = subprocess.run(["ls"], shell=True)
    completed =subprocess.run(["FlexPDE6","-S",FlexFileName],timeout=20)#,shell=True)
    print('returned: ', completed.returncode)
    
    with open("test.txt") as f:
        data=sp.loadtxt(f, skiprows=7)

        acc.append(data[-1:,1])
      

plt.plot(iRange,acc)
plt.title('Railgun Accuracy for Various Currents')
plt.legend(iRange)
plt.show()

found = False
i = 0
highest = 0
while found==False:
    if acc[i+1]>acc[i]:
        i = i+1
    if acc[i+1]<acc[i]:
        if acc[i-1]>acc[i]:
            i = i-1
        if acc[i-1]<acc[i]:
            highest = acc[i]
            found = True
            break
highestLin = 0
for i in range(len(acc)):
    if acc[i]>highestLin:
        highestLin = acc[i]
        bestILin = iRange[i]
print(acc)
bestI = iRange[i]
print(bestI)
print(bestILin)