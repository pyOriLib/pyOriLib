# -*- coding: utf-8 -*-
"""
Created on Wed May 26 08:56:21 2021

@author: lecomte5
"""
import matplotlib.pyplot as plt
import numpy as np
degree = np.pi/180

def norm(x,*args):
    """
 norm - Normalization of a vector 
 Usage : norm(v)
         norm(x,y,z)
         
    """
    if len(args)==0:
        if isinstance(x,vector3d):
            x,y,z=x.x,x.y,x.z
        else:
            x,y,z=float(x[0]),float(x[1]),float(x[2])  
    elif len(args)==2:
        x,y,z=x,args[0],args[1]
    else:
        x,y,z=0,0,0
    return np.sqrt(x*x+y*y+z*z)

def dot(x,y):
    if isinstance(x,vector3d):
        a1,a2,a3=x.x,x.y,x.z
    else:
        a1,a2,a3=x[0],x[1],x[2]
    if isinstance(y,vector3d):
        b1,b2,b3=y.x,y.y,y.z
    else:
        b1,b2,b3=y[0],y[1],y[2]
    return a1*b1+a2*b2+a3*b3
    
def vect_vect_theta( v1, v2 ):
    # calcul de l'angle entre deux vecteurs
    dummy=norm(v1)*norm(v2);
    if (dummy !=0):
        return np.arccos(dot(v1,v2)/dummy)
    
class vector3d:
    """
    
 The class vector3d describes three dimensional vectors, given by
 their coordinates x, y, z and allows to calculate with them as
 comfortable as with real numbers.

 Syntax
   v = vector3d(x,y,z)
   v = vector3d(x,y,z,antipodal)
   v = vector3d.byPolar(theta,rho)

 Input
  x,y,z - cart. coordinates

 Output
  v - @vector3d

 Flags
  antipodal - <VectorsAxes.html consider vector as an axis>

 Class Properties
  x, y, z      - cart. coordinates
  isNormalized - whether the vector is a direction
  antipodal    - <VectorsAxes.html whether the vector is an axis>

 Dependent Class Properties
  theta      - polar angle in radiant
  rho        - azimuthal angle in radiant
  resolution - mean distance between the points on the sphere
  xyz        - cart. coordinates as matrix
    """
    
    def __init__(self,x=0,y=0,z=0,antipodal = False,isNormalized=False):
        self.x=x
        self.y=y
        self.z=z
        self.antipodal=antipodal
        self.isNormalized=isNormalized
        
    def as_vector3d(self):
        return (self.x, self.y, self.z)     
    
    # def __repr__(self):
    #     return '[%.3f %.3f %.3f]' %(self.x,self.y,self.z)
    def __str__(self):
        return '[{}, {}, {}]'.format(*self.as_vector3d())    
    def __neg__(self):
        return vector3d(-self.x,-self.y,-self.z)
    def __mul__(self,operand):
        return np.array([[self.x*operand.x, self.x*operand.y,self.x*operand.z],
                         [self.y*operand.x, self.y*operand.y,self.y*operand.z],
                         [self.z*operand.x, self.z*operand.y,self.z*operand.z],
                         ])
    def byPolar(self,polarAngle,azimuthAngle):
        x = np.sin(polarAngle) * np.cos(azimuthAngle)
        y = np.sin(polarAngle) * np.sin(azimuthAngle)
        z = np.cos(polarAngle) 
        return vector3d(x, y, z)
    def rand(self):
        x=np.random.random()
        y=np.random.random()
        z=np.random.random()
        return vector3d(x,y,z)
    def __add__(self, other):
        return vector3d(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return vector3d(self.x - other.x, self.y - other.y, self.z - other.z)


	# Scalar Multiplication

    def __mul__(self, number):
        return vector3d(self.x * number, self.y * number, self.z * number)

    def __rmul__(self, number):
        return self.__mul__(number)

 	# Cross product
 	# cross = a ** b
    def __pow__(self, operand):
        return vector3d(self.y*operand.z - self.z*operand.y, 
 			                self.z*operand.x - self.x*operand.z, 
 			                self.x*operand.y - self.y*operand.x)

 	# Dot Project
 	# dp = a & b
    def __and__(self, operand):
        return (self.x * operand.x) + \
 		       (self.y * operand.y) + \
 		       (self.z * operand.z)
                
    # def __norm__(self):
    #     return norm(self.x,self.y,self.z)
    # Operations
    def Normalize(self):
        if self.isNormalized:
            return
        else:    
            dummy=norm(self)
            
            self.isNormalized=True
        return vector3d(self.x/dummy,self.y/dummy,self.z/dummy) 
def axisangle2mat(r, theta):
    if isinstance(r,vector3d):
        r[0],r[1],r[2]=r.x,r.y,r.z
    nr=norm(r)
    eps=1e-6
    for i in range(3):
        r[i]=r[i]/nr
    g= np.zeros((3,3))
    C = np.cos(theta);
    S = np.sin(theta);
 
    g[0,0] = r[0] * r[0] * (1- C) + C;
    g[0,1] = r[0] * r[1] * (1- C) + r[2] * S;
    g[0,2] = r[0] * r[2] * (1- C) - r[1] * S;
 
    g[1,0] = r[1] * r[0] * (1- C) - r[2] * S;
    g[1,1] = r[1] * r[1] * (1- C) + C;
    g[1,2] = r[1] * r[2] * (1- C) + r[0] * S;
 
    g[2,0] = r[2] * r[0] * (1- C) + r[1] * S;
    g[2,1] = r[2] * r[1] * (1- C) - r[0] * S;
    g[2,2] = r[2] * r[2] * (1- C) + C;
    for i in range(3):
        for j in range(3):
            if (abs(g[i,j])<eps):
                g[i,j]=0
    return np.array(g,float)

def eul2mat( e ):
    g= np.zeros((3,3))
    e[0]=e[0]*degree
    e[1]=e[1]*degree
    e[2]=e[2]*degree
    g[0,0] = np.cos (e[0]) * np.cos (e[2]) - np.sin (e[0]) * np.sin (e[2]) * np.cos (e[1]);
    g[0,1] = np.sin (e[0]) * np.cos (e[2]) + np.cos (e[0]) * np.sin (e[2]) * np.cos (e[1]);
    g[0,2] = np.sin (e[2]) * np.sin (e[1]);
 
    g[1,0] = -np.cos (e[0]) * np.sin (e[2]) - np.sin (e[0]) * np.cos (e[2]) * np.cos (e[1]);
    g[1,1] = -np.sin (e[0]) * np.sin (e[2]) + np.cos (e[0]) * np.cos (e[2]) *np. cos (e[1]);
    g[1,2] = np.cos (e[2]) * np.sin (e[1]);
 
    g[2,0] = np.sin (e[0]) * np.sin (e[1]);
    g[2,1] = -np.cos (e[0]) * np.sin (e[1]);
    g[2,2] = np.cos (e[1]);
    return np.array(g,float)

def mat2eul( g ):
    eps = 1E-6;
    e=np.zeros(3)
    e[1] = np.arccos (g[2][2]);
    if (abs (e[1]) > eps) and (abs (e[1] - np.pi) > eps):
        e[0] = np.arctan2 (g[2][0], -g[2][1]);
        e[2] = np.arctan2 (g[0][2],  g[1][2]);
    else:
        e[0] = np.arctan2 (g[0][1], g[0][0]);
        e[2] = 0;
    return e

def miller2mat ( hkl,uvw ):
    M =norm(hkl);
    N =norm(uvw);
    g= np.zeros((3,3))
    g[0,0] =  uvw[0] / N;
    g[1,0] =  uvw[1] / N;
    g[2,0] =  uvw[2] / N;
 
    g[0,1] =( hkl[1] *  uvw[2] -  hkl[2] *  uvw[1]) / (M * N);
    g[1,1] =( hkl[2] *  uvw[0] -  hkl[0] *  uvw[2]) / (M * N);
    g[2,1] =( hkl[0] *  uvw[1] -  hkl[1] *  uvw[0]) / (M * N);
 
    g[0,2] =  hkl[0] / M;
    g[1,2] =  hkl[1] / M;
    g[2,2] =  hkl[2] / M;
    
    return np.array(g,float)

def axis2indices( Axis, MaxIndices ):
    maX=(2*MaxIndices+1)*(2*MaxIndices+1)*(2*MaxIndices+1)-1;    
    q=0;AxisIndices=np.zeros((3,maX));AngleOff=np.zeros(maX);
    
    for i in range(-MaxIndices,MaxIndices+1):
        for j in range(-MaxIndices,MaxIndices+1):
            for k in range(-MaxIndices,MaxIndices+1):
                if not((i==0)and(j==0)and(k==0)):                
                    AxisIndices[:,q]=[i, j, k];                    
                    AngleOff[q]=vect_vect_theta(Axis,AxisIndices[:,q]);                    
                    q=q+1;
    #BubbleSort
    n = len(AngleOff);
    
    for i in range(n):
        for j in range(0, n-i-1):
            if AngleOff[j]>AngleOff[j+1]:
                AngleOff[j], AngleOff[j+1] = AngleOff[j+1], AngleOff[j]
                AxisIndices[:,j], AxisIndices[:,j+1]=AxisIndices[:,j+1], AxisIndices[:,j]     
    return (AxisIndices[:,0],AngleOff[0])

def maxi(liste):
    maxi = liste[0]
    ind=0
    longueur=len(liste)
    for i in range(longueur):
        if liste[i] >= maxi:
            maxi = liste[i]
            ind = i
    return maxi,ind

def mat2axisangle ( g ):
    eps = 1e-6
    r=np.zeros(3)
    ptheta = np.arccos((g[0][0] + g[1][1] + g[2][2] - 1) / 2);
    if ((ptheta) < eps):
        r[0] = 0;
        r[1] = 0;
        r[2] = 0;
    elif (ptheta < (1 - eps) * np.pi):
        r[0] = (g[1][2] - g[2][1]) / (2 * np.sin (ptheta));
        r[1] = (g[2][0] - g[0][2]) / (2 * np.sin (ptheta));
        r[2] = (g[0][1] - g[1][0]) / (2 * np.sin (ptheta));    
    else:    
        r[0] = np.sqrt ((g[0][0] + 1) / 2);
        r[1] = np.sqrt ((g[1][1] + 1) / 2);
        r[2] = np.sqrt ((g[2][2] + 1) / 2);
 
    C,m = maxi(r);
    for i in range(3):
        if (i != m): 
            r[i] = r[i]*np.sign (g[i][m]);
            
    return np.array(r,float),ptheta

def mat2rodrigues ( g ):
    r,theta= mat2axisangle(g);
    R=np.zeros(3)
    for i in range(3):
        R[i] = r[i] * np.tan (theta/ 2);
    return np.array(R,float) 

def rodrigues2mat (R):
    nr=norm(R);g= np.zeros((3,3))
    r=np.zeros(3)
    theta = 2 * np.arctan (nr);
    # Calculation of the angle/axis pair from R 
    for i in range(3):
        r[i] = R[i] / nr;
 
    g=axisangle2mat(r, theta)  
    return np.array(g,float)   

def axisangle2quaternion (r, theta):
    q=np.zeros(4)
    q[0] = np.cos (theta / 2);
    for i in range(3):
        q[i + 1] = r[i] * np.sin (theta / 2);
    return np.array(q,float)

def quaternion2axisangle (q):
    eps = 1e-6
    
    ptheta = 2*np.arccos (q[0]);
    if (ptheta < eps):
        r[0] = 1;
        r[1] = 0;
        r[2] = 0;
    else:
        for i in range(3):
            r[i] = q[i + 1] / np.sin (ptheta / 2);
    return np.array(r,float),ptheta

def mat2stproj(g,axe):
    eps = 1e-6
    p=g.dot(axe)
    # if abs(p)<eps:
    #     alpha=0;beta=0
    # else:
    p[0]=p[0]/norm(p)
    p[1]=p[1]/norm(p)
    p[2]=p[2]/norm(p)
    alpha=np.arccos(p[2])
    if abs(alpha)<eps :
        beta=0
    else:
        beta=np.arctan2(p[1]/np.sin(alpha),p[0]/np.sin(alpha))
        if (alpha > np.pi/2) :
            alpha=np.pi-alpha
            beta=beta+np.pi
    return np.tan(alpha/2),beta

def mat2eaproj( g, axe ):
    p1,p2 = mat2stproj(g, axe);    
    return np.sqrt(2)*np.sin( np.arctan(p1)),p2

def mat2eaprojxy( g, axe ):
    p1,p2 = mat2eaproj(g, axe);    
    return p1*np.cos(p2),p1*np.sin(p2)

def mat2stprojxy(g, axe):
    p1,p2=mat2stproj(g,axe);
    return p1*np.cos(p2),p1*np.sin(p2)
class crystalSymmetry:
    def __init__(self,name='iron',comment='Fe',LaueId=11,Laue="m3m",a=1,b=1,c=1,alpha=90*degree,beta=90*degree,gamma=90*degree,**kwargs):
        self.name=name
        self.comment=comment
        if Laue!="":
            self.LaueId=Laue2LaueId(Laue)
        self.a=a
        self.b=b
        self.c=c
        self.alpha=alpha      # angle between b and c
        self.beta=beta        # angle between c and a
        self.gamma=gamma      # angle between a and b
        self.MetricTensorDuel=metric_tensor_duel(a,b,c,alpha,beta,gamma)
        self.MetricTensor=metric_tensor_direct(a,b,c,alpha,beta,gamma)
        self.DirMatrix=cristal_matrix(a,b,c,alpha,beta,gamma,**kwargs)
        self.DirInvMat=np.linalg.inv(self.DirMatrix)
        self.NorMatrix=np.transpose(self.DirInvMat)
        if LaueId!="" :
            self.Laue=Space2Laue(LaueId)
        self.nbOp,self.rSym=symGenerator(LaueId);
        #self.Space=Space
    # def __repr__(self):
    #     return '#s <#s #s #s (#.2f° #.2f° #.2f° ) #s>'  %(
    #         self.name, self.a,self.b,self.c,self.alpha,self.beta,self.gamma,
    #         self.Laue)
        
def metric_tensor_duel(a,b,c,alpha,beta,gamma) :
    g= np.zeros((3,3))
    V=a*a*b*b*c*c*(1+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)-np.cos(alpha)*np.cos(alpha)-np.cos(beta)*np.cos(beta)-np.cos(gamma)*np.cos(gamma))
    g[0][0] = b*b*c*c*np.sin(alpha)*np.sin(alpha)/V
    g[0][1] = a*b*c*c*(np.cos(alpha)*np.cos(beta)-np.cos(gamma))/V
    g[0][2] = a*b*b*c*(np.cos(alpha)*np.cos(gamma)-np.cos(beta))/V
    
    g[1][0] = a*b*c*c*(np.cos(alpha)*np.cos(beta)-np.cos(gamma))/V
    g[1][1] = a*a*c*c*np.sin(beta)*np.sin(beta)/V
    g[1][2] = a*a*b*c*(np.cos(beta)*np.cos(gamma)-np.cos(alpha))/V
    
    g[2][0] = a*b*b*c*(np.cos(alpha)*np.cos(gamma)-np.cos(beta))/V
    g[2][1] = a*a*b*c*(np.cos(beta)*np.cos(gamma)-np.cos(alpha))/V
    g[2][2] = a*a*b*b*np.sin(gamma)*np.sin(gamma)/V   
    return np.array(g,float)

def metric_tensor_direct(a,b,c,alpha,beta,gamma) :
    g= np.zeros((3,3))
    g[0][0] = a*a
    g[0][1] = a*b*np.cos(gamma)
    g[0][2] = a*c*np.cos(beta)
   
    g[1][0] = a*b*np.cos(gamma)
    g[1][1] = b*b
    g[1][2] = b*c*np.cos(alpha)
   
    g[2][0] = a*c*np.cos(beta)
    g[2][1] = b*c*np.cos(alpha)
    g[2][2] = c*c
    return np.array(g,float)

def cristal_matrix(a,b,c,alpha,beta,gamma,convention='X//a',**kwargs):
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    cg = np.cos(gamma)
    sg=np.sin(gamma)
    dd=ca*ca+cb*cb+cg*cg   
    # print(convention)
    # print(convention=='Y//b')
    ax,ay,az,bx,by,bz,cx,cy,cz=0,0,0,0,0,0,0,0,0
    
    if (convention=='X//a'):
        #print('X//a')
        ax = a
        ay = 0 
        az = 0 
        
        bx = b * cg
        by = b * sg
        bz = 0
        
        cx = c * cb
        cy = c * (ca-cb*cg)/sg    
        cz = c * np.sqrt(1.0+2.0*ca*cb*cg-dd)/sg 
    
    if (convention=='Y//b'):
        ax = a * cg
        ay = a * sg
        az = 0 
        
        bx = 0
        by = b
        bz = 0
        
        cx = c * (cb-ca*cg)/sg 
        cy = c * ca 
        cz = c * np.sqrt(1.0+2.0*ca*cb*cg-dd)/sg         
        
    if (convention=='Z//c'):
        
        ax = a * np.sqrt(1.0+2.0*ca*cb*cg-dd)/sa
        ay = a * (cg-ca*cb)/sa
        az = a * cb
        
        bx = 0
        by = b * sa
        bz = b * ca
        
        cx = 0
        cy = 0
        cz = c

        
    StructureMat= np.array([[ax , bx, cx ],
                            [ay , by, cy ],
                            [az , bz, cz]])
        
    return StructureMat

def symOperators(LaueId):
    table = {
            1: np.array(
                [0]),
            2: np.array(
                [1,2,0,1,0]),
            3: np.array(
                [2,2,0,0,1,2,1,0,0]),
            4: np.array(
                [1,4,0,0,1]),
            5: np.array(
                [2,4,0,0,1,2,1,0,0]),
            6: np.array(
                [1,3,0,0,1]), 
            7: np.array(
                [2,3,0,0,1,2,1,0,0]),
            8: np.array(
                [1,6,0,0,1]),
            9: np.array(
                [2,6,0,0,1,2,0,1,0]),
            10: np.array(
                [3,2,0,0,1,2,1,0,0,3,1,1,1]),
            11: np.array(
                [3,4,0,0,1,2,1,1,0,3,1,1,1])
            }
    return table[LaueId]

def symGenerator(LaueId):    
    g=symOperators(LaueId)
    #print(symOperators(LaueId))
    nOper=0
    Oper={}
    q=0
    #print(g)
    if g[0]==0:
        q=1
        Oper=np.eye(3,3)
    elif g[0]==1: #un seul operateur
        angle=g[1]
        axis=[g[2],g[3],g[4]]
        for i in range(angle):
            ang=(2*np.pi*i)/angle            
            Oper[i]=axisangle2mat(axis, ang)
    elif g[0]==2: #deux operateurs
        angle1=g[1]
        angle2=g[5]
        axis1=[g[2],g[3],g[4]]
        axis2=[g[6],g[7],g[8]]
        q=0
        for i in range(g[1]):
            ang=(2*np.pi*i)/angle1 
            mat1=axisangle2mat(axis1, ang)
            for j in range(g[5]):
                ang=(2*np.pi*i)/angle2 
                mat2=axisangle2mat(axis2, ang)
                mat=mat1.dot(mat2)
                Oper[q]=mat
                q+=1
    elif g[0]==3: #trois operateurs
        angle1=g[1]
        angle2=g[5]
        angle3=g[9]
        axis1=[g[2],g[3],g[4]]
        axis2=[g[6],g[7],g[8]]
        axis3=[g[10],g[11],g[12]]
        q=0
        for i in range(g[1]):
            ang=(2*np.pi*i)/angle1 
            mat1=axisangle2mat(axis1, ang)
            for j in range(g[5]):
                ang=(2*np.pi*j)/angle2 
                mat2=axisangle2mat(axis2, ang)
                mat=mat1.dot(mat2)
                for k in range(g[9]):
                    ang=(2*np.pi*k)/angle3 
                    mat2=axisangle2mat(axis3, ang)
                    mat=mat.dot(mat2)
                    Oper[q]=mat
                    q+=1
    nOper=q
    return nOper,Oper
def SymSetMat( CS,SS, g ):    
    n = len(CS.rSym)
    # print(n)
    # print(CS.rSym)
    m,r = symGenerator(SS.LaueId)
    #print(m,r)
    q=np.zeros((n*m,3,3), dtype=float)
    w=0
    for i in range(n):
        #print(CS.rSym[i])
        q[i]=g.dot(CS.rSym[i])
        #print(q)
        w+=1      
    return q
def Space2Laue( space ):
    dico=('-1','2/m', 'mmm', '4/m', '4/mmm', '-3', '-3m','6/m', '6/mmm', 'm3', 'm3m')      
         
    return dico[space-1]  
def Laue2LaueId( Laue ):
    for i in range(11):
        if Space2Laue(i+1)==Laue:
            LaueId=i+1
    return LaueId  
# Generation Orientation
def Gen_vect (n1,n2):
    if n1<0 or n1>1 or n2<0 or n2 > 0:
        return
    alpha=np.arccos(2*n1-1)
    beta=2*np.pi*n2
    r=np.zeros(3)
    r[0]=np.sin(alpha)*np.cos(beta)
    r[1]=np.sin(alpha)*np.sin(beta)
    r[2]=np.cos(beta)
    return r

def Gen_OI_nb_eul(phi1,phi,phi2,dev,nb,cs):
    eul=np.zeros(3)
    ss=crystalSymmetry(LaueId=1)# triclinic
    eul[0]=phi*degree;eul[1]=phi*degree;eul[2]=phi2*degree
    g=eul2mat(eul)
    rmat=SymSetMat(cs,ss,g)    
    m=len(rmat) 
    v=np.zeros(3)
    gGen = np.zeros((3,3))
    eulGen=np.zeros((nb,3))
    choix=np.random.randint(m, size=nb) # chiffre 
    for j in range(nb):
        for i in range(m):
            r,thetamax=mat2axisangle(rmat[i])  
            n1,n2,n3 = np.random.rand(3)
            sg1,sg2,sg3=np.random.randint(2,size=3)
            
            if sg1==0:
                sign1=-1
            else:
                sign1=1
            if sg2==0:
                sign2=-1
            else:
                sign2=1 
            if sg3==0:
                sign3=-1
            else:
                sign3=1    
        # n2 = np.random.rand()
        # n3 = np.random.rand()
            beta=np.arccos(r[2])
            alpha=np.arccos(r[0])*np.sin(beta)
            betaG=beta+sign1*dev*degree*n1
            alphaG=alpha+sign2*dev*degree*n2
            theta=thetamax+sign3*dev*degree*n3
            v[0]=np.cos(alphaG)*np.sin(betaG)
            v[1]=np.sin(alphaG)*np.sin(betaG)
            v[2]=np.cos(betaG)
            gGen=axisangle2mat(v,theta)
            if choix[j]==i :
                eulGen[j]=mat2eul(gGen)
                break
    # r,thetamax=mat2axisangle(g)
    # v=np.zeros((nb,3))
    # n1 = np.random.rand(nb)
    # n2 = np.random.rand(nb)
    # beta=np.arccos(r[2])
    # alpha=np.arccos(r[0])*np.sin(beta)
    # betaG=beta*dev*degree*n1
    # alphaG=alpha*dev*degree*n2
    # for i in range(nb):
    #     v[i][0]=np.cos(alphaG[i])*np.sin(betaG[i])
    #     v[i][1]=np.sin(alphaG[i])*np.sin(betaG[i])
    #     v[i][2]=np.cos(betaG[i])
    # #gGen = np.zeros((nb,3,3))
    # gGen=np.array((3,3))
    # eulGen=np.zeros((nb,3))
    # n3 = np.random.rand(nb)
    # #theta=thetamax+dev*degree*np.power(dummy,1/3)
    # theta=thetamax+dev*degree*n3
    # #print(theta)
    # for i in range(nb):
    #     #gGen[i]=axisangle2mat(r,theta[i])
    #     gGen=axisangle2mat(v[i],theta[i])
    #     eulGen[i]=mat2eul(gGen)
    return eulGen
        
# Graphique
def cercle(ax,Xc,Yc,R,coul='k',linewidth=2):
    theta = np.linspace(0, 2*np.pi, 100)
    x1 = R*np.cos(theta)
    x2 = R*np.sin(theta)
    ax.plot(x1, x2, coul,linewidth)
    
def EquiArea(ax,R=200,Xc=0,Yc=0,step=30):
    #Rw=np.tan(15*np.pi/180);
    Rw=0;
    cercle(ax,Xc,Yc,R);
    for i in range(step,90,step):
        rh=0.5*i*degree;
        xi=(R*np.tan(rh));
        cercle(ax,Xc,Yc,xi,coul='#C0C0C0',linewidth=1);
    for i in range(0,190,step):
        ph=i*degree;
        x=(R*np.cos(ph));	y=R*np.sin(ph);
        xi=x;		yi=y;
        xf=(x*Rw);	yf=(y*Rw);        
        ax.plot([xi, xf], [yi, yf],'#C0C0C0') 
        ax.plot([-xi, -xf], [-yi, -yf],'#C0C0C0')
def wulff(ax,R=200,Xc=0,Yc=0,step=15):    
	#grands cercles
    xi=0;yi=0;
    for i in range(-85, 5,step):
        delt=i*degree;
        cosD=np.cos(delt);	sinD=np.sin(delt);
        for j in range(10,175,5):
            epsi=j*degree;
            sinEp=np.sin(epsi);
            rho=np.arccos(cosD*sinEp);
            TgRhos2=np.tan(rho/2);
            if (j==90):
                phi=np.pi/2; 
            else:
                phi=np.arctan(sinEp*sinD/np.cos(epsi));
            if (phi<0): 
                phi+=np.pi;
            xf=R*np.cos(phi)*TgRhos2;
            yf=-R*np.sin(phi)*TgRhos2;
            if (j==10):
                xi=xf;	yi=yf;
            if (i%10==0) :
                ax.plot([xi, xf], [yi, yf],'tab:gray') 
                ax.plot([-xi, -xf], [-yi, -yf],'tab:gray')
            else:
                ax.plot([xi, xf], [yi, yf],'#D3D3D3') 
                ax.plot([-xi, -xf], [-yi, -yf],'#D3D3D3') 
            xi=xf;yi=yf;
    #petits cercles
    for i in range(10, 175,step):
        delt=i*degree;
        tanD=np.tan(delt);	sinD=np.sin(delt);
        for j in range(0,180+step,step):
            epsi=j*degree;
            sinEp=np.sin(epsi);
            rho=np.arccos(sinD*sinEp);
            TgRhos2=np.tan(rho/2);
            if (j==90): 
                phi=np.pi/2; 
            else: 
                phi=np.arctan(1/(tanD*np.cos(epsi)));
            if (phi<0): phi+=np.pi;
            yf=R*np.cos(phi)*TgRhos2;
            xf=-R*np.sin(phi)*TgRhos2;
            if (j==0) :xi=xf;	yi=yf;
            if (i%10==0):
                ax.plot([xi, xf], [yi, yf],'tab:gray') 
                ax.plot([-xi, -xf], [-yi, -yf],'tab:gray')
            else:    
                ax.plot([xi, xf], [yi, yf],'#D3D3D3') 
                ax.plot([-xi, -xf], [-yi, -yf],'#D3D3D3') 
            xi=xf;	yi=yf
    cercle(ax,Xc,Yc,R)                
#if __name__ == "__main__":
    # cs=crystalSymmetry('Copper')
    # phi1=90;phi=35;phi2=45;nb=20
    # print(Gen_OI_nb_eul(phi1,phi,phi2,nb,cs))
    # for i in range(11):
    #     print(Space2Laue(i+1))
    
    # r=[1,1,1];theta=np.pi/6
    # print(axisangle2mat(r, theta))
    # r=vector3d(1, 2, 3)
    # rnorm=r.Normalize()
    # print(rnorm)
    # v=vector3d(2,3,4)
    # print(r+v)
    # print(r-v)
    # print(r*v)
    # print(2 *v*3)
    # polar_angle = 90*degree;
    # azimuth_angle = 45*degree;
    # p=vector3d().byPolar(polar_angle,azimuth_angle)
    # print(p)  
    # q=r.normalize()
    # print(q)
    # g=[[ 0.11774846,0.96075772,  0.25115713],
    #   [-0.8042109,  -0.0561144 ,  0.59168911],
    #   [ 0.58256342, -0.27165378,  0.76604444]]
    # eul=mat2eul( g )
    # print('(%.3f° %.3f° %.3f°)' %(eul[0]/degree,eul[1]/degree,eul[2]/degree))
    # hkl=[1,2,3];uvw=[1,1,0]
    # print(miller2mat ( hkl,uvw ))
    # Axis=[0.8, 0.7,0.1]; MaxIndices=5;
    # AxisInd,AngleOff=axis2indices( Axis, MaxIndices )# 
    # print('[%.3f , %.3f, %.3f] =[%s %s %s] à %.3f° près' %(Axis[0],Axis[1],Axis[2],AxisInd[0],AxisInd[1],AxisInd[2],AngleOff/degree)) 
    #print(mat2axisangle ( g ))
    # R=mat2rodrigues ( g )
    # print(R)
    # g=rodrigues2mat (R)
    # print(g)
    # r=[1,1,1];theta=45*degree
    # q=axisangle2quaternion (r, theta)
    # print(q)
    # r, theta=quaternion2axisangle(q)
    # print('[%s %s %s] %.3f°' %(r[0],r[1],r[2],theta/degree))
    #print(symOperators(9))
    #print(symGenerator(6))
    # cs=cristal_matrix(1,2,3,90*degree,90*degree,120*degree,convention='Z//c')
    # print(cs)
    # fig=plt.figure(figsize=(8,8))
    # ax = fig.add_subplot()
    # ax.xaxis.set_visible(False)
    # ax.yaxis.set_visible(False)
    # ax.set_frame_on(False)
    # cs=crystalSymmetry('Copper')
    # # #print(cs)
    # ss=crystalSymmetry(LaueId=1)
    # #print(ss)
    # AxisHKL=[1,1,1]
    # plt.title(str(AxisHKL),fontsize=20)#[%s %s %s]' %(AxisHKL[0],AxisHKL[1],AxisHKL[2])
    # axe=cs.DirInvMat.dot(AxisHKL)
    # axe=axe[:]/norm(axe) 
    #phi1=90;phi=35;phi2=45 # Texture Cuivre
    #phi1=59;phi=37;phi2=63 # Texture S
    #phi1=0;phi=45;phi2=90  # Texture Goss G
    #phi1=35;phi=45;phi2=90 # Laiton
    # phi1=0;phi=0;phi2=0    # Cube
    # nb=200
    
    # eul=Gen_OI_nb_eul(phi1,phi,phi2,5,nb,cs)
    # print(eul)
    #print('%s %s %s' %(eul[:,0]/degree, eul[:,1]/degree,eul[:,2]/degree))
    #eul=[90*degree,35*degree,45*degree]
    # print(eul)
    # R=200
    # for i in range(nb):
    #     g = eul2mat(eul[i])
    # # print(g)
    #     rmat=SymSetMat(cs,ss,g)
    # # print(rmat)
    #     m=len(rmat)        
    #     for j in range(m):
    # #       #x,y=mat2stprojxy(rmat[i],axe)
    #         x,y=mat2eaprojxy( rmat[j],axe) 
    #         ax.scatter(x*R, y*R, edgecolor = 'none',s = 20,marker = 'o',c = 'red',zorder=10)

    # EquiArea(ax)
    # ax.axis('equal')
    
        # Qian
        
from numpy.linalg import inv
def IPFcercle(ax,Xc=0,Yc=0,R=200,coul='k',linewidth=2,step=30):
    theta = np.linspace(0, np.pi/6, 30)
    x1 = R*np.cos(theta)
    x2 = R*np.sin(theta)
    ax.plot(x1, x2, coul,linewidth)
    
    for i in range(0,31,step):
        ph=i*degree;
        x=(R*np.cos(ph));	y=R*np.sin(ph);

        ax.plot([x, 0], [y, 0],'k') 
        
def mat2IPFproj(p):
    eps = 1e-6
    
    # if abs(p)<eps:
    #     alpha=0;beta=0
    # else:
    p[0]=p[0]/norm(p)
    p[1]=p[1]/norm(p)
    p[2]=p[2]/norm(p)
    alpha=np.arccos(p[2])
    if abs(alpha)<eps :
        beta=0
    else:
        beta=np.arctan2(p[1]/np.sin(alpha),p[0]/np.sin(alpha))
        if (alpha > np.pi/2) :
            alpha=np.pi-alpha
            beta=beta+np.pi
    return np.tan(alpha/2),beta

def mat2IPFprojxy(p):
    p1,p2=mat2IPFproj(p);
    return p1*np.cos(p2),p1*np.sin(p2)

def mat2IPFeaproj(p):
    p1,p2 = mat2IPFproj(p);    
    return np.sqrt(2)*np.sin( np.arctan(p1)),p2

def mat2IPFeaprojxy(p):
    p1,p2 = mat2IPFeaproj(p);    
    return p1*np.cos(p2),p1*np.sin(p2)

def eul2mat( e ):
    g= np.zeros((3,3))
    e[0]=e[0]*degree
    e[1]=e[1]*degree
    e[2]=e[2]*degree
    g[0,0] = np.cos (e[0]) * np.cos (e[2]) - np.sin (e[0]) * np.sin (e[2]) * np.cos (e[1]);
    g[0,1] = np.sin (e[0]) * np.cos (e[2]) + np.cos (e[0]) * np.sin (e[2]) * np.cos (e[1]);
    g[0,2] = np.sin (e[2]) * np.sin (e[1]);
 
    g[1,0] = -np.cos (e[0]) * np.sin (e[2]) - np.sin (e[0]) * np.cos (e[2]) * np.cos (e[1]);
    g[1,1] = -np.sin (e[0]) * np.sin (e[2]) + np.cos (e[0]) * np.cos (e[2]) *np. cos (e[1]);
    g[1,2] = np.cos (e[2]) * np.sin (e[1]);
 
    g[2,0] = np.sin (e[0]) * np.sin (e[1]);
    g[2,1] = -np.cos (e[0]) * np.sin (e[1]);
    g[2,2] = np.cos (e[1]);
    return np.array(g,float)
