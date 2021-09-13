# -*- coding: utf-8 -*-
"""
@author: Qian WANG
First, change the folder in line 45
"""
# plot IPDF
import numpy as np
import pyOriLib.pyOriLib as pyo
import matplotlib.pyplot as plt
import matplotlib

degree = np.pi/180

# draw IPF triangle
      


fig=plt.figure()
ax = fig.add_subplot(1,1,1)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

ax.set_frame_on(False)
pyo.IPFcercle(ax)
ax.axis('equal')

# orientation data
def LoadFileEul( fname ):
    eul= [];hardness= [];
    try:
        fid = open(fname, "r")
        lines = fid.readlines() 
       
        for line in lines:
            line = line.strip()
            dum1,dum2, dum3, dum4= line.split('\t')  
            
            eul+=[[float(dum1),float(dum2),float(dum3)]]
            hardness+=[[float(dum4)]]
    except:
        print("Le fichier", fname, "est introuvable")    
    return np.array(eul,float), np.array(hardness,float)


data_folder ="E:\\PHD\\specilist\\python\\Groupe_Travail_Jeudi 24 juin\\"
file_to_open = data_folder + "qqti.eul"#"TEST.eul"
eul,hardness=LoadFileEul( file_to_open )

# draw orientation in IPF




cs=pyo.crystalSymmetry(name='titane',comment='Ti',Laue='6/mmm',LaueId=9,
                       a=2.95,b=2.95,c=4.68,alpha=90*degree,beta=90*degree,gamma=120*degree,
                       convention='Y//b')
ss=pyo.crystalSymmetry(name='ss',comment='ss',LaueId=1,Laue='-1')
ss.LauId=1
stereo=0

from numpy.linalg import inv

ND=np.array([[0],[0],[1]])

axe=ND[:]/pyo.norm(ND) 
N=len(eul)


rmat=np.empty((12,),dtype=object)

rmat[0]=np.array([[1,0,0],[0,1,0],[0,0,1]])
rmat[1]=np.array([[0.5,np.sqrt(3)/2,0],[-np.sqrt(3)/2,0.5,0],[0,0,1]])
rmat[2]=np.array([[-0.5,np.sqrt(3)/2,0],[-np.sqrt(3)/2,-0.5,0],[0,0,1]])
rmat[3]=np.array([[-1,0,0],[0,-1,0],[0,0,1]])
rmat[4]=np.array([[-0.5,-np.sqrt(3)/2,0],[np.sqrt(3)/2,-0.5,0],[0,0,1]])
rmat[5]=np.array([[0.5,-np.sqrt(3)/2,0],[np.sqrt(3)/2,0.5,0],[0,0,1]])
rmat[6]=np.array([[1,0,0],[0,-1,0],[0,0,-1]])
rmat[7]=np.array([[0.5,np.sqrt(3)/2,0],[np.sqrt(3)/2,-0.5,0],[0,0,-1]])
rmat[8]=np.array([[-0.5,np.sqrt(3)/2,0],[np.sqrt(3)/2,0.5,0],[0,0,-1]])
rmat[9]=np.array([[-1,0,0],[0,1,0],[0,0,-1]])
rmat[10]=np.array([[-0.5,-np.sqrt(3)/2,0],[-np.sqrt(3)/2,0.5,0],[0,0,-1]])
rmat[11]=np.array([[0.5,-np.sqrt(3)/2,0],[-np.sqrt(3)/2,-0.5,0],[0,0,-1]])


g=np.empty((45,),dtype=object)
for j in range(N):
    g[j] = pyo.eul2mat(eul[j]) 

    
    
    m=len(rmat)   
    x=np.zeros(N);y=np.zeros(N);k=np.zeros(m);a=np.zeros(m);b=np.zeros(m)     
    for i in range(12):
        G=np.empty((12,),dtype=object)
        p=rmat[i].dot(g[j]).dot(axe)

        a[i],b[i]=pyo.mat2IPFprojxy(p)
       
        k[i]=np.arctan2(b[i],a[i])/degree
        if 0<=k[i]<=30:
            y[j]=b[i]
            x[j]=a[i]
        # else:
        #     b[i]=np.NaN
        #     a[i]=np.NaN



    
    R=200
    x=x[:]*R
    y=y[:]*R
    h=hardness[:]  
    plt.scatter(x, y, edgecolor = 'none',s = 100,marker = 'o',c =h,zorder=10)
# plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.family"] = "normal"
plt.rcParams.update({'font.size': 15})
plt.colorbar(label='Hardness')

plt.text(-20,-15,'[0001]')
plt.text(180,-15,r"[10$\overline{1}$0]")
plt.text(160,110,r"[1$\overline{2}$10]")


fig.savefig('IPF_hardness.tif', format='tif', dpi=1200)
   



