# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 08:32:41 2021

@author: lecomte5
"""
from tkinter import *
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from tkinter.constants import DISABLED, NORMAL
import numpy as np
import pyOriLib.pyOriLib as pyo
from tkinter.colorchooser import askcolor 
#-------------------------------------------------------
#      DECLARATION VARIABLES GLOBALES
#-------------------------------------------------------
degree = np.pi/180
R=200;Xc=0;Yc=0;
#cs=pyo.crystalSymmetry('Copper')
cs=pyo.crystalSymmetry(name='titane',comment='Ti',Laue='6/mmm',LaueId=9,
                       a=2.95,b=2.95,c=5.68,alpha=90*degree,beta=90*degree,gamma=120*degree,
                       convention='Y//b')
ss=pyo.crystalSymmetry(LaueId=1)
ss.LaueId=1
stereo=0
eul1=[46*degree,56*degree,89*degree]
eul2=[0*degree,0*degree,0*degree]
AxisHKL=[0,0,1]
coul1="red" 
coul2="green"
marker1='o'
marker2='o'
size1=10
size2=10
shape=[ ".",",","o","v","^","<",">","1","2","3","4","8","s","p","P","*","h","H",
       "+","x","X","D","d"]
sizeList=["10","20","30","50"]
# boite de couleur
def ChoixCouleur1():
    global coul1
    result = askcolor(title = "Tkinter Color Chooser")
    coul1=result[1]
def ChoixCouleur2():
    global coul2
    result = askcolor(title = "Tkinter Color Chooser")
    coul2=result[1]
#-------------------------------------------------------
#      CREATION D'UNE FENTRE SECONDAIRE : FEN2 : CS
#-------------------------------------------------------
def createCS(): 
    #--------------------------------
    # Création de quelques fonctions
    #--------------------------------
    def RecupData():
        # Cette fonction permet de récupérer les paramètres cristallins et sauve dans cs
        global cs

        Laue=combo.get()
        LaueId=pyo.Laue2LaueId(Laue)
        cs=pyo.crystalSymmetry(name=varName.get(), comment=varComment.get(),
                           LaueId=LaueId,Laue=Laue, a=varA.get(),b=varB.get(),c=varC.get(),
                           alpha=varAlpha.get()*degree,beta=varBeta.get()*degree,gamma=varGamma.get()*degree)
        # manque encore les conventions X//a, Y//b ou Z//c
        
    def clic_combo(event): 
        # Fonction pour la gestion de la ComboBox (liste des groupes de Laue)
        Laue=combo.get()
        if Laue=='-1':
            EntA["state"]=NORMAL;EntB["state"]=NORMAL;EntC["state"]=NORMAL
            EntAlpha['state'] = NORMAL;EntBeta['state'] = NORMAL;EntGamma['state'] = NORMAL
        if Laue=='2/m':
            EntA["state"]=NORMAL;EntB["state"]=NORMAL;EntC["state"]=NORMAL
            EntAlpha['state'] = DISABLED;EntBeta['state'] = NORMAL;EntGamma['state'] = DISABLED
            varAlpha.set("90");varGamma.set("90");        
        if Laue=='mmm':
            EntA["state"]=NORMAL;EntB["state"]=NORMAL;EntC["state"]=NORMAL
            EntAlpha['state'] = DISABLED;EntBeta['state'] = DISABLED;EntGamma['state'] = DISABLED
            varAlpha.set("90");varBeta.set("90");varGamma.set("90");
        if Laue=='4/mmm' or Laue=='4/m':
            EntA["state"]=NORMAL;EntB["state"]=DISABLED;EntC["state"]=NORMAL
            EntAlpha['state'] = DISABLED;EntBeta['state'] = DISABLED;EntGamma['state'] = DISABLED
            varAlpha.set("90");varBeta.set("90");varGamma.set("90");
        if Laue=='-3' or Laue=='-3m':
            EntA["state"]=NORMAL;EntB["state"]=DISABLED;EntC["state"]=NORMAL
            EntAlpha['state'] = DISABLED;EntBeta['state'] = DISABLED;EntGamma['state'] = DISABLED
            varAlpha.set("90");varBeta.set("90");varGamma.set("120");
        if Laue=='6/m' or Laue=='6/mmm':
            EntA["state"]=NORMAL;EntB["state"]=DISABLED;EntC["state"]=NORMAL
            EntAlpha['state'] = DISABLED;EntBeta['state'] = DISABLED;EntGamma['state'] = DISABLED
            varAlpha.set("90");varBeta.set("90");varGamma.set("120");
        if Laue=='m3' or Laue=='m3m':
            EntA["state"]=NORMAL;EntB["state"]=DISABLED;EntC["state"]=DISABLED
            EntAlpha['state'] = DISABLED;EntBeta['state'] = DISABLED;EntGamma['state'] = DISABLED
            varAlpha.set("90.0");varBeta.set("90.0");varGamma.set("90.0");
            
    def EnterVar(event):
        # Cette fonction permet de "recopier" ce qui est écrit sur Entry : EntA
        Laue=combo.get()
        if Laue=='m3' or Laue=='m3m' :
            result = EntA.get()
            varB.set(result)
            varC.set(result)
        if Laue=='4/mmm' or Laue=='4/m' or Laue=='6/m' or Laue=='6/mmm':
            result = EntA.get()
            varB.set(result)
               
    
    #------------------------------------         
    #         Instanciation de fen2   
    #------------------------------------     
    fen2=Toplevel(fen1)
    fen2.resizable(0,0) 
    
    varName=StringVar();varName.set(cs.name)
    varComment=StringVar();varComment.set(cs.comment)
    varA=DoubleVar();varB=DoubleVar();varC=DoubleVar()
    varA.set(cs.a);varB.set(cs.b);varC.set(cs.c);
    varAlpha=DoubleVar();varBeta=DoubleVar();varGamma=DoubleVar()
    varAlpha.set(cs.alpha/degree);varBeta.set(cs.beta/degree);varGamma.set(cs.gamma/degree);
    fen2.title("Crystal Symmetry")
    fen2.geometry("350x280")
    fen2.minsize(250,200) 
    fen2.iconbitmap("LOGO LEM3.ico")
    fen2.config(bg="#1AEEC1")
    fram1 = Frame(fen2,bg="#1AEEC1")
    Label(fram1,text="Crystal Reference",bg="#1AEEC1",font=("Arial",20),fg="black").pack(side=TOP)
    Label(fram1,text="Crystal Symmetry",bg="#1AEEC1", font=("Arial,25"),fg="white").pack(side=BOTTOM)
    fram1.grid(row=0,column=0,sticky="n")
    
    fram2 = LabelFrame(fen2,text="Mineral",bg="#1AEEC1",fg="white")
    Label(fram2,text="mineral name",bg="#1AEEC1").grid(row=0,column=0)
    Entry(fram2,textvariable=varName).grid(row=0,column=1)
    Label(fram2,text="comment",bg="#1AEEC1").grid(row=1,column=0)
    Entry(fram2,textvariable=varComment).grid(row=1,column=1)
    fram2.grid(row=1,column=0)
    
    fram3=LabelFrame(fen2,text="Crystal Coordinate System",bg="#1AEEC1",fg="white")
    Label(fram3,text="Point Group",bg="#1AEEC1").grid(row=0,column=0)
    vlist = ['-1','2/m','mmm','4/m','4/mmm','-3','-3m','6/m','6/mmm','m3','m3m',]    
    combo=ttk.Combobox(fram3, values = vlist,width=10, state= 'readonly')
    combo.set(cs.Laue)    
    combo.grid(row=0,column=1, columnspan=3)
    Label(fram3,text='Axis length',bg="#1AEEC1").grid(row=1,column=0)
    Label(fram3,text='a',bg="#1AEEC1").grid(row=1,column=2)
    EntA=Entry(fram3,textvariable=varA,width=5)
    EntA.bind("<KeyRelease>",EnterVar)
    EntA.grid(row=1,column=3)
    Label(fram3,text='b',bg="#1AEEC1").grid(row=1,column=4)
    EntB=Entry(fram3,textvariable=varB,width=5)
    #EntB.bind("<KeyRelease>",EnterVar)
    EntB.grid(row=1,column=5)
    Label(fram3,text='c',bg="#1AEEC1").grid(row=1,column=6)
    EntC=Entry(fram3,textvariable=varC,width=5)
    EntC.grid(row=1,column=7)
    Label(fram3,text='Axis Angle',bg="#1AEEC1").grid(row=2,column=0)
    Label(fram3,text='alpha',bg="#1AEEC1").grid(row=2,column=2)
    EntAlpha=Entry(fram3,textvariable=varAlpha,width=5,state = NORMAL)
    EntAlpha.grid(row=2,column=3)
    Label(fram3,text='beta',bg="#1AEEC1").grid(row=2,column=4)
    EntBeta=Entry(fram3,textvariable=varBeta,width=5)
    EntBeta.grid(row=2,column=5)
    Label(fram3,text='gamma',bg="#1AEEC1").grid(row=2,column=6)
    EntGamma=Entry(fram3,textvariable=varGamma,width=5,state='disabled')
    EntGamma.grid(row=2,column=7)
    
    fram3.grid(row=2,column=0)
    btnOK=Button(fen2,text="Save CS before Quit", command = RecupData).grid(row=3,column=0) 
    Button(fen2,text="Quit",command=fen2.destroy).grid(row=4,column=0)
    combo.bind("<<ComboboxSelected>>", clic_combo)
#-------------------------------------------------------
#      PROGRAMME PRINCIPAL --> FEN1 : FENETRE PRINCIPAL
#-------------------------------------------------------    
fen1 = Tk()
fen1.title("OriPlot")
fen1.geometry("800x400")
fen1.minsize(R,R) # au mini le rayon
fen1.iconbitmap("LOGO LEM3.ico")
# Déclaration des variables globales
varPhi11=DoubleVar()
varPhi11.set(str(eul1[0]/degree))
varPhi1=DoubleVar()
varPhi1.set(str(eul1[1]/degree))
varPhi21=DoubleVar()
varPhi21.set(str(eul1[2]/degree))
varPhi12=DoubleVar()
varPhi12.set(str(eul2[0]/degree))
varPhi2=DoubleVar()
varPhi2.set(str(eul2[1]/degree))
varPhi22=DoubleVar()
varPhi22.set(str(eul2[2]/degree))
varH = IntVar()
varH.set(str(AxisHKL[0]))
varK = IntVar()
varK.set(str(AxisHKL[1]))
varL = IntVar()
varL.set(str(AxisHKL[2]))

#FRAME 1
frame1 = LabelFrame(fen1, text="Ploting Parameter")
frame1.pack(fill='both', expand=False, side=LEFT)
# fonction de contrôle --> cs
def calculDG():
    
    eul1[0] = varPhi11.get()*degree
    eul1[1] = varPhi1.get()*degree
    eul1[2] = varPhi21.get()*degree
    
    eul2[0] = varPhi12.get()*degree
    eul2[1] = varPhi2.get()*degree
    eul2[2] = varPhi22.get()*degree
    
    g1 = pyo.eul2mat(eul1) 
    print(g1)
    g2 = pyo.eul2mat(eul2)
    print(g2)
    # dg=np.linalg.inv(g1).dot(g2)
    # print(dg)
    rmat1=pyo.SymSetMat(cs,ss,g1) 
    rmat2=pyo.SymSetMat(cs,ss,g2)
    m1=len(rmat1)
    m2=len(rmat2)
    # for i in range(m1):
    #     print(np.linalg.det(rmat1[i]))

    dg=np.zeros((m1*m2,3,3),dtype=float)
    print(type(dg))
    print('m1 = %s et m2 = %s' %(m1,m2))
    count=0
    for i in range(m1):
        for j in range(m2):
            dg[count]=np.linalg.inv(rmat1[i]).dot(rmat2[j])
            count+=1
    omegamin=1000
    rmin=np.zeros(3)
    ind=0
    dgmin=np.zeros((3,3))
    for i in range(count):
        r,omega=pyo.mat2axisangle(dg[i])
        if omega<omegamin:
            omegamin=omega
            rmin=r
            ind=i
            dgmin=dg[i]
    print(dgmin)
    print('[%.2f %.2f %.2f] - %.1f°' %(*rmin,omegamin/degree))
    print(rmin[0])
    miller,ang=pyo.axis2indices(rmin,4)
    print('[%2d %2d %2d] = [%.2f %.2f %.2f] à %.1f°' %(*miller,*rmin,ang/degree))


          
def printCS():
    print(cs.name)
    print(cs.comment)
    print("%s %s %s -- %s %s %s " %(cs.a,cs.b,cs.c,cs.alpha/degree,cs.beta/degree,cs.gamma/degree))    
# Fonction principale : Trace la figure de Pôle
def PlotPDF():  
    global eul
    global AxisHKL
    global ss
    global stereo
    global marker1
    global marker2
    global size1
    global size2
    
    eul1[0] = varPhi11.get()*degree
    eul1[1] = varPhi1.get()*degree
    eul1[2] = varPhi21.get()*degree
    
    eul2[0] = varPhi12.get()*degree
    eul2[1] = varPhi2.get()*degree
    eul2[2] = varPhi22.get()*degree
    
#     r=[1,1,1]
#     dg = pyo.axisangle2mat(r,np.pi)
    
    
    AxisHKL[0] = varH.get()
    AxisHKL[1] = varK.get()
    AxisHKL[2] = varL.get()  
    
    stereo=varPrj.get()
    
    ss.LaueId=varSS.get()
    #print(ss.LaueId)
    marker1=comboShape1.get()
    marker2=comboShape2.get()
    
    size1=int(comboSize1.get())*10
    size2=int(comboSize2.get())*10
    
    ax.clear()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    ax.set_frame_on(False)
    ax.set_title(str(AxisHKL))
    axe=cs.DirInvMat.dot(AxisHKL)
    axe=axe[:]/pyo.norm(axe) 
    g1 = pyo.eul2mat(eul1) 
    g2 = pyo.eul2mat(eul2)
    #g2=g1.dot(dg)
    rmat1=pyo.SymSetMat(cs,ss,g1) 
    rmat2=pyo.SymSetMat(cs,ss,g2)
    m1=len(rmat1)   
    x=np.zeros(m1);y=np.zeros(m1)     
    for i in range(m1):
        if stereo==1 :
            x[i],y[i]=pyo.mat2stprojxy(rmat1[i],axe)
        else:
            x[i],y[i]=pyo.mat2eaprojxy(rmat1[i],axe)
    x=Xc+x[:]*R
    y=Yc+y[:]*R
    ax.scatter(x, y, edgecolor = 'none',s = size1,marker = marker1,c = coul1,zorder=10)
    m2=len(rmat2)   
    x=np.zeros(m2);y=np.zeros(m2)     
    for i in range(m2):
        if stereo==1 :
            x[i],y[i]=pyo.mat2stprojxy(rmat2[i],axe)
        else:
            x[i],y[i]=pyo.mat2eaprojxy(rmat2[i],axe)
    x=Xc+x[:]*R
    y=Yc+y[:]*R
    ax.scatter(x, y, edgecolor = 'none',s = size2,marker = marker2,c = coul2,zorder=10)
    if stereo==1:       
        pyo.wulff(ax,step=5)
    else:
        pyo.EquiArea(ax,step=5)
    ax.axis('equal')
    
    fig.canvas.draw()
    fig.canvas.flush_events()
        
    
sframe1= LabelFrame(frame1, text="Euler Angle 1")
sframe1.grid(row=0, column=0,sticky='n')
Label(sframe1, text="phi1").grid(row=0, column=0)
Entry(sframe1, textvariable=varPhi11,width=5).grid(row=0, column=1)
Label(sframe1, text="phi").grid(row=0, column=2)
Entry(sframe1, textvariable=varPhi1,width=5).grid(row=0, column=3)
Label(sframe1, text="phi2").grid(row=0, column=4)
Entry(sframe1, textvariable=varPhi21,width=5).grid(row=0, column=5)
Button(sframe1, text="Plot", command=PlotPDF).grid(row=1, column=3, columnspan=3)
bt1=Button(sframe1, text="Color", command=ChoixCouleur1).grid(row=1, column=0)
comboShape1=ttk.Combobox(sframe1, values = shape,width=3, state= 'readonly')
comboShape1.set('o') 
comboShape1.grid(row=1,column=1)
comboSize1=ttk.Combobox(sframe1, values = sizeList,width=2, state= 'readonly')
comboSize1.set('10') 
comboSize1.grid(row=1,column=2)

sframe1bis= LabelFrame(frame1, text="Euler Angle 2")
sframe1bis.grid(row=1, column=0,sticky='n')
Label(sframe1bis, text="phi1").grid(row=0, column=0)
Entry(sframe1bis, textvariable=varPhi12,width=5).grid(row=0, column=1)
Label(sframe1bis, text="phi").grid(row=0, column=2)
Entry(sframe1bis, textvariable=varPhi2,width=5).grid(row=0, column=3)
Label(sframe1bis, text="phi2").grid(row=0, column=4)
Entry(sframe1bis, textvariable=varPhi22,width=5).grid(row=0, column=5)
Button(sframe1bis, text="Plot", command=PlotPDF).grid(row=1, column=3, columnspan=3)
bt2=Button(sframe1bis, text="Color", command=ChoixCouleur2).grid(row=1, column=0)
comboShape2=ttk.Combobox(sframe1bis, values = shape,width=3, state= 'readonly')
comboShape2.set('o') 
comboShape2.grid(row=1,column=1)
comboSize2=ttk.Combobox(sframe1bis, values = sizeList,width=3, state= 'readonly')
comboSize2.set('10') 
comboSize2.grid(row=1,column=2)

sframe2= LabelFrame(frame1, text="Axis")
sframe2.grid(row=2, column=0)
Label(sframe2, text="h").grid(row=0, column=0)
Spinbox(sframe2, from_=-10, to=10,width=5, textvariable=varH).grid(row=0, column=1)
Label(sframe2, text="k").grid(row=0, column=2)
Spinbox(sframe2, from_=-10, to=10,width=5, textvariable=varK).grid(row=0, column=3)
Label(sframe2, text="l").grid(row=0, column=4)
Spinbox(sframe2, from_=-10, to=10,width=5, textvariable=varL).grid(row=0, column=5)
Button(sframe2, text="Plot", command=PlotPDF).grid(row=1, column=2, columnspan=3)

sframe3= LabelFrame(frame1, text="Symmetry")
sframe3.grid(row=3, column=0, sticky="nsew")
Button(sframe3, text="CS",command=createCS).grid(row=0, column=0)
valSS=[1,3];etiqs=['triclinic','orthorhombic'];varSS=IntVar();varSS.set(ss.LaueId)
Radiobutton(sframe3,variable=varSS,text=etiqs[0],value=valSS[0]).grid(row=0, column=1)
Radiobutton(sframe3,variable=varSS,text=etiqs[1],value=valSS[1]).grid(row=0, column=2)
valProj=[1,0];etiqs=['Stereo','equiArea'];varPrj=IntVar();varPrj.set(stereo)
Radiobutton(sframe3,variable=varPrj,text=etiqs[0],value=valProj[0]).grid(row=1,column=0)
Radiobutton(sframe3,variable=varPrj,text=etiqs[1],value=valProj[1]).grid(row=1,column=1)
Button(sframe3, text="Print CS", command=printCS).grid(row=2, column=0, columnspan=2)
Button(sframe3, text="Calcul Dg", command=calculDG).grid(row=2, column=2, columnspan=2)

# FRAME 2
frame2 = LabelFrame(fen1, text="Pole Figure")
frame2.pack(fill='both', expand=True, side=RIGHT)

fig=plt.figure()
ax = fig.add_subplot(1,1,1)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)
if stereo==1:
    pyo.wulff(ax)
else:
    pyo.EquiArea(ax)
ax.axis('equal')

can= FigureCanvasTkAgg(fig, master=frame2)  # ici le lien entre fig et fen1
can.draw()
can.get_tk_widget().pack(side="top",fill='both',expand=True)

fen1.mainloop()
