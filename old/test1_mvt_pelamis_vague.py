#essai d'obtention du mouvement du pelamis dans la vague
#manque pas mal d'information pour finir comme
# l'expression de la hauteur d'eau \eta(M,t)
# l'expression du champ de vitesse des particules d'eau
import numpy as np
import matplotlib.pyplot as plt

def passag_cyl_xyz(pos,theta,alpha):
    """ prend un vecteur position dans la base cylindrique du pelamis et renvoie la position en cartersien dans 
    le référentiel du laboratoire.
    les deux angles sont aussi passés en argument"""
    
    P=np.array([[-np.sin(theta)*np.sin(alpha),-np.cos(theta)*np.sin(alpha),np.cos(alpha)],[np.cos(theta),-np.sin(theta),0],[np.sin(theta)*np.cos(alpha),np.cos(theta)*np.cos(alpha),np.sin(alpha)]])
    return np.dot(P,pos)
    
def kvag(w):
    """ fonction qui renvoie la valeur du vecteur d'onde pour une pulsation donnée
    Dépend du modèle choisi d'eau profonde ou non """
    #premier modele TROP simple juste pour tester
    c=1
    return w/c

def eta(XYZ,t):
    """ fonction qui renvoie l'altitude de la vague par rapport au niveau de mer plate 
    prend en entrée la position sous forme de vecteur dans la base cartésienne du labo et le temps"""
    return h0*np.cos(w*t-kvag(w)*XYZ[0])

def v_mer(XYZ,t):
    """ fonction qui renvoie la vitesse des particules d'eau en fonction de leur position XYZ et du temps"""
    k=kvag(w)
    vx=h0*w*np.exp(k*z)*np.cos(w*t-k*x)
    vz=-h0*w*np.exp(k*z)*np.sin(w*t-k*x)
    return [vx,0,vz]


def Pression(XYZ,t):
    """ fonction qui renvoie la valeur de la surpression due aux vagues au niveau de la 
    position XYZ dans le référentiel du labo et à l'instant t"""
    eta_vague=eta(XYZ,t)
    if XYZ[2] > eta_vague:
        """ on est au-dessus du niveau d'eau """
        return 0
    else:
        return rho_e*g*(eta_vague-XYZ[2])+rho_e/2*(XYZ[0]**2+XYZ[1]**2+XYZ[2]**2)**2
        
def Force_Moment(pos,t,theta,alpha,zOc):
    """ fonction qui calcule et renvoie la valeur de la force de pression (poussée d'archimède) suivant l'axe vertical vers le haut et le moment de la force de pression suivant l'axe ey horizontal perpendiculaire aux vagues au niveau du point M de position pos dans le référentiel du cylindre
    Attention a la translation suivant Oz de zOc"""

    # print("Force_Moment")
    
    XYZ=passag_cyl_xyz(pos,theta,alpha)
    # print(XYZ)
    XYZ[2]=XYZ[2]+zOc
    sP=Pression(XYZ,t)
    dangl=2*np.pi/Ntheta
    dL=L/NL
    dF=sP*R*dangl*dL*passag_cyl_xyz([1,0,0],theta,alpha)[2]
    # passag_cyl_xyz([1,0,0],theta,alpha)[2] est le produit scalaire er_cylindre * ez
    dM=pos[2]*(-1)*sP*R*dangl*dL*passag_cyl_xyz([1,0,0],theta,alpha)[1]
    return dF,dM

def Force_Moment_total(t,alpha,zOc):
    """ renvoie la force et le moment total sur le cylindre"""
    # print("Force_Moment_total")
    Ftot=0
    Mtot=0
    for i in range(NL):
        zc=(i+0.5)*L/NL
        for j in range(Ntheta):
            theta=(j+0.5)*2*np.pi/Ntheta
            pos=np.array([R,0,zc])
            dF,dM=Force_Moment(pos,t,theta,alpha,zOc)
            Ftot=Ftot+dF
            Mtot=Mtot+dM
    return Ftot,Mtot
    

            
            

    
### definition des parametres utiles
# pour les vagues
h0=1#m amplitude des vagues
f=0.5 #Hz  fréquence des vagues
w=2*np.pi*f #pulation des vagues
rho_e=1e3 #kg/m**3 masse volumique de l'eau
g=9.81 

#pour le pelamis
R=0.5 #m le rayon du pelamis
L=10 #m la longueur du pelamis
Ntheta=10 #nbre de points utilisés pour la discrétisation suivant theta
NL=10 #nbre de points utilisés pour la discrétisation suivant L
J=25 #moment d'inertie par rapport à l'axe Oy de rotation
m=100 #kg masse du pelamis
###
    
    

Nt=1e2
T=5
#initialisation des vecteurs de stockages
tt=np.zeros(int(Nt))
vz=np.zeros(int(Nt))
dalpha=np.zeros(int(Nt))
zOc=np.zeros(int(Nt))
alpha=np.zeros(int(Nt))

print("on débute")
for i in range(1,int(Nt)):
    
    dt=T/Nt
    tt[i]=tt[i-1]+dt
    Ft,Mt=Force_Moment_total(i*dt,alpha[i-1],zOc[i-1])
    
    vz[i]=vz[i-1]+dt*(Ft/m-g)
    dalpha[i]=dalpha[i-1]+Mt/J*dt
    zOc[i]=zOc[i-1]+vz[i]*dt
    alpha[i]=alpha[i-1]+dalpha[i]*dt







    
    
    
    