#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import matplotlib
import pylab
import numpy as np
import pandas as pd
from random import seed
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm


# In[2]:


def rotation(x_comp, y_comp, z_comp, axis, ang):
    item = np.array([x_comp,y_comp,z_comp])
    if axis == 'x':
        Rx = np.matrix([[1.,0.,0.], [0.,np.cos(ang),-1*np.sin(ang)], [0.,np.sin(ang),np.cos(ang)]])
        return Rx.dot(item)
    elif axis == 'y':
        Ry = np.matrix([[np.cos(ang),0.,np.sin(ang)], [0.,1.,0.], [-1.*np.sin(ang),0,np.cos(ang)]])
        return Ry.dot(item)
    else:
        Rz = np.matrix([[np.cos(ang),-1*np.sin(ang),0.], [np.sin(ang),np.cos(ang),0.], [0.,0.,1.]])
        return Rz.dot(item)


# In[3]:


def generator(lamda, n, k0, mom_jet, dist):
    #np.random.seed(8)#for test with same random data
    ngen = np.random.poisson(lamda)#number of hadrons
    print(ngen,"particles")
    theta = np.pi*2.*np.random.uniform(size=ngen)#angle
    z_comp = mom_jet
    r = np.random.uniform(size=ngen)
    x = 1-np.power(r,1/(n+1))#momentum fraction without normalized
    #print("x: ", x)
    kt_1 = np.zeros(ngen)
    
    for i in range(ngen):
        kt_1[i]=20.0
        while(kt_1[i]>10.0):
            r = np.random.uniform()
            kt_1[i]=k0*np.tan(r*np.pi/2)

    
    x_sum = np.sum(x)
    x = x/x_sum #normalized momentum fraction
    z_comp = z_comp*x

    #print("trans: ", kt_1)
    #print("theta: ", np.arctan2(kt_1,x*300.), "\n")
    phi = np.pi*2*np.random.uniform(size=ngen)
    x_comp = kt_1*np.cos(phi)
    y_comp = kt_1*np.sin(phi)
    sum_x = np.sum(x_comp)
    sum_y = np.sum(y_comp)

    #normalization
    x_comp = x_comp-(sum_x*x)
    y_comp = y_comp-(sum_y*x)
    
    #test
    #print(np.sum(x_comp))
    #print(np.sum(y_comp))
    
    #theta_x = np.arctan(x_comp/z_comp)
    #theta_y = np.arctan(y_comp/z_comp)
    #print(theta_x, "\n", theta_y)
    #x_position = dist*np.tan(theta_x)
    #y_position = dist*np.tan(theta_y)
    
    #print("gen-test")
    
    return z_comp, y_comp, x_comp


# In[18]:


def _3dplot(x_, y_, z_):
    x_array = x_#theta
    y_array = y_#phi
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    hist, xedges, yedges = np.histogram2d(x_array, y_array, bins=(48,48), range=[[0.,np.pi],[-1*np.pi,+np.pi]], weights=z_)
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1])
    #print(hist,"\n","\n")
    #print(xedges,"\n","\n")
    #xpos = xpos.T
    #ypos = ypos.T
    xpos = xpos.flatten("C")
    ypos = ypos.flatten("C")
    zpos = np.zeros_like (xpos)
    print(xpos,"\n","\n",ypos,"\n","\n")
    #ypos = ypos.T
    #print(ypos,"\n","\n")
    hist = hist.T
    #based on size of read out number of bins should be fixed in the plot
    dx = 0.05*np.pi
    dy = 0.05*np.pi*2.
    dz = hist.flatten()
    print(dz)

    cmap = cm.get_cmap('jet') # Get desired colormap - you can change this!
    max_height = np.max(dz)   # get range of colorbars so we can normalize
    min_height = np.min(dz)
    # scale each z to [0,1], and get their rgb values
    rgba = [cmap((k-min_height)/max_height) for k in dz] 
    
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')
    #print(x_array)
    plt.title("Momentum Magetitude")
    plt.xlabel("Theta Angle")
    plt.ylabel("Phi Angle")
    
    plt.show()
    return 0


# In[19]:


def _2dplot(x_, y_, z_):
    plt.hist2d(x_, y_, bins=(20, 20), weights=z_, range=([0.,np.pi],[-1*np.pi,+np.pi]),cmap=plt.cm.Reds)
    plt.xlabel("Theta")
    plt.ylabel("Phi")
    plt.show()
    return 0


# In[20]:


def main():
    loop = 100#number of Jets get collision
    lamda = 7#average number of Hadrons in one Jet
    n = 3.0#quark number/2.0 for poin/3.0 for proton  physics parameter
    k0 = 1.0#momentum p0 = 1.0 GeV
    mom_jet = 300.#unit GeV
    dist = 4.72#unit meter
    ang = np.pi/3.#rotate angle
    axis = 'x'
    theta = np.array([0])#[0,pi]
    phi = np.array([0])#[-pi,pi]
    mag_mom = np.array([0])
    
    for index in range(loop):
        #generate jets
        z_comp, y_comp, x_comp = generator(lamda, n, k0, mom_jet, dist)
        #print(len(x_comp), len(y_comp), len(z_comp))
        #print(x_comp)#test
        
        #rotate jets
        item = rotation(x_comp, y_comp, z_comp, axis, ang)
        x_new = item[0].flatten()
        y_new = item[1].flatten()
        z_new = item[2].flatten()
        #print(len(x_new[0]), len(y_new), len(z_new))
        #print(np.sum(x_new),'\n',np.sum(y_new),'\n',np.sum(z_new))#test
        
        mag_mom = np.append(mag_mom, np.sqrt(np.multiply(x_new,x_new)+np.multiply(y_new,y_new)+np.multiply(z_new,z_new)))
        phi_ele = np.arctan2(y_new, x_new)
        theta_ele = np.arctan2(np.sqrt(np.multiply(x_new,x_new) + np.multiply(y_new,y_new)), z_new)
        phi = np.append(phi, phi_ele)
        theta = np.append(theta, theta_ele)
        
        
    phi = phi.flatten('C')
    theta = theta.flatten('C')
    mag_mom = mag_mom.flatten('C')
    phi = np.delete(phi, 0)
    theta = np.delete(theta, 0)   
    mag_mom = np.delete(mag_mom, 0)
    #print("Theta: ", theta)
    #print("Phi: ", phi)
    #theta = theta/np.pi
    #phi = phi/np.pi
    #_2dplot(theta, phi, mag_mom)
    _3dplot(theta, phi, mag_mom)
    #values(x_, y_, z_)
    
    return 0


# In[21]:


main()


# In[ ]:





# In[ ]:




