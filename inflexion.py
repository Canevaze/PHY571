# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:11:05 2023

@author: Léo
"""
import matplotlib.pyplot as plt
import numpy as np

def trouver_point_inflexion(x, y):
    # Calculer la dérivée seconde
    dy2 = np.gradient(np.gradient(y))

    # Trouver les indices où la dérivée seconde change de signe
    indices_inflexion = np.where(np.diff(np.sign(dy2)))[0]

    # Sélectionner le point d'inflexion comme le point médian entre ces indices
    if len(indices_inflexion) > 0:
        point_inflexion_index = indices_inflexion[len(indices_inflexion) // 2]
        point_inflexion = (x[point_inflexion_index], y[point_inflexion_index])
        return point_inflexion
    else:
        return None
v=np.load('v_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')    
eta=np.load('eta_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')    

v2=np.load('v_it_200_N_40_V_0.03_L_3.1_R_1_avg_50_nb_pts_50.npy')
eta2=np.load('eta_it_200_N_40_V_0.03_L_3.1_R_1_avg_50_nb_pts_50.npy')

v3=np.load('v_it_200_N_400_V_0.03_L_10_R_1_avg_3_nb_pts_50.npy')
eta3=np.load('eta_it_200_N_400_V_0.03_L_10_R_1_avg_3_nb_pts_50.npy')

plt.figure(figsize=(6, 6))
plt.scatter(eta2,v2,color='r',label='N=40',marker='s',edgecolors='black')
plt.scatter(eta,v,color='b',label='N=100', marker ='x')
plt.scatter(eta3,v3,color='g',label='N=400', marker ='+')
plt.axis([0,5.0,0,1.03])
plt.legend()
plt.grid()
plt.xlabel("$\eta$", fontsize=20)
ylabel=plt.ylabel("$v_{a}$", fontsize=20)
ylabel.set_rotation(0)
plt.gca().yaxis.set_label_coords(-0.15, 0.45)
title=plt.title("Evolution of the order parameter with respect to $\eta$", fontsize = 20)
title.set_position([0.5, 1.2])
# print("Point inflexion pour N=100 :",trouver_point_inflexion(eta, v))
# print("Point inflexion pour N=40 :",trouver_point_inflexion(eta2, v2))