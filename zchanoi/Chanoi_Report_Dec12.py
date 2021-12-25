#!/usr/bin/env python
# coding: utf-8

# <div align="center">
# 
# ## The Need for Multi-Material Selective Laser Sintering Computational Simulations using Finite Volume Method - Review and Prelminary Calculations <br>
# 
# <div align="center">
#     
# ### Zachary Aidan Chanoi 
# <div align="center">
# Department of Mechanical Engineering, The University of Texas El Paso <br><br><br>
# 
# <div style="text-align: justify"> 
# 
# ### Abstract  
# Selective Laser Sintering is an Additive Manufacturing technique that sinters material powders together using a laser in a heated enclosure. Typically, SLS printing of thermoplastic polymers is done with one material. Technology advances have introduced the possibility of multi-material (M2) SLS printing involving more than one different thermoplastics being sintered in one print. However, limited understanding of the physcial processes behind M2 sintering keeps the technology from being practicle. This beckons a better understanding of the thermal mechanisms behind M2 SLS printing. This report presents reviews previous attempts to model SLS and introduces the preliminary calculations and modeling of a 2-dimensional (2D) M2 SLS print using finite volume method (FVM).  
# 
# ### Introduction 
# 
# Selective Laser Sintering (SLS) is an Additive Manufacturing (AM) technique that sinters polymer or metal powders together using a laser in a heated enclosure to produce functional parts. Typically, SLS printing of thermoplastic polymers is done with one material. Technology advances have introduced the possibility of multi-material (M2) SLS printing involving more than one thermoplastic to be sintered in one part. M2 AM is an attractive domain of AM because of the potential to create highly-customizable mechanically-tailored structures for many industries.  However, differences in thermal properties between materals result in frequent print errors which create weakness in printed structures, poor dimensional stability, and aborted prints (Figure 1). This beckons a better understanding of the thermal mechanisms of the M2 SLS printing process. This report presents the preliminary calculations and modeling of a 2-dimensional (2D) M2 SLS print using finite volume method (FVM). 
#     
#   
# <img src="FailedM2prints.png" alt="Drawing" style="width: 500px;"/>
# <div align="center">
# <b> Figure 1 - Preliminary M2 SLS Prints with Poor Dimensional Stability </b></figcaption>   
#     
# <div style="text-align: justify">
#     
# ### Background 
# 
# By taking a computer-aided design (CAD) and slicing it into layers of data, AM can turn the CAD design into a physcial part. AM is not held back by the same limitations of subtractive manufacturing such as metal machining. This enables the capability to produce complex geometries previously unachievable. Today, AM has applications in many industries and is a rapidly growing manufacturing industry. In biomedical engineering, AM is especially promising for prosthetics where parts manufactured specifically for a patient could mean preventing a lifetime of debilitating pain. 
# 
# Recent research and developent of AM has focused on introducing multi-material (M2) capabilities into existing AM technologies [2, 3]. Currently, most AM parts are made up of only a single material. The inherent layer-by-layer nature of AM makes it so switching materials mid-layer is challenging. Some AM technologies have the ability to stop prints to swap materials. However, because AM is a heat-sensitive process, disrupting prints can lead to catastropic print failures which wastes material and time (Figure 2).
# 
# <img src="CommonSLSprintfails.png" alt="Drawing" style="width: 500px;"/>
# <div align="center">
# <b> Figure 2 - Common SLS Printing Failures </b></figcaption>   
#     
# <div style="text-align: justify">
# 
# 
# Cheng et al. (2004) [4] modeled the SLS of mixed metal powders using FVM to study the effects that laser intensity, scanning veloicty, and the amount of previously sintered layers have on the sintering process. The SLS process was described as a "three-dimensional transient melting and resolidification" process, the model was reduced to a two-dimensional, two component metal powder simulation with an infinite horizontal domain. The laser was described to be Gaussian that only interacts with the top surface of the loose powder. The sintering depth was set to just exceed the interface between the loose powder layer and the next sintered layer underneath. A diagram of their physical model can be seen in Figure 3. 
# 
# <p align="center">
# <img src="Cheng2004.png" alt="Drawing" style="width: 400px;"/>
# <div align="center">
# <b> Figure 3 - Physical Model of the SLS Printing Process [Cheng et al. (2004)] </b></figcaption>   
#     
# <div style="text-align: justify">
# 
#     
# A temperature-transforming model was used and the dimensinless governening equation used and discretized by FVM was:
# 
# $$ -U\frac{\partial (CT)}{\partial X} + W \frac{\partial (CT)}{\partial Z} = \frac{\partial}{\partial X}(K\frac{\partial T}{\partial Z}) - (-U \frac{\partial S}{\partial X} + W \frac{\partial S}{\partial Z}) \tag 1  $$
# 
# Where $ W $ is the shrinkage velocity of the powder, $ C $ is the heat capacity, $ S $ is the source term, and $ K $ is the thermal conductivity. All these terms were made dependent on their location and were defined by piecewise functions.
# 
# A simulation was done with 0, 1, 3, 5, 10, and 50 previously sintered layers below. As expected, it was found that the sintering depth increased substantially with increasing laser intensity, and decreases with increasing scanning velocity. It was also determined that the size and shape of the melting pool of powder is also dependent on the laser intensity and velocity. 
#     
# Wang et al. (2020) [5] highlights the urgent need for a computationally efficient, yet accurate, method for the SLS process which cannot be achieved using traditional numberical simulation techniques. They developed a 3D FVM model to descrive the printing of Ti6Al4V powder. They used the follow equation to govern the SLS process: 
# 
# $$ Q(x,y,z,t) - \rho C_{P}\frac{dT}{dt} + \nabla \cdot \kappa \nabla T = 0 \tag 2 $$
# 
# Where $ Q(x,y,z,t) $ is the body heat source of laser beam, the spatial coordinates, and time, $ T $ is the temperature, $ \rho $, $ C_{P} $, and $ \kappa $ are the temperature-dependent physical parameters for density, specific heat, and thermal conductivity, respectively. The boundary heat conditions were definied by convection and radiation:
# 
# $$ q_{s}(x,t) = h(T - T_{\infty}) + \varepsilon \sigma(T^{4} - T^{4}_{\infty}) \tag 3 $$ 
#     
# Where $ T $ is the temperature, $ T_{\infty} $ is the ambient temperature, $ h $ is the convection coefficient $ \varepsilon $ is the emissivity, and $ \sigma $ is the Stefan-Boltzmann constant.
#     
# Wang et al. utalized 3D FVM to discretize the governing partial diffential equation into algebraic equatiosns:  
# 
# $$ a_{P}T_{P} = a_{W}T_{W} + a_{E}T_{E} + a_{S}T_{S} + a_{N}T_{N} + a_{B}T_{B} + a_{T}T_{T} + a^{0}_{P}T^{0}_{P} + S_{u} \tag 4 $$
#  
# $$ a_{P} = a_{W} + a_{E} + a_{S} + a_{N} + a_{B} + a_{T} + a^{0}_{P} - a_{S} \tag 5 $$
# 
# $$ a^{0}_{P} = \frac{\rho C_{P} \Delta A}{\Delta t} \tag 6 $$
# 
# $$ a_{i} = \frac{k_{i}A_{i}}{d_{i}} \tag 7 $$
#     
# The computational efficiency of the FVM was compared with an finite element method (FEM) model. For the same simulations, the FVM model was about 4.5 times faster than the FEM model - taking 174 minutes, 679 minutes, and 1,630 minutes, compared to 767 minutes, 3,041 minutes, and 7,045 minutes, respectively, with similar results. Thus, the superiority of FVM over FEM for SLS modeling was demonstrated. 
#     
# Mokrane et al. [6] also used FVM to model SLS. Previously, Mokran used Discrete Element Method to model laser-powder interactions, however, it was decided it was too computationally intense to be implemented in a larger scale simualtion. FVM was chosen because of it's simplicity to consider powder volume shrinkage and deposition of fresh powder layers. The powder leyer was considered to be a homogeneous medium without any flow due to high polymer viscosity, and the governing equation used was: 
#     
# $$ \rho C_{P}\frac{dT}{dt} = + \nabla \cdot (\gamma \nabla T) + Q + S_{f} + S_{c} \tag 8 $$    
#     
# Where $ Q $ is the volumetric heat source, $ T $ is the temperature, $ \rho $, $ C_{P} $, and $ \gamma $ are the density, specific heat, and thermal conductivity, of the powder bed, S_{f} is the melting latent heat, and S{c} is the cooling crystallization enthalpy.   
#     
# The SLS of PA12 was investigated. A bird's eye view of the the temperature distribution following a laser scan can be seen in Figure 4. A plot of the thermal history for the scanning, cooling, new layer application, and additional laser scan can be seen in Figure 5. Finally, a cross-section of the printing bed shows the temperature distribution after five layers have been sintered can be seen in Figure 6. 
# 
# <p align="center">
# <img src="Mokrane.png" alt="Drawing" style="width: 300px;"/>
# <div align="center">
# <b> Figure 4 - Bird's eye view of the the temperature distribution [Mokrane et al. (2018)] </b></figcaption>   
#     
# <p align="center">
# <img src="Mokrane2.png" alt="Drawing" style="width: 300px;"/>
# <div align="center">
# <b> Figure 5 - Thermal History of the Scanning, Cooling, New Layer Application, and Additional Laser Scan [Mokrane et al. (2018)] </b></figcaption>     
#     
# <p align="center">
# <img src="Mokrane3.png" alt="Drawing" style="width: 300px;"/>
# <div align="center">
# <b> Figure 6 - Cross-section of the Printing Bed Showing the Temperature Distribution after Five Layers have been Sintered     [Mokrane et al. (2018)] </b></figcaption>   
# <div style="text-align: justify">    
#  
#     
# ### Problem Statment and Methodology 
#     
# The heat transfer of a homogenous bed of PA12 powder is investigated using FVM. The governing equation used to model the change in temperature with respect to time is the same as the one used in Wang (2020) [5] only it has been simplified to two-dimensions: 
#     
# $$ Q - \rho C_{P}\frac{dT}{dt} + \nabla \cdot \kappa \nabla T = 0 $$
#     
# Which is reduced to: 
#     
# 
# $$ \dfrac{\partial}{\partial x}(\kappa \dfrac{\partial T}{\partial x}) + \dfrac{\partial}{\partial y}(\kappa \dfrac{\partial T}{\partial y}) + S_{\phi} = 0 \tag 9 $$ 
# 
# Where:
#     
# $$ S_{\phi} = q - \rho C_{P} \frac{dT}{dt}  $$ 
#     
# The equation is discreteized into a two-dimensional finite volume domain: the central node P is sourrounded and bounded by the N, S, E, W, nodes (Figure 7):
#     
# $$ \biggl(\kappa_e A_e \biggl(\dfrac{\partial T}{\partial x}\biggr)_e-\kappa_w A_w \biggl(\dfrac{\partial T}{\partial x}\biggr)_w \biggr)+\biggl(\kappa_n A_n \biggl(\dfrac{\partial T}{\partial y}\biggr)_n-\kappa_s A_s \biggl(\dfrac{\partial T}{\partial y}\biggr)_s\biggl)+ S_{\phi}\Delta V = 0 $$    
#   
# <p align="center">     
# <img src="Grid.png" alt="Drawing" style="width: 250px;"/>
# <div align="center">
# <b> Figure 7 - Standard FVM Grid with N, S, E, W Nodes </b></figcaption>  > 
#     
# <div style="text-align: justify">     
#     
# Where $ Q $ is the heat from the laser, $ T $ is the temperature, $ \rho $, $ C_{P} $, and $ \kappa $ are the density, specific heat, and thermal conductivity, respectively, of the PA12 powder.  
# 
# This report makes the following assumptions:
# 1. The PA12 powder is homogeneous 
# 2. The $ \rho $, $C_{P} $, and $ \kappa $, of the PA12 powder is constant
# 3. There is no heat flux across any boundary except boundary N - where the heat flux from the laser is considered 
# 4. The temperature of each boundary starts at 163$^{\circ}$C
#     
# Figure 8 shows a visual interpitation of the system being modeled:
#  
# <p align="center">     
# <img src="PD.png" alt="Drawing" style="width: 200px;"/>
# <div align="center">
# <b> Figure 8 - Visual Problem Definiton of SLS Printing of PA12 Powder </b></figcaption>  > 
#     
# <div style="text-align: justify">      
# 
# ### Results and Discussion 
# 
# The resulting temperture distribution for a 5W laser that is heating the PA12 for 30, 60, 120 seconds is seen in Figure 9, Figure 10, and Figure 11, respectivly. PA12 sinters at 200$^{\circ}$C, thus, 30 seconds of energy from the 5W laser is not sufficient to sinter the top layer of powder. As expected, the max temperature for each test was in the center of the north boundary where the heat flux from the laser was applied. The temperature distribution is in the shape of a U, with higher temperatures penetrating into the powder than laterally across. This could be a significant issue in the context of SLS printing because previosuly layers of powder would be right below the top layer being sintered. Thus, layers that had just been sintered and then covered may burn once a fresh layer is sintered on top of it. 
#     
# The results of this FVM of the SLS process is far from accurately representing the SLS process in the physical-world. The density of PA12, its specific heat, and its thermal conductivity, would not be constant because of changes in temperature and phase changes. Additionally, the laser's position would be dynamic and follow a predetermined path. Introducing a way to model fresh (cool) powder on top of a previously sintered layer would be desired such as what was seen in Mokrane et al. 
#     
# <p align="center">
# <img src="R30.png" alt="Drawing" style="width: 100;"/>
# <div align="center">
# <b> Figure 9 - Thermal Distribution of PA12 Powder after 30 Seconds of 5W Laser </b></figcaption>  > 
#     
# <p align="center">
# <img src="R60.png" alt="Drawing" style="width: 100;"/>
# <div align="center">
# <b> Figure 10 - Thermal Distribution of PA12 Powder after 60 Seconds of 5W Laser </b></figcaption>  > 
#     
# <p align="center">
# <img src="R120.png" alt="Drawing" style="width: 100;"/>
# <div align="center">
# <b> Figure 11 - Thermal Distribution of PA12 Powder after 120 Seconds of 5W Laser </b></figcaption>  > 
#     
# <div style="text-align: justify">     
#        
# ### Future Work 
# 
# Correcting the limitations discussed in the Results must be done for a better simulation of SLS printing. Furthermore, to begin studying multi-material printing, a second material must be introduced into the simulation. This can be done in many ways depending on the multi-material printing attempt. For example, more elastic thermoplastics, such as TPE, have been used as a base for PA12 to help minimize thermal warping (Figure 12). A piecewise function that defines what material is where spatially should be introduced. If multiple powders are desired to be mixed together, a rule of mixtures approach could be utalized to define new thermal characteristics which could then be introduced into the simulation. 
#     
# <p align="center">
# <img src="Problem%20Definition%20Visual.png" alt="Drawing" style="width: 100px;"/>
# <div align="center">
# <b> Figure 12 - Concept Problem Definition of SLS Print with TPE Base to Minimize Warpping </b></figcaption>  > 
# <div style="text-align: justify"> 
#     
# ### Conclusion 
# 
# This report discussed the need for sophisticated simulations of multi-material SLS printing. M2 AM is a rapidly growing field with applications in many industries. However, a lack of understanding of the thermal kinetics of M2 SLS leads to frequent and catastrophic print errors which prevent M2 printing to be consistent enough for practical use. Computationally modeling single material SLS printing has been attempted many times. Because of improved computational time and accuracy, FVM has become the popular way to model SLS. This report showed a very simplified attempt to model the heat transfer through a bed of PA12 powder using FVM. Future work was discussed that presents what must be done to accuratly model M2 SLS using FVM.  
#     
# ### References 
# 
# [1] Kumar., V. Computational Fluid Dynamics. The University of Texas at El Paso. (2021) 
#     
# [2] Green, J., Slager, J., Lopez, M., Chanoi, Z., Stewart, C., Gonzalez, G. Local control of mechanical properties through multi-material blending in fused filament fabrication with an active-mixing hotend. Additive Manufacturing. (2021). 
#     
# [3] Stichel, T., Laumer, T., Raths, M., & Roth, S. (2018). Multi-material deposition of polymer powders with vibrating nozzles for a new approach of laser sintering. Journal of Laser Micro Nanoengineering. 
#     
# [4] Chen, T., & Zhang, Y. (2004). Numerical simulation of two-dimensional melting and resolidification of a two-component metal powder layer in selective laser sintering process. Numerical Heat Transfer. 
# 
# [5] Wang, J., Wang, Y., Shi, J. On efficiency and effectiveness of finite volume method for thermal analysis of selective laser melting (2020). \
# 
# [6] Mokrane, A., Boutaous, M., & Xin, S. (2018). Process of selective laser sintering of polymer powders: Modeling, simulation, and validation.
#     
# [7] Zhao, M., Drummer, D., Wudy, K., & Drexler, M. (2015). Sintering Study of Polyamide 12 Particles for Selective Laser Melting. International Journal of Recent Contributions from Engineering. 
#     
# [8] Introduction to Computational Fluid Dynamics. Publisher. (Author)    
# 

# In[90]:


import matplotlib.pyplot as plt
import numpy as np
import json

Lx=.110; Ly=.09; n = 50; m = 50; nn=n*m; # Geom descretization/mesh
Ta=162; Tb=162; Tc=162; Td=None; # BCs
k=0.3; Cp = 3; rho = 930;  # Properties
LP = 5; Lt = 120;

dx = Lx/n;  dy = Ly/m;

x = np.linspace(dx/2,Lx-dx/2,n); y = np.linspace(dy/2,Ly-dy/2,m); [X, Y]=np.meshgrid(x,y);

A=np.zeros([nn,nn]); b=np.zeros([nn]); d2=np.zeros([m,n])

dz=0.01;  G=k; qa=0; qb=0; qc=0; qd=LP*Lt;

for j in range(m):
     for i in range(n):
        P = j*n+i; W = P-1; E = P+1; N=P+n; S=P-n;
        Aw = dy*dz; Ae=Aw; An=dx*dz; As = An; Su=0; Sp=0; Sua=qa*Aw; Sub=qb*Ae; Suc=qc*As; Sud=qd*An;
        aW = G*Aw/dx; aE = G*Ae/dx; aN = G*An/dy; aS = G*As/dy;
        if(i>0):
            A[P,W]=-aW;
        else:
            if Ta:
                aW=0; Sp=Sp- rho*Cp*Aw/dx; Su=Su + rho*Cp*Aw/dx*Ta;  
            else:
                aW=0; Sp=Sp; Su=Sua + Su;
     
        if(i<n-1): 
            A[P,E]=-aE; 
        else:
            if Tb:
                aE=0; Sp=Sp - rho*Cp*Ae/dx; Su=Su + rho*Cp*Ae/dx*Tb;
            else:
                aE=0; Sp=Sp; Su=Sub + Su;
        if(j>0): 
            A[P,S]=-aS; 
        else:
            if Tc:
                aS=0; Sp=Sp- rho*Cp*As/dy; Su=Su + rho*Cp*As/dy*Tc;
            else:
                aS=0; Sp=Sp; Su=Suc + Su;           
        if(j<m-1): 
            A[P,N]=-aN; 
        else:
            if Td:
                aN=0; Sp=Sp- rho*Cp*An/dy; Su=Su + rho*Cp*An/dy*Td;
            else:
                aN=0; Sp=Sp; Su=Sud + Su;           
        aP = aW + aE + aS + aN - Sp; 
        A[P,P] = aP; 
        #print(A[P,P])
        #print(aP)
        b[P]=Su; 

d=np.linalg.solve(A,b);  

for j in range(m):
    for i in range(n):
        IN=j*n+i; d2[j,i]=d[IN];
        
fig = plt.figure()
plt.contourf(X, Y, d2)
plt.colorbar()
plt.show()

print("Max Temp =",max(d))
print(d)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




