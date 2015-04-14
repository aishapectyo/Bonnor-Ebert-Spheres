#!/usr/bin/env python
from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np
from scipy import integrate


"""Info: Code solves lane-emden differential equation for following boundary conditions:
        psi(0) = 0; dpsi/dzi(zi=0) = 0. 

        The external pressure is determine from the following expression for the
        critical mass: 1.18*sigma^4/sqrt(G^3*P_o).
        A temperature of 20K is assume for the medium.
        The sound of speed is found from: sqrt(RT/mu); where R is the gas constant and mu is the
        molecular weight.

        After solving the ODE, the code will solve for dimensionless mass 'm'
        given an M and P_o. The value of rho_c is estimated. 
        Then, the code solves for zi_o for a given rho_c until the desired 
        TOTAL critical mass is found. This will give us the value of the BE boundary.   

By Aisha Mahmoud-Perez
Star Formation
"""

#----Functions-----#
def get_ext_pressure(T, m_crit):
        k = 1.38E-23
        p = 79.4 * T**4 * k / m_crit**2 #in SI units.
        p = p *10 #convert to bar.
        return p

def get_dmnless_mass(m_crit, p_out, c_sound):
        G = 6.67E-11
        p_out = p_out/10.0 #Pascals
        c_sound = c_sound / 100.0 #SI units
        m = (np.sqrt(p_out) * G**(1.5) * m_crit) / c_sound**4
        return m

def get_zi(m_dmnless,rho_out, rho_in, psi_diff):
        con = rho_out / (4 * np.pi * rho_in)
        zi_squared = m_dmnless / (con**(0.5) * psi_diff)
        return zi_squared

def get_mass(rho_in, c_sound, zi_sq, psi_diff):
        G = 6.67E-8 #cgs units
        con1 = 4 * np.pi * rho_in
        con2 = c_sound**2 / (4 * np.pi * G *rho_in)
        m = con1 * con2**(1.5) * zi_sq * psi_diff
        m = m / 1000.0 #change to kg
        m_s = m / 1.98E30 #chanege to solar masses
        return m_s
        
#----Constants----#
mu = 1.5 #any value between 1-2
T = 10 #temperature of gas in K
R = 8.31E7 # gas constant in cgs units
cs = np.sqrt(R*T/mu) #cgs units
P_o = get_ext_pressure(T, 1)
rho_o = P_o/(cs**2)

#----Solve Lane-Emden----#
y_init=[0, 0] #initial conditions
t = sp.linspace(0.0001, 8, 5000) #create array of values. Equidistant integration steps.
solution = sp.integrate.odeint(f_lane_emden, y_init, t)
psi = solution[:,0] #dump psi values here
dpsi = solution[:,1] #dump dpsi values here
rho_frac = sp.exp(-psi) #rho(r)/rho_o = exp(-psi)

#----Plot basic Lane-Emden----#
font = {'family' : 'ubuntu',
                        'size'   : 13}

plt.rc('font', **font)
plt.figure(1)
plt.plot(t, psi, color='LightSeaGreen', linewidth = 2, label='$\psi$')
plt.plot(t, rho_frac, color='Plum', linewidth = 2,label='$\\rho$/$\\rho_c$')
plt.xlabel('Nondimentsional radius $\\xi$')
plt.legend()

plt.rc('font', **font)
plt.figure(2)
plt.plot(-0.5, 0, color='b', linewidth = 2,label='log($\\rho$/$\\rho_c$)')
plt.xlim(-0.3,4.3)
plt.ylim(-1.2, 0)
plt.xlabel('Radius (AU)')
plt.legend(loc='lower left')


plt.show()

#----Find zi_o----#
rho_c = np.linspace(1.5E-20, 2.9E-18, 2000) #find the array of values of rho_c. Guesstimate!
rho_c_len = len(rho_c)
dpsi_len = len(dpsi)
# Loop trough each value of rho_c and find a value of zi^2. Use that zi^2 to find
# the total enclosed mass. If the mass desired is reached, save the
# rho_c and zi_sq (and their indexes)  and exit the loop.  
# Convert the zi_sq into a physical radius.
dimensionless_mass = get_dmnless_mass(0.5*1.98E30, P_o, cs) #do not use solar units
for i in xrange(rho_c_len):
        final_mass_computed = 0.48
        for n in xrange(dpsi_len):
                dimensionless_radius_squared = get_zi(dimensionless_mass,rho_o, rho_c[i], dpsi[n])
                if dimensionless_radius_squared < 6.4:
                        mass_cloud = get_mass(rho_c[i], cs, dimensionless_radius_squared, dpsi[n])
                        if mass_cloud > 0.4999 and mass_cloud < 0.5001:
                                if mass_cloud > final_mass_computed:
                                        final_mass_computed=mass_cloud
                                        total_mass = mass_cloud
                                        index_rho_c = i
                                        index_dpsi = n
                                        z_real = np.sqrt(dimensionless_radius_squared)
                                        rho_c_real = rho_c[i]

#----Get physical values for plotting----#
radius_au = physical_r(z_real, rho_c_real, cs)
radius_cutoff = []
rho_frac_cutoff = []
rho_c_m = rho_c_real*1E-4 #convert to kg/m^3
for m in xrange(len(t)):
        radius_list = physical_r(t, rho_c_real, cs)
        if radius_list[m] < radius_au:
                rho_frac_cutoff.append(rho_c_m*rho_frac[m])
                radius_cutoff.append(radius_list[m])


#----Plot results----#
rho_frac_cutoff = np.array(rho_frac_cutoff)
radius_cutoff = np.array(radius_cutoff)
font = {'family' : 'ubuntu',
                                        'size'   : 13}

plt.rc('font', **font)
plt.figure(3)
plt.plot(radius_cutoff, rho_frac_cutoff, color='b', linewidth = 2,label='$\\rho(r)$')
plt.xlabel('Radius (AU)')
plt.ylabel('$\\rho$(r) $kg/m^3$')
plt.legend()

plt.rc('font', **font)
plt.figure(4)
plt.plot(np.log10(radius_cutoff), np.log10(rho_frac_cutoff/rho_c_m), color='b', linewidth = 2,label='log($\\rho$/$\\rho_c$)')
plt.xlabel('Radius (AU)')
plt.legend(loc='lower left')
plt.show()

