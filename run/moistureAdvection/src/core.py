import numpy as np
import math
import sys
import os
import time

#Rotate x and z by theta.
def rotate(x0,z0,theta):
    x = x0*np.cos(theta) - z0*np.sin(theta)
    z = x0*np.sin(theta) + z0*np.cos(theta)
    return x,z
    
#Find converging amplitude for damping source term.
def finalAmplitude(L,w):
    T = 2*np.pi/w 
    frac = T/L
    A = (1. - np.exp(-0.5*frac))/(1. - np.exp(-frac))
    return A
    
#Find converging amplitude for velocity source term.
def finalAmplitude2(L,r):
    frac = 4.*r/L
    A = (1. - np.exp(-0.5*frac))/(1. - np.exp(-frac))
    return A

#Produce Schaer profile about dist_center_x,z with certain half_width.
def schaerProfile(x,z,dist_center_x,dist_center_z,half_width):
    x = x - dist_center_x
    z = z - dist_center_z 
    r = np.sqrt(x*x + z*z)*1./half_width
    if r < 1:
        return np.cos(r*np.pi/2)**2
    else:
        return 0.

#Produce analytic solution of schaer profile after time t.
#Source terms are damping and are based on central schaer profile point only.        
def schaerProfileDamped(x,z,dist_center_x,dist_center_z,half_width,theta,w,t,L):
    #Find "orginal" coord by rotating back to initial non-rotated condition.
    x0,z0 = rotate(x,z,-theta)
    
    T = 2*np.pi/w
    #How many rotations have occured?
    cycle_phase = (t/T)%T
    #How many half-rotations have occured?
    cycle_phase2 = (2*t/T)%T
    
    #If point within the Schaer radius, calculate the cells value.
    x = x - dist_center_x
    z = z - dist_center_z 
    r = np.sqrt(x*x + z*z)*1./half_width
    if r < 1:
        expon = np.exp(-T/(2.*L))
        factor = 1.
        #If in first half of cycle.
        if int(cycle_phase2)%2 == 0:
            for i in range(int(cycle_phase)):
                factor *= expon
                factor = 1 - (1-factor)*expon
            factor *= np.exp(-(t%T)/(1.*L))
        #If in second half of cycle.
        else:
            for i in range(int(cycle_phase)+1):
                factor *= expon
                if i != int(cycle_phase):
                    factor = 1 - (1-factor)*expon
            factor = 1 - (1-factor)*np.exp(-(t%(T/2))/(1.*L))
            
        return factor*np.cos(r*np.pi/2)**2
    else:
        return 0.
        
#Produce analytic solution of schaer profile after time t.
#Source terms are damping and are accurate for each data point (not dependent on central point).
def schaerProfileDamped2(x,z,dist_center_x,dist_center_z,half_width,theta,w,t,L):
    #Find "orginal" coord by rotating back to initial non-rotated condition.
    x0,z0 = rotate(x,z,-theta)
    
    #If point within the Schaer radius, calculate the cells value.
    x = x - dist_center_x
    z = z - dist_center_z 
    r = np.sqrt(x*x + z*z)*1./half_width
    if r < 1:
        #Angle perturbation relative to central Schaer profile point.
        theta_dash = math.atan2(z0,x0) - np.pi/2
        T = 2*np.pi/w
        expon = np.exp(-T/(2.*L))
        factor = 1. 

        #How many rotations and half-rotations have occured?
        cycle_phase = ((t+theta_dash/w)/T)%T
        cycle_phase2 = (2*(t+theta_dash/w)/T)%T
        
        #Point is initially to left of central point.
        if theta_dash >= 0.:

            #If in first half of cycle.
            if int(cycle_phase2)%2 == 0:
                
                for i in range(int(cycle_phase)):
                    #Special case for first half-cycle.
                    if i == 0:
                        factor *= np.exp(-(0.5*T-theta_dash/w)/L)
                    else:
                        factor *= expon
                    factor = 1 - (1-factor)*expon
                #Special case for first half-cycle.
                if int(cycle_phase) == 0:
                    factor *= np.exp(-(t%T)/(1.*L))
                else:
                    factor *= np.exp(-((t+theta_dash/w)%T)/(1.*L))
            #If in second half of cycle.
            else:
                
                for i in range(int(cycle_phase)+1):
                    #Special case for first half-cycle
                    if i == 0:
                        factor *= np.exp(-(0.5*T-theta_dash/w)/L)
                    else:
                        factor *= expon
                    if i != int(cycle_phase):
                        factor = 1 - (1-factor)*expon
                factor = 1 - (1-factor)*np.exp(-((t+theta_dash/w)%(T/2))/(1.*L))
        #Point is initially to right of central point.
        else:
            #As rho_2 initially 0, no source term until point has crossed boundary.
            if t < abs(theta_dash)/w:
                pass
            #First half of cycle.
            elif int(cycle_phase2)%2 == 0:

                for i in range(int(cycle_phase)):
                    factor *= expon
                    factor = 1 - (1-factor)*expon
                factor *= np.exp(-((t+theta_dash/w)%T)/(1.*L))
            #Second half of cycle.
            else:
                
                for i in range(int(cycle_phase)+1):
                    factor *= expon
                    if i != int(cycle_phase):
                        factor = 1 - (1-factor)*expon
                factor = 1 - (1-factor)*np.exp(-((t+theta_dash/w)%(T/2))/(1.*L))
            
        return factor*np.cos(r*np.pi/2)**2
    else:
        return 0.

#Produce analytic solution of schaer profile after time t.
#Source terms are velocity-dependent 
#and are accurate for each data point (not dependent on central point).        
def schaerProfileVelocity(x,z,dist_center_x,dist_center_z,half_width,theta,w,t,L):
    #Find "orginal" coord by rotating back to initial non-rotated condition.
    x0,z0 = rotate(x,z,-theta)

    #If point within the Schaer radius, calculate the cells value.    
    r = np.sqrt(x*x + z*z)
    x = x - dist_center_x
    z = z - dist_center_z 
    r0 = np.sqrt(x*x + z*z)*1./half_width
    if r0 < 1:
        #Angle perturbation relative to central Schaer profile point.
        theta_dash = math.atan2(z0,x0) - np.pi/2
        T = 2*np.pi/w
        expon = np.exp(-2*r/L)
        factor = 1.

        #How many rotations and half-rotations have occured?
        cycle_phase = ((t+theta_dash/w)/T)%T
        cycle_phase2 = (2*(t+theta_dash/w)/T)%T
        
        #Point is initially to left of central point.
        if theta_dash >= 0.:
            #If in first half of cycle.
            if int(cycle_phase2)%2 == 0:
                for i in range(int(cycle_phase)):
                    #Special case for first half-cycle.
                    if i == 0:
                        factor *= np.exp( (-1. - np.sin(theta_dash + np.pi/2.)) * r/L )
                    else:
                        factor *= expon
                    factor = 1 - (1-factor)*expon
                #Special case for first half-cycle.
                if int(cycle_phase) == 0:
                    factor *= np.exp( ( np.sin( w*(t%T) + theta_dash + np.pi/2. ) - np.sin(theta_dash + np.pi/2.) ) * r/L )
                else:
                    factor *= np.exp( ( np.sin( w*((t+theta_dash/w)%T) + np.pi/2. ) - 1 ) * r/L )
            #If in second half of cycle.
            else:
                for i in range(int(cycle_phase)+1):
                    if i == 0:
                        factor *= np.exp( (-1. - np.sin(theta_dash + np.pi/2.)) * r/L )
                    else:
                        factor *= expon
                    if i != int(cycle_phase):
                        factor = 1 - (1-factor)*expon
                factor = 1 - (1-factor) * np.exp( - (np.sin( ((t+theta_dash/w)%(T/2))*w + 3*np.pi/2.) + 1) * r/L )
   
        #Point is initially to right of central point.
        else:
            #As rho_2 initially 0, no source term until point has crossed boundary.
            if t < abs(theta_dash)/w:
                pass
            #If in first half of cycle.
            elif int(cycle_phase2)%2 == 0:
                for i in range(int(cycle_phase)):
                    factor *= expon
                    factor = 1 - (1-factor)*expon
                factor *= np.exp( ( np.sin( w*((t+theta_dash/w)%T) + np.pi/2. ) - 1 ) * r/L )
            #If in second half of cycle.
            else:
                for i in range(int(cycle_phase)+1):
                    factor *= expon
                    if i != int(cycle_phase):
                        factor = 1 - (1-factor)*expon
                factor = 1 - (1-factor) * np.exp( - (np.sin( ((t+theta_dash/w)%(T/2))*w + 3*np.pi/2.) + 1) * r/L )
            
        return factor*np.cos(r0*np.pi/2)**2
    else:
        return 0.
        
#Read log file and extract values of T1 fraction.
def readLogFile(filename,size):
    file = open(os.path.join(sys.path[0],filename),'r+')
    string = file.read()
    file.close()
    
    values = np.zeros(size)
    temp = 0
    # string_to_find = "T1 fraction: (sum(T1)|sum(T)) [0 0 0 0 0 0 0] "
    string_to_find = "T1 fraction: (sum((q*T))|sum(T)) [0 0 0 0 0 0 0] "
    
    for n in xrange(size):
        start = string.find(string_to_find,temp) + len(string_to_find)
        end = string.find("\n",start)
        values[n] = float(string[start:end])
        temp = end
        
    return values
    
#Read log file and extract values of T1 fraction.
def readBoundary(filename,size):
    file = open(os.path.join(sys.path[0],filename),'r+')
    string = file.read()
    file.close()
    
    values = np.zeros(size)
    values2 = np.zeros(size)
    temp = 0
    # string_to_find = "T1 fraction: (sum(T1)|sum(T)) [0 0 0 0 0 0 0] "
    string_to_find = "T goes from "
    string_to_find2 = "to "
    
    for n in xrange(size):
        start = string.find(string_to_find,temp) + len(string_to_find)
        end = string.find(" ",start)
        values[n] = float(string[start:end])
        temp = end
        
        start = string.find(string_to_find2,temp) + len(string_to_find2)
        end = string.find("\n",start)
        values2[n] = float(string[start:end])
        temp = end
        
    return values,values2
    
#Read log file and extract values of T1 fraction.
def readConservation(filename,size):
    file = open(os.path.join(sys.path[0],filename),'r+')
    string = file.read()
    file.close()
    
    values = np.zeros(size)
    temp = 0
    # string_to_find = "T1 fraction: (sum(T1)|sum(T)) [0 0 0 0 0 0 0] "
    string_to_find = "Total T in system: sum(T) [0 0 0 0 0 0 0] "
    
    for n in xrange(size):
        start = string.find(string_to_find,temp) + len(string_to_find)
        end = string.find("\n",start)
        values[n] = float(string[start:end])
        temp = end
        
    return values
    
#Read field file and extract field array.
def readField(filename,size):
    file = open(os.path.join(sys.path[0],filename),'r+')
    string = file.read()
    file.close()
    
    field = np.zeros((size,size))
    
    temp = string.find("internalField",0)
    temp = string.find("(",temp)
    temp = string.find("\n",temp) + 1

    for k in xrange(size):
        for j in xrange(size):
            start = temp
            end = string.find("\n",temp)
            field[k][j] = float( string[start:end] )
            temp = end + 1
    
    return field