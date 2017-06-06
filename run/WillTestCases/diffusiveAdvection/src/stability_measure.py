import sys 
import os
import numpy as np
import matplotlib.pyplot as plt

analysis = 1
plots = 1

dt = np.arange(1.9,2.25,0.05)
#dt = np.array([3.039,3.0392,3.0393,3.0394,3.0395,3.0397,3.0398])
#dt = np.array([3.])

if analysis == 1:

    controlDictString = '''FoamFile {
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      controlDict;
    }

    startFrom       startTime;
    startTime       0;
    stopAt          endTime;
    endTime         610;
    deltaT          %s;
    writeControl    runTime;
    writeInterval   %s;
    purgeWrite      0;
    writeFormat     ascii;
    writePrecision  6;
    writeCompression uncompressed;
    timeFormat      general;
    timePrecision   6;
    runTimeModifiable no;
    adjustTimeStep no;
    maxCo          1;
    libs           ("libfiniteVolumeAtmos.so");
    '''

    id = "cubic_exp"
    dir = sys.path[0]
    base_dir = os.path.join(dir,"..")
    base_base_dir = os.path.join(base_dir,"..")
    system = os.path.join(base_dir,"system")
    controlDict = os.path.join(system,"controlDict")
    for n in xrange(len(dt)):
        print n
        writeInterval = (int(600/dt[n]) + 1)*dt[n]
        file = open(controlDict,"wb")
        print controlDictString % (dt[n],writeInterval)
        file.write(controlDictString % (dt[n],writeInterval))
        file.close()
        
        os.system( os.path.join(base_dir,"Allclean") )
        os.system( os.path.join(base_dir,"Allrun") )
        os.system( "mv " + os.path.join(base_dir,"{}/rho".format(writeInterval)) + " " + os.path.join(base_dir,"rho_{}_{}".format(id,n)) )
        os.system( "mv " + os.path.join(base_dir,"{}/rho".format(int(round(writeInterval)))) + " " + os.path.join(base_dir,"rho_{}_{}".format(id,n)) )
    
'''
POSTPROCESSING
'''
if plots == 1:
    execfile(os.path.join(sys.path[0],"core.py"))

    def l2ErrorNorm(phi, phiExact):
	    "Calculates the l2 error norm (RMS error) of phi in comparison to"
	    "phiExact, ignoring the boundaries"
	
	    phi_err = phi - phiExact
	    return np.sqrt(np.sum(phi_err**2)/np.sum(phiExact**2)) 

    length = 50

    errors = np.zeros_like(dt)
    profile = np.zeros((length,length))
    data = np.zeros((length,length))

    w = 0.01*2*np.pi/6
    L = 1000.        #Damping.
    t = (int(600/dt[n]) + 1)*dt[n]
    theta = w*t

    x = np.linspace(-1000,1000,length)
    z = np.linspace(0,2000,length)
    x_center = 0.
    z_center = 1000.
    x = x - x_center
    z = z - z_center
    dist_center_x0 = 0.
    dist_center_z0 = 500.
    half_width = 300.
    dist_center_x,dist_center_z = rotate(dist_center_x0,dist_center_z0,theta)
    for j in xrange(len(x)):
        for k in xrange(len(z)):        
            profile[k][j] = schaerProfileVelocity(x[j],z[k],dist_center_x,dist_center_z,half_width,theta,w,t,L)
            
    for n in xrange(len(dt)):
        t = writeInterval
        theta = w*t
        dist_center_x,dist_center_z = rotate(dist_center_x0,dist_center_z0,theta)
        for j in xrange(len(x)):
            for k in xrange(len(z)):        
                profile[k][j] = 0.5 + 0.5*schaerProfileVelocity(x[j],z[k],dist_center_x,dist_center_z,half_width,theta,w,t,L)
        data = readField(os.path.join("..","rho_cubic_exp_{}".format(n)),length)
        errors[n] = l2ErrorNorm(data,profile)
        
    plt.clf()
    plt.plot(dt,errors,'k')
    plt.plot(dt,errors,'ko')
    plt.yscale("log")
    plt.xlabel("$dt$")
    plt.ylabel("l2 Error Norm")
    plt.title("Distribution error variation with time step.")
    plt.savefig(os.path.join(sys.path[0],"stability.png"))
