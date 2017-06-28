import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def l2ErrorNorm(phi, phiExact):
	"Calculates the l2 error norm (RMS error) of phi in comparison to"
	"phiExact, ignoring the boundaries"
	
	phi_err = phi - phiExact
	return np.sqrt(np.sum(phi_err**2)/np.sum(phiExact**2)) 

def rvs(z):
    theta0 = 300
    P0 = 100000
    g = 9.81
    c_p = 1004
    muv = 0.018
    mud = 0.029
    
    z = z + 1000
    T = theta0*np.exp(-g*z/(c_p*theta0))
    P = P0*np.exp(-g*z/(c_p*theta0))
    es = 611.2*np.exp(17.67*(T-273.15)/(T-29.65))
    return muv/mud * es/(P-es)
    
execfile(os.path.join(sys.path[0],"core.py"))

####### SETTINGS
directory = sys.path[0]
velocity = 1
cmap = 0
q = 1
show_negative = 0

length = 100
x = np.linspace(-1000,1000,length)
z = np.linspace(0,2000,length)
x = np.linspace( (x[0] + x[1])/2., (x[-2] + x[-1])/2., length )
z = np.linspace( (z[0] + z[1])/2., (z[-2] + z[-1])/2., length )

rho1_analytic = np.zeros((length,length))
rho2_analytic = np.zeros((length,length))
rho12_analytic = np.zeros((length,length))
rho_analytic = np.zeros((length,length))

error = np.zeros((length,length))
error2 = np.zeros((length,length))
error3 = np.zeros((length,length))

x_center = 0.
z_center = 1000.

x = x - x_center
z = z - z_center
        
dist_center_x0 = 0.
dist_center_z0 = 500.
half_width = 300.

w = 0.01*2*np.pi/6
t = 0
L = 100000.        #Damping.
rhoAir = 1.
rho0 = 0.05

        
dir = sys.path[0]
dir = os.path.join(dir,"solidBodyRotationTemperature")

console = "magick convert -delay 20 "
folders = range(0,1818,18)
# folders = range(0,606,6)
# folders = range(0,101,1)
# folders = np.array([0])

l2err_rho1 = np.zeros(len(folders))
l2err_rho2 = np.zeros(len(folders))
l2err_rho12 = np.zeros(len(folders))

abserr_rho1 = np.zeros(len(folders))
abserr_rho2 = np.zeros(len(folders))
abserr_rho12 = np.zeros(len(folders))

total_rho1 = np.zeros(len(folders))
total_rho2 = np.zeros(len(folders))
total_rho12 = np.zeros(len(folders))

total_sa_rho1 = np.zeros(len(folders))
total_sa_rho2 = np.zeros(len(folders))
total_sa_rho12 = np.zeros(len(folders))
total_a_rho1 = np.zeros(len(folders))
total_a_rho2 = np.zeros(len(folders))
total_a_rho12 = np.zeros(len(folders))

min_rho1 = np.zeros(len(folders))
min_rho2 = np.zeros(len(folders))
min_rho12 = np.zeros(len(folders))

max_rho1 = np.zeros(len(folders))
max_rho2 = np.zeros(len(folders))
max_rho12 = np.zeros(len(folders))

for i2 in range(len(folders)):
    i = folders[i2]
    print i
    t = i
    directory = os.path.join(dir,str(int(i)))
    console += os.path.join(dir,"rho1_{0}.png ".format(i))
    
    theta = w*t
    dist_center_x,dist_center_z = rotate(dist_center_x0,dist_center_z0,theta)

    rho_numeric = readField(os.path.join(directory,"rho"),length)
    rho1_numeric = readField(os.path.join(directory,"q1"),length)
    rho2_numeric = rho_numeric - rho1_numeric
    rho2_numeric = readField(os.path.join(directory,"q2"),length)
    S_numeric = readField(os.path.join(directory,"S"),length)
    
    rho1_semi_analytic = readField(os.path.join(directory,"q1_analytic"),length)
    rho2_semi_analytic = readField(os.path.join(directory,"q2_analytic"),length)
    
    for j in xrange(len(x)):
        for k in xrange(len(z)):        
            rho12_analytic[k][j] = rho0*schaerProfile(x[j],z[k],dist_center_x,dist_center_z,half_width)
            rho2_analytic[k][j] = min(rho12_analytic[k][j],rvs(z[k]))
            rho1_analytic[k][j] = rho12_analytic[k][j] - rho2_analytic[k][j]
    
    total_rho1[i2] = np.sum(rho1_numeric)
    total_rho2[i2] = np.sum(rho2_numeric)
    # total_rho12[i2] = np.sum(rho_numeric-rhoAir)
    total_rho12[i2] = np.sum(rho1_numeric+rho2_numeric)
    total_sa_rho1[i2] = np.sum(rho1_semi_analytic)
    total_sa_rho2[i2] = np.sum(rho2_semi_analytic)
    total_sa_rho12[i2] = np.sum(rho_numeric-rhoAir)
    total_a_rho1[i2] = np.sum(rho1_analytic)
    total_a_rho2[i2] = np.sum(rho2_analytic)
    total_a_rho12[i2] = np.sum(rho12_analytic)
    
    min_rho1[i2] = np.min(rho1_numeric)
    min_rho2[i2] = np.min(rho2_numeric)
    min_rho12[i2] = np.min(rho1_numeric+rho2_numeric)
    
    max_rho1[i2] = np.max(rho1_numeric)
    max_rho2[i2] = np.max(rho2_numeric)
    max_rho12[i2] = np.max(rho1_numeric+rho2_numeric)
    
    if show_negative == 1:
        for j in xrange(len(x)):
            for k in xrange(len(z)):
                if rho1_numeric[k][j] > rho_numeric[k][j]:
                    rho1_numeric[k][j] = -rho0
                    # rho1_numeric[k][j] = 0.
                # if rho2_numeric[k][j] > 0.:
                    # rho2_numeric[k][j] = 0.
                if rho2_numeric[k][j] < 0.:
                    rho2_numeric[k][j] = -rho0
    
    plt.figure(figsize=(16,14))
    # plt.figure(figsize=(16,6))

    if cmap == 0:
        my_cmap = plt.cm.get_cmap('hot_r')
    elif cmap == 1:
        my_cmap = plt.cm.get_cmap('seismic')
    
    plt.subplot(221)
    data1 = rho1_analytic
    # data1 = rho_numeric
    CS = plt.pcolor(x, z, data1, cmap=my_cmap, vmin=0, vmax=rho0)
    # CS = plt.pcolor(x, z, data1, cmap=my_cmap, vmin=rhoAir, vmax=rhoAir+rho0)
    # CS = plt.pcolor(x, z, data3, cmap=my_cmap)
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel('Moisture Content')
    plt.axis([min(x),max(x),min(z),max(z)])
    plt.xlabel("$x$ (m)")
    plt.ylabel("$z$ (m)")
    # plt.title("$\\rho$, Min={:.4f}, Max={:.4f}".format(np.min(data1),np.max(data1)))
    plt.title("$r_l$ Analytic, Min={:.4f}, Max={:.4f}".format(np.min(data1),np.max(data1)))

    plt.subplot(222)
    data2 = rho2_analytic
    # data2 = S_numeric
    CS = plt.pcolor(x, z, data2, cmap=my_cmap, vmin=0, vmax=rho0)
    # CS = plt.pcolor(x, z, data2, cmap=plt.cm.get_cmap('seismic'), vmin=-0.001, vmax=0.001)
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel('Moisture Content')
    plt.axis([min(x),max(x),min(z),max(z)])
    plt.xlabel("$x$ (m)")
    plt.ylabel("$z$ (m)")
    plt.title("$r_v$ Analytic, Min={:.4f}, Max={:.4f}".format(np.min(data2),np.max(data2)))
    
    plt.subplot(223)
    data3 = rho1_numeric
    CS = plt.pcolor(x, z, data3, cmap=my_cmap, vmin=0., vmax=rho0)
    # CS = plt.pcolor(x, z, data3, cmap=my_cmap)
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel('Moisture Content')
    plt.axis([min(x),max(x),min(z),max(z)])
    plt.xlabel("$x$ (m)")
    plt.ylabel("$z$ (m)")
    plt.title("$r_l$ Numeric, Min={:.4f}, Max={:.4f}".format(np.min(data3),np.max(data3)))

    plt.subplot(224)
    data4 = rho2_numeric
    CS = plt.pcolor(x, z, data4, cmap=my_cmap, vmin=0., vmax=rho0)
    # CS = plt.pcolor(x, z, data4, cmap=my_cmap)
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel('Moisture Content')
    plt.axis([min(x),max(x),min(z),max(z)])
    plt.xlabel("$x$ (m)")
    plt.ylabel("$z$ (m)")
    plt.title("$r_v$ Numeric, Min={:.4f}, Max={:.4f}".format(np.min(data4),np.max(data4)))

    if q == 1:
        plt.suptitle("Schaer Profile: $t=${} s, L={}, Mean Error: {:.5f}, L2 Error: {:.5f}".format(t,L,abserr_rho1[i2],l2err_rho1[i2]),size=18)
    elif q == 2:
        plt.suptitle("Schaer Profile: $t=${} s, L={}, Mean Error: {:.5f}, L2 Error: {:.5f}".format(t,L,abserr_rho2[i2],l2err_rho2[i2]),size=18)

    plt.savefig(os.path.join(dir,"rho1_{0}.png".format(i)))
    plt.close()

# plt.figure()
# plt.subplot(211)
# plt.plot(folders,l2err_rho1,'b',label="$\\rho_1$")
# plt.plot(folders,l2err_rho2,'r',label="$\\rho_2$")
# plt.plot(folders,l2err_rho12,'k',label="$\\rho_1+ \\rho_2$")
# plt.ylabel("l2 Error Norm.")
# plt.title("Error Variation Relative to Analytic Function")
# plt.subplot(212)
# plt.plot(folders,abserr_rho1,'b',label="$\\rho_1$")
# plt.plot(folders,abserr_rho2,'r',label="$\\rho_2$")
# plt.plot(folders,abserr_rho12,'k',label="$\\rho_1+ \\rho_2$")
# plt.xlabel("Time ($t$)")
# plt.ylabel("Mean Abs Error")
# plt.legend()
# plt.savefig(os.path.join(sys.path[0],"output_errors.png"))
# plt.close()

plt.figure()
plt.plot(folders,total_a_rho1,'b:')
plt.plot(folders,total_a_rho2,'r:')
plt.plot(folders,total_a_rho12,'k:')
plt.plot(folders,total_sa_rho1,'b--',linewidth=1.5)
plt.plot(folders,total_sa_rho2,'r--',linewidth=1.5)
plt.plot(folders,total_sa_rho12,'k--',linewidth=1.5)
plt.plot(folders,total_rho1,'b',linewidth=2.5,label="Total $r_l$")
plt.plot(folders,total_rho2,'r',linewidth=2.5,label="Total $r_v$")
plt.plot(folders,total_rho12,'k',linewidth=2.5,label="Total $r_l+r_v$")
plt.ylabel("Total Tracer")
plt.title("Total Tracer Variation")
plt.xlabel("Time ($t$)")
plt.legend()
plt.savefig(os.path.join(sys.path[0],"output_total.png"))
plt.close()

plt.figure()
plt.subplot(211)
plt.plot(folders,min_rho1,'b',label="$r_l$")
plt.plot(folders,min_rho2,'r',label="$r_v$")
plt.plot(folders,min_rho12,'k',label="$r_l + r_v$")
plt.ylabel("Minimum Value")
plt.legend()
plt.title("Minimum and Maximum Variation with Time")
plt.subplot(212)
plt.plot(folders,max_rho1,'b',label="$r_l$")
plt.plot(folders,max_rho2,'r',label="$r_v$")
plt.plot(folders,max_rho12,'k',label="$\\rho$")
plt.xlabel("Time ($t$)")
plt.ylabel("Maximum Value")
plt.savefig(os.path.join(sys.path[0],"output_minmax.png"))
plt.close()

# plt.figure()
# plt.plot(x,rho_numeric[38],'b',label="Numeric")
# plt.plot(x,rho_analytic[38],'k',label="Analytic")
# plt.xlabel("$x$")
# plt.legend()
# plt.savefig(os.path.join(sys.path[0],"output_zzz1.png"))
# plt.close()

# plt.figure()
# plt.plot(z,rho_numeric[:,25],'b',label="Numeric")
# plt.plot(z,rho_analytic[:,25],'k',label="Analytic")
# plt.xlabel("$z$")
# plt.legend(loc=3)
# plt.savefig(os.path.join(sys.path[0],"output_zzz2.png"))
# plt.close()
    
if show_negative == 0:
    console += "-loop 0 {}".format( os.path.join(sys.path[0],"output_q%s.gif" % (q)) )
elif show_negative == 1:
    console += "-loop 0 {}".format( os.path.join(sys.path[0],"output_q%s_enhanced.gif" % (q)) )
# print console
os.system(console)
