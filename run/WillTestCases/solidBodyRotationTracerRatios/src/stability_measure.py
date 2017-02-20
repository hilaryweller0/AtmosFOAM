import sys 
import os
import numpy as np

dt = np.linspace(1.0,1.8,3)
# dt = np.array([0.9])

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
    file.write(controlDictString % (dt[n],writeInterval))
    file.close()
    
    os.system( os.path.join(base_dir,"Allclean") )
    os.system( os.path.join(base_dir,"Allrun") )
    os.system( "mv " + os.path.join(base_dir,"{}/T".format(writeInterval)) + " " + os.path.join(base_base_dir,"T_{}_{}".format(id,n)) )
    