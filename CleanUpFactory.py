import os, os.path

for root, _, files in os.walk("DetectorFactory/"):
    for f in files:
        fullpath = os.path.join(root, f)
# remove amap, dmap and gdml, not the builder ! 
        if ( fullpath.find(".amap") != -1 or fullpath.find(".dmap") != -1 or fullpath.find(".gdml") != -1  )  :
           if ( fullpath.find("MyStore") == -1 or fullpath.find("Imports") == -1 ) :
                print '- erasing ' + fullpath
                os.remove(fullpath)
#
