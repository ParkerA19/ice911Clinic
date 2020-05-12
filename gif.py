# ICE911 Clinic Team
# Gif Generator 

import imageio
import os

def compatible(fname):
    return fname[0] != '.' # This is so we do not include any hidden system files

def gif(basepath, gifname):
    """
    This works when you have a directory of all images that you want to concatenate into a gif
    basebath is the path to the folder with the images
    gifname is the name of the gif you would to output
    """
    listdir = os.listdir(basepath)
    listdir = [l.lower() for l in listdir]
    listdir = sorted(listdir)

    # set that output file
    gifpath = basepath + gifname + '.gif'
    with imageio.get_writer(gifpath, mod='I', duration=0.25) as writer:   # Cnat have an existing gif already for this to work
        for fname in listdir:
            if compatible(fname):
                path = basepath + '/' + fname       # Set that image path
                image = imageio.imread(path)        # Read that image
                writer.append_data(image)           # Append that image
