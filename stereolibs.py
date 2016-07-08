#!/usr/bin/python

import pdb
import sys
import ConfigParser
from math import sqrt
from termcolor import colored

###############################################################
#
# Pixel-based - the naive algorithm -- ALG 0
#
###############################################################

# pixel based
def pixelbased(ref_image, tar_image, out_image, settings_file):

    # read settings
    print colored("pixelbased> ", "blue", attrs=["bold"]) + "Reading algorithm settings"
    config = ConfigParser.ConfigParser()
    config.read(settings_file)
    settings = {}
    settings["policy"] = config.get("pixelbased", "policy")
    settings["disp_range"] = config.getint("pixelbased", "disp_range")
    try:
        settings["threshold"] = config.getint("pixelbased", "threshold")
    except:
        pass

    # get height and width
    print colored("pixelbased> ", "blue", attrs=["bold"]) + "Reading image properties"
    width = max(out_image.shape[1], out_image.shape[1])
    height = max(out_image.shape[0], out_image.shape[0])

    # iterate over the pixels
    print colored("pixelbased> ", "blue", attrs=["bold"]) + "Building disparity map"
    y = 0
    for pixel in xrange(height * width):
        
        # determine x and y
        y = pixel / width
        x = pixel % width

        # get a pixel from the reference image
        ref_pix = ref_image[y,x]
           
        # per-pixel initialization
        disparity = 0
        min_distance = sys.maxint

        # calculate the output value
        for xx in xrange(max(x - settings["disp_range"], 0), min(x + settings["disp_range"], width)):
            tar_pix = tar_image[y, xx]

            # matching cost
            if settings["policy"] == "ABS_DIF":
                d = abs(sqrt((int(tar_pix[2]) - int(ref_pix[2]))**2 + (int(tar_pix[1]) - int(ref_pix[1]))**2 + (int(tar_pix[0]) - int(ref_pix[0]))**2))
            elif settings["policy"] == "SQR_DIF":
                d = (int(tar_pix[2]) - int(ref_pix[2]))**2 + (int(tar_pix[1]) - int(ref_pix[1]))**2 + (int(tar_pix[0]) - int(ref_pix[0]))**2
            elif settings["policy"] == "TRA_DIF":
                d = min(abs(sqrt((int(tar_pix[2]) - int(ref_pix[2]))**2 + (int(tar_pix[1]) - int(ref_pix[1]))**2 + (int(tar_pix[0]) - int(ref_pix[0]))**2)), settings["threshold"])

            if d < min_distance:
                min_distance = d
                disparity = x - xx

        # determine the pixel value for the output image
        pv = int(float(255 * abs(disparity)) / (settings["disp_range"]))
        out_image.itemset((y, x, 0), pv)

    # return
    return out_image


###############################################################
#
# Integral Images - optimization function
#
###############################################################

def build_integral_image_matrix(img):

    """This function is used to calculate the integral image
    for a given picture"""

    # get height and width
    width = img.shape[1]
    height = img.shape[0]

    # initialize the ii matrix
    ii = []

    # algorithm
    for y in xrange(height):
        
        # add a line
        ii.append([])

        for x in xrange(width):

            ii[y].append(None)
            norm = None

            if (y > 0) and (x > 0):
                norm = sqrt(img[y][x][0]**2 + img[y][x][1]**2 + img[y][x][2]**2)
                ii[y][x] = ii[y-1][x] + ii[y][x-1] - ii[y-1][x-1] + norm

            elif (y > 0) and (x == 0):
                norm = sqrt(img[y][x][0]**2 + img[y][x][1]**2 + img[y][x][2]**2)
                ii[y][x] = ii[y-1][x] + norm

            elif (y == 0) and (x > 0):
                norm = sqrt(img[y][x][0]**2 + img[y][x][1]**2 + img[y][x][2]**2)
                ii[y][x] = ii[y][x-1] + norm
                
            else:
                norm = sqrt(img[y][x][0]**2 + img[y][x][1]**2 + img[y][x][2]**2)
                ii[y][x] = norm

    # return
    return ii


def get_integral(x, y, width, height, ii, window_size):

    """This function is used to get the result of the integral
    for a given window"""

    d = ii[min(y + window_size, height-1)][min(x + window_size, width - 1)]
    c = ii[min(y + window_size, height - 1)][max(0, x - window_size)]
    b = ii[max(0, y - window_size)][min(x + window_size, width - 1)]
    a = ii[max(0, y - window_size)][max(0, x - window_size)]

    return d - c - b + a


###############################################################
#
# Fixed Window -- ALG 1
#
###############################################################

def fixedwindow(ref_image, tar_image, out_image, settings_file):

    # read settings
    print colored("fixedwindow> ", "blue", attrs=["bold"]) + "Reading algorithm settings"
    config = ConfigParser.ConfigParser()
    config.read(settings_file)
    settings = {}
    settings["disp_range"] = config.getint("fixedwindow", "disp_range")
    settings["window_size"] = config.getint("fixedwindow", "window_size")
    settings["policy"] = config.get("fixedwindow", "policy")
    try:
        settings["threshold"] = config.getint("fixedwindow", "threshold")
    except:
        pass

    # get height and width
    print colored("fixedwindow> ", "blue", attrs=["bold"]) + "Reading image properties"
    width = max(out_image.shape[1], out_image.shape[1])
    height = max(out_image.shape[0], out_image.shape[0])

    # build the integral images matrices
    print colored("fixedwindow> ", "blue", attrs=["bold"]) + "Building integral images matrices"
    ref_ii = build_integral_image_matrix(ref_image)
    tar_ii = build_integral_image_matrix(tar_image)

    # iterate over the pixels
    print colored("fixedwindow> ", "blue", attrs=["bold"]) + "Building disparity map"
    y = 0
    for pixel in xrange(height * width):

        # determine x and y
        y = pixel / width
        x = pixel % width

        # initialize disparity and distance
        min_distance = sys.maxint
        disparity = 0
        
        # aggregation for the reference pixel
        ref_sum = get_integral(x, y, width, height, ref_ii, settings["window_size"])

        # iterate over the pixel of the target image
        # between x-d and x+d 
        for xx in xrange(max(x-settings["disp_range"], 0), min(x+settings["disp_range"], width)):
            
            d = 0
            tar_sum = get_integral(xx, y, width, height, tar_ii, settings["window_size"])                    

            # matching cost
            if settings["policy"] == "ABS_DIF":
                d = abs(tar_sum - ref_sum)
            elif settings["policy"] == "SQR_DIF":
                d = (tar_sum - ref_sum) ** 2
            elif settings["policy"] == "TRA_DIF":
                d = min(abs(tar_sum - ref_sum), settings["threshold"])

            if d < min_distance:
                min_distance = d
                disparity = x - xx
        
        # determine the pixel value for the output image
        pixel_value = int(float(255 * abs(disparity)) / (settings["disp_range"]))
        out_image.itemset((y, x, 0), pixel_value)

    # return
    return out_image


###############################################################
#
# Segmentation
#
###############################################################

def segment(img):
    
    """This function is used to segment an image"""
    
    # segmentation
    print colored("main> ", "blue", attrs=["bold"]) + "Performing segmentation on image"
    seg_img = cv2.pyrMeanShiftFiltering(img, 30, 30)

    # return
    return seg_img
