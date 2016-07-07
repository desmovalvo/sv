#!/usr/bin/python

# requirements
from scipy.spatial import distance
from termcolor import colored
from random import randint
import numpy as np
import getopt
import time
import sys
import cv2

# main
if __name__ == "__main__":

    # read command line parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], "vl:r:", ["verbose", "leftimg=", "rightimg="])
    except getopt.GetoptError as err:
        sys.exit(2)

    verbose = False
    leftimage = None
    rightimage = None
    for opt, arg in opts:
        if opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-l", "--leftimg"):
            leftimage = arg
        elif opt in ("-l", "--rightimg"):
            rightimage = arg
        else:
            assert False, "unhandled option"

    # algorithm parameters
    disp_range = 10

    # open images
    print colored("main> ", "blue", attrs=["bold"]) + "Opening image %s" % leftimage
    l_img = cv2.imread(leftimage)
    print colored("main> ", "blue", attrs=["bold"]) + "Image %s has resolution %sx%s" % (leftimage, l_img.shape[1], l_img.shape[0])

    print colored("main> ", "blue", attrs=["bold"]) + "Opening image %s" % rightimage
    r_img = cv2.imread(rightimage)
    print colored("main> ", "blue", attrs=["bold"]) + "Image %s has resolution %sx%s" % (rightimage, r_img.shape[1], r_img.shape[0])

    # initialization of outputimage
    print colored("main> ", "blue", attrs=["bold"]) + "Initializing output image"
    outimage = "out.png"
    o_width = max(r_img.shape[1], l_img.shape[1])
    o_height = max(r_img.shape[0], l_img.shape[0])
    size = (o_height, o_width, 1)
    o_img = np.zeros(size, dtype=np.uint8)

    # segmentation
    # print colored("main> ", "blue", attrs=["bold"]) + "Performing segmentation on image %s" % leftimage,
    # start_time = time.clock() * 1000
    # l_s_img = cv2.pyrMeanShiftFiltering(l_img, 30, 30)
    # end_time = time.clock() * 1000
    # print "%s ms" % round(end_time - start_time, 3)

    # print colored("main> ", "blue", attrs=["bold"]) + "Performing segmentation on image %s" % rightimage,
    # start_time = time.clock() * 1000
    # r_s_img = cv2.pyrMeanShiftFiltering(r_img, 30, 1)
    # end_time = time.clock() * 1000
    # print "%s ms" % round(end_time - start_time, 3)

    y = 0
    for pixel in xrange(o_height * o_width):
        
        # determine x and y
        newy = pixel / o_width
        if y != newy:
            y = newy
            print "RIGA %s" % y
        x = pixel % o_width
        
        # for the given pixel we should scan the entire
        # line of the target image and calculate the residual
            
        # the starting pixel is at (x,y) in the reference image
        # let's check every pixel in the target image at (x+d,y)            
        min_distance = 100000
        disparity = 0
        ref_pix = l_img[y,x]
        for xx in xrange(x-disp_range, x+disp_range):
            try:
                tar_pix = r_img[y, xx]
                d = distance.euclidean(ref_pix, tar_pix) 
                if d < min_distance:
                    min_distance = d
                    disparity = x - xx
            except:
                pass
                
        # pixel value
        pv = (float(255 * abs(disparity)) / (2 * disp_range))
        o_img.itemset((y, x, 0), pv)

    # display the output image
    cv2.imshow("output image", o_img)
    cv2.waitKey(0)

    # write the image to file
    cv2.imwrite(outimage, o_img)
