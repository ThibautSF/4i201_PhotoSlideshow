#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main program

Author : SIMON-FINE & DE RYCKE
"""

import photo
import math


# Photos dicts : key,value → photoID,Photo
allPhotos = {}
hPhotos = {}
vPhotos = {}


def read_instance_file(filename, percent=100):
    file = open(filename, 'r')

    if percent < 0:
        percent = 0
    elif percent > 100:
        percent = 100

    percent = percent/100

    n = int(file.readline())

    lmax = math.floor(percent*n)

    i=0
    for line in file:
        if lmax == i:
            break

        lst = line.split(" ")

        orient = lst[0]
        nb_words = int(lst[1])
        words = lst[2:2+nb_words]

        one_photo = photo.Photo(str(i),orient,words)

        allPhotos[str(i)] = one_photo
        if orient == "H":
            hPhotos[str(i)] = one_photo
        elif orient == "V":
            vPhotos[str(i)] = one_photo

        i += 1


def simple_slideshow():
    slide_show = []  # list of lists (one element if H / two elements if V)

    for image_id in hPhotos:
        slide_show.append([image_id])

    vlst = []
    for image_id in vPhotos:
        vlst.append(image_id)
        if len(vlst)>=2:
            slide_show.append(vlst)
            vlst = []

    return slide_show


def write_slideshow(slideshow,filename="output.txt"):
    file = open(filename,'w')
    file.write(str(len(slideshow))+'\n')

    for image in slideshow:
        file.write(" ".join(image)+'\n')


if __name__ == '__main__':
    read_instance_file("ressources/a_example.txt",100)
    print(allPhotos)
    print(hPhotos)
    print(vPhotos)

    simple_ss = simple_slideshow()
    print(simple_ss)

    write_slideshow(simple_ss)
