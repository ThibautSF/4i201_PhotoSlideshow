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

    l_max = math.floor(percent*n)

    i=0
    for line in file:
        if l_max == i:
            break

        line = line.rstrip()

        lst = line.split(" ")

        orient = lst[0]
        nb_words = int(lst[1])
        words = set(lst[2:2+nb_words])

        one_photo = photo.Photo(str(i),orient,words)

        allPhotos[str(i)] = one_photo
        if orient == "H":
            hPhotos[str(i)] = one_photo
        elif orient == "V":
            vPhotos[str(i)] = one_photo

        i += 1


def simple_slide_show():
    slide_show = []  # list of lists (one element if H / two elements if V)

    for image_id in hPhotos:
        slide_show.append([image_id])

    v_lst = []
    for image_id in vPhotos:
        v_lst.append(image_id)
        if len(v_lst) >= 2:
            slide_show.append(v_lst)
            v_lst = []

    return slide_show


def write_slide_show(slide_show,filename="output.txt"):
    file = open(filename,'w')
    file.write(str(len(slide_show))+'\n')

    for image in slide_show:
        file.write(" ".join(image)+'\n')


def calc_score(slide_show):
    score = 0
    i = 0
    while i < len(slide_show)-1:
        pos_slide = slide_show[i]
        next_slide = slide_show[i+1]

        pos_slide_words = set()
        for photo_id in pos_slide:
            pos_slide_words |= allPhotos[photo_id].keyWords

        next_slide_words = set()
        for photo_id in next_slide:
            next_slide_words |= allPhotos[photo_id].keyWords

        score += min(len(pos_slide_words & next_slide_words), len(pos_slide_words - next_slide_words), len(next_slide_words - pos_slide_words))

        i += 1

    return score


if __name__ == '__main__':
    read_instance_file("ressources/a_example.txt",100)
    print(allPhotos)
    print(hPhotos)
    print(vPhotos)

    simple_ss = simple_slide_show()
    print(simple_ss)

    write_slide_show(simple_ss)

    print(calc_score(simple_ss))
