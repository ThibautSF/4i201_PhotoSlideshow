#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main program

Author : SIMON-FINE & DE RYCKE
"""

import photo
import math
import random
import collections
import copy
import time
import constraint as cs
import pulp as plp

# Photos dicts : key,value → photoID,Photo
allPhotos = collections.OrderedDict()
hPhotos = collections.OrderedDict()
vPhotos = collections.OrderedDict()


def read_instance_file(filename, percent=100):
    file = open(filename, 'r')

    if percent < 0:
        percent = 0
    elif percent > 100:
        percent = 100

    percent = percent / 100

    n = int(file.readline())

    l_max = math.floor(percent * n)

    i = 0
    for line in file:
        if l_max == i:
            break

        line = line.rstrip()

        lst = line.split(" ")

        orient = lst[0]
        nb_words = int(lst[1])
        words = set(lst[2:2 + nb_words])

        one_photo = photo.Photo(str(i), orient, words)

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


def write_slide_show(slide_show, filename="output.txt"):
    file = open(filename, 'w')
    file.write(str(len(slide_show)) + '\n')

    for image in slide_show:
        file.write(" ".join(image) + '\n')


def calc_score(slide_show):
    score = 0
    i = 0
    while i < len(slide_show) - 1:
        score += score_transition(slide_show[i], slide_show[i + 1])

        i += 1

    return score


def score_transition(slide_a, slide_b):
    score = 0

    a_words = set()
    for photo_id in slide_a:
        a_words |= allPhotos[photo_id].keyWords

    b_words = set()
    for photo_id in slide_b:
        b_words |= allPhotos[photo_id].keyWords

    score += min(len(a_words & b_words), len(a_words - b_words),
                 len(b_words - a_words))

    return score


def compute_slides():
    slides = []  # list of lists (one element if H / two elements if V)

    for image_id in hPhotos:
        slides.append([image_id])

    lst_v_keys = list(vPhotos.keys())
    while len(lst_v_keys) > 1:
        image1 = lst_v_keys[0]
        lst_v_keys.remove(image1)
        words = set()
        words |= allPhotos[image1].keyWords

        best_image = lst_v_keys[0]
        best_size = len(words | allPhotos[best_image].keyWords)
        j = 1
        while j < len(lst_v_keys):
            image = lst_v_keys[j]
            words2 = words | allPhotos[image].keyWords
            if len(words2) > best_size:
                best_image = image
                best_size = len(words2)
            j += 1

        slides.append([image1, best_image])
        lst_v_keys.remove(best_image)

    return slides


def glouton_slide_show(isfast=True):
    slide_show = []  # list of lists (one element if H / two elements if V)
    loc_h_photos = hPhotos.copy()
    loc_v_photos = vPhotos.copy()

    # Begin init slideshow
    a_photo_id, a_photo = random.choice(list(allPhotos.items()))

    first_slide = list()
    first_slide.append(a_photo_id)

    if a_photo.orient == "H":
        loc_h_photos.pop(a_photo_id)
    elif a_photo.orient == 'V':
        loc_v_photos.pop(a_photo_id)
        b_photo_id = random.choice(list(loc_v_photos.keys()))
        first_slide.append(b_photo_id)
        loc_v_photos.pop(b_photo_id)

    slide_show.append(first_slide)
    # End init slideshow

    while len(loc_h_photos) > 0 or len(loc_v_photos) > 1:
        if not isfast:
            slide_show, loc_h_photos, loc_v_photos = glouton_add_slide(slide_show, loc_h_photos, loc_v_photos)
        else:
            slide_show, loc_h_photos, loc_v_photos = glouton_addfast_slide(slide_show, loc_h_photos, loc_v_photos)

    return slide_show


def glouton_add_slide(slide_show, h_photos, v_photos):
    best_transition = list()
    score = -1

    # Check H photos
    for h_photo_id in h_photos:
        score_t = score_transition(slide_show[-1], [h_photo_id])
        if score_t > score:
            score = score_t
            best_transition = [h_photo_id]

    # Check V photos
    lst_v_keys = list(v_photos.keys())
    i = 0
    while i < len(lst_v_keys) - 1:
        j = i + 1
        while j < len(lst_v_keys):
            score_t = score_transition(slide_show[-1], [lst_v_keys[i], lst_v_keys[j]])
            if score_t > score:
                score = score_t
                best_transition = [lst_v_keys[i], lst_v_keys[j]]
            j += 1
        i += 1

    slide_show.append(best_transition)

    # Remove choosed photo from choice
    if len(best_transition) == 1:
        h_photos.pop(best_transition[0])
    elif len(best_transition) == 2:
        v_photos.pop(best_transition[0])
        v_photos.pop(best_transition[1])

    return slide_show, h_photos, v_photos


def glouton_addfast_slide(slide_show, h_photos, v_photos):
    best_transition = list()
    score = -1

    # Check H photos
    for h_photo_id in h_photos:
        score_t = score_transition(slide_show[-1], [h_photo_id])
        if score_t > score:
            score = score_t
            best_transition = [h_photo_id]
            if score > 0:
                break

    # Check V photos
    lst_v_keys = list(v_photos.keys())
    brk = False
    i = 0
    while i < len(lst_v_keys) - 1 and not brk:
        j = i + 1
        while j < len(lst_v_keys) and not brk:
            score_t = score_transition(slide_show[-1], [lst_v_keys[i], lst_v_keys[j]])
            if score_t > score:
                score = score_t
                best_transition = [lst_v_keys[i], lst_v_keys[j]]
                if score > 0:
                    brk = True
            j += 1
        i += 1

    slide_show.append(best_transition)

    # Remove choosed photo from choice
    if len(best_transition) == 1:
        h_photos.pop(best_transition[0])
    elif len(best_transition) == 2:
        v_photos.pop(best_transition[0])
        v_photos.pop(best_transition[1])

    return slide_show, h_photos, v_photos


def descente_stochastique(slide_show, taille):
    score = calc_score(slide_show)
    for i in range(taille):
        slide_show_res = copy.deepcopy(slide_show)
        index1 = random.randint(0, len(slide_show) - 1)
        first_image = slide_show[index1]
        if allPhotos[first_image[0]].orient == 'V':
            index2 = random.randint(0, len(slide_show) - 1)
            second_image = slide_show[index2]
            while allPhotos[second_image[0]].orient != 'V':
                index2 = random.randint(0, len(slide_show) - 1)
                second_image = slide_show[index2]
            first = random.randint(0, 1)
            second = random.randint(0, 1)
            slide_show_res[index1][first] = second_image[second]
            slide_show_res[index2][second] = first_image[first]
        else:
            index2 = random.randint(0, len(slide_show) - 1)
            second_image = slide_show[index2]
            while allPhotos[second_image[0]].orient != 'H':
                index2 = random.randint(0, len(slide_show) - 1)
                second_image = slide_show[index2]

            slide_show_res[index1][0] = second_image[0]
            slide_show_res[index2][0] = first_image[0]

        ancienne_transition = 0
        if index1 < len(slide_show) - 1:  # Vérifie si index1 est le dernier élément de la liste
            ancienne_transition += score_transition(slide_show[index1], slide_show[index1 + 1])
        if index1 > 0:
            ancienne_transition += score_transition(slide_show[index1], slide_show[index1 - 1])
        if index2 < len(slide_show) - 1:
            ancienne_transition += score_transition(slide_show[index2], slide_show[index2 + 1])
        if index2 > 0:
            ancienne_transition += score_transition(slide_show[index2], slide_show[index2 - 1])

        nouvelle_transition = 0
        if index1 < len(slide_show) - 1:
            nouvelle_transition += score_transition(slide_show_res[index1], slide_show_res[index1 + 1])
        if index1 > 0:
            nouvelle_transition += score_transition(slide_show_res[index1], slide_show_res[index1 - 1])
        if index2 < len(slide_show) - 1:
            nouvelle_transition += score_transition(slide_show_res[index2], slide_show_res[index2 + 1])
        if index2 > 0:
            nouvelle_transition += score_transition(slide_show_res[index2], slide_show_res[index2 - 1])

        nouveau_score = score - ancienne_transition + nouvelle_transition

        if score < nouveau_score:
            slide_show = copy.deepcopy(slide_show_res)
            score = nouveau_score

    return slide_show


def plne_slide_show(lst_slides):
    # Get all slide couple
    problem = cs.Problem()
    indices_lst = range(len(lst_slides))
    problem.addVariables(["i","j"], indices_lst)
    problem.addConstraint(cs.AllDifferentConstraint())
    couples = problem.getSolutions()

    del problem

    model, s_vals = plne_problem(lst_slides, couples)

    cycle_solution = dict()
    for variable in model.variables():
        if variable.varValue > 0 and 'X' in variable.name:
            # print("{} = {}".format(variable.name, variable.varValue))
            lst = variable.name.split("_")
            i = int(lst[1])
            j = int(lst[2])
            score = s_vals[frozenset({i, j})]

            # print("{};{};{};{}".format(lst[0],lst[1],lst[2],variable.varValue))

            if i not in cycle_solution:
                cycle_solution[i] = (j, variable.varValue, score)
            else:
                if variable.varValue > cycle_solution[i][1]:
                    cycle_solution[i] = (j, variable.varValue, score)

    print(cycle_solution)

    score_cycle = 0
    min_score = -1
    worst_transition = None
    for i in cycle_solution:
        j = cycle_solution[i][0]
        score = cycle_solution[i][2]

        score_cycle += score

        if min_score == -1 or score < min_score:
            min_score = score
            worst_transition = (i, j)

    print("Score plne = {0} / Score cycle = {1}".format(plp.value(model.objective), score_cycle))
    print("Transition du cycle à briser {0} de score {1}".format(worst_transition, min_score))

    slide_show_sol = list()
    i = worst_transition[1]
    while i != worst_transition[0]:
        slide_show_sol.append(lst_slides[i])
        i = cycle_solution[i][0]
    slide_show_sol.append(lst_slides[i])

    return slide_show_sol


def plne_problem(lst_slides, couples):
    # Prepare slideshow problem
    model = plp.LpProblem(name="Slideshow")

    # k_vals = range(1, len(lst_slides))
    v = len(lst_slides) - 1
    s_vals = dict()  # Toutes constantes S(i,j) → score transition i à j
    x_vars = dict()  # Toutes variables X(i,j) → 1 si slide i et j se suivent / 0 sinon
    z_vars = dict()  # Toutes les variables Z(i,j) rassemblées par i,j
    transition_out_i = dict()  # Toutes variables X(i,j) rassemblées par i (pour contrainte 1 slide après)
    transition_in_j = dict()  # Toutes variables X(i,j) rassemblées par j (pour contrainte 1 slide avant)

    transition_cycle = dict()  # Toutes variables X1(i,j) X2(j,i) rassemblées par i-j

    t1 = time.clock()
    for one_couple in couples:
        i = one_couple["i"]
        j = one_couple["j"]

        var_s_name = "S_{}_{}".format(i, j)
        var_x_name = "X_{}_{}".format(i, j)
        var_z_name = "Z_{}_{}".format(i, j)
        t_name = frozenset({i, j})  # Pour pas qu'il y ait de doublon

        # Constante S(i,j) → score transition de slide position i à position j
        if t_name not in s_vals:
            s_vals[t_name] = score_transition(lst_slides[i], lst_slides[j])

        # Variable X(i,j) → 1 si slide i et j se suivent / 0 sinon
        x_var = plp.LpVariable(cat=plp.LpBinary, name=var_x_name)

        if i not in x_vars:
            x_vars[i] = dict()
        x_vars[i][j] = x_var

        z_var = plp.LpVariable(cat=plp.LpContinuous, name=var_z_name, lowBound=0)
        if i not in z_vars:
            z_vars[i] = dict()
        z_vars[i][j] = z_var

        if i not in transition_out_i:
            transition_out_i[i] = set()

        transition_out_i[i].add(x_var)
        if len(transition_out_i[i]) == v:
            c = plp.LpConstraint(e=plp.lpSum(var for var in transition_out_i[i]),
                                 sense=plp.LpConstraintEQ,
                                 rhs=1,
                                 name="constraint_after_{0}".format(i))
            model.addConstraint(c)
            transition_out_i.pop(i)

        if j not in transition_in_j:
            transition_in_j[j] = set()

        transition_in_j[j].add(x_var)
        if len(transition_in_j[j]) == v:
            c = plp.LpConstraint(e=plp.lpSum(var for var in transition_in_j[j]),
                                 sense=plp.LpConstraintEQ,
                                 rhs=1,
                                 name="constraint_before_{0}".format(j))
            model.addConstraint(c)
            transition_in_j.pop(j)

        if t_name not in transition_cycle:
            transition_cycle[t_name] = set()
        transition_cycle[t_name].add(x_var)
        if len(transition_cycle[t_name]) == 2:
            c = plp.LpConstraint(e=plp.lpSum(var for var in transition_cycle[t_name]),
                                 sense=plp.LpConstraintLE,
                                 rhs=1,
                                 name="constraint_cycle_{}-{}".format(i, j))
            model.addConstraint(c)

    t2 = time.clock()

    print("Var computing in : " + str(t2 - t1) + "s")

    t1 = time.clock()
    # Contrainte : somme sur j privé de 1 de Z(1,j) = |V|-1
    c = plp.LpConstraint(e=plp.lpSum(z_vars[0][j] for j in z_vars[0]),
                         sense=plp.LpConstraintEQ,
                         rhs=v,
                         name="constraint_z1")
    model.addConstraint(c)

    # Contraintes : somme sur j privé de {1 et i } de Z(i,j,k) = La somme sur j privé de {i} de Z(j,i,k) ( avec k
    # != 1 et i différent de 1 et k) print("k="+str(k))
    for i in z_vars:
        if i == 0:
            continue

        # print("i="+str(i))
        left_c = plp.lpSum(z_vars[i][j] for j in z_vars[i] if (j != 0 and j != i)) + 1
        # print(left_c)
        right_c = plp.lpSum(z_vars[j][i] for j in z_vars[i] if j != i)
        # print(right_c)
        c = left_c == right_c
        model.addConstraint(c, "constraint_z({0},j)_z(j,{0})".format(i))
        # print(str(left_c) + " == " + str(right_c))

    # Contraintes Z(i,j) + Z(j,i) <= (|V|-1)(X(i,j) + X(j,i))
    for i in range(len(lst_slides)):
        for j in range(len(lst_slides)):
            if j == 0 or j == i:
                continue

            left_c = z_vars[i][j] + z_vars[j][i]
            transition = frozenset({i, j})
            right_c = v * (plp.lpSum(var for var in transition_cycle[transition]))

            c = left_c <= right_c
            model.addConstraint(c, "constraint_zinfx_{}_{}".format(i, j))

    t2 = time.clock()

    print("Constraint computing in : " + str(t2 - t1) + "s")

    t1 = time.clock()
    objective = plp.lpSum((-s_vals[frozenset({i, j})] * x_vars[i][j] for j in x_vars[i]) for i in x_vars)
    model.sense = plp.LpMinimize
    model.setObjective(objective)
    t2 = time.clock()
    print("Objective computing in : " + str(t2 - t1) + "s")

    t1 = time.clock()
    if plp.GLPK().available():
        model.solve(plp.GLPK(options=[]))
    else:
        model.solve()

    t2 = time.clock()
    print("Problem computing in : " + str(t2 - t1) + "s")

    return model, s_vals


def heuristique_arrondi(photo,hphoto,vphoto):

    liste_base = []
    for image_id in hPhotos:
        liste_base.append([image_id])

    v_lst = []
    for image_id in vPhotos:
        v_lst.append(image_id)
        if len(v_lst) >= 2:
            liste_base.append(v_lst)
            v_lst = []

    y = {}

    for i in range(len(liste_base)-1):
        for j in range(i+1, len(liste_base)):
            y[str(i) + ";" + str(j)] = score_transition(liste_base[i],liste_base[j])
    y = sorted(y.items(), key=lambda t: t[1])
    y.reverse()

    liste_couple_interdit = []
    liste_nb_pris=[]
    liste_couple_select = []
    for i in y:
        petit = i[0][0:i[0].find(';')]
        grand = i[0][i[0].find(';')+1:len(i[0])]
        if int(petit) not in liste_nb_pris and int(grand) not in liste_nb_pris and i[0] not in liste_couple_interdit:
            if liste_couple_interdit != []:
                jpetit = -1
                jgrand = -1
                valp = ""
                valg = ""
                for j in range(len(liste_couple_interdit)):
                    if liste_couple_interdit[j][0:liste_couple_interdit[j].find(";")] == petit:
                        valp =  liste_couple_interdit[j][liste_couple_interdit[j].find(";")+1:len(liste_couple_interdit[j])]
                        jpetit = j
                        liste_nb_pris.append(int(petit))
                    elif liste_couple_interdit[j][liste_couple_interdit[j].find(";")+1:len(liste_couple_interdit[j])] == petit:
                        valp = liste_couple_interdit[j][0:liste_couple_interdit[j].find(";")]
                        jpetit = j
                        liste_nb_pris.append(int(petit))
                    elif liste_couple_interdit[j][0:liste_couple_interdit[j].find(";")] == grand:
                        valg =  liste_couple_interdit[j][liste_couple_interdit[j].find(";")+1:len(liste_couple_interdit[j])]
                        jgrand = j
                        liste_nb_pris.append(int(grand))

                    elif liste_couple_interdit[j][liste_couple_interdit[j].find(";")+1:len(liste_couple_interdit[j])] == grand:
                        valg =  liste_couple_interdit[j][0:liste_couple_interdit[j].find(";")]
                        jgrand = j
                        liste_nb_pris.append(int(grand))
                if jgrand == -1 and jpetit == -1:
                    liste_couple_interdit.append(i[0])
                elif jgrand != -1 and jpetit != -1:
                    liste_couple_interdit.pop(jgrand)
                    if jgrand < jpetit:
                        liste_couple_interdit.pop(jpetit-1)
                    else:
                        liste_couple_interdit.pop(jpetit)
                    liste_couple_interdit.append(str(min(int(valp), int(valg))) + ";" + str(max(int(valp), int(valg))))
                elif jgrand != -1:
                    liste_couple_interdit.pop(jgrand)
                    liste_couple_interdit.append(str(min(int(valg), int(petit))) + ";" + str(max(int(valg), int(petit))))
                elif jpetit != -1:
                    liste_couple_interdit.pop(jpetit)
                    liste_couple_interdit.append(str(min(int(valp), int(grand))) + ";" + str(max(int(valp), int(grand))))


            else:
                liste_couple_interdit.append(i[0])
            liste_couple_select.append(i[0])
        if len(liste_nb_pris) == len(photo) - 2:
            break

    slide_show = []
    true_slide_show =[]
    compare = range(len(liste_nb_pris)+1)
    slide_show.append([(list(set(compare) - set(liste_nb_pris)))[0]])
    true_slide_show.append(liste_base[(list(set(compare) - set(liste_nb_pris)))[0]])
    while len(liste_couple_select) != 0:
        for i in range(len(liste_couple_select)):
            if liste_couple_select[i][0:liste_couple_select[i].find(';')] == str(slide_show[len(slide_show) - 1][0]):
                slide_show.append(
                    [liste_couple_select[i][liste_couple_select[i].find(';') + 1:len(liste_couple_select[i])]])
                true_slide_show.append(liste_base[int(slide_show[len(slide_show) - 1][0])])
                liste_couple_select.pop(i)
                break
            elif liste_couple_select[i][liste_couple_select[i].find(';') + 1:len(liste_couple_select[i])] == str(
                    slide_show[len(slide_show) - 1][0]):
                slide_show.append([liste_couple_select[i][0:liste_couple_select[i].find(';')]])
                true_slide_show.append(liste_base[int(slide_show[len(slide_show) - 1][0])])
                liste_couple_select.pop(i)
                break
    # print("Le slide show obtenu : " + str(true_slide_show))
    # print("Le score de ce slide show : " + str(calc_score(true_slide_show)))
    return true_slide_show


if __name__ == '__main__':
    debug = True

    # read_instance_file("ressources/a_example.txt",100)
    # read_instance_file("ressources/b_lovely_landscapes.txt", 1)
    read_instance_file("ressources/c_memorable_moments.txt", 5)

    if debug:
        print('Datas :')
        print(allPhotos)
        print(hPhotos)
        print(vPhotos)
        print('')

    simple_ss = simple_slide_show()
    write_slide_show(simple_ss)
    print(simple_ss)
    print('Score simple : ' + str(calc_score(simple_ss)))

    print()

    t1 = time.clock()
    glouton_ss = glouton_slide_show(False)
    t2 = time.clock()
    write_slide_show(glouton_ss, "output_glouton.txt")
    print(glouton_ss)
    print('Score glouton : ' + str(calc_score(glouton_ss)) + " en " + str(t2-t1) + "s")

    t1 = time.clock()
    glouton_ssf = glouton_slide_show()
    t2 = time.clock()
    write_slide_show(glouton_ssf, "output_glouton_fast.txt")
    print(glouton_ssf)
    print('Score glouton fast : ' + str(calc_score(glouton_ssf)) + " en " + str(t2 - t1) + "s")

    new_ss = descente_stochastique(glouton_ssf, 1000)
    print("Score Stochastique : " + str(calc_score(new_ss)) + "\nListe Stochastique" + str(new_ss))

    slides = simple_ss
    t1 = time.clock()
    plne_ss = plne_slide_show(slides)
    t2 = time.clock()
    write_slide_show(plne_ss, "output_plne.txt")
    print(plne_ss)
    print('Score plne : ' + str(calc_score(plne_ss)) + " en " + str(t2 - t1) + "s")

    t1 = time.clock()
    arrondi_ss = heuristique_arrondi(allPhotos,hPhotos,vPhotos)
    t2 = time.clock()
    write_slide_show(arrondi_ss, "output_arrondi.txt")
    print(arrondi_ss)
    print('Score arrondi : ' + str(calc_score(arrondi_ss)) + " en " + str(t2 - t1) + "s")
