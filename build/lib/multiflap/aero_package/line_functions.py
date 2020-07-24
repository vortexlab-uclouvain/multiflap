'''
file line_functions.py
@author Victor Colognesi
@copyright Copyright © UCLouvain 2020

multiflap is a Python tool for finding periodic orbits and assess their stability via the Floquet multipliers.

Copyright <2020> <Université catholique de Louvain (UCLouvain), Belgique>

List of the contributors to the development of multiflap, Description and complete License: see LICENSE and NOTICE files.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
'''
import numpy as np
from .distances_functions import distance_point2line, distance_line2line

def quarterchord_naive(leadingedge, trailingedge):
    line = .75*leadingedge + .25*trailingedge
    chord_leadingedge = leadingedge
    chord_trailingedge = trailingedge
    return line, chord_leadingedge, chord_trailingedge


def smooth_line(line, npt):

    smooth_radius = 0.05

    # Getting the length of the curve
    notok = 1
    l = np.zeros((3, npt))
    whilecicle = 0
    while notok == 1:
        whilecicle = whilecicle + 1
        n = np.size(line[0])
        lcurve = 0
        for i in range(n-1):
            lcurve = lcurve + np.sqrt(np.sum((line[:,i+1]-line[:,i])**2))

        ds = lcurve/(npt-1)

        l[:,0] = line [:,0]

        j = 0  #Until here is OK

        for i in range(1,npt):
            dist = np.sqrt(np.sum((line[:,j+1]-l[:,i-1])**2))
            while dist < ds:
                j = j+1
                if j < (np.size(line[0])-1):
                    dist = dist + np.sqrt(np.sum((line[:,j+1]-line[:,j])**2))
                else:
                    dist = ds + 1

            if j < (np.size(line[0])-1):
                fac = dist - ds
                den = np.sqrt(np.sum((line[:,j+1]-line[:,j])**2))
                fac = fac/den
                l[:,i] = (1-fac)*line[:,j+1]+fac*line[:,j]
            else:
                l[:,i] = line [:,-1]

        l2 = l

        notok = 0

        for i in range(1,npt-1):
            vec1 = l[:,i] - l[:,i-1]
            vec2 = l[:,i+1] - l[:,i]
            coss = np.sum(vec1*(vec2))/((np.linalg.norm(vec1))*(np.linalg.norm(vec2)))
            # =============================================================================
            # since truncation errors might lead to having coss slightly greater than 1
            # =============================================================================
            if coss > (1):
                angle = np.arccos(1.0)
            else:
                angle = np.arccos(coss)
            d = np.linalg.norm(vec1)
            smooth_angle = d/smooth_radius
            while angle > smooth_angle:
                notok = 1
                point = (l[:,i+1] + l[:,i-1])/2
                l2[:,i] = (l2[:,i] + point)/2
                vec1 = l2[:,i] - l[:, i-1]
                vec2 = l[:,i+1] - l2[:,i]
                coss = np.sum(vec1*(vec2))/((np.linalg.norm(vec1))*(np.linalg.norm(vec2)))
                angle = np.arccos(coss)
                d = np.linalg.norm(vec1)
                smooth_angle = d/smooth_radius

        line = l2

    return line


def smooth_quarterline(oldline,oldchordle,oldchordte):

    npt = np.size((oldline[0]))
    lcurve = 0
    for i in range(npt-1):
        lcurve = lcurve + np.sqrt(np.sum((oldline[:,i+1]-oldline[:,i])**2))

    ds = lcurve/(npt-1)

    line = np.zeros_like(oldline)
    chord_leadingedge = np.zeros_like((oldline))
    chord_trailingedge = np.zeros_like((oldline))
    line[:,0] = 1*oldline[:,0]
    chord_leadingedge[:,0] = oldchordle[:,0]
    chord_trailingedge[:,0] = oldchordte[:,0]
    j = 0

    for i in range(1,npt):
        dist = np.sqrt(np.sum((oldline[:,j+1]-line[:,i-1])**2))
        while dist < ds:
            j = j+1
            if j < (np.size(oldline[0])-1):
                dist = dist + np.sqrt(np.sum((oldline[:,j+1]-oldline[:,j])**2))
            else:
                dist = ds + 1

        if j < (np.size(line[0])-1):
            fac = dist - ds
            den = np.sqrt(np.sum((oldline[:,j+1]-oldline[:,j])**2))
            fac = fac/den
            line[:,i] = (1-fac)*oldline[:,j+1]+fac*oldline[:,j]
            chord_leadingedge[:,i] = (1-fac)*oldchordle[:,j+1]+fac*oldchordle[:,j]
            chord_trailingedge[:,i] = (1-fac)*oldchordte[:,j+1]+fac*oldchordte[:,j]
        else:
            line[:,i] = 1*oldline [:,-1]
            chord_leadingedge[:,i] = oldchordle[:,-1]
            chord_trailingedge[:,i] = oldchordte[:,-1]

    return line, chord_leadingedge, chord_trailingedge


def improveline_new(line, chord_leadingedge, chord_trailingedge, leadingedge, trailingedge):
    n = np.size(line[0])
    nle = np.size(leadingedge[0])
    nte = np.size(trailingedge[0])

    newline = np.zeros_like(line)
    newchordle = np.zeros_like(chord_leadingedge)
    newchordte = np.zeros_like(chord_trailingedge)
    newline[:,0] = line[:,0]
    newline[:,-1] = line[:,-1]
    newchordle[:,0] = chord_leadingedge[:,0]
    newchordle[:,-1] = chord_leadingedge[:,-1]
    newchordte[:,0] = chord_trailingedge[:,0]
    newchordte[:,-1] = chord_trailingedge[:,-1]

    # External Iteration

    for i in range(1,n-1):
        linevec = line[:,i+1] - line[:,i-1]
        linevec = linevec/np.linalg.norm(linevec)

        # orthogonalizing chordvec wrt linevec
        chordvec = chord_trailingedge[:,i] - chord_leadingedge[:,i]
        chordvec = chordvec - np.dot((np.sum(linevec*chordvec)),linevec)    # BUG HERE
        chordvec = chordvec/np.linalg.norm(chordvec)

        """
        getting points in le closest to chord line
        take the first two points of le as first references
        """
        jref1 = 0
        jref2 = 1

        [_ , _, dref1] = distance_point2line(leadingedge[:,jref1], line[:,i], chordvec) # If it does not work DP2L[][2] LAST ELEMENT OF OUTPUT
        [_, _, dref2] = distance_point2line(leadingedge[:,jref2], line[:,i], chordvec)


        """
        First IF Routine, checked (Below)
        """
        if (dref2<dref1):
            a = jref2
            jref2 = jref1
            jref1 = a
            a = dref2
            dref2 = dref1
            dref1 = a

        for j in range(2, nle):
            [_,_,dtest] = distance_point2line(leadingedge[:,j], line[:,i], chordvec)
            if (dtest < dref2):
                dref2 = dtest
                jref2 = j
            if (dref2 < dref1):
                a = jref2
                jref2 = jref1
                jref1 = a
                a = dref2
                dref2 = dref1
                dref1 = a

        levec = leadingedge[:,jref2] - leadingedge[:,jref1]
        levec = levec/np.linalg.norm(levec)
        [_, _, newchordle[:,i], _, _] = distance_line2line(leadingedge[:,jref1],levec,line[:,i],chordvec)

        # Repeat this for the trailing edge

        jref1 = 0
        jref2 = 1

        [_, _, dref1] = distance_point2line(trailingedge[:,jref1], line[:,i], chordvec) # If it does not work DP2L[][2] LAST ELEMENT OF OUTPUT
        [_, _, dref2] = distance_point2line(trailingedge[:,jref2], line[:,i], chordvec)

        if dref2 < dref1:
            a = jref2
            jref2 = jref1
            jref1 = a
            a = dref2
            dref2 = dref1
            dref1 = a

        for j in range(2, nte):

            [_,_,dtest] = distance_point2line(trailingedge[:,j], line[:,i], chordvec)
            if (dtest < dref2):
                dref2 = dtest
                jref2 = j
            if (dref2 < dref1):
                a = jref2
                jref2 = jref1
                jref1 = a
                a = dref2
                dref2 = dref1
                dref1 = a

        tevec = trailingedge[:,jref2] - trailingedge[:,jref1]
        tevec = tevec/np.linalg.norm(tevec)
        [_, _, newchordte[:,i], _, _] = distance_line2line(trailingedge[:,jref1],tevec,line[:,i],chordvec)
        newline[:,i] = .75*newchordle[:,i] + .25*newchordte[:,i]

    return(newline, newchordle, newchordte)

def improveline_iteration(line, chord_leadingedge, chord_trailingedge, leadingedge, trailingedge, it):

    relfactor = .25

    [line, chord_leadingedge, chord_trailingedge] = smooth_quarterline(line,chord_leadingedge,chord_trailingedge)
    [newline, newchordle, newchordte] = improveline_new(line, chord_leadingedge, chord_trailingedge,leadingedge,trailingedge)

    chgt = 0

    for i in range(np.size(line[0]) - 1):

        local_max = np.maximum(0, i-1)
        local_min = np.minimum(np.size(line[0]), i + 1)
        norm1 = np.linalg.norm(newline[:,i] - line[:,i])
        norm2 = np.linalg.norm(line[:,local_max]-line[:,local_min])
        chgt = chgt + relfactor*norm1/norm2

    chgt = chgt/np.size(line[0])

    line = relfactor*newline + (1-relfactor)*line
    chord_leadingedge = relfactor*newchordle + (1-relfactor)*chord_leadingedge
    chord_trailingedge = relfactor*newchordte + (1-relfactor)*chord_trailingedge

    return line, chord_leadingedge, chord_trailingedge, chgt
