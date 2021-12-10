#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 22:14:33 2021

@author: jguterl
"""

f = open('/Users/jguterl/Dropbox/pyFACE/face_original/FACE/src/modFACE_header.f90','r')
lines = f.readlines()
out = []

lines = [l.strip() for l in lines]
for l in lines:
    l=l.strip()
    if l.count('::')>0:
        fval =''
        ftype = l.split('::')[0].split(',')[0]
        if l.split('::')[0].count(','):
            fopt =  l.split('::')[0].split(',')[1:]
        if ftype ==' real(DP)':
            ftype='real'
        
        r = l.split('::')[1]
        if r.count('!'):
             fcom = '#' + ''.join(r.split('!')[1:])
        else:
            fcom = '#'
        rr = r.split('!')[0]    
        if rr.count('='):
            fvar = rr.split('=')[0]
            fval = rr.split('=')[1]
        else:
            fvar = rr
            fval =''
               
        out.append('{} {} /{}/ {}'.format(fvar, ftype, fval,fcom))
    else:
        out.append(l)    
    
with open('parse_header.txt','w') as f:
    for l in out:
        f.write('{}\n'.format(l))