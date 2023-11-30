#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  unbenannt.py
#
#  Copyright 2013 Kevin Li <kevin.shing.bruce.li@cern.ch>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#


from . import f90sussix


class F90Sussix():
    '''
        @class: Sussix
        Class to manage sussix calls
    '''
    def __init__(self):
        self.tunex = []
        self.tuney = []
        self.tunez = []
        self.amplitude = []
        self.phase = []
        self.ox = []
        self.ax = []
        self.oy = []
        self.ay = []
        self.oz = []
        self.az = []

    def sussix_inp(self, nt1=1, nt2=128, idam=2, ir=1,
                   tunex=0.2, tuney=0.2, istun=None, nline=0):
        si = open("sussix.inp", "w")
        si.write('C\nC INPUT FOR SUSSIX_V4 ---17/09/1997---\n')
        si.write('C DETAILS ARE IN THE MAIN PROGRAM SUSSIX_V4.F\nC\n\n')
        si.write('ISIX  = 0\n')
        si.write('NTOT  = 1\n')
        si.write('IANA  = 1\n')
        si.write('ICONV = 0\n')
        si.write('TURNS = '+str(nt1)+'  '+str(nt2)+'\n')
        si.write('NARM  = 200\n')
        if not istun:
            si.write('ISTUN = 1  0.05  0.05  0.05\n')
        else:
            si.write('ISTUN = {0}  {1}  {2}  {3}\n'.format(*istun))
        si.write('TUNES = '+str(tunex)+'  '+str(tuney)+' .001\n')
        si.write('NSUS  = 0\n')
        si.write('IDAM  = '+str(idam)+'\n')
        si.write('NTWIX = 1\n')
        si.write('IR    = '+str(ir) + '\n')
        si.write('IMETH = 2\n')
        si.write('NRC   = 6\n')
        si.write('EPS   = 1D-3\n')
        si.write('NLINE = {}\n'.format(nline))
        si.write('L,M,K =  \n')
        si.write('IDAMX = 1\n')
        si.write('NFIN  = 500\n')
        si.write('ISME  = 1\n')
        si.write('IUSME = 200\n')
        si.write('INV   = 0\n')
        si.write('IINV  = 250\n')
        si.write('ICF   = 0\n')
        si.write('IICF  = 350\n')
        si.close()


    def sussix(self, x, xp, y, yp, z, zp):
        z = f90sussix.sussixnoo(x, xp, y, yp, z, zp, len(x))
        self.tunex = z[0]
        self.tuney = z[1]
        self.tunez = z[2]
        self.amplitude = z[3]
        self.phase = z[4]
        self.ox = z[5]
        self.ax = z[6]
        self.oy = z[7]
        self.ay = z[8]
        self.oz = z[9]
        self.az = z[10]

        return z
