       SUBROUTINE dwrpnorm(DTHETA, DMU, DSD, NSIZE, NMU, IK, dd) 
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     functions for dwrpnorm
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Venezia
C             30125 VENEZIA
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     March, 17, 2003
C
C     Version: 0.1-1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2003 Claudio Agostinelli
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C   This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PARAMETER:
C     NAME:     I/O:    TYPE:  DIMENSIONS:   DESCRIPTIONS:
C     DTHETA    input    D      NSIZE        vector of the data
C     DMU       input    D      1            mean
C     DSD       input    D      1            sd
C     NSIZE     input    I      1            length of the data
C     NMU       inout    I      1            length of mu 
C     IK        input    I      1            number of index to sum
C     dd        output   D      NSIZE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

      dimension dtheta(nsize), dmu(nmu), dd(nmu, nsize)

      do 10 i=1, nsize
          do 20 j=1, nmu
          dd(j, i) = dzero
 20   continue
 10   continue

      do 30 i=1, nsize
         do 40 j=1, nmu   
        dd(j,i) = dexp(-((dtheta(i) - dmu(j))**ddue)/(ddue * dsd**ddue))
            do 50 k=1,ik
              dd(j,i) = dd(j,i) 
     &        + dexp(-((dtheta(i) - dmu(j) + ddue*dpi*k)**ddue)
     &        /(ddue * dsd**ddue))
     &        + dexp(-((dtheta(i) - dmu(j) - ddue*dpi*k)**ddue)
     &        /(ddue * dsd**ddue))
 50         continue
 40      continue
 30   continue

      return
      end




















