       SUBROUTINE mlewrpno(DTHETA, DMU, DSD, NSIZE, IK, IM, IR, 
     &  dw, dwk, dwm) 
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     functions for mlewrpno
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Venezia
C             30125 VENEZIA
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     March, 04, 2003
C
C     Version: 0.1
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
C     IK        input    I      1            number of index to sum
C     IM        input    I      1            if 1 you are estimating the mean
C     IR        input    I      1            if 1 you are estimating the sd
C     dw        output   D      NSIZE
C     dwk       output   D      NSIZE      
C     dwm       output   D      NSIZE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

      dimension dtheta(nsize), dw(nsize), dwk(nsize), dwm(nsize)

      do 10 i=1, nsize
          dw(i) = dzero
          dwk(i) = dzero
          dwm(i) = dzero
 10   continue

      do 20 i=1, nsize   
          dw(i) = dexp(-((dtheta(i) - dmu)**ddue)/(ddue * dsd**ddue))
          dwm(i) = (dtheta(i) - dmu)**ddue * 
     &             dexp(-((dtheta(i) - dmu)**ddue)/(ddue * dsd**ddue))
          do 30 j=1,ik
              dw(i) = dw(i) 
     &        + dexp(-((dtheta(i) - dmu + ddue*dpi*j)**ddue)
     &        /(ddue * dsd**ddue))
     &        + dexp(-((dtheta(i) - dmu - ddue*dpi*j)**ddue)
     &        /(ddue * dsd**ddue))

              if (im.eq.1) then
                 dwk(i) = dwk(i) 
     &           + j*dexp(-((dtheta(i) - dmu + ddue*dpi*j)**ddue)
     &           /(ddue * dsd**ddue))
     &           - j*dexp(-((dtheta(i) - dmu - ddue*dpi*j)**ddue)
     &           /(ddue * dsd**ddue))
              endif

              if (ir.eq.1) then
                 dwm(i) = dwm(i)
     &           + (dtheta(i) - dmu + ddue*j*dpi)**ddue
     &           * dexp(-((dtheta(i) - dmu + ddue*dpi*j)**ddue)
     &           /(ddue * dsd**ddue))
     &           + (dtheta(i) - dmu - ddue*j*dpi)**ddue
     &           * dexp(-((dtheta(i) - dmu - ddue*dpi*j)**ddue)
     &           /(ddue * dsd**ddue))
              endif
 30       continue
 20   continue

      return
      end




















