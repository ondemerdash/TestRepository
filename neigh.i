c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  neigh.i  --  pairwise neighbor list indices and storage  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     lbuffer     width of the neighbor list buffer region
c     lbuf2       square of half the neighbor list buffer width
c     vbuf2       square of vdw cutoff plus neighbor list buffer
c     cbuf2       square of charge cutoff plus neighbor list buffer
c     mbuf2       square of multipole cutoff plus neighbor list buffer
c     vbufx       square of vdw cutoff plus twice the list buffer
c     cbufx       square of charge cutoff plus twice the list buffer
c     mbufx       square of multipole cutoff plus twice the list buffer
c     xvold       x-coordinate at last vdw neighbor list update
c     yvold       y-coordinate at last vdw neighbor list update
c     zvold       z-coordinate at last vdw neighbor list update
c     xcold       x-coordinate at last charge neighbor list update
c     ycold       y-coordinate at last charge neighbor list update
c     zcold       z-coordinate at last charge neighbor list update
c     xmold       x-coordinate at last multipole neighbor list update
c     ymold       y-coordinate at last multipole neighbor list update
c     zmold       z-coordinate at last multipole neighbor list update
c     nvlst       number of sites in list for each vdw site
c     vlst        site numbers in neighbor list of each vdw site
c     nelst       number of sites in list for each electrostatic site
c     elst        site numbers in list of each electrostatic site
c     dovlst      logical flag to rebuild vdw neighbor list
c     doclst      logical flag to rebuild charge neighbor list
c     domlst      logical flag to rebuild multipole neighbor list
c     domollst2bod
c     domollst3bod

      integer, pointer :: nvlst(:)
      integer, pointer :: vlst(:,:)
      integer, pointer :: nelst(:)
      integer, pointer :: elst(:,:)
      integer, pointer :: nmollst(:)
      integer, pointer :: mollst(:,:)
      integer, pointer :: nmollst3(:,:)
      integer, pointer :: mollst3(:,:,:)
      integer, pointer :: nmollst3mod(:)
      integer, pointer :: mollst3mod(:,:)
      real*8 lbuffer,lbuf2,molbuf2
      real*8 vbuf2,cbuf2,mbuf2,molbufx
      real*8 vbufx,cbufx,mbufx
      real*8, pointer :: xvold(:)
      real*8, pointer :: yvold(:)
      real*8, pointer :: zvold(:)
      real*8, pointer :: xcold(:)
      real*8, pointer :: ycold(:)
      real*8, pointer :: zcold(:)
      real*8, pointer :: xmold(:)
      real*8, pointer :: ymold(:)
      real*8, pointer :: zmold(:)
      real*8, pointer :: xmolold(:)
      real*8, pointer :: ymolold(:)
      real*8, pointer :: zmolold(:)
      real*8, pointer :: xmolold3(:)
      real*8, pointer :: ymolold3(:)
      real*8, pointer :: zmolold3(:)
      logical dovlst,doclst,domlst
      logical domollst2bod,domollst3bod
      common /neigh/ lbuffer,lbuf2,vbuf2,cbuf2,mbuf2,vbufx,cbufx,mbufx,
     &               xvold,yvold,zvold,xcold,ycold,zcold,xmold,ymold,
     &               zmold,nvlst,vlst,nelst,elst,dovlst,doclst,domlst,
     &               domollst2bod,domollst3bod,nmollst,mollst,xmolold,
     &               ymolold,zmolold,xmolold3,ymolold3,zmolold3,molbuf2,
     &               molbufx,nmollst3,mollst3,nmollst3mod,mollst3mod
