c
c
c     ###############################################################
c     ##  COPYRIGHT (C) 2006 by David Gohara & Jay William Ponder  ##
c     ##                    All Rights Reserved                    ##
c     ###############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nblist  --  maintain pairwise neighbor lists  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nblist" constructs and maintains nonbonded pair neighbor lists
c     for vdw and electrostatic interactions
c
c
      subroutine nblist
      implicit none
      include 'cutoff.i'
      include 'potent.i'
c
c
c     update the vdw and electrostatic neighbor lists
c
      if (use_vdw .and. use_vlist)  call vlist
      if (use_charge .and. use_clist)  call clist
      if ((use_mpole.or.use_polar) .and. use_mlist)  call mlist
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine vlist  --  build van der Waals neighbor lists  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "vlist" performs an update or a complete rebuild of the
c     van der Waals neighbor list
c
c
      subroutine vlist
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'neigh.i'
      include 'vdw.i'
      integer i,j,k
      integer ii,iv
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius
      real*8 rdn,r2
      real*8, allocatable :: xred(:)
      real*8, allocatable :: yred(:)
      real*8, allocatable :: zred(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xred(n))
      allocate (yred(n))
      allocate (zred(n))
c
c     apply reduction factors to find coordinates for each site
c
      do i = 1, nvdw
         ii = ivdw(i)
         iv = ired(ii)
         rdn = kred(ii)
         xred(i) = rdn*(x(ii)-x(iv)) + x(iv)
         yred(i) = rdn*(y(ii)-y(iv)) + y(iv)
         zred(i) = rdn*(z(ii)-z(iv)) + z(iv)
      end do
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(vbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' VLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (dovlst) then
         dovlst = .false.
         if (octahedron) then
            do i = 1, nvdw
               call vbuild (i,xred,yred,zred)
            end do
         else
            call vfull (xred,yred,zred)
         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, nvdw
         xi = xred(i)
         yi = yred(i)
         zi = zred(i)
         xr = xi - xvold(i)
         yr = yi - yvold(i)
         zr = zi - zvold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .ge. lbuf2) then
c
c     store coordinates to reflect update of this site
c
            xvold(i) = xi
            yvold(i) = yi
            zvold(i) = zi
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, nvdw
               xr = xi - xvold(k)
               yr = yi - yvold(k)
               zr = zi - zvold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. vbuf2) then
                  j = j + 1
                  vlst(j,i) = k
               end if
            end do
            nvlst(i) = j
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xi - xvold(k)
               yr = yi - yvold(k)
               zr = zi - zvold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. vbuf2) then
                  do j = 1, nvlst(k)
                     if (vlst(j,k) .eq. i)  goto 20
                  end do
                  nvlst(k) = nvlst(k) + 1
                  vlst(nvlst(k),k) = i
   20             continue
               else if (r2 .le. vbufx) then
                  do j = 1, nvlst(k)
                     if (vlst(j,k) .eq. i) then
                        vlst(j,k) = vlst(nvlst(k),k)
                        nvlst(k) = nvlst(k) - 1
                        goto 30
                     end if
                  end do
   30             continue
               end if
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xred)
      deallocate (yred)
      deallocate (zred)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine vbuild  --  make vdw pair list for one site  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "vbuild" performs a complete rebuild of the van der Waals
c     pair neighbor list for a single site
c
c
      subroutine vbuild (i,xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'iounit.i'
      include 'neigh.i'
      include 'vdw.i'
      integer i,j,k
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
      real*8 xred(*)
      real*8 yred(*)
      real*8 zred(*)
c
c
c     store coordinates to reflect update of the site
c
      xi = xred(i)
      yi = yred(i)
      zi = zred(i)
      xvold(i) = xi
      yvold(i) = yi
      zvold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
      j = 0
      do k = i+1, nvdw
         xr = xi - xred(k)
         yr = yi - yred(k)
         zr = zi - zred(k)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .le. vbuf2) then
            j = j + 1
            vlst(j,i) = k
         end if
      end do
      nvlst(i) = j
c
c     check to see if the neighbor list is too long
c
      if (nvlst(i) .ge. maxvlst) then
         write (iout,10)
   10    format (/,' VBUILD  --  Too many Neighbors;',
     &              ' Increase MAXVLST')
         call fatal
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine vfull  --  make vdw pair list for all sites  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "vfull" performs a complete rebuild of the van der Waals
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine vfull (xred,yred,zred)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'neigh.i'
      include 'vdw.i'
      integer i,j,k
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8 xred(*)
      real*8 yred(*)
      real*8 zred(*)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      nlight = (ncell+1) * nvdw
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nvdw
         nvlst(i) = 0
         xvold(i) = xred(i)
         yvold(i) = yred(i)
         zvold(i) = zred(i)
         xsort(i) = xred(i)
         ysort(i) = yred(i)
         zsort(i) = zred(i)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(vbuf2)
      call lights (off,nvdw,xsort,ysort,zsort)
c
c     loop over all atoms computing the interactions
c
      do i = 1, nvdw
         xi = xsort(rgx(i))
         yi = ysort(rgy(i))
         zi = zsort(rgz(i))
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - xsort(j)
            yr = yi - ysort(kgy)
            zr = zi - zsort(kgz)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. vbuf2) then
               if (i .lt. k) then
                  nvlst(i) = nvlst(i) + 1
                  vlst(nvlst(i),i) = k
               else
                  nvlst(k) = nvlst(k) + 1
                  vlst(nvlst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nvdw
         if (nvlst(i) .ge. maxvlst) then
            write (iout,30)
   30       format (/,' VFULL  --  Too many Neighbors;',
     &                 ' Increase MAXVLST')
            call fatal
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine clist  --  build partial charge neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "clist" performs an update or a complete rebuild of the
c     electrostatic neighbor lists for partial charges
c
c
      subroutine clist
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'charge.i'
      include 'iounit.i'
      include 'neigh.i'
      integer i,j,k,ii
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
c
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(cbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' CLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (doclst) then
         doclst = .false.
         if (octahedron) then
            do i = 1, nion
               call cbuild (i)
            end do
         else
            call cfull
         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, nion
         ii = kion(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xcold(i)
         yr = yi - ycold(i)
         zr = zi - zcold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .ge. lbuf2) then
c
c     store coordinates to reflect update of this site
c
            xcold(i) = xi
            ycold(i) = yi
            zcold(i) = zi
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, nion
               xr = xi - xcold(k)
               yr = yi - ycold(k)
               zr = zi - zcold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cbuf2) then
                  j = j + 1
                  elst(j,i) = k
               end if
            end do
            nelst(i) = j
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xi - xcold(k)
               yr = yi - ycold(k)
               zr = zi - zcold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. cbuf2) then
                  do j = 1, nelst(k)
                     if (elst(j,k) .eq. i)  goto 20
                  end do
                  nelst(k) = nelst(k) + 1
                  elst(nelst(k),k) = i
   20             continue
               else if (r2 .le. cbufx) then
                  do j = 1, nelst(k)
                     if (elst(j,k) .eq. i) then
                        elst(j,k) = elst(nelst(k),k)
                        nelst(k) = nelst(k) - 1
                        goto 30
                     end if
                  end do
   30             continue
               end if
            end do
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine cbuild  --  make charge pair list for one site  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "cbuild" performs a complete rebuild of the partial charge
c     electrostatic neighbor list for a single site
c
c
      subroutine cbuild (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'charge.i'
      include 'iounit.i'
      include 'neigh.i'
      integer i,j,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     store new coordinates to reflect update of the site
c
      ii = kion(i)
      xi = x(ii)
      yi = y(ii)
      zi = z(ii)
      xcold(i) = xi
      ycold(i) = yi
      zcold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
      j = 0
      do k = i+1, nion
         kk = kion(k)
         xr = xi - x(kk)
         yr = yi - y(kk)
         zr = zi - z(kk)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .le. cbuf2) then
            j = j + 1
            elst(j,i) = k
         end if
      end do
      nelst(i) = j
c
c     check to see if the neighbor list is too long
c
      if (nelst(i) .ge. maxelst) then
         write (iout,10)
   10    format (/,' CBUILD  --  Too many Neighbors;',
     &              ' Increase MAXELST')
         call fatal
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine cfull  --  make charge pair list for all sites  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "cfull" performs a complete rebuild of the partial charge
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine cfull
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'charge.i'
      include 'iounit.i'
      include 'light.i'
      include 'neigh.i'
      integer i,j,k,ii
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      nlight = (ncell+1) * nion
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nion
         nelst(i) = 0
         ii = kion(i)
         xcold(i) = x(ii)
         ycold(i) = y(ii)
         zcold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(cbuf2)
      call lights (off,nion,xsort,ysort,zsort)
c
c     loop over all atoms computing the interactions
c
      do i = 1, nion
         xi = xsort(rgx(i))
         yi = ysort(rgy(i))
         zi = zsort(rgz(i))
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - xsort(j)
            yr = yi - ysort(kgy)
            zr = zi - zsort(kgz)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. cbuf2) then
               if (i .lt. k) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
               else
                  nelst(k) = nelst(k) + 1
                  elst(nelst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nion
         if (nelst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' CFULL  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mlist  --  build atom multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mlist" performs an update or a complete rebuild of the
c     electrostatic neighbor lists for atomic multipoles
c
c
      subroutine mlist
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      integer i,j,k,ii
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 radius,r2
c
c
c     neighbor list cannot be used with the replicates method
c
      radius = sqrt(mbuf2)
      call replica (radius)
      if (use_replica) then
         write (iout,10)
   10    format (/,' MLIST  --  Pairwise Neighbor List cannot',
     &              ' be used with Replicas')
         call fatal
      end if
c
c     perform a complete list build instead of an update
c
      if (domlst) then
         domlst = .false.
         if (octahedron) then
            do i = 1, npole
               call mbuild (i)
            end do
         else
            call mfull
         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, npole
         ii = ipole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
         xr = xi - xmold(i)
         yr = yi - ymold(i)
         zr = zi - zmold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .ge. lbuf2) then
c
c     store coordinates to reflect update of this site
c
            xmold(i) = xi
            ymold(i) = yi
            zmold(i) = zi
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, npole
               xr = xi - xmold(k)
               yr = yi - ymold(k)
               zr = zi - zmold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. mbuf2) then
                  j = j + 1
                  elst(j,i) = k
               end if
            end do
            nelst(i) = j
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xi - xmold(k)
               yr = yi - ymold(k)
               zr = zi - zmold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. mbuf2) then
                  do j = 1, nelst(k)
                     if (elst(j,k) .eq. i)  goto 20
                  end do
                  nelst(k) = nelst(k) + 1
                  elst(nelst(k),k) = i
   20             continue
               else if (r2 .le. mbufx) then
                  do j = 1, nelst(k)
                     if (elst(j,k) .eq. i) then
                        elst(j,k) = elst(nelst(k),k)
                        nelst(k) = nelst(k) - 1
                        goto 30
                     end if
                  end do
   30             continue
               end if
            end do
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mbuild  --  make mpole pair list for one site  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mbuild" performs a complete rebuild of the atomic multipole
c     electrostatic neighbor list for a single site
c
c
      subroutine mbuild (i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      integer i,j,k,ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
c
c
c     store new coordinates to reflect update of the site
c
      ii = ipole(i)
      xi = x(ii)
      yi = y(ii)
      zi = z(ii)
      xmold(i) = xi
      ymold(i) = yi
      zmold(i) = zi
c
c     generate all neighbors for the site being rebuilt
c
      j = 0
      do k = i+1, npole
         kk = ipole(k)
         xr = xi - x(kk)
         yr = yi - y(kk)
         zr = zi - z(kk)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .le. mbuf2) then
            j = j + 1
            elst(j,i) = k
         end if
      end do
      nelst(i) = j
c
c     check to see if the neighbor list is too long
c
      if (nelst(i) .ge. maxelst) then
         write (iout,10)
   10    format (/,' MBUILD  --  Too many Neighbors;',
     &              ' Increase MAXELST')
         call fatal
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mfull  --  make mpole pair list for all sites  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mfull" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine mfull
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      integer i,j,k,ii
      integer kgy,kgz
      integer start,stop
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r2,off
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      nlight = (ncell+1) * npole
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, npole
         nelst(i) = 0
         ii = ipole(i)
         xmold(i) = x(ii)
         ymold(i) = y(ii)
         zmold(i) = z(ii)
         xsort(i) = x(ii)
         ysort(i) = y(ii)
         zsort(i) = z(ii)
      end do
c
c     use the method of lights to generate neighbors
c
      off = sqrt(mbuf2)
      call lights (off,npole,xsort,ysort,zsort)
c
c     loop over all atoms computing the interactions
c
      do i = 1, npole
         xi = xsort(rgx(i))
         yi = ysort(rgy(i))
         zi = zsort(rgz(i))
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - xsort(j)
            yr = yi - ysort(kgy)
            zr = zi - zsort(kgz)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. mbuf2) then
               if (i .lt. k) then
                  nelst(i) = nelst(i) + 1
                  elst(nelst(i),i) = k
               else
                  nelst(k) = nelst(k) + 1
                  elst(nelst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
      do i = 1, npole
         if (nelst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' MFULL  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine imagen  --  neighbor minimum image distance  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "imagen" takes the components of pairwise distance between
c     two points and converts to the components of the minimum
c     image distance; fast version for neighbor list generation
c
c
      subroutine imagen (xr,yr,zr)
      implicit none
      include 'boxes.i'
      real*8 xr,yr,zr
c
c
c     for orthogonal lattice, find the desired image directly;
c     to save time, this only returns the correct magnitudes
c
      if (orthogonal) then
         xr = abs(xr)
         yr = abs(yr)
         zr = abs(zr)
         if (xr .gt. xbox2)  xr = xr - xbox
         if (yr .gt. ybox2)  yr = yr - ybox
         if (zr .gt. zbox2)  zr = zr - zbox
c
c     for monoclinic lattice, convert "xr" and "zr" specially
c
      else if (monoclinic) then
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, use general conversion equations
c
      else if (triclinic) then
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, remove the corner pieces
c
      else if (octahedron) then
         if (abs(xr) .gt. xbox2)  xr = xr - sign(xbox,xr)
         if (abs(yr) .gt. ybox2)  yr = yr - sign(ybox,yr)
         if (abs(zr) .gt. zbox2)  zr = zr - sign(zbox,zr)
         if (abs(xr)+abs(yr)+abs(zr) .gt. box34) then
            xr = xr - sign(xbox2,xr)
            yr = yr - sign(ybox2,yr)
            zr = zr - sign(zbox2,zr)
         end if
      end if
      return
      end


c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull2body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mfull" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull2body
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii
      integer kgy,kgz,l1,lenmol
      integer start,stop
      real*8 xi,yi,zi,M1
      real*8 xr,yr,zr
      real*8 r2,off,xcm1,ycm1,zcm1
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
c      print*, "Before allocation Molfull2body"

      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
      domollst2bod = .false.
c      print*, "After allocation Molfull2body"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nmol
         nmollst(i) = 0
c         ii = ipole(i)
         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
c             j = pnum(l1)
            k=l1-1
            j=imol(1,i)+k                
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1
         
         xmolold(i) = xcm1
         ymolold(i) = ycm1
         zmolold(i) = zcm1
         xsort(i) = xcm1
         ysort(i) = ycm1
         zsort(i) = zcm1
      end do
c
c     use the method of lights to generate neighbors
c
c      print*,"Before lights Call in Molfull2body"
      molbuf2=25.0d0
c      off = sqrt(molbuf2)
      off=molbuf2
      off=10.0d0
      call lights (off,nmol,xsort,ysort,zsort)

c      print*,"After lights Call in Molfull2body"
c
c     loop over all atoms computing the interactions
c
      do i = 1, nmol
         xi = xsort(rgx(i))
         yi = ysort(rgy(i))
         zi = zsort(rgz(i))
         if (kbx(i) .le. kex(i)) then
            repeat = .false.
            start = kbx(i) + 1
            stop = kex(i)
         else
            repeat = .true.
            start = 1
            stop = kex(i)
         end if
   10    continue
         do j = start, stop
            k = locx(j)
            kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if
            xr = xi - xsort(j)
            yr = yi - ysort(kgy)
            zr = zi - zsort(kgz)
            call imagen (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. molbuf2) then
               if (i .lt. k) then
                  nmollst(i) = nmollst(i) + 1
                  mollst(nmollst(i),i) = k
               else
                  nmollst(k) = nmollst(k) + 1
                  mollst(nmollst(k),k) = i
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(i) + 1
            stop = nlight
            goto 10
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
c      do i = 1, nmol
c         if (nmollst(i) .ge. maxelst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')
c            call fatal
c         end if
c      end do
      return
      end


c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull3body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "molfull3body" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull3body
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol
      integer start,stop
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,xcm1,ycm1,zcm1
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      print*, "Before allocation Molfull3body"

      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))

      print*, "After allocation Molfull3body"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nmol-1
         do j=i+1,nmol
         nmollst3(i,j) = 0
         end do
      end do 
c         ii = ipole(i)
      do i = 1, nmol
         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
c             j = pnum(l1)
            k=l1-1
            j=imol(1,i)+k                
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1
         
         xmolold(i) = xcm1
         ymolold(i) = ycm1
         zmolold(i) = zcm1
         xsort(i) = xcm1
         ysort(i) = ycm1
         zsort(i) = zcm1
      end do
c
c     use the method of lights to generate neighbors
      print*,"Before lights Call in Molfull3body"
c      molbuf2=4.25+3
      molbuf2=5.5+2
      off = sqrt(molbuf2)
      call lights (off,nmol,xsort,ysort,zsort)

      print*,"After lights Call in Molfull3body"
c
c     loop over all atoms computing the interactions
c
      do i = 1, nmol-1
         print*, "In loop i=",i
         do l1=i+1,nmol
            xi = xsort(rgx(i))
            yi = ysort(rgy(i))
            zi = zsort(rgz(i))
      
            xl1 = xsort(rgx(l1))
            yl1 = ysort(rgy(l1))
            zl1 = zsort(rgz(l1)) 
        
            if (kbx(l1) .le. kex(l1)) then
               repeat = .false.
               start = kbx(l1) + 1
c               start = kbx(i) + 2
               stop = kex(l1)
            else
               repeat = .true.
               start = 1
c                start = 2
              stop = kex(l1)
            end if
            print*, "In inner loop i=",i," l1=",l1
   10    continue

            do j = start, stop
            if(j .ne. i) then
               k = locx(j)
               kgy = rgy(k)
            if (kby(l1) .le. key(l1)) then
               if (kgy.lt.kby(l1) .or. kgy.gt.key(l1))  goto 20
            else
               if (kgy.lt.kby(l1) .and. kgy.gt.key(l1))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(l1) .le. kez(l1)) then
               if (kgz.lt.kbz(l1) .or. kgz.gt.kez(l1))  goto 20
            else
               if (kgz.lt.kbz(l1) .and. kgz.gt.kez(l1))  goto 20
            end if

               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
           
               call imagen (xr,yr,zr)

               xr1=xr
               yr1=yr
               zr1=zr
             
               xr = xl1 - xsort(j)
               yr = yl1 - ysort(kgy)
               zr = zl1 - zsort(kgz)
               
               call imagen (xr,yr,zr)
          
               xr2 = xr
               yr2 = yr
               zr2 = zr

               xr = xi - xl1
               yr = yi - yl1
               zr = zi - zl1
               
               call imagen(xr,yr,zr)
              
               xr3=xr
               yr3=yr
               zr3=zr

               r2 = (xr1*xr1 + yr1*yr1 + zr1*zr1 +
     &            xr2*xr2 + yr2*yr2 + zr2*zr2 +
     &            xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
               if (r2 .le. molbuf2) then
                  if (l1 .lt. k) then
                     nmollst3(i,l1) = nmollst3(i,l1) + 1
                     mollst3(nmollst3(i,l1),i,l1) = k
                  else
                     nmollst3(i,k) = nmollst3(i,k) + 1
                     mollst3(nmollst3(i,k),i,k) = l1
                  end if
               end if
c   20       continue
            else
             goto 20
            end if
   20       continue
            print*, "In inner loop i=",i," l1=",l1," j=",j
            end do

            if (repeat) then
               repeat = .false.
               start = kbx(l1) + 1
c               start = kbx(i) +2
               stop = nlight
               goto 10
            end if

         end do
      end do

c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nmol
         if (nmollst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' MFULL  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
      return
      end


c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull3body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "molfull3body" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull3bodymod2
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol
      integer start,stop,start2,stop2
      integer k1,kgy1,kgz1,j1,done(nmol,nmol,nmol)
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,xcm1,ycm1,zcm1
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat,repeat2
c
c
c     perform dynamic allocation of some local arrays
c
      print*, "Before allocation Molfull3bodymod2"

      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))

      print*, "After allocation Molfull3bodymod2"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nmol
         do j=1,nmol
         nmollst3(i,j) = 0
         end do
      end do 

      do i=1,nmol
         do j=1,nmol
            do k=1,nmol
               done(i,j,k)=0
            end do
         end do 
      end do 

      do i = 1, nmol
         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k=l1-1
            j=imol(1,i)+k                
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1
         
         xmolold(i) = xcm1
         ymolold(i) = ycm1
         zmolold(i) = zcm1
         xsort(i) = xcm1
         ysort(i) = ycm1
         zsort(i) = zcm1
      end do
c
c     use the method of lights to generate neighbors
      print*,"Before lights Call in Molfull3bodymod2"
c      molbuf2=4.25+3
      molbuf2=5.5+2
      off = sqrt(molbuf2)
      call lights (off,nmol,xsort,ysort,zsort)

      print*,"After lights Call in Molfull3bodymod2"
c
c     loop over all atoms computing the interactions
c
      do i = 1, nmol-1
         print*, "In loop i=",i
         do l1=i+1,nmol
            xi = xsort(rgx(i))
            yi = ysort(rgy(i))
            zi = zsort(rgz(i))
      
            xl1 = xsort(rgx(l1))
            yl1 = ysort(rgy(l1))
            zl1 = zsort(rgz(l1)) 
        
            if (kbx(i) .le. kex(i)) then
               repeat = .false.
               start = kbx(i) + 1
               stop = kex(i)
            else
               repeat = .true.
               start = 1
              stop = kex(i)
            end if
            print*, "In inner loop i=",i," l1=",l1
   10    continue

            do j = start, stop
c            if(j .ne. i) then
               k = locx(j)
               kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if

               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
           
               call imagen (xr,yr,zr)

               xr1=xr
               yr1=yr
               zr1=zr

               xr = xl1 - xsort(j)
               yr = yl1 - ysort(kgy)
               zr = zl1 - zsort(kgz)

               call imagen (xr,yr,zr)

               xr2 = xr
               yr2 = yr
               zr2 = zr

               xr = xi - xl1
               yr = yi - yl1
               zr = zi - zl1

               call imagen(xr,yr,zr)

               xr3=xr
               yr3=yr
               zr3=zr

               r2 = (xr1*xr1 + yr1*yr1 + zr1*zr1 +
     &            xr2*xr2 + yr2*yr2 + zr2*zr2 +
     &            xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
               if (r2 .le. molbuf2 .and. done(i,l1,k).eq.0
     &            .and.done(k,l1,i).eq.0 .and. k.ne.l1)then
                  if (i.lt.k) then
                     nmollst3(i,l1) = nmollst3(i,l1) + 1
                     mollst3(nmollst3(i,l1),i,l1) = k
                     done(i,l1,k)=1
                     done(k,l1,i)=1
                  else
                     nmollst3(k,l1) = nmollst3(k,l1) + 1
                     mollst3(nmollst3(k,l1),k,l1) = i
                     done(i,l1,k)=1
                     done(k,l1,i)=1 
                  end if
               end if

               if (kbx(l1) .le. kex(l1)) then
                  repeat2 = .false.
                  start2 = kbx(l1) + 1
                  stop2 = kex(l1)
               else
                  repeat2 = .true.
                  start2 = 1
                  stop2 = kex(l1)
               end if

   30    continue

               do j1 = start2, stop2
                  
                  k1=locx(j1)
                  kgy1 =rgy(k1)
                  if (kby(l1) .le. key(l1)) then
                   if(kgy1.lt.kby(l1)
     &                .or. kgy1.gt.key(l1)) goto 40
                  else
                   if (kgy1.lt.kby(l1) 
     &               .and. kgy1.gt.key(l1)) goto 40                    
                  end if
                  kgz1=rgz(k1)
                  if (kbz(l1) .le. kez(l1)) then
                   if (kgz1.lt.kbz(l1)
     &                 .or. kgz1.gt.kez(l1))  goto 40
                  else
                   if (kgz1.lt.kbz(l1) 
     &                 .and. kgz1.gt.kez(l1))  goto 40
                  end if
                  
                  xr = xi - xsort(j1)
                  yr = yi - ysort(kgy1)
                  zr = zi - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr1=xr
                  yr1=yr
                  zr1=zr

                  xr = xl1 - xsort(j1)
                  yr = yl1 - ysort(kgy1)
                  zr = zl1 - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr2 = xr
                  yr2 = yr
                  zr2 = zr
                      

                  r2 = (xr1*xr1 + yr1*yr1 + zr1*zr1 +
     &                 xr2*xr2 + yr2*yr2 + zr2*zr2 +
     &                xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
               if (r2 .le. molbuf2 .and. done(i,l1,k1).eq.0
     &            .and.done(i,k1,l1).eq.0 .and. k1.ne.i) then
                  if (l1 .lt. k1) then
                     nmollst3(i,l1) = nmollst3(i,l1) + 1
                     mollst3(nmollst3(i,l1),i,l1) = k1
                     done(i,l1,k1)=1
                     done(i,k1,l1)=1
                  else
                     nmollst3(i,k1) = nmollst3(i,k1) + 1
                     mollst3(nmollst3(i,k1),i,k1) = l1
                     done(i,l1,k1)=1
                     done(i,k1,l1)=1
                  end if
               end if
   40       continue
               end do
               if (repeat2) then
                  repeat2 = .false.
                  start2 = kbx(l1) + 1
                  stop2 = nlight
                  goto 30
               end if

c            else
c             goto 20
c            end if

   20       continue
            print*, "In inner loop i=",i," l1=",l1," j=",j
            end do
            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nlight
               goto 10
            end if

         end do
      end do

c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nmol
         if (nmollst(i) .ge. maxelst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
      return
      end



c     ##  subroutine molfull3body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "molfull3body" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull3bodymod
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol
      integer start,stop
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,xcm1,ycm1,zcm1
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat
c
c
c     perform dynamic allocation of some local arrays
c
      print*, "Before allocation Molfull3body"

      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))

      print*, "After allocation Molfull3body"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nmol-1
         do j=i+1,nmol
         nmollst3(i,j) = 0
         end do
      end do 
c         ii = ipole(i)
      do i = 1, nmol
         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
c             j = pnum(l1)
            k=l1-1
            j=imol(1,i)+k                
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1
         
         xmolold(i) = xcm1
         ymolold(i) = ycm1
         zmolold(i) = zcm1
         xsort(i) = xcm1
         ysort(i) = ycm1
         zsort(i) = zcm1
      end do
c
c     use the method of lights to generate neighbors
      print*,"Before lights Call in Molfull3body"
c      molbuf2=4.25+3
      molbuf2=6.5+2
      off = sqrt(molbuf2)
      call lights (off,nmol,xsort,ysort,zsort)

      print*,"After lights Call in Molfull3body"
c
c     loop over all atoms computing the interactions
c
      do i = 1, nmol
         print*, "In loop i=",i
         do a=1,nmollst(i)
            xi = xsort(rgx(i))
            yi = ysort(rgy(i))
            zi = zsort(rgz(i))
      
            l1=mollst(a,i)
            xl1 = xsort(rgx(l1))
            yl1 = ysort(rgy(l1))
            zl1 = zsort(rgz(l1)) 
        
            if (kbx(l1) .le. kex(l1)) then
               repeat = .false.
               start = kbx(l1) + 1
c               start = kbx(i) + 2
               stop = kex(l1)
            else
               repeat = .true.
               start = 1
c                start = 2
              stop = kex(l1)
            end if

   10    continue

            do j = start, stop
            if((j .ne. i) .and. (j .ne. l1)) then
               k = locx(j)
               kgy = rgy(k)
            if (kby(l1) .le. key(l1)) then
               if (kgy.lt.kby(l1) .or. kgy.gt.key(l1))  goto 20
            else
               if (kgy.lt.kby(l1) .and. kgy.gt.key(l1))  goto 20
            end if
            kgz = rgz(k)
            if (kbz(l1) .le. kez(l1)) then
               if (kgz.lt.kbz(l1) .or. kgz.gt.kez(l1))  goto 20
            else
               if (kgz.lt.kbz(l1) .and. kgz.gt.kez(l1))  goto 20
            end if

               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
           
               call imagen (xr,yr,zr)

               xr1=xr
               yr1=yr
               zr1=zr
             
               xr = xl1 - xsort(j)
               yr = yl1 - ysort(kgy)
               zr = zl1 - zsort(kgz)
               
               call imagen (xr,yr,zr)
          
               xr2 = xr
               yr2 = yr
               zr2 = zr

               xr = xi - xl1
               yr = yi - yl1
               zr = zi - zl1
               
               call imagen(xr,yr,zr)
              
               xr3=xr
               yr3=yr
               zr3=zr

               r2 = (xr1*xr1 + yr1*yr1 + zr1*zr1 +
     &            xr2*xr2 + yr2*yr2 + zr2*zr2 +
     &            xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
               print*, "Before nb list update j=",j,"r2=",r2
               if (r2 .le. molbuf2*2) then
                  if (l1 .lt. k) then
                     nmollst3(i,l1) = nmollst3(i,l1) + 1
                     mollst3(nmollst3(i,l1),i,l1) = k
                  else
                     nmollst3(i,k) = nmollst3(i,k) + 1
                     mollst3(nmollst3(i,k),i,k) = l1
                  end if
               end if
               print*, "After nb list update j=",j,"r2=",r2
            end if
   20       continue
c            end if
            end do

            if (repeat) then
               repeat = .false.
               start = kbx(l1) + 1
c               start = kbx(i) +2
               stop = nlight
               goto 10
            end if

         end do
      end do

c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nmol
         if (nmollst(i) .ge. maxelst) then
            write (iout,30)
   30       format (/,' MFULL  --  Too many Neighbors;',
     &                 ' Increase MAXELST')
            call fatal
         end if
      end do
      return
      end

c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull3body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "molfull3body" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull3bodymod3
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol
      integer start,stop,start2,stop2
      integer k1,kgy1,kgz1,j1,done(nmol,nmol,nmol)
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,xcm1,ycm1,zcm1
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat,repeat2
c
c
c     perform dynamic allocation of some local arrays
c
      print*, "Before allocation Molfull3bodymod2"

      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))

      print*, "After allocation Molfull3bodymod2"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nmol
         do j=1,nmol
         nmollst3(i,j) = 0
         end do
      end do

      do i=1,nmol
         do j=1,nmol
            do k=1,nmol
               done(i,j,k)=0
            end do
         end do
      end do

      do i = 1, nmol
         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k=l1-1
            j=imol(1,i)+k
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1

         xmolold(i) = xcm1
         ymolold(i) = ycm1
         zmolold(i) = zcm1
         xsort(i) = xcm1
         ysort(i) = ycm1
         zsort(i) = zcm1
      end do
c
c     use the method of lights to generate neighbors
      print*,"Before lights Call in Molfull3bodymod3"
c      molbuf2=4.25+3
      molbuf2=5.5+2
      off = sqrt(molbuf2)
      call lights (off,nmol,xsort,ysort,zsort)

      print*,"After lights Call in Molfull3bodymod3"
c
c     loop over all atoms computing the interactions
c
      do i = 1, nmol
         print*, "In loop i=",i

            xi = xsort(rgx(i))
            yi = ysort(rgy(i))
            zi = zsort(rgz(i))

            if (kbx(i) .le. kex(i)) then
               repeat = .false.
               start = kbx(i) + 1
               stop = kex(i)
            else
               repeat = .true.
               start = 1
              stop = kex(i)
            end if

            print*, "In inner loop i=",i

   10    continue

            do j = start, stop
c            if(j .ne. i) then
               k = locx(j)
               kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
               kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if

               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)

               call imagen (xr,yr,zr)

               xr1=xr
               yr1=yr
               zr1=zr

               if (kbx(k) .le. kex(k)) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = kex(k)
               else
                  repeat2 = .true.
                  start2 = 1
                  stop2 = kex(k)
               end if

   30    continue


               do j1 = start2, stop2

                  k1=locx(j1)
                  kgy1 =rgy(k1)
                  if (kby(k) .le. key(k)) then
                   if(kgy1.lt.kby(k)
     &                .or. kgy1.gt.key(k)) goto 40
                  else
                   if (kgy1.lt.kby(k)
     &               .and. kgy1.gt.key(k)) goto 40
                  end if
                  kgz1=rgz(k1)
                  if (kbz(k) .le. kez(k)) then
                   if (kgz1.lt.kbz(k)
     &                 .or. kgz1.gt.kez(k))  goto 40
                  else
                   if (kgz1.lt.kbz(k)
     &                 .and. kgz1.gt.kez(k))  goto 40
                  end if


                  xr = xi - xsort(j1)
                  yr = yi - ysort(kgy1)
                  zr = zi - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr2=xr
                  yr2=yr
                  zr2=zr

                  xr = xsort(j) - xsort(j1)
                  yr = ysort(kgy) - ysort(kgy1)
                  zr = zsort(kgz) - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr3 = xr
                  yr3 = yr
                  zr3 = zr


                  r2 = (xr1*xr1 + yr1*yr1 + zr1*zr1 +
     &                 xr2*xr2 + yr2*yr2 + zr2*zr2 +
     &                xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
               if (r2 .le. molbuf2 .and. done(i,k,k1).eq.0
     &           .and. done(i,k1,k).eq.0 .and. done(k,i,k1).eq.0 
     &           .and. done(k,k1,i).eq.0 .and. done(k1,i,k).eq.0
     &           .and. done(k1,k,i).eq.0) then
                  if (i.lt.k1 .and. k1.lt.k) then
                     nmollst3(i,k1) = nmollst3(i,k1) + 1
                     mollst3(nmollst3(i,k1),i,k1) = k
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
                  if (i.lt.k .and. k.lt.k1) then
                     nmollst3(i,k) = nmollst3(i,k)+1
                     mollst3(nmollst3(i,k),i,k) =k1
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
                  if (k.lt.i .and. i.lt.k1) then
                     nmollst3(k,i) = nmollst3(k,i)+1
                     mollst3(nmollst3(k,i),k,i)=k1
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1 
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
                  if (k.lt.k1 .and. k1.lt.i) then
                     nmollst3(k,k1) = nmollst3(k,k1)+1
                     mollst3(nmollst3(k,k1),k,k1)=i
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1 
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if                     
                  if (k1.lt.i .and. i.lt.k) then
                     nmollst3(k1,i) = nmollst3(k1,i)+1
                     mollst3(nmollst3(k1,i),k1,i)=k
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
                  if (k1.lt.k .and. k.lt.i) then
                     nmollst3(k1,k) = nmollst3(k1,k)+1
                     mollst3(nmollst3(k1,k),k1,k)=i
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
               end if
   40       continue
               end do
               if (repeat2) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = nlight
                  goto 30
               end if


   20       continue
            print*, "In inner loop i=",i," j=",j
            end do
            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nlight
               goto 10
            end if

      end do

c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c
c     check to see if the neighbor lists are too long
c
      do i = 1, nmol
         if (nmollst(i) .ge. maxelst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')            
            call fatal
         end if
      end do
      return
      end



c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull3body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "molfull3body" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull3bodymod4
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol
      integer start,stop,start2,stop2
c      integer k1,kgy1,kgz1,j1,done(nmol,nmol,nmol)
      integer k1,kgy1,kgz1,j1
      integer, allocatable :: done(:,:,:)
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,xcm1,ycm1,zcm1
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat,repeat2
c
c
c     perform dynamic allocation of some local arrays
c
c      print*, "Before allocation Molfull3bodymod4"
      domollst3bod = .false.
      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c      print*, "After xyzsort(nlight) alloction"
c      allocate (done(nmol,nmol,nmol))
c      print*, "After done allocation Molfull3bodymod4"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nmol
       nmollst3mod(i)=0
      end do

c      do i=1,nmol
c         do j=1,nmol
c            do k=1,nmol
c               done(i,j,k)=0
c            end do
c         end do
c      end do

      do i = 1, nmol
         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k=l1-1
            j=imol(1,i)+k
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1

         xmolold(i) = xcm1
         ymolold(i) = ycm1
         zmolold(i) = zcm1
         xsort(i) = xcm1
         ysort(i) = ycm1
         zsort(i) = zcm1
      end do

c      allocate (done(nmol,nmol,nmol))
c
c     use the method of lights to generate neighbors
c      print*,"Before lights Call in Molfull3bodymod4"
c      molbuf2=4.25+3
c      molbuf2=8.5+2
      molbuf2=25.0d0
      off = sqrt(molbuf2)
      off =10.0d0
c      off=molbuf2
      call lights (off,nmol,xsort,ysort,zsort)

c      print*,"After lights Call in Molfull3bodymod4"
c
c     loop over all atoms computing the interactions
c
      do i = 1, nmol
c         print*, "In loop i=",i

            xi = xsort(rgx(i))
            yi = ysort(rgy(i))
            zi = zsort(rgz(i))

            if (kbx(i) .le. kex(i)) then
               repeat = .false.
               start = kbx(i) + 1
               stop = kex(i)
            else
               repeat = .true.
               start = 1
              stop = kex(i)
            end if

c            print*, "In inner loop i=",i

   10    continue

            do j = start, stop
c            if(j .ne. i) then
               k = locx(j)
               kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
               kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if

               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)

               call imagen (xr,yr,zr)

               xr1=xr
               yr1=yr
               zr1=zr

               if (kbx(k) .le. kex(k)) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = kex(k)
               else
                  repeat2 = .true.
                  start2 = 1
                  stop2 = kex(k)
               end if

   30    continue


               do j1 = start2, stop2

                  k1=locx(j1)
                  kgy1 =rgy(k1)
                  if (kby(k) .le. key(k)) then
                   if(kgy1.lt.kby(k)
     &                .or. kgy1.gt.key(k)) goto 40
                  else
                   if (kgy1.lt.kby(k)
     &               .and. kgy1.gt.key(k)) goto 40
                  end if
                  kgz1=rgz(k1)
                  if (kbz(k) .le. kez(k)) then
                   if (kgz1.lt.kbz(k)
     &                 .or. kgz1.gt.kez(k))  goto 40
                  else
                   if (kgz1.lt.kbz(k)
     &                 .and. kgz1.gt.kez(k))  goto 40
                  end if


                  xr = xi - xsort(j1)
                  yr = yi - ysort(kgy1)
                  zr = zi - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr2=xr
                  yr2=yr
                  zr2=zr

                  xr = xsort(j) - xsort(j1)
                  yr = ysort(kgy) - ysort(kgy1)
                  zr = zsort(kgz) - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr3 = xr
                  yr3 = yr
                  zr3 = zr


                  r2 = (xr1*xr1 + yr1*yr1 + zr1*zr1 +
     &                 xr2*xr2 + yr2*yr2 + zr2*zr2 +
     &                xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
c               if (r2 .le. molbuf2 .and. done(i,k,k1).eq.0
c     &           .and. done(i,k1,k).eq.0 .and. done(k,i,k1).eq.0 
c     &           .and. done(k,k1,i).eq.0 .and. done(k1,i,k).eq.0
c     &           .and. done(k1,k,i).eq.0) then
               if (r2 .le. molbuf2) then
                  if (i.lt.k1 .and. k1.lt.k) then
c                     nmollst3(i,k1) = nmollst3(i,k1) + 1
c                     mollst3(nmollst3(i,k1),i,k1) = k
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k1
c                     print*,"Nmollst3mod(i)=",nmollst3mod(i),"i=",i,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(i),i)
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k                     
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(i)=",nmollst3mod(i),"i=",i,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(i),i)
                  end if
                  if (i.lt.k .and. k.lt.k1) then
c                     nmollst3(i,k) = nmolst3(i,k)+1
c                     mollst3(nmollst3(i,k),i,k) =k1
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k
c                     print*,"Nmollst3mod(i)=",nmollst3mod(i),"i=",i,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(i),i)
                     nmollst3mod(i)=nmollst3mod(i)+1                    
                     mollst3mod(nmollst3mod(i),i)=k1
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(i)=",nmollst3mod(i),"i=",i,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(i),i)
                  end if
                  if (k.lt.i .and. i.lt.k1) then
c                     nmollst3(k,i) = nmolst(k,i)+1
c                     mollst3(nmollst3(k,i),k,i)=k1
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=i
c                     print*,"Nmollst3mod(k)=",nmollst3mod(k),"k=",k,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k),k)
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=k1
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1 
c                    done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(k)=",nmollst3mod(k),"k=",k,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k),k)
                  end if
                  if (k.lt.k1 .and. k1.lt.i) then
c                     nmollst3(k,k1) = nmollst3(k,k1)+1
c                     mollst3(nmollst3(k,k1),k,k1)=i
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=k1
c                     print*,"Nmollst3mod(k)=",nmollst3mod(k),"k=",k,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k),k)
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=i
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1 
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(k)=",nmollst3mod(k),"k=",k,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k),k)
                  end if                     
                  if (k1.lt.i .and. i.lt.k) then
c                     nmollst3(k1,i) = nmollst3(k1,i)+1
c                     mollst3(nmollst3(k1,i),k1,i)=k
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=i
c                     print*,"Nmollst3mod(k1)=",nmollst3mod(k1),"k1=",k1,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k1),k1)
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=k
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(k1)=",nmollst3mod(k1),"k1=",k1,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k1),k1)                  
                  end if
                  if (k1.lt.k .and. k.lt.i) then
c                     nmollst3(k1,k) = nmollst3(k1,k)+1
c                     mollst3(nmollst3(k1,k),k1,k)=i
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=k
c                     print*,"Nmollst3mod(k1)=",nmollst3mod(k1),"k1=",k1,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k1),k1)
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=i
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(k1)=",nmollst3mod(k1),"k1=",k1,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k1),k1)                  
                  end if
               end if
   40       continue
          
               end do
               if (repeat2) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = nlight
                  goto 30
               end if

c            else
c             goto 20
c            end if

   20       continue
c            print*, "In inner loop i=",i," j=",j
            end do
            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nlight
               goto 10
            end if

      end do

c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c      deallocate (done)
c
c     check to see if the neighbor lists are too long
c
c      do i = 1, nmol
c         if (nmollst3mod(i) .ge. maxelst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')            
c            call fatal
c            print*, "Too many in nmollst3mod!!"
c         end if
c      end do
      return
      end



c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull3body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "molfull3body" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull3bodymod4nominimage
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol
      integer start,stop,start2,stop2
c      integer k1,kgy1,kgz1,j1,done(nmol,nmol,nmol)
      integer k1,kgy1,kgz1,j1
      integer, allocatable :: done(:,:,:)
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,xcm1,ycm1,zcm1
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical repeat,repeat2
c
c
c     perform dynamic allocation of some local arrays
c
      print*, "Before allocation Molfull3bodymod4nominimage"

      nlight = (ncell+1) * nmol
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
      allocate (done(nmol,nmol,nmol))

      print*, "After allocation Molfull3bodymod44nominimage"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nmol
       nmollst3mod(i)=0
      end do

      do i=1,nmol
         do j=1,nmol
            do k=1,nmol
               done(i,j,k)=0
            end do
         end do
      end do

      do i = 1, nmol
         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k=l1-1
            j=imol(1,i)+k
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1

         xmolold(i) = xcm1
         ymolold(i) = ycm1
         zmolold(i) = zcm1
         xsort(i) = xcm1
         ysort(i) = ycm1
         zsort(i) = zcm1
      end do
c
c     use the method of lights to generate neighbors
      print*,"Before lights Call in Molfull3bodymod4nominimage"
c      molbuf2=4.25+3
      molbuf2=10
c      off = sqrt(molbuf2)
      call lights (off,nmol,xsort,ysort,zsort)

      molbuf2=15
      print*,"After lights Call in Molfull3bodymod4nominimage"
c
c     loop over all atoms computing the interactions
c
      do i = 1, nmol
c         print*, "In loop i=",i

            xi = xsort(rgx(i))
            yi = ysort(rgy(i))
            zi = zsort(rgz(i))

            if (kbx(i) .le. kex(i)) then
               repeat = .false.
               start = kbx(i) + 1
               stop = kex(i)
            else
               repeat = .true.
               start = 1
              stop = kex(i)
            end if

c            print*, "In inner loop i=",i

   10    continue

            do j = start, stop
c            if(j .ne. i) then
               k = locx(j)
               kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
               kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if

               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)

c               call imagen (xr,yr,zr)

               xr1=xr
               yr1=yr
               zr1=zr

               if (kbx(k) .le. kex(k)) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = kex(k)
               else
                  repeat2 = .true.
                  start2 = 1
                  stop2 = kex(k)
               end if

   30    continue


               do j1 = start2, stop2

                  k1=locx(j1)
                  kgy1 =rgy(k1)
                  if (kby(k) .le. key(k)) then
                   if(kgy1.lt.kby(k)
     &                .or. kgy1.gt.key(k)) goto 40
                  else
                   if (kgy1.lt.kby(k)
     &               .and. kgy1.gt.key(k)) goto 40
                  end if
                  kgz1=rgz(k1)
                  if (kbz(k) .le. kez(k)) then
                   if (kgz1.lt.kbz(k)
     &                 .or. kgz1.gt.kez(k))  goto 40
                  else
                   if (kgz1.lt.kbz(k)
     &                 .and. kgz1.gt.kez(k))  goto 40
                  end if


                  xr = xi - xsort(j1)
                  yr = yi - ysort(kgy1)
                  zr = zi - zsort(kgz1)

c                  call imagen (xr,yr,zr)

                  xr2=xr
                  yr2=yr
                  zr2=zr

                  xr = xsort(j) - xsort(j1)
                  yr = ysort(kgy) - ysort(kgy1)
                  zr = zsort(kgz) - zsort(kgz1)

c                  call imagen (xr,yr,zr)

                  xr3 = xr
                  yr3 = yr
                  zr3 = zr


                  r2 = (xr1*xr1 + yr1*yr1 + zr1*zr1 +
     &                 xr2*xr2 + yr2*yr2 + zr2*zr2 +
     &                xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
               if (r2 .le. molbuf2 .and. done(i,k,k1).eq.0
     &           .and. done(i,k1,k).eq.0 .and. done(k,i,k1).eq.0 
     &           .and. done(k,k1,i).eq.0 .and. done(k1,i,k).eq.0
     &           .and. done(k1,k,i).eq.0) then
                  if (i.lt.k1 .and. k1.lt.k) then
c                     nmollst3(i,k1) = nmollst3(i,k1) + 1
c                     mollst3(nmollst3(i,k1),i,k1) = k
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k1
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k                     
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
                  if (i.lt.k .and. k.lt.k1) then
c                     nmollst3(i,k) = nmolst3(i,k)+1
c                     mollst3(nmollst3(i,k),i,k) =k1
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k
                     nmollst3mod(i)=nmollst3mod(i)+1                    
                     mollst3mod(nmollst3mod(i),i)=k1
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
                  if (k.lt.i .and. i.lt.k1) then
c                     nmollst3(k,i) = nmolst(k,i)+1
c                     mollst3(nmollst3(k,i),k,i)=k1
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=i
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=k1
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1 
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
                  if (k.lt.k1 .and. k1.lt.i) then
c                     nmollst3(k,k1) = nmollst3(k,k1)+1
c                     mollst3(nmollst3(k,k1),k,k1)=i
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=k1
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=i
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1 
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if                     
                  if (k1.lt.i .and. i.lt.k) then
c                     nmollst3(k1,i) = nmollst3(k1,i)+1
c                     mollst3(nmollst3(k1,i),k1,i)=k
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=i
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=k
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
                  if (k1.lt.k .and. k.lt.i) then
c                     nmollst3(k1,k) = nmollst3(k1,k)+1
c                     mollst3(nmollst3(k1,k),k1,k)=i
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=k
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=i
                     done(i,k1,k)=1
                     done(i,k,k1)=1
                     done(k,i,k1)=1
                     done(k,k1,i)=1
                     done(k1,i,k)=1
                     done(k1,k,i)=1
                  end if
               end if
   40       continue
               end do
               if (repeat2) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = nlight
                  goto 30
               end if

c            else
c             goto 20
c            end if

   20       continue
c            print*, "In inner loop i=",i," j=",j
            end do
            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nlight
               goto 10
            end if

      end do

c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      deallocate (done)
c
c     check to see if the neighbor lists are too long
c
c      do i = 1, nmol
c         if (nmollst(i) .ge. maxelst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')            
c            call fatal
c         end if
c      end do
      return
      end



c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine molfull3body--make 2body pair list for all mol ##
c     ##                                                            ##
c     ################################################################
c
c
c     "molfull3body" performs a complete rebuild of the atomic multipole
c     pair neighbor list for all sites using the method of lights
c
c
      subroutine molfull3bodycobar
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cell.i'
      include 'iounit.i'
      include 'light.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,a
      integer kgy,kgz,l1,lenmol
      integer start,stop,start2,stop2
c      integer k1,kgy1,kgz1,j1,done(nmol,nmol,nmol)
      integer k1,kgy1,kgz1,j1
      integer, allocatable :: done(:,:,:)
      real*8 xi,yi,zi,M1,xl1,yl1,zl1
      real*8 xr,yr,zr,xr1,yr1,zr1
      real*8 xr2,yr2,zr2,xr3,yr3,zr3
      real*8 r2,off,x1,y1,z1,r1,r3
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      real*8  molbuf_cobar,shellsum
      integer shell1,shell2,shell3,shellsum_mod
      logical repeat,repeat2
c
c
c     perform dynamic allocation of some local arrays
c
c      print*, "Before allocation Molfull3bodymod4"
      domollst3bod = .false.
      nlight = (ncell+1) * nmol
      print*,"nlight in molfull3bodycobar",nlight 
      allocate (xsort(nlight))
      allocate (ysort(nlight))
      allocate (zsort(nlight))
c      print*, "After xyzsort(nlight) alloction"
c      allocate (done(nmol,nmol,nmol))
c      print*, "After done allocation Molfull3bodymod4"
c
c     transfer interaction site coordinates to sorting arrays
c
      do i = 1, nmol
       nmollst3mod(i)=0
      end do

c      do i=1,nmol
c         do j=1,nmol
c            do k=1,nmol
c               done(i,j,k)=0
c            end do
c         end do
c      end do

      do i = 1, nmol
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
            k=l1-1
            j=imol(1,i)+k
           if(name(j).eq.'O') then
            x1 = x(j)
            y1 = y(j)
            z1 = z(j)
           end if
         end do


         xmolold3(i) = x1
         ymolold3(i) = y1
         zmolold3(i) = z1
         xsort(i) = x1
         ysort(i) = y1
         zsort(i) = z1
      end do

c      allocate (done(nmol,nmol,nmol))
c
c     use the method of lights to generate neighbors
c      print*,"Before lights Call in Molfull3bodymod4"
c      molbuf2=4.25+3
c      molbuf2=8.5+2
c      molbuf2=4
c      off = sqrt(molbuf2)
      molbuf_cobar=8.0d0
      off =10.0d0
c      off=molbuf2
      call lights (off,nmol,xsort,ysort,zsort)

c      print*,"After lights Call in Molfull3bodymod4"
c
c     loop over all atoms computing the interactions
c
      do i = 1, nmol
c         print*, "In loop i=",i

            xi = xsort(rgx(i))
            yi = ysort(rgy(i))
            zi = zsort(rgz(i))

            if (kbx(i) .le. kex(i)) then
               repeat = .false.
               start = kbx(i) + 1
               stop = kex(i)
            else
               repeat = .true.
               start = 1
              stop = kex(i)
            end if

c            print*, "In inner loop i=",i

   10    continue

            do j = start, stop
c            if(j .ne. i) then
               k = locx(j)
               kgy = rgy(k)
            if (kby(i) .le. key(i)) then
               if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 20
            else
               if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 20
            end if
               kgz = rgz(k)
            if (kbz(i) .le. kez(i)) then
               if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 20
            else
               if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 20
            end if

               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)

               call imagen (xr,yr,zr)

               xr1=xr
               yr1=yr
               zr1=zr

               if (kbx(k) .le. kex(k)) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = kex(k)
               else
                  repeat2 = .true.
                  start2 = 1
                  stop2 = kex(k)
               end if

   30    continue


               do j1 = start2, stop2

                  k1=locx(j1)
                  kgy1 =rgy(k1)
                  if (kby(k) .le. key(k)) then
                   if(kgy1.lt.kby(k)
     &                .or. kgy1.gt.key(k)) goto 40
                  else
                   if (kgy1.lt.kby(k)
     &               .and. kgy1.gt.key(k)) goto 40
                  end if
                  kgz1=rgz(k1)
                  if (kbz(k) .le. kez(k)) then
                   if (kgz1.lt.kbz(k)
     &                 .or. kgz1.gt.kez(k))  goto 40
                  else
                   if (kgz1.lt.kbz(k)
     &                 .and. kgz1.gt.kez(k))  goto 40
                  end if


                  xr = xi - xsort(j1)
                  yr = yi - ysort(kgy1)
                  zr = zi - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr2=xr
                  yr2=yr
                  zr2=zr

                  xr = xsort(j) - xsort(j1)
                  yr = ysort(kgy) - ysort(kgy1)
                  zr = zsort(kgz) - zsort(kgz1)

                  call imagen (xr,yr,zr)

                  xr3 = xr
                  yr3 = yr
                  zr3 = zr
                  
            r1=sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)
            r2=sqrt(xr2*xr2 + yr2*yr2 + zr2*zr2)
            r3=sqrt(xr3*xr3 + yr3*yr3 + zr3*zr3)

            if( (r1.lt.r3).and.(r2.lt.r3) ) then
c              shell1=int(r1/3.1d0)+1
c              shell2=int(r2/3.1d0)+1
c              shell3=int(r3/3.1d0)+1
c              shellsum=shell1+shell2
              shellsum=r1+r2
c              shellsum_mod=shell1+shell2+shell3
            else if ( (r1.lt.r2).and.(r3.lt.r2)) then
c              shell1=int(r1/3.1d0)+1
c              shell2=int(r3/3.1d0)+1
c              shell3=int(r2/3.1d0)+1
c              shellsum=shell1+shell2
              shellsum=r1+r3
c              shellsum_mod=shell1+shell2+shell3
            else if ( (r2.lt.r1).and.(r3.lt.r1)) then
c              shell1=int(r2/3.1d0)+1
c              shell2=int(r3/3.1d0)+1
c              shell3=int(r1/3.1d0)+1
c              shellsum=shell1+shell2
              shellsum=r2+r3
c              shellsum_mod=shell1+shell2+shell3
            end if


c               if (r2 .le. molbuf2 .and. done(i,k,k1).eq.0
c     &           .and. done(i,k1,k).eq.0 .and. done(k,i,k1).eq.0 
c     &           .and. done(k,k1,i).eq.0 .and. done(k1,i,k).eq.0
c     &           .and. done(k1,k,i).eq.0) then
               if (shellsum .le. molbuf_cobar) then
                  if (i.lt.k1 .and. k1.lt.k) then
c                     nmollst3(i,k1) = nmollst3(i,k1) + 1
c                     mollst3(nmollst3(i,k1),i,k1) = k
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k1
c                     print*,"Nmollst3mod(i)=",nmollst3mod(i),"i=",i,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(i),i)
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k                     
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(i)=",nmollst3mod(i),"i=",i,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(i),i)
                  end if
                  if (i.lt.k .and. k.lt.k1) then
c                     nmollst3(i,k) = nmolst3(i,k)+1
c                     mollst3(nmollst3(i,k),i,k) =k1
                     nmollst3mod(i)=nmollst3mod(i)+1
                     mollst3mod(nmollst3mod(i),i)=k
c                     print*,"Nmollst3mod(i)=",nmollst3mod(i),"i=",i,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(i),i)
                     nmollst3mod(i)=nmollst3mod(i)+1                    
                     mollst3mod(nmollst3mod(i),i)=k1
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(i)=",nmollst3mod(i),"i=",i,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(i),i)
                  end if
                  if (k.lt.i .and. i.lt.k1) then
c                     nmollst3(k,i) = nmolst(k,i)+1
c                     mollst3(nmollst3(k,i),k,i)=k1
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=i
c                     print*,"Nmollst3mod(k)=",nmollst3mod(k),"k=",k,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k),k)
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=k1
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1 
c                    done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(k)=",nmollst3mod(k),"k=",k,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k),k)
                  end if
                  if (k.lt.k1 .and. k1.lt.i) then
c                     nmollst3(k,k1) = nmollst3(k,k1)+1
c                     mollst3(nmollst3(k,k1),k,k1)=i
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=k1
c                     print*,"Nmollst3mod(k)=",nmollst3mod(k),"k=",k,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k),k)
                     nmollst3mod(k)=nmollst3mod(k)+1
                     mollst3mod(nmollst3mod(k),k)=i
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1 
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(k)=",nmollst3mod(k),"k=",k,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k),k)
                  end if                     
                  if (k1.lt.i .and. i.lt.k) then
c                     nmollst3(k1,i) = nmollst3(k1,i)+1
c                     mollst3(nmollst3(k1,i),k1,i)=k
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=i
c                     print*,"Nmollst3mod(k1)=",nmollst3mod(k1),"k1=",k1,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k1),k1)
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=k
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(k1)=",nmollst3mod(k1),"k1=",k1,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k1),k1)                  
                  end if
                  if (k1.lt.k .and. k.lt.i) then
c                     nmollst3(k1,k) = nmollst3(k1,k)+1
c                     mollst3(nmollst3(k1,k),k1,k)=i
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=k
c                     print*,"Nmollst3mod(k1)=",nmollst3mod(k1),"k1=",k1,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k1),k1)
                     nmollst3mod(k1)=nmollst3mod(k1)+1
                     mollst3mod(nmollst3mod(k1),k1)=i
c                     done(i,k1,k)=1
c                     done(i,k,k1)=1
c                     done(k,i,k1)=1
c                     done(k,k1,i)=1
c                     done(k1,i,k)=1
c                     done(k1,k,i)=1
c                     print*,"Nmollst3mod(k1)=",nmollst3mod(k1),"k1=",k1,
c     &               "mollst3mod=",mollst3mod(nmollst3mod(k1),k1)                  
                  end if
               end if
   40       continue
          
               end do
               if (repeat2) then
                  repeat2 = .false.
                  start2 = kbx(k) + 1
                  stop2 = nlight
                  goto 30
               end if

c            else
c             goto 20
c            end if

   20       continue
c            print*, "In inner loop i=",i," j=",j
            end do
            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nlight
               goto 10
            end if

      end do

c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
c      deallocate (done)
c
c     check to see if the neighbor lists are too long
c
c      do i = 1, nmol
c         if (nmollst3mod(i) .ge. maxelst) then
c            write (iout,30)
c   30       format (/,' MFULL  --  Too many Neighbors;',
c     &                 ' Increase MAXELST')            
c            call fatal
c            print*, "Too many in nmollst3mod!!"
c         end if
c      end do
      return
      end



c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mollist -- build atom multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mollist" performs an update or a complete rebuild of the
c     body/molecule neighbor lists for 3-body approx of polariz. energy
c
c
      subroutine mollist2body
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,l1,lenmol
      real*8 xi,yi,zi,xcm1,ycm1,zcm1
      real*8 xr,yr,zr,M1
      real*8 radius,r2,mol_lbuffer,mol_lbuf2
c
c
c     neighbor list cannot be used with the replicates method
c
c      radius = sqrt(mbuf2)
c      call replica (radius)
c      if (use_replica) then
c         write (iout,10)
c   10    format (/,' MLIST  --  Pairwise Neighbor List cannot',
c     &              ' be used with Replicas')
c         call fatal
c      end if
c
c     perform a complete list build instead of an update
c
c      print*, "Hello from Mollist2body"
c         mbuf2 = (mpolecut+lbuffer)**2
c         mbufx = (mpolecut+2.0d0*lbuffer)**2

      molbuf2=25.0d0
      mol_lbuffer=4.0d0
      mol_lbuf2=(0.5d0*mol_lbuffer)**2
      molbufx=81.0d0

      if (domollst2bod) then
         domollst2bod = .false.
c         if (octahedron) then
c            do i = 1, npole
c               call mbuild (i)
c            end do
c         else
c          print*, "Within Mollist2body if domollst2bod"
            call molfull2body
c         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, nmol
c         ii = ipole(i)

         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
c             j = pnum(l1)
            k=l1-1
            j=imol(1,i)+k
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1

         xr = xcm1 - xmolold(i)
         yr = ycm1 - ymolold(i)
         zr = zcm1 - zmolold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .ge. mol_lbuf2) then
c
c     store coordinates to reflect update of this site
c
            xmolold(i) = xcm1
            ymolold(i) = ycm1
            zmolold(i) = zcm1
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, nmol
               xr = xcm1 - xmolold(k)
               yr = ycm1 - ymolold(k)
               zr = zcm1 - zmolold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. molbuf2) then
                  j = j + 1
                  mollst(j,i) = k
               end if
            end do
            nmollst(i) = j
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xcm1 - xmolold(k)
               yr = ycm1 - ymolold(k)
               zr = zcm1 - zmolold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. molbuf2) then
                  do j = 1, nmollst(k)
                     if (mollst(j,k) .eq. i)  goto 20
                  end do
                  nmollst(k) = nmollst(k) + 1
                  mollst(nmollst(k),k) = i
   20             continue
               else if (r2 .le. molbufx) then
                  do j = 1, nmollst(k)
                     if (mollst(j,k) .eq. i) then
                        mollst(j,k) = mollst(nmollst(k),k)
                        nmollst(k) = nmollst(k) - 1
                        goto 30
                     end if
                  end do
   30             continue
               end if
            end do
         end if
      end do
      return
      end

c     #################################################################
c     ##                                                             ##
c     ##  subroutine mollist -- build atom multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mollist" performs an update or a complete rebuild of the
c     body/molecule neighbor lists for 3-body approx of polariz. energy
c
c
      subroutine mollist3body
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
c      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,k1,ii,l1,lenmol,listtemp(800),tempcount
      real*8 xi,yi,zi,x1,y1,z1
      real*8 xr,yr,zr,xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3,M1
      real*8 radius,r1,r3,r2,r2outer,molbuf_cobar
      real*8 mol_lbuffer,mol_lbuf2,shellsum
c
c
c     neighbor list cannot be used with the replicates method
c
c      radius = sqrt(mbuf2)
c      call replica (radius)
c      if (use_replica) then
c         write (iout,10)
c   10    format (/,' MLIST  --  Pairwise Neighbor List cannot',
c     &              ' be used with Replicas')
c         call fatal
c      end if
c
c     perform a complete list build instead of an update
c
c      print*, "Hello from Mollist3body"

c         mbuf2 = (mpolecut+lbuffer)**2
c         mbufx = (mpolecut+2.0d0*lbuffer)**2


      molbuf_cobar=8.0d0
c      molbuf_cobar2=16.0d0
      mol_lbuffer=4.0d0
      mol_lbuf2=(0.5d0*mol_lbuffer)**2
c      molbuf_cobarx=256.0d0

      if (domollst3bod) then
         domollst3bod = .false.
c         if (octahedron) then
c            do i = 1, npole
c               call mbuild (i)
c            end do
c         else
c          print*, "Within Mollist3body if domollst3bod"
c            call molfull3body
c            call molfull3bodymod2
c            call molfull3bodymod3
            call molfull3bodycobar
c            call molfull3bodymod4nominimage
c         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, nmol
c         ii = ipole(i)
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
c             j = pnum(l1)
            k=l1-1
            j=imol(1,i)+k
           if(name(j).eq.'O') then
            x1 = x(j)
            y1 = y(j)
            z1 = z(j)
           end if
         end do

         xr = x1 - xmolold3(i)
         yr = y1 - ymolold3(i)
         zr = z1 - zmolold3(i)
         call imagen (xr,yr,zr)
         r2outer = xr*xr + yr*yr + zr*zr
         if (r2outer .ge. mol_lbuf2) then
c
c     store coordinates to reflect update of this site
c
            xmolold3(i) = x1
            ymolold3(i) = y1
            zmolold3(i) = z1
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, nmol-1
               xr = x1 - xmolold3(k)
               yr = y1 - ymolold3(k)
               zr = z1 - zmolold3(k)
               call imagen (xr,yr,zr)
               xr1=xr
               yr1=yr
               zr1=zr
               do k1=k+1, nmol
                   xr = x1 - xmolold3(k1)
                   yr = y1 - ymolold3(k1)
                   zr = z1 - zmolold3(k1)
                   call imagen (xr,yr,zr)
                   xr2=xr
                   yr2=yr
                   zr2=zr
                   xr = xmolold3(k) - xmolold3(k1)
                   yr = ymolold3(k) - ymolold3(k1)
                   zr = zmolold3(k) - zmolold3(k1)
                   call imagen (xr,yr,zr)
                   xr3=xr
                   yr3=yr
                   zr3=zr

                  r1=sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)
                  r2=sqrt(xr2*xr2 + yr2*yr2 + zr2*zr2)
                  r3=sqrt(xr3*xr3 + yr3*yr3 + zr3*zr3)

                  if( (r1.lt.r3).and.(r2.lt.r3) ) then
                     shellsum=r1+r2
                  else if ( (r1.lt.r2).and.(r3.lt.r2)) then
                     shellsum=r1+r3
                  else if ( (r2.lt.r1).and.(r3.lt.r1)) then
                     shellsum=r2+r3
                  end if

                  if (shellsum .le. molbuf_cobar) then
                       j =j+1
                       mollst3mod(j,i)=k
                       j =j+1
                       mollst3mod(j,i)=k1
                  end if
               end do  
            end do
            nmollst3mod(i) = j
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               do k1=k+1,nmol
                 if (k1.ne.i) then
                   xr = x1 - xmolold3(k)
                   yr = y1 - ymolold3(k)
                   zr = z1 - zmolold3(k)
                   call imagen (xr,yr,zr)
                   xr1=xr
                   yr1=yr
                   zr1=zr
c                   r2 = xr*xr + yr*yr + zr*zr
                   xr = x1 - xmolold3(k1)
                   yr = y1 - ymolold3(k1)
                   zr = z1 - zmolold3(k1)
                   call imagen (xr,yr,zr)
                   xr2=xr
                   yr2=yr
                   zr2=zr

                   xr=xmolold3(k)-xmolold3(k1)
                   yr=ymolold3(k)-ymolold3(k1)
                   zr=zmolold3(k)-zmolold3(k1)
                   call imagen (xr,yr,zr)
                   xr3=xr
                   yr3=yr
                   zr3=zr
                
                  r1=sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)
                  r2=sqrt(xr2*xr2 + yr2*yr2 + zr2*zr2)
                  r3=sqrt(xr3*xr3 + yr3*yr3 + zr3*zr3)

                  if( (r1.lt.r3).and.(r2.lt.r3) ) then
                     shellsum=r1+r2
                  else if ( (r1.lt.r2).and.(r3.lt.r2)) then
                     shellsum=r1+r3
                  else if ( (r2.lt.r1).and.(r3.lt.r1)) then
                     shellsum=r2+r3
                  end if


c      LEFT OFF HERE
c                   if (r2 .le. molbuf2) then
                   if(shellsum .le. molbuf_cobar) then
                     do j = 1, nmollst3mod(k),2
                       if(((mollst3mod(j-1,k).eq.k1).and.
     &(mollst3mod(j,k).eq.i)).or.((mollst3mod(j-1,k).eq.i).and.
     &                  (mollst3mod(j,k).eq.k1))) goto 20
                     end do
                     nmollst3mod(k) = nmollst3mod(k) + 2
                     if(k1.lt.i) then
                       mollst3mod(nmollst3mod(k)-1,k)=k1
                       mollst3mod(nmollst3mod(k),k)=i
                     else
                       mollst3mod(nmollst3mod(k)-1,k)=i
                       mollst3mod(nmollst3mod(k),k)=k1
                     end if
   20                continue
c                   else if (r2 .le. molbufx) then
                   else
                     do j = 1, nmollst3mod(k),2
                       if(((mollst3mod(j-1,k).eq.k1).and.
     &(mollst3mod(j,k).eq.i)).or.((mollst3mod(j-1,k).eq.i).and.
     &(mollst3mod(j,k).eq.k1))) then
                          mollst3mod(j-1,k)=0
                          mollst3mod(j,k)=0
c                         mollst(j,k) = mollst(nmollst(k),k)
c                         nmollst(k) = nmollst(k) - 1
                         goto 30
                       end if
                     end do
   30                continue
                   end if
                 end if
               end do

               do j=1,800
                 listtemp(j)=0
               end do 
               tempcount=0  
               do j=1,nmollst3mod(k)
                  if(mollst3mod(j,k).ne.0) then
                     tempcount=tempcount+1
                    listtemp(tempcount)=mollst3mod(j,k)
                  end if
               end do  
               nmollst3mod(k)=0
               do j=1,tempcount
                 nmollst3mod(k)=nmollst3mod(k)+1
                 mollst3mod(nmollst3mod(k),k)=listtemp(j)
               end do 
            end do
         end if
      end do
      return
      end

