c Returns the coordinates of grid points that lie within a TOLERANCE,
c on electron density isosurface of value ISOVALUE
c
c ISOVALUE == 0.002 ~ Vdw radius
c
c Mark J. Williamson Summer 2015
c
      subroutine get_points_from_iso(rtdb, basis, geom, isovalue,
     > tolerance, sub_nprp, i_sub_prp_c)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
      integer rtdb                ! [Input] rtdb        
      integer basis               ! [Input] Basis set
      integer geom                ! [Input] Geometry
      double precision isovalue   ! [Input]
      double precision tolerance  ! [Input]
      integer sub_nprp            ! [Output]
      integer i_sub_prp_c         ! [Output]

      logical util_get_io_unit
      integer i,j,n,nat,iat
      integer nprp
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      character*16 element, at_tag
      integer i_prp_c
      integer h_sub_prp_c
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt
      integer l_eden, k_eden
      double precision xp,yp,zp

      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elpmap: geom_ncent',911,GEOM_ERR)

      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_elpmap: ma failed',911,MA_ERR)

      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_elpmap: ma failed',911,MA_ERR)

c Copy atoms over to k_xyzpt array
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_elpmap: geom_cent_get',911,GEOM_ERR)
   30 continue


c Set the grid up for the initial electron density calculation
      call prop_grid_initialize(rtdb,nat,dbl_mb(k_xyzpt))
c Get the number of points and their coordinate point index, i_prp_c
      call prop_grid_get_r_ptr(nprp,i_prp_c)

c allocate memory and an index, k_eden, to the density values
      if (.not. ma_push_get(mt_dbl,nprp,'eden',l_eden,k_eden))
     &    call errquit('hnd_elpmap: eden ma failed',911,MA_ERR)

c Calculates density matrix
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)

c Calculates the electronic contribution of the
c electronic integral defined by nder for a given density at 
c the grid points defined in points. 
      call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(i_prp_c),nprp, 
     & dbl_mb(k_eden),-1)

c Scan ahead and work out how many points (sub_nprp) lie on the isosurface
      sub_nprp = 0
      do i=1, nprp
        if ( abs(dbl_mb(k_eden + i) - isovalue) < tolerance ) then
          sub_nprp = sub_nprp + 1
        end if
      end do

c Allocate array for storting the location of these points
      if(.not.ma_alloc_get(MT_DBL, 3*sub_nprp, 'grid array',
     &      h_sub_prp_c, i_sub_prp_c) ) call errquit(
     &      'grid_initialize: unable to allocate heap space',
     &      0, MA_ERR)

      sub_nprp = 0
      do i=1, nprp
        if ( abs(dbl_mb(k_eden + i) - isovalue) < tolerance ) then
c Coordinates of grid point
           xp = dbl_mb(i_prp_c  +3*(i-1))
           yp = dbl_mb(i_prp_c+1+3*(i-1))
           zp = dbl_mb(i_prp_c+2+3*(i-1))

c Copy matching coordinates to new array
           dbl_mb(i_sub_prp_c  +3*(sub_nprp)) = xp
           dbl_mb(i_sub_prp_c+1+3*(sub_nprp)) = yp
           dbl_mb(i_sub_prp_c+2+3*(sub_nprp)) = zp

           sub_nprp = sub_nprp + 1
        end if
      end do
 

c Clean up phase
      call prop_grid_destroy()

      if (.not.ma_pop_stack(l_eden)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_eden failed',911,MA_ERR)

      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)

      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
    
1000  format(4F13.8)
      end subroutine








c Utility method for writing ESP points out to a file
      subroutine write_ESP_to_grid_file(sub_nprp, i_sub_prp_c,
     & k_tepot, esp_name)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
      integer un,iat 
      integer sub_nprp, i_sub_prp_c
      integer k_tepot
      double precision xp, yp, zp
      character*255 esp_name

c output file
      if(.not.util_get_io_unit(un))
     >    call errquit("no free file unit",0,0)

      open(unit=un,file=esp_name,
     >       form='formatted',status='unknown')

      write(un, '(I6X80A)') sub_nprp,"ESP value dump"
    
      do iat=1,sub_nprp
         xp = dbl_mb(i_sub_prp_c  +3*(iat-1))
         yp = dbl_mb(i_sub_prp_c+1+3*(iat-1))
         zp = dbl_mb(i_sub_prp_c+2+3*(iat-1))

         write(un,'(3F10.5,F15.6)') xp,yp,zp,
     &          dbl_mb(k_tepot+iat-1)

      end do

      close(un)

      end subroutine
