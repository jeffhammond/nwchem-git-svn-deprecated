      subroutine hnd_elpiso(rtdb,basis,geom)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "stdio.fh"
#include "bas.fh"
C
C     This routine calculates the electrostatic potential
C     on a given electron density isosurface.
C
      integer rtdb           ! [Input] rtdb
      integer basis          ! [Input] Basis set
      integer geom           ! [Input] Geometry


      double precision iso   ! Electron density iso surface to work on
      double precision tol   ! tolerance for above

      integer nat            ! Number of atoms in system
      character*16 at_tag
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      character*30 theory
      integer nbf

      integer l_xyzpt     ! Memhandle to atoms coordinate array
      integer k_xyzpt     ! Index to atoms coordinate array
      integer l_zanpt
      integer k_zanpt

      integer nprp           ! Number of points on cube
      integer h_prp_c        ! Memhandle to points on cube
      integer i_prp_c        ! Index to points on cube

      integer sub_nprp           ! Number of points on isosurface
      integer sub_h_prp_c        ! Memhandle to points on isosuface
      integer sub_i_prp_c        ! Index to points on isosuface

      integer l_epot         ! Memhandle to TODO
      integer k_epot         ! Index to TODO

      integer l_tepot        ! Memhandle to total electric field array
      integer k_tepot        ! Index to total electric field array

      integer iat, i
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision elpotn
      double precision rr

      character*255 esp_name ! File name for output cube

C set inputs
      if(.not.rtdb_get(rtdb,'prop:espiso:iso',MT_DBL,1,iso)) then
         iso = 0.002D0
      end if
      if(.not.rtdb_get(rtdb,'prop:espiso:tol',MT_DBL,1,tol)) then
         tol = 0.00003D0
      end if

c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Get density matrix
c
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     + call errquit('task: no task input for theory?',0, RTDB_ERR)

      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)

C     Determine how many atoms we have
      if (.not.geom_ncent(geom,nat))
     &    call errquit('hnd_elpmap: geom_ncent',911,GEOM_ERR)

      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_elpmap: ma failed',911,MA_ERR)

      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_elpmap: ma failed',911,MA_ERR)

c    Place atoms
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_elpmap: geom_cent_get',911,GEOM_ERR)
   30 continue

C     Calculate the electron density and extract a subset of points
      call get_points_from_iso(rtdb,basis,geom,iso,tol,
     & sub_nprp,sub_i_prp_c)

C     Allocate total electric field array
      if (.not. ma_push_get(mt_dbl,sub_nprp,'tot epot',l_tepot,k_tepot))
     &    call errquit('hnd_elpiso: ma failed',911,MA_ERR)

C     TODO
      if (.not. ma_push_get(mt_dbl,sub_nprp,'epot pnt',l_epot,k_epot))
     &    call errquit('hnd_elpiso: ma failed',911,MA_ERR)


c     This routine calculates the electronic contribution of the
c     electronic integral defined by nder for a given density at
c     the grid points defined in points.
      call hnd_elfcon(basis,geom,g_dens(ndens),
     &       dbl_mb(sub_i_prp_c),sub_nprp,dbl_mb(k_epot),0)


c Main loop
      if (ga_nodeid().gt.0) goto 300
c
c     ----- calculate electrostatic potential -----
c
      do 230  iat=1,sub_nprp
         xp = dbl_mb(sub_i_prp_c  +3*(iat-1))
         yp = dbl_mb(sub_i_prp_c+1+3*(iat-1))
         zp = dbl_mb(sub_i_prp_c+2+3*(iat-1))

         elpotn = -dbl_mb(k_epot+iat-1)
c
c     ----- add nuclear contribution -----
c
         do 210 i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)

            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if (rr.lt.1.0d-3) go to 210
            elpotn = elpotn + zan/rr
  210    continue

         dbl_mb(k_tepot+iat-1) = elpotn

  230    continue ! Assembling and printing next atom



  300 call ga_sync()


C     Write to a cube file
      if(ga_nodeid().eq.0) then

C  This is commented out since it causes any multiproc (np > 1)
C  to hang
C
        if (.not. rtdb_cget(rtdb, "prop:grid:output",1,esp_name))
     >     call util_file_prefix("esp.grid",esp_name)

        call prop_grid_write_cube_unformatted(geom,sub_nprp,
     >      dbl_mb(k_tepot), sub_i_prp_c, esp_name)
      end if


C
C Clean up phase
C
      if (.not.ma_pop_stack(l_epot)) call errquit
     &   ('hnd_elpiso, ma_pop_stack of l_epot failed',911,MA_ERR)

      if (.not.ma_pop_stack(l_tepot)) call errquit
     &   ('hnd_elpiso, ma_pop_stack of l_tepot failed',911,MA_ERR)

      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)

      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_elpiso, ma_pop_stack of l_xyzpt failed',911,MA_ERR)

C
C     Terminate integrals
C
      call schwarz_tidy()
      call int_terminate()

      end subroutine













C Returns the coordinates of grid points that lie within a TOLERANCE,
C on electron density isosurface of value ISOVALUE
C
C ISOVALUE == 0.002 ~ Vdw radius
C
C Mark J. Williamson Summer 2015
C
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
     &    errquit('hnd_elpiso: geom_ncent',911,GEOM_ERR)

      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_elpiso: ma failed',911,MA_ERR)

      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_elpiso: ma failed',911,MA_ERR)

C Copy the atoms' coordinates over to the k_xyzpt array
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_elpiso: geom_cent_get',911,GEOM_ERR)
   30 continue


C Set the grid up for the initial electron density calculation
      call prop_grid_initialize(rtdb,nat,dbl_mb(k_xyzpt))
C Get the number of points and their coordinate point index, i_prp_c
      call prop_grid_get_r_ptr(nprp,i_prp_c)

C allocate memory and an index, k_eden, to the density values
      if (.not. ma_push_get(mt_dbl,nprp,'eden',l_eden,k_eden))
     &    call errquit('hnd_elpiso: eden ma failed',911,MA_ERR)

C Calculates density matrix
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)

C Calculates the electronic contribution of the
C electronic integral defined by nder for a given density at
C the grid points defined in points.
C nder=-1; compute the 1-electron densit
      call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(i_prp_c),nprp, 
     & dbl_mb(k_eden),-1)

C Scan ahead and work out how many points (sub_nprp) lie on the isosurface
      sub_nprp = 0
      do i=1, nprp
        if ( abs(dbl_mb(k_eden + (i-1)) - isovalue) < tolerance ) then
          sub_nprp = sub_nprp + 1
        end if
      end do

C Allocate array for storing the location of these points
      if(.not.ma_alloc_get(MT_DBL, 3*sub_nprp, 'grid array',
     &      h_sub_prp_c, i_sub_prp_c) ) call errquit(
     &      'grid_initialize: unable to allocate heap space',
     &      0, MA_ERR)

      sub_nprp = 0
      do i=1, nprp
        if ( abs(dbl_mb(k_eden + (i-1)) - isovalue) < tolerance ) then
C Coordinates of grid point
           xp = dbl_mb(i_prp_c  +3*(i-1))
           yp = dbl_mb(i_prp_c+1+3*(i-1))
           zp = dbl_mb(i_prp_c+2+3*(i-1))

C Copy matching coordinates to new array
           dbl_mb(i_sub_prp_c  +3*(sub_nprp)) = xp
           dbl_mb(i_sub_prp_c+1+3*(sub_nprp)) = yp
           dbl_mb(i_sub_prp_c+2+3*(sub_nprp)) = zp

           sub_nprp = sub_nprp + 1
        end if
      end do
 

C Clean up phase
      call prop_grid_destroy()

      if (.not.ma_pop_stack(l_eden)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_eden failed',911,MA_ERR)

      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)

      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_elpmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
    
1000  format(4F13.8)
      end subroutine

