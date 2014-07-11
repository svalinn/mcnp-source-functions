!+ $Id: source.F90,v 1.2 2006/10/03 02:10:16 mashnik Exp $
! Copyright LANS/LANL/DOE - see file COPYRIGHT_INFO

subroutine source
  ! dummy subroutine.  aborts job if source subroutine is missing.
  ! if nsr=0, subroutine source must be furnished by the user.
  ! at entrance, a random set of uuu,vvv,www has been defined.  the
  ! following variables must be defined within the subroutine:
  ! xxx,yyy,zzz,icl,jsu,erg,wgt,tme and possibly ipt,uuu,vvv,www.
  ! subroutine srcdx may also be needed.
  use mcnp_global
  use mcnp_debug

  implicit real(dknd) (a-h,o-z)
  logical, save :: first_run = .true.
  real(dknd), save :: bin_width
  real(dknd) :: radial_sample ! sampled radial position
  integer :: radial_index ! the index of the sampled position
  real(dknd) :: r,z ! the r & z coords
  real(dknd) :: angle
  integer,save :: icl_tmp
  integer :: j
  
  if(first_run.eqv..true.)then
  ! sets up the plasma source
    call setup(rdum(1),rdum(2),rdum(3),rdum(4),rdum(5),rdum(6), &
               rdum(7),rdum(8),rdum(9),rdum(10),rdum(11),rdum(12), &
	       rdum(13),rdum(14)," ",idum(1),idum(2))

   ! set the grid spacing	   
   bin_width = rdum(10)/real(idum(2))	       	      

   do m=1,mxa
     if(idum(3).eq.ncl(m))then
        icl_tmp = m
        cycle
     endif
   enddo

  endif	       

! resample source point
200 continue     
  ! sample source pos
  call sample_source(idum(2),bin_width,radial_sample,radial_index,rang(),rang())

  ! convert to cylindrical
  ! minor rad = 10, major = 11, elongation = 12, triangularity = 13
  ! shaf =14
  call convert_rad_to_rz(rdum(11),rdum(10),rdum(12),rdum(13),rdum(14), &
                         radial_sample,rang(),r,z)
  ! convert to xyz
  call convert_rz_to_xyz(r,rang(),xxx,yyy,rdum(15),rdum(16),angle)
  xxx=100.*xxx
  yyy=100.*yyy
  zzz=100.*z

  call chkcel(icl_tmp,2,j)
  if(j.ne.0)then
    goto 200
  endif

  call sample_energy(radial_index,rang(),rang(),erg)

  wgt=1.0
  jsu=0
  tme=0.0

  icl=icl_tmp

  ! 
  return
end subroutine source
