! Class-safe access to unchanged, otherwise unconstructible mixture kernels.
! This probe does not exercise the public constructor on the kernel path.
program meos_ideal_mixture_reference
  use eoslibinit, only: init_multiparameter
  use eos_parameters, only: meos_idealmix, single_eos_alloc
  use multiparameter_idealmix
  implicit none
  type(meos_idealmix) :: mix
  character(len=16), parameter :: names(5) = &
       [character(len=16) :: 'C3','R134A','N-H2','O-H2','P-H2']
  character(len=32) :: mode
  integer :: nc, phase, i, found, k
  integer, allocatable :: ids(:)
  real :: t, p, z, zt, zp, h, ht, hp, s, st, sp, g, gt, gp
  real :: gi, git, gip, physical_gp, rho, actual_p, pv_rho, cv
  real, allocatable :: x(:), zn(:), hn(:), sn(:), phi(:), phit(:), phip(:)
  real, allocatable :: phin(:,:), volumes(:), gs(:), checks(:,:), values(:)
  call get_command_argument(1,mode)
  if (trim(mode)=='public') then
    write(*,'(A)') 'CONSTRUCTOR_BEGIN'
    flush(6)
    call init_multiparameter('C3,R134A','NIST_MEOS_MIX','DEFAULT')
    write(*,'(A)') 'CONSTRUCTOR_RETURNED'
    stop
  endif
  read(*,*) nc,t,p,phase
  if (nc<1 .or. nc>64) error stop 'invalid constituent count'
  allocate(ids(nc),x(nc),zn(nc),hn(nc),sn(nc),phi(nc),phit(nc),phip(nc))
  allocate(phin(nc,nc),volumes(nc),gs(nc),checks(4,nc))
  read(*,*) ids
  read(*,*) x
  if (any(ids<1) .or. any(ids>5)) error stop 'invalid constituent id'
  ! Establish an ordinary active model for the pure init methods. Each
  ! internal constituent is then allocated with the supported pure label.
  call init_multiparameter('C3','NIST_MEOS','DEFAULT')
  allocate(mix%nist(nc))
  do i=1,nc
    call single_eos_alloc(trim(names(ids(i))),mix%nist(i)%meos,'NIST_MEOS')
    call mix%nist(i)%meos%init()
  enddo
  call calc_multiparameter_idealmix_zfac(nc,mix,t,p,x,phase,z,zt,zp,zn)
  call calc_multiparameter_idealmix_enthalpy(nc,mix,t,p,x,phase,h,ht,hp,hn)
  call calc_multiparameter_idealmix_entropy(nc,mix,t,p,x,phase,s,st,sp,sn)
  call calc_multiparameter_idealmix_Gres(nc,mix,t,p,x,phase,g,gt,gp)
  call calc_multiparameter_idealmix_fugacity(nc,mix,t,p,x,phase,phi,phit,phip,phin)
  physical_gp=0
  do i=1,nc
    call mix%nist(i)%meos%densitySolver(t,p,phase,rho,found)
    volumes(i)=1/rho
    call mix%nist(i)%meos%calc_resgibbs(t,p,[1.0],phase,gi,git,gip)
    gs(i)=gi
    physical_gp=physical_gp+x(i)*gip
    call mix%nist(i)%meos%mp_pressure(rho,t,actual_p,pv_rho)
    cv=mix%nist(i)%meos%cv(t,volumes(i))
    checks(:,i)=[actual_p,-rho*rho*pv_rho,cv,real(found)]
  enddo
  allocate(values(13+8*nc+nc*nc))
  k=1
  values(k:k+2)=[z,zt,zp]; k=k+3
  values(k:k+nc-1)=zn; k=k+nc
  values(k:k+2)=[h,ht,hp]; k=k+3
  values(k:k+nc-1)=hn; k=k+nc
  values(k:k+2)=[s,st,sp]; k=k+3
  values(k:k+nc-1)=sn; k=k+nc
  values(k:k+3)=[g,gt,gp,physical_gp]; k=k+4
  values(k:k+nc-1)=gs; k=k+nc
  values(k:k+nc-1)=phi; k=k+nc
  values(k:k+nc-1)=phit; k=k+nc
  values(k:k+nc-1)=phip; k=k+nc
  values(k:k+nc*nc-1)=reshape(phin,[nc*nc]); k=k+nc*nc
  values(k:k+nc-1)=volumes
  write(*,'(A)',advance='no') 'VALUES='
  write(*,'(*(ES26.17E3,1X))') values
  write(*,'(A)',advance='no') 'CHECKS='
  write(*,'(*(ES26.17E3,1X))') checks
  do i=1,nc
    deallocate(mix%nist(i)%meos)
  enddo
  deallocate(mix%nist)
end program
