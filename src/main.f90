program lasso_fortran
  use wrappers, only: lasso_for_acps
  implicit none

  integer, parameter :: lu = 10
  integer*8 :: natoms, nexp, nrows, ncols, naddsub, nadd, addmaxl, nsub
  integer*8 :: nmaxcoef
  integer*1, allocatable :: lmax(:)
  character*2, allocatable :: atoms(:)
  character(:), allocatable :: yaddnames(:)
  real*8, allocatable :: explist(:)
  real*8, allocatable :: w(:), x(:,:), yref(:), yempty(:), yadd(:,:), ynofit(:,:)
  real*8, allocatable :: maxcoef(:), wsqrt(:), y(:), beta(:)
  real*8 :: wrms

  integer*8 :: i
  integer :: ios

  open(lu,file="octavedump.dat",status="old",form="unformatted",access="stream",iostat=ios)
  if (ios /= 0) then
     write (*,*) "error opening octavedump.dat"
     stop 1
  end if

  ! integers
  read (lu) natoms, nexp, nrows, ncols, naddsub, nadd, addmaxl
  nsub = naddsub - nadd
  write (*,*) "## Reading from binary file:"
  write (*,*) "# atoms: ",natoms
  write (*,*) "# exponents: ",nexp
  write (*,*) "# rows: ",nrows
  write (*,*) "# columns: ",ncols
  write (*,*) "# additional method evaluations: ",naddsub

  ! atom names
  allocate(atoms(natoms))
  do i = 1, natoms
     read (lu) atoms(i)
  end do

  ! additional method names
  allocate(character(len=addmaxl) :: yaddnames(nadd))
  do i = 1, nadd
     read (lu) yaddnames(i)
  end do

  ! small data arrays
  allocate(lmax(natoms),explist(nexp))
  read (lu) lmax, explist

  ! large data arrays
  allocate(w(nrows),x(nrows,ncols),yref(nrows),yempty(nrows),yadd(nrows,nadd),&
     ynofit(nrows,nsub))
  read (lu) w, x, yref, yempty
  if (nadd > 0) then
     read (lu) yadd
  else
     yadd = 0d0
  end if
  if (nsub > 0) then
     read (lu) ynofit
  else
     ynofit =0d0
  endif

  ! maxcoef
  read (lu) nmaxcoef
  allocate(maxcoef(nmaxcoef))
  if (nmaxcoef > 0) then
     read (lu) maxcoef
     if (nmaxcoef /= ncols) then
        maxcoef = 0d0
     end if
  else
     maxcoef = 0d0
  end if

  write (*,*) "# maximum coefficients: ",nmaxcoef

  ! apply the weights and transform the matrices for the fit
  allocate(wsqrt(nrows),y(nrows))
  wsqrt = sqrt(w)
  do i = 1, nrows
     x(i,:) = x(i,:) * wsqrt(i)
     yadd(i,:) = yadd(i,:) * wsqrt(i)
  end do
  y = yref - yempty
  do i = 1, nsub
     y = y - ynofit(:,i)
  end do
  y = y * wsqrt

  write (*,*) "## DONE"
  close(lu)

  beta = lasso_for_acps(int(nrows,4),int(ncols,4),x,y,10d0,maxcoef=maxcoef)
  wrms = sqrt(sum((y - matmul(x,beta))**2))
  write (*,*) "sum(abs(beta)) = ", sum(abs(beta))
  write (*,*) "sum(beta^2) = ", sum(beta**2)
  write (*,*) "max(abs(beta)) = ", maxval(abs(beta))
  write (*,*) "wrms = ", wrms
  write (*,*) "non-zero terms = ", count(abs(beta) > 1d-9)

end program lasso_fortran

