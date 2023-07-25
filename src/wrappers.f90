module wrappers
  implicit none

  private
  public :: lasso_for_acps
  public :: lasso_cv

contains

  ! Run LASSO for ACP development. The Fortran version of runlasso.m
  ! in acpdb.
  function lasso_for_acps(rows,cols,x,y,t,maxcoef,yadd,verbose) result(w)
    use lassofun, only: lasso
    integer, intent(in) :: rows, cols
    real*8, intent(in) :: x(rows,cols), y(rows)
    real*8, intent(in) :: t
    real*8, intent(in), optional :: maxcoef(cols)
    real*8, intent(in), optional :: yadd(rows,*)
    logical, intent(in), optional :: verbose
    real*8 :: w(cols)

    real*8, allocatable :: factor(:), xtilde(:,:)
    logical, allocatable :: idx(:)
    logical :: verbose_, maxcoef_
    real*8 :: ratio, wrms
    integer :: i

    verbose_ = .false.
    if (present(verbose)) verbose_ = verbose
    maxcoef_ = .false.
    if (present(maxcoef)) maxcoef_ = .true.

    if (present(yadd)) then
       write (*,*) "not implemented yet!"
       stop 1
    end if
    if (maxcoef_) then
       ! initialize
       allocate(factor(cols),xtilde(rows,cols))
       factor = 1d0
       xtilde = 0d0

       do while (.true.)
          do i = 1, cols
             xtilde(:,i) = x(:,i) * factor(i)
          end do
          w = lasso(xtilde,y,t)
          do i = 1, cols
             w(i) = w(i) * factor(i)
          end do
          wrms = sqrt(sum((y - matmul(x,w))**2))
          idx = (abs(w) > maxcoef)
          if (verbose_) then
             write (*,'("norm-1 = ",F12.4,", wrms = ",F12.4,", maxcoef: ok = ",I4," / violations = ",I4)') &
                sum(abs(w)),wrms,cols-count(idx),count(idx)
          end if
          if (count(idx) == 0) exit

          do i = 1, cols
             if (abs(w(i)) > maxcoef(i)) then
                ratio = min(maxcoef(i) / abs(w(i)),0.99)
                factor(i) = factor(i) * ratio**(1.1d0)
             end if
          end do
       end do
       deallocate(xtilde,factor)
    else
       w = lasso(x,y,t)
    end if

  end function lasso_for_acps

  ! Determine the constraint by shuffle-split cross-validation using
  ! bracketing and golden section search. Parameters are hard-wired
  ! for now.
  subroutine lasso_cv(rows,cols,x,y,t,w,verbose)
    use lassofun, only: lasso
    integer, intent(in) :: rows, cols
    real*8, intent(in) :: x(rows,cols), y(rows)
    real*8, intent(out) :: t
    real*8, intent(out) :: w(cols)
    logical, intent(in) :: verbose

    real*8, parameter :: train_fraction = 0.80
    integer, parameter :: nsplit = 8
    real*8, parameter :: gr = (1d0 + sqrt(5d0)) / 2d0
    real*8, parameter :: etthres = 0.01d0

    real*8 :: urand, tt, wrms, wrmsn(nsplit)
    integer, allocatable :: s(:), r(:), seed(:)
    integer :: i, j, rows_t, n, et, nstep
    integer, allocatable :: usetrain(:,:), usevalid(:,:)
    real*8, allocatable :: xt(:,:), yt(:), xv(:,:), yv(:)
    integer, parameter :: et_min = -8
    integer, parameter :: et_max = 5
    logical, allocatable :: aux(:)
    real*8 :: tsave(3), wsave(3), etg(4), fetg(4), etgbest, fetgbest, tbest

    ! initialize random seed
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    write (*,'("* Started LASSO-CV")')
    write (*,*) "* Seed = ", seed

    ! generate training/validation data
    rows_t = int(train_fraction * rows)
    allocate(usetrain(rows_t,nsplit),usevalid(rows-rows_t,nsplit),s(rows),r(rows_t),aux(rows))
    do n = 1, nsplit
       do i = 1, rows
          s(i) = i
       end do
       r = s(1:rows_t)
       do i = rows_t+1, rows
          call random_number(urand)
          j = floor(i*urand) + 1
          if (j <= rows_t) r(j) = s(i)
       end do
       usetrain(:,n) = r
       aux = .true.
       do i = 1, rows_t
          aux(r(i)) = .false.
       end do
       usevalid(:,n) = pack(s,aux)
    end do
    deallocate(s,r,aux)

    ! bracket the t
    tsave = -1d0
    wsave = 0d0
    write (*,'("* Cross-validation: bracketing the constraint")')
    write (*,'("#exp  lambda     mean-wrms   ------                 individual wrms           ------")')
    allocate(xt(rows_t,cols),yt(rows_t),xv(rows_t,cols),yv(rows_t))
    do et = et_min, et_max
       tt = 10d0**et
       do n = 1, nsplit
          call calculate_wrms(n,wrmsn(n),real(et,8),(et==et_min .and. n==1))
       end do
       wrms = sum(wrmsn) / real(nsplit,8)
       write (*,'(I2,X,1p,E10.2,X,F10.2,"    ",999(F8.2,X))') et, tt, wrms, wrmsn

       tsave(1:2) = tsave(2:3)
       tsave(3) = tt
       wsave(1:2) = wsave(2:3)
       wsave(3) = wrms
       if (wsave(2) < wsave(1) .and. wsave(2) < wsave(3)) then
          if (any(tsave < 0)) then
             write (*,*) "error: minimum found too early in bracketing"
             stop 1
          end if
          exit
       end if
    end do

    ! golden section search
    write (*,'("* Done with bracketing")')
    etg(1) = log10(tsave(1))
    etg(4) = log10(tsave(3))
    fetg(1) = wsave(1)
    fetg(4) = wsave(3)
    write (*,'("* Starting golden section search in exponent range: ",2(F10.2,X))') etg(1), etg(4)
    write (*,'("* Required steps to convergence: ",I5)') ceiling(log(etthres / (etg(4)-etg(1))) / log(1d0/gr))
    etg(2) = etg(4) - (etg(4) - etg(1)) / gr
    call calculate_mean_wrms(fetg(2),etg(2),.true.)
    etg(3) = etg(1) + (etg(4) - etg(1)) / gr
    call calculate_mean_wrms(fetg(3),etg(3),.false.)
    write (*,'("#nstep     exp-1      wrms-1        exp-2      wrms-2        exp-3      wrms-3        exp-4      wrms-4")')
    nstep = 0
    write (*,'(I3,X,4(2X,F10.4,X,F10.4,2X))') nstep, etg(1), fetg(1),&
       etg(2), fetg(2), etg(3), fetg(3), etg(4), fetg(4)
    do while (abs(etg(4) - etg(1)) > etthres)
       if (fetg(2) < fetg(3)) then
          etg(4) = etg(3)
          fetg(4) = fetg(3)
          etg(3) = etg(2)
          fetg(3) = fetg(2)
          etg(2) = etg(4) - (etg(4) - etg(1)) / gr
          call calculate_mean_wrms(fetg(2),etg(2),.false.)
       else
          etg(1) = etg(2)
          fetg(1) = fetg(2)
          etg(2) = etg(3)
          fetg(2) = fetg(3)
          etg(3) = etg(1) + (etg(4) - etg(1)) / gr
          call calculate_mean_wrms(fetg(3),etg(3),.false.)
       end if
       nstep = nstep + 1
       write (*,'(I3,X,4(2X,F10.4,X,F10.4,2X))') nstep, etg(1), fetg(1),&
          etg(2), fetg(2), etg(3), fetg(3), etg(4), fetg(4)
    end do
    if (fetg(2) < fetg(3)) then
       etgbest = etg(2)
       fetgbest = fetg(2)
    else
       etgbest = etg(3)
       fetgbest = fetg(3)
    end if
    tbest = 10d0**(etgbest)
    write (*,'("* Golden section search converged with: ")')
    write (*,'(" exponent = ",F10.4)') etgbest
    write (*,'(" constraint = ",1p,E14.4)') tbest
    write (*,'(" wrms = ",F10.4)') fetgbest

    ! final fit
    write (*,'("* Final fit with the whole training set: ")')
    t = tbest
    w = lasso(x,y,t,w0=w)
    wrms = sqrt(sum((yv - matmul(xv,w))**2))
    write (*,'("* FINAL WRMS = ",F10.4/)') wrms

  contains
    subroutine calculate_wrms(n,wrms,et,cold)
      integer, intent(in) :: n
      real*8, intent(out) :: wrms
      real*8, intent(in) :: et
      logical, intent(in) :: cold

      real*8 :: tt

      tt = 10d0**et
      xt = x(usetrain(:,n),:)
      yt = y(usetrain(:,n))
      xv = x(usevalid(:,n),:)
      yv = y(usevalid(:,n))
      if (cold) then
         w = lasso(xt,yt,tt)
      else
         w = lasso(xt,yt,tt,w0=w)
      end if
      wrms = sqrt(sum((yv - matmul(xv,w))**2))

    end subroutine calculate_wrms

    subroutine calculate_mean_wrms(wrms,et,cold)
      real*8, intent(out) :: wrms
      real*8, intent(in) :: et
      logical, intent(in) :: cold

      real*8 :: wrmsn(nsplit)

      do n = 1, nsplit
         call calculate_wrms(n,wrmsn(n),et,cold)
      end do
      wrms = sum(wrmsn) / real(nsplit,8)

    end subroutine calculate_mean_wrms

  end subroutine lasso_cv

end module wrappers
