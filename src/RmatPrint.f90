module RmatPrint_m

use MatranUtil_m
use Rmat_m
implicit none

#ifdef OVERVIEW

   RmatPrint_m contains programs for pretty printing a Rmat or its
   array.  The generic function is print.  RmPrint prints a Rmat,
   including its nonnumeric components.  RmPrintArray prints only the
   array.

   Author: Pete Stewart
   Jun 17 2003

#endif

   interface Print
      module procedure RmPrintArray, RmPrint
   end interface Print

contains

   ! RmPrintArray prints the mxn subarray in the northwest 
   ! corner of A.  

   subroutine RmPrintArray(A, m, n, w, d, e, lw, nbl)
      real(wp), intent(in) :: A(:,:)
      integer,  intent(in) :: m  ! Number of rows to print
      integer,  intent(in) :: n  ! Number of columns to print
      integer,  intent(in) :: w  ! Print in e<w>.<d>e<e> format
      integer,  intent(in) :: d
      integer,  intent(in),&
                optional   :: e
      integer,  intent(in),&
                optional   :: lw ! Line width
      logical,  intent(in),&
                optional   :: nbl! If present and true, do not print
                                 ! a blank line before printing
                                 ! the array.


      integer :: lnwdth, i, j, jl, ju, nil, mm, rlw
      integer, parameter :: ROWMAXEXP = 5, LINEWIDTH = 80
      character(50) :: colhead, line, line1
      if (m<=0 .or. n<=0)&
           call MatranError('RmPrintArray in RmatPrint:&
                           & Illegal dimensions.')

      ! Determine line width (lnwdth)

      if (present(lw)) then
         lnwdth = lw
      else
         lnwdth = LINEWIDTH  
      end if


      ! Determine row label width (rlw)

      mm = m
      do i=1,ROWMAXEXP
         mm = mm/10
         if (mm == 0) exit
      end do
      if (i > ROWMAXEXP)&
         call MatranError('RmPrintArray in RmatPrint:&
                           & Row too long.')
      rlw = i+1

      ! Determine number of items in a line (nil)

      nil = (lnwdth-rlw-1)/w


      ! Compute the output formats

      write(colhead, "('(',i4,'x,',i4,'i',i4,')')") rlw+1, nil, w

      if (present(e)) then
         write(line1, "('(1p,i',i4,',1x,',i4,'e',i4,'.',i4,'e',i4,')')") &
                         rlw, nil, w, d, e
         write(line, "('(1p,',i4,'x,',i4,'e',i4,'.',i4,'e',i4,')')") &
                         rlw+1, nil, w, d, e
      else
         write(line1, "('(1p,i',i4,',1x,',i4,'e',i4,'.',i4,'e',i4,')')") &
                         rlw, nil, w, d, 3
         write(line, "('(1p,',i4,'x,',i4,'e',i4,'.',i4,'e',i4,')')") &
                         rlw+1, nil, w, d, 3
      end if


      if (.not.present(nbl)) then
         print *, ' '
      else
         if (.not.nbl) print *, ' '
      end if

      do i=1,m
         jl = 1
         do while (jl <= n)
            ju = min(jl+nil-1, n)
            write (*, colhead) (j, j=jl,ju)
            if (jl == 1) then
               write(*, line1) i, (A(i,j), j=jl,ju)
            else
               write(*, line) (A(i,j), j=jl,ju)
            end if
            jl = jl + nil
         end do
      end do

   end subroutine RmPrintArray

   ! RmPrint prints a Rmat, including its nonnumeric components.

   subroutine RmPrint(A, w, d, note, e, lw)
      type(Rmat), intent(in) :: A
      integer, intent(in)    :: w  ! Print in e<w>.<d>e<e> format
      integer, intent(in)    :: d
      character(*), intent(in), &
                 optional    :: note ! Annotation
      integer, intent(in), &
               optional :: e
      integer, intent(in), &
               optional :: lw ! Line width


      call GuardTemp(A)
      print *, ' '
      if (present(note)) then
         print *, note
      end if

      if (associated(A%temporary)) then
          print '(a,i0,a,i0,a,i0,a,i0,a,a2,a,l1,a,i0 )', &
               ' ', A%nrow, ' ', A%ncol, ' ', A%narow, ' ', A%nacol,&
               ' ', A%tag, ' ', A%adjustable, ' ', A%temporary
      else
          print '(a,i0,a,i0,a,i0,a,i0,a,a2,a,l1,a,i0 )', &
               ' ', A%nrow, ' ', A%ncol, ' ', A%narow, ' ', A%nacol,&
               ' ', A%tag, ' ', A%adjustable, ' ', 0
      end if

      if (.not.associated(A%a)) then
         print *, "Rmat array not associated."
      else
         if (A%nrow>0 .and. A%ncol>0) &
            call Print(A%a, A%nrow, A%ncol, w, d, e, lw, .true.)
      end if
      call CleanTemp(A)

   end subroutine RmPrint

end module RmatPrint_m



