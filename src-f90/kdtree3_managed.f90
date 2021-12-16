!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.
!
! module kdtree2_precision_module
  
!   integer, parameter :: sp = kind(0.0)
!   integer, parameter :: dp = kind(0.0d0)

!   private :: sp, dp

!   !
!   ! You must comment out exactly one
!   ! of the two lines.  If you comment
!   ! out kdkind = sp then you get single precision
!   ! and if you comment out kdkind = dp 
!   ! you get double precision.
!   !

!   ! integer, parameter :: kdkind = sp  
!   integer, parameter :: kdkind = dp  
!   public :: kdkind

! end module kdtree2_precision_module


module kdtree3_module

#ifdef OP2_WITH_CUDAFOR
  use CUDAFOR
#endif

  ! use kdtree2_precision_module
  ! use kdtree2_priority_queue_module
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric is supported. 
  !
  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module. 
  !
  !-------------DATA TYPE, CREATION, DELETION---------------------
  ! public :: kdkind
  public :: kdtree2, kdtree2_result, tree_node, kdtree2_create


  integer, parameter :: bucket_size = 12
  ! The maximum number of points to keep in a terminal node.

  type interval
      real(8), allocatable, managed :: lower(:), upper(:)
  end type interval

  type :: tree_node
      ! an internal tree node
      private
      integer, allocatable, managed :: cut_dim(:)
      ! the dimension to cut
      real(8), allocatable, managed  :: cut_val(:)
      ! where to cut the dimension
      real(8), allocatable, managed  :: cut_val_left(:), cut_val_right(:) 
      ! improved cutoffs knowing the spread in child boxes.
      integer, allocatable, managed   :: l(:), u(:)
      type (tree_node), pointer, managed :: left(:), right(:)
      type(interval), pointer, managed :: box(:) => null()
      ! child pointers
      ! Points included in this node are indexes[k] with k \in [l,u] 


  end type tree_node

  type :: kdtree2
      ! Global information about the tree, one per tree
      integer, allocatable, managed  :: dimen(:)=0, n(:)=0
      ! dimensionality and total # of points
      real(8), pointer, managed :: the_data(:,:) => null()
      ! pointer to the actual data array 
      ! 
      !  IMPORTANT NOTE:  IT IS DIMENSIONED   the_data(1:d,1:N)
      !  which may be opposite of what may be conventional.
      !  This is, because in Fortran, the memory layout is such that
      !  the first dimension is in sequential order.  Hence, with
      !  (1:d,1:N), all components of the vector will be in consecutive
      !  memory locations.  The search time is dominated by the
      !  evaluation of distances in the terminal nodes.  Putting all
      !  vector components in consecutive memory location improves
      !  memory cache locality, and hence search speed, and may enable 
      !  vectorization on some processors and compilers. 

      integer, pointer, managed :: ind(:) => null()
      ! permuted index into the data, so that indexes[l..u] of some
      ! bucket represent the indexes of the actual points in that
      ! bucket.
      logical, allocatable, managed        :: sort(:) = .false.
      ! do we always sort output results?
      logical, allocatable, managed        :: rearrange(:) = .false. 
      real(8), pointer, managed :: rearranged_data(:,:) => null()
      ! if (rearrange .eqv. .true.) then rearranged_data has been
      ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
      ! permitting search to use more cache-friendly rearranged_data, at
      ! some initial computation and storage cost.
      type (tree_node), pointer, managed :: root(:) => null()
      ! root pointer of the tree
  end type kdtree2


  private
  ! everything else is private.

contains

  function kdtree2_create(input_data,dim,sort,rearrange) result (mr)
    !
    ! create the actual tree structure, given an input array of data.
    !
    ! Note, input data is input_data(1:d,1:N), NOT the other way around.
    ! THIS IS THE REVERSE OF THE PREVIOUS VERSION OF THIS MODULE.
    ! The reason for it is cache friendliness, improving performance.
    !
    ! Optional arguments:  If 'dim' is specified, then the tree
    !                      will only search the first 'dim' components
    !                      of input_data, otherwise, dim is inferred
    !                      from SIZE(input_data,1).
    !
    !                      if sort .eqv. .true. then output results
    !                      will be sorted by increasing distance.
    !                      default=.false., as it is faster to not sort.
    !                      
    !                      if rearrange .eqv. .true. then an internal
    !                      copy of the data, rearranged by terminal node,
    !                      will be made for cache friendliness. 
    !                      default=.true., as it speeds searches, but
    !                      building takes longer, and extra memory is used.
    !
    ! .. Function Return Cut_value ..
    type (kdtree2), pointer  :: mr
    integer, intent(in), optional      :: dim
    logical, intent(in), optional      :: sort
    logical, intent(in), optional      :: rearrange
    ! ..
    ! .. Array Arguments ..
    real(8), target :: input_data(:,:)
    !
    integer :: i
    ! ..
    allocate (mr)
    mr%the_data => input_data
    ! pointer assignment

    if (present(dim)) then
       mr%dimen = dim
    else
       mr%dimen = size(input_data,1)
    end if
    mr%n = size(input_data,2)

    if (mr%dimen > mr%n) then
       !  unlikely to be correct
       write (*,*) 'KD_TREE_TRANS: likely user error.'
       write (*,*) 'KD_TREE_TRANS: You passed in matrix with D=',mr%dimen
       write (*,*) 'KD_TREE_TRANS: and N=',mr%n
       write (*,*) 'KD_TREE_TRANS: note, that new format is data(1:D,1:N)'
       write (*,*) 'KD_TREE_TRANS: with usually N >> D.   If N =approx= D, then a k-d tree'
       write (*,*) 'KD_TREE_TRANS: is not an appropriate data structure.'
       stop
    end if

    call build_tree(mr)

    if (present(sort)) then
       mr%sort = sort
    else
       mr%sort = .false.
    endif

    if (present(rearrange)) then
       mr%rearrange = rearrange
    else
       mr%rearrange = .true.
    endif

    if (mr%rearrange) then
       allocate(mr%rearranged_data(mr%dimen,mr%n))
       do i=1,mr%n
          mr%rearranged_data(:,i) = mr%the_data(:, &
           mr%ind(i))
       enddo
    else
       nullify(mr%rearranged_data)
    endif

  end function kdtree2_create

    subroutine build_tree(tp)
      type (kdtree2), pointer :: tp
      ! ..
      integer :: j
      type(tree_node), pointer :: dummy => null()
      ! ..
      allocate (tp%ind(tp%n))
      forall (j=1:tp%n)
         tp%ind(j) = j
      end forall
      tp%root => build_tree_for_range(tp,1,tp%n, dummy)
    end subroutine build_tree

    recursive function build_tree_for_range(tp,l,u,parent) result (res)
      ! .. Function Return Cut_value ..
      type (tree_node), pointer :: res
      ! ..
      ! .. Structure Arguments ..
      type (kdtree2), pointer :: tp
      type (tree_node),pointer           :: parent
      ! ..
      ! .. Scalar Arguments ..
      integer, intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      integer :: i, c, m, dimen
      logical :: recompute
      real(8)    :: average

!!$      If (.False.) Then 
!!$         If ((l .Lt. 1) .Or. (l .Gt. tp%n)) Then
!!$            Stop 'illegal L value in build_tree_for_range'
!!$         End If
!!$         If ((u .Lt. 1) .Or. (u .Gt. tp%n)) Then
!!$            Stop 'illegal u value in build_tree_for_range'
!!$         End If
!!$         If (u .Lt. l) Then
!!$            Stop 'U is less than L, thats illegal.'
!!$         End If
!!$      Endif
!!$      
      ! first compute min and max
      dimen = tp%dimen
      allocate (res)
      allocate(res%box(dimen))

      ! First, compute an APPROXIMATE bounding box of all points associated with this node.
      if ( u < l ) then
         ! no points in this box
         nullify(res)
         return
      end if

      if ((u-l)<=bucket_size) then
         !
         ! always compute true bounding box for terminal nodes.
         !
         do i=1,dimen
            call spread_in_coordinate(tp,i,l,u,res%box(i))
         end do
         res%cut_dim = 0
         res%cut_val = 0.0
         res%l = l
         res%u = u
         res%left =>null()
         res%right => null() 
      else
         ! 
         ! modify approximate bounding box.  This will be an
         ! overestimate of the true bounding box, as we are only recomputing 
         ! the bounding box for the dimension that the parent split on.
         !
         ! Going to a true bounding box computation would significantly
         ! increase the time necessary to build the tree, and usually
         ! has only a very small difference.  This box is not used
         ! for searching but only for deciding which coordinate to split on.
         !
         do i=1,dimen
            recompute=.true.
            if (associated(parent)) then
               if (i .ne. parent%cut_dim) then
                  recompute=.false.
               end if
            endif
            if (recompute) then
               call spread_in_coordinate(tp,i,l,u,res%box(i))
            else
               res%box(i) = parent%box(i)
            endif
         end do
         

         c = maxloc(res%box(1:dimen)%upper-res%box(1:dimen)%lower,1)
         !
         ! c is the identity of which coordinate has the greatest spread.
         !
         
         if (.false.) then
            ! select exact median to have fully balanced tree.
            m = (l+u)/2
            call select_on_coordinate(tp%the_data,tp%ind,c,m,l,u)
         else
            !
            ! select point halfway between min and max, as per A. Moore,
            ! who says this helps in some degenerate cases, or 
            ! actual arithmetic average. 
            !
            if (.true.) then
               ! actually compute average
               average = sum(tp%the_data(c,tp%ind(l:u))) / real(u-l+1,8)
            else
               average = (res%box(c)%upper + res%box(c)%lower)/2.0
            endif
               
            res%cut_val = average
            m = select_on_coordinate_value(tp%the_data,tp%ind,c,average,l,u)
         endif
            
         ! moves indexes around
         res%cut_dim = c
         res%l = l
         res%u = u
!         res%cut_val = tp%the_data(c,tp%ind(m))

         res%left => build_tree_for_range(tp,l,m,res)
         res%right => build_tree_for_range(tp,m+1,u,res)

         if (associated(res%right) .eqv. .false.) then
            res%box = res%left%box
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = res%cut_val_left
         elseif (associated(res%left) .eqv. .false.) then
            res%box = res%right%box
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val = res%cut_val_right
         else
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = (res%cut_val_left + res%cut_val_right)/2


            ! now remake the true bounding box for self.  
            ! Since we are taking unions (in effect) of a tree structure,
            ! this is much faster than doing an exhaustive
            ! search over all points
            res%box%upper = max(res%left%box%upper,res%right%box%upper)
            res%box%lower = min(res%left%box%lower,res%right%box%lower) 
         endif
      end if
    end function build_tree_for_range

    integer function select_on_coordinate_value(v,ind,c,alpha,li,ui) &
     result(res)
      ! Move elts of ind around between l and u, so that all points
      ! <= than alpha (in c cooordinate) are first, and then
      ! all points > alpha are second. 

      !
      ! Algorithm (matt kennel). 
      !
      ! Consider the list as having three parts: on the left,
      ! the points known to be <= alpha.  On the right, the points
      ! known to be > alpha, and in the middle, the currently unknown
      ! points.   The algorithm is to scan the unknown points, starting
      ! from the left, and swapping them so that they are added to
      ! the left stack or the right stack, as appropriate.
      ! 
      ! The algorithm finishes when the unknown stack is empty. 
      !
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, li, ui
      real(8), intent(in) :: alpha
      ! ..
      real(8) :: v(1:,1:)
      integer :: ind(1:)
      integer :: tmp  
      ! ..
      integer :: lb, rb
      !
      ! The points known to be <= alpha are in
      ! [l,lb-1]
      !
      ! The points known to be > alpha are in
      ! [rb+1,u].  
      !
      ! Therefore we add new points into lb or
      ! rb as appropriate.  When lb=rb
      ! we are done.  We return the location of the last point <= alpha.
      !
      ! 
      lb = li; rb = ui

      do while (lb < rb)
         if ( v(c,ind(lb)) <= alpha ) then
            ! it is good where it is.
            lb = lb+1
         else
            ! swap it with rb.
            tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
            rb = rb-1
         endif
      end do
      
      ! now lb .eq. ub 
      if (v(c,ind(lb)) <= alpha) then
         res = lb
      else
         res = lb-1
      endif
      
    end function select_on_coordinate_value

    subroutine select_on_coordinate(v,ind,c,k,li,ui)
      ! Move elts of ind around between l and u, so that the kth
      ! element
      ! is >= those below, <= those above, in the coordinate c.
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, k, li, ui
      ! ..
      integer :: i, l, m, s, t, u
      ! ..
      real(8) :: v(:,:)
      integer :: ind(:)
      ! ..
      l = li
      u = ui
      do while (l<u)
         t = ind(l)
         m = l
         do i = l + 1, u
            if (v(c,ind(i))<v(c,t)) then
               m = m + 1
               s = ind(m)
               ind(m) = ind(i)
               ind(i) = s
            end if
         end do
         s = ind(l)
         ind(l) = ind(m)
         ind(m) = s
         if (m<=k) l = m + 1
         if (m>=k) u = m - 1
      end do
    end subroutine select_on_coordinate

   subroutine spread_in_coordinate(tp,c,l,u,interv) 
      ! the spread in coordinate 'c', between l and u. 
      !
      ! Return lower bound in 'smin', and upper in 'smax', 
      ! ..
      ! .. Structure Arguments ..
      type (kdtree2), pointer :: tp
      type(interval), intent(out) :: interv
      ! ..
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, l, u
      ! ..
      ! .. Local Scalars ..
      real(8) :: last, lmax, lmin, t, smin,smax
      integer :: i, ulocal
      ! ..
      ! .. Local Arrays ..
      real(8), pointer :: v(:,:)
      integer, pointer :: ind(:)
      ! ..
      v => tp%the_data(1:,1:)
      ind => tp%ind(1:)
      smin = v(c,ind(l))
      smax = smin

      ulocal = u

      do i = l + 2, ulocal, 2
         lmin = v(c,ind(i-1))
         lmax = v(c,ind(i))
         if (lmin>lmax) then
            t = lmin
            lmin = lmax
            lmax = t
         end if
         if (smin>lmin) smin = lmin
         if (smax<lmax) smax = lmax
      end do
      if (i==ulocal+1) then
         last = v(c,ind(ulocal))
         if (smin>last) smin = last
         if (smax<last) smax = last
      end if

      interv%lower = smin
      interv%upper = smax

    end subroutine spread_in_coordinate



end module kdtree3_module

