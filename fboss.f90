module boss

    implicit none

    ! public boss_kernel

contains

pure function calc_angle(a, b, c) result(angle)

    implicit none

    double precision, intent(in), dimension(3) :: a
    double precision, intent(in), dimension(3) :: b
    double precision, intent(in), dimension(3) :: c

    double precision, dimension(3) :: v1
    double precision, dimension(3) :: v2

    double precision :: cos_angle
    double precision :: angle

    v1 = a - b
    v2 = c - b

    v1 = v1 / norm2(v1)
    v2 = v2 / norm2(v2)

    cos_angle = dot_product(v1,v2)

    ! Clipping
    if (cos_angle > 1.0d0) cos_angle = 1.0d0
    if (cos_angle < -1.0d0) cos_angle = -1.0d0

    angle = acos(cos_angle)

end function calc_angle

pure function calc_cos_angle(a, b, c) result(cos_angle)

    implicit none

    double precision, intent(in), dimension(3) :: a
    double precision, intent(in), dimension(3) :: b
    double precision, intent(in), dimension(3) :: c

    double precision, dimension(3) :: v1
    double precision, dimension(3) :: v2

    double precision :: cos_angle

    v1 = a - b
    v2 = c - b

    v1 = v1 / norm2(v1)
    v2 = v2 / norm2(v2)

    cos_angle = dot_product(v1,v2)

end function calc_cos_angle

function distance_matrix(coords) result(dmat)

    double precision, dimension(:,:), intent(in) :: coords

    double precision, allocatable, dimension(:,:) :: dmat

    integer :: n, i, j

    n = size(coords, dim=1)

    allocate(dmat(n,n))

    do i = 1, n
        do j = i, n

            dmat(i,j) = norm2(coords(i,:) - coords(j,:))
            dmat(j,i) = dmat(i,j)

        enddo
    enddo

end function distance_matrix

subroutine ang_cube(x, angmat)

    double precision, dimension(:,:), intent(in) :: x 

    double precision, dimension(:,:,:) :: angmat

    double precision :: theta

    integer :: n, i, j, k

    n = size(x, dim=1)

    angmat = 0.0d0

    do i = 1, n
        do j = 1, n
            if (i == j) cycle

            do k = j+1, n
                if (i == k) cycle

                theta = calc_angle(x(j, :), x(i, :), x(k, :))

                ! write (*,*) i,j,k,theta
                angmat(i,j,k) = theta
                angmat(i,k,j) = theta

            enddo
        enddo
    enddo

end subroutine ang_cube


subroutine initialize_two_body(coords, charges, dmat, distance, scaling, counts, max_id, power)

    implicit none
    
    double precision, dimension(:,:), intent(in) :: coords
    integer, dimension(:), intent(in) :: charges 
    double precision, dimension(:,:), intent(in) :: dmat

    double precision, dimension(:,:), intent(out) :: distance
    double precision, dimension(:,:), intent(out) :: scaling
    integer, dimension(:), intent(out) :: counts
    
    integer, intent(in) :: max_id
    double precision, intent(in) :: power

    integer :: i, j, n, key
    integer :: qi, qj
    
    n = size(charges)
    ! max_id = maxval(charges)

    distance = 0.0d0
    scaling  = 0.0d0
    counts   = 0


    do i = 1, n
        do j = i+1, n

            qi = charges(i)
            qj = charges(j)
            key = (min(qi, qj) - 1)*max_id + max(qi, qj)

            ! write (*,*) min(qi, qj), max(qi, qj), key, size(counts), size(distance, dim=1), size(distance, dim=2)

            counts(key) = counts(key) + 1
            ! write (*,*) key, counts(key), dmat(i,j)
            
            distance(key,counts(key)) = dmat(i,j)
            scaling(key,counts(key)) = dmat(i,j)**(-power)


        enddo
    enddo       
   
    ! write (*,*) counts 

    ! stop
end subroutine initialize_two_body

function two_body_kernel(two_body_distance1, two_body_scaling1, two_body_counts1, &
    two_body_distance2, two_body_scaling2, two_body_counts2, width) result(k)

    implicit none
    
    double precision, dimension(:,:) :: two_body_distance1
    double precision, dimension(:,:) :: two_body_scaling1
    integer, dimension(:) :: two_body_counts1

    double precision, dimension(:,:) :: two_body_distance2
    double precision, dimension(:,:) :: two_body_scaling2
    integer, dimension(:) :: two_body_counts2

   
    double precision :: width, cut
    double precision :: inv_width
    double precision :: k, mu, scaling, contrib

    double precision :: s11, s22, s12

    integer i, j, a, b
        
    cut = 4.0d0*width
    k = 0.0d0

    inv_width = -1.0d0 / (4.0d0 * width**2)

    ! j = 0

    do i = 1, size(two_body_counts1)

        if (two_body_counts1(i) < 1) cycle
        if (two_body_counts2(i) < 1) cycle

        do a = 1, two_body_counts1(i)
            do b = 1, two_body_counts2(i)

                mu = (two_body_distance1(i,a) - two_body_distance2(i,b))**2

                if (mu > cut) cycle
                scaling = two_body_scaling1(i,a) * two_body_scaling2(i,b)

                ! contrib = two_body_scaling1(i,a)**2 + two_body_scaling2(i,b)**2 - &
                !                 & 2.0d0 * exp(mu * inv_width) * scaling
                
                contrib = exp(mu * inv_width) * scaling
                ! s11 = lal

                k = k + contrib

                ! write (*,*), i, k, mu, scaling, contrib
                ! j = j + 1
            end do
        end do
    end do

    ! write (*,*) "TOTAL", j
end function two_body_kernel

subroutine initialize_three_body(coords, charges, dmat, angcube, distance, scaling, counts, max_id )

    implicit none
    
    double precision, dimension(:,:), intent(in) :: coords
    integer, dimension(:), intent(in) :: charges 
    double precision, dimension(:,:), intent(in) :: dmat
    double precision, dimension(:,:,:), intent(in) :: angcube

    double precision, dimension(:,:), intent(out) :: distance
    double precision, dimension(:,:), intent(out) :: scaling
    integer, dimension(:), intent(out) :: counts

    integer :: i, j, n, max_id, key, k
    integer :: qi, qj, qk

    double precision :: factor
    
    n = size(charges)

    distance = 0.0d0
    scaling  = 0.0d0
    counts   = 0


    do i = 1, n
        do j = 1, n
        if (i == j) cycle
            do k = j+1, n
                if (i == k) cycle


                qi = charges(i)
                qj = charges(j)
                qk = charges(k)
                key = three2key(qi, min(qj, qk), max(qj, qk), max_id)

                ! write (*,*) qi, min(qj, qk), max(qj, qk), key, size(counts), size(distance, dim=1), size(distance, dim=2)

                counts(key) = counts(key) + 1
                ! write (*,*) key, counts(key), angcube(i,j,k)
            ! 
                distance(key,counts(key)) = angcube(i,j,k)


                factor = (1.0d0 + 3.0d0 * cos(angcube(i,j,k))  * cos(angcube(k,i,j))  &
                       & * cos(angcube(j,k,i))) / (dmat(i,j) * dmat(j,k) * dmat(k,i))**3
                scaling(key,counts(key)) = factor


            enddo
        enddo
    enddo       
   
    ! write (*,*) counts 

    ! stop
end subroutine initialize_three_body

function three2key(a, b, c, max_id) result(key)

    implicit none

    integer, intent(in) :: a, b, c, max_id
    integer :: key

    key = (a-1)*max_id**2 + (b-1)*max_id + c

end function three2key


!function unique_two_body(nuclear_charges) result(unique)

!end function unique_two_body


end module boss


subroutine boss_kernel(X1, X2, nm1, nm2, sigma, d_width, d_power, k)

    use boss, only: distance_matrix, ang_cube, initialize_two_body, two_body_kernel, &
        & initialize_three_body

    implicit none

    double precision, dimension(:,:,:), intent(in) :: X1
    double precision, dimension(:,:,:), intent(in) :: X2
    
    integer, intent(in) :: nm1, nm2

    double precision, dimension(:,:,:), allocatable :: coordinates1
    double precision, dimension(:,:,:), allocatable :: coordinates2

    integer, dimension(:,:), allocatable :: nuclear_charges1 
    integer, dimension(:,:), allocatable :: nuclear_charges2
    
    double precision, dimension(nm1,nm2), intent(out) :: k

    double precision, dimension(:,:), allocatable :: dmat1 
    double precision, dimension(:,:), allocatable :: dmat2
    
    double precision, allocatable, dimension(:,:,:) :: two_body_distance1
    double precision, allocatable, dimension(:,:,:) :: two_body_scaling1
    integer, allocatable, dimension(:,:) :: two_body_counts1

    double precision, allocatable, dimension(:,:,:) :: two_body_distance2
    double precision, allocatable, dimension(:,:,:) :: two_body_scaling2
    integer, allocatable, dimension(:,:) :: two_body_counts2
    
    integer :: max_id
    integer :: max_size
    integer :: max_two_body
    integer :: max_three_body

    double precision, parameter :: pi = 4.d0 * atan(1.d0)

    ! double precision, parameter :: d_width = 0.1d0
    ! double precision, parameter :: d_power = 6.0d0
    ! double precision, parameter :: sigma = 100.0d0

    double precision, intent(in) :: d_width
    double precision, intent(in) :: d_power
    double precision, intent(in) :: sigma

    double precision, parameter :: t_width = 0.2d0
   
    logical, parameter :: three = .false.
    ! logical, parameter :: three = .true.

    double precision, parameter :: t_weight = 0.0d0


    double precision, dimension(nm1) :: s11
    double precision, dimension(nm2) :: s22
    double precision :: s12
    double precision :: l2

    integer, dimension(:), allocatable :: n1, n2

    integer :: na1, na2

    integer :: a, b, i
    
    double precision, dimension(:,:,:), allocatable :: acube1 
    double precision, dimension(:,:,:), allocatable :: acube2
    
    double precision, allocatable, dimension(:,:,:) :: three_body_distance1
    double precision, allocatable, dimension(:,:,:) :: three_body_scaling1
    integer, allocatable, dimension(:,:) :: three_body_counts1

    double precision, allocatable, dimension(:,:,:) :: three_body_distance2
    double precision, allocatable, dimension(:,:,:) :: three_body_scaling2
    integer, allocatable, dimension(:,:) :: three_body_counts2
    
    max_size = max(size(X1, dim=3),size(X2, dim=3))
    
    allocate(nuclear_charges1(max_size,nm1))
    allocate(nuclear_charges2(max_size,nm2))
   
    allocate(n1(nm1))
    allocate(n2(nm2))

    n1 = 0
    n2 = 0
    nuclear_charges1 = -1
    nuclear_charges2 = -1

    max_id = 0
    do a = 1, nm1
        do i = 1, max_size
            if (X1(a,1,i) > 0) then
                nuclear_charges1(i,a) = int(X1(a,1,i))
                n1(a) = n1(a) + 1
                max_id = max(max_id, nuclear_charges1(i,a))
            endif
        enddo
    enddo
    
    do a = 1, nm2
        do i = 1, max_size
            if (X2(a,1,i) > 0) then
                nuclear_charges2(i,a) = int(X2(a,1,i))
                n2(a) = n2(a) + 1
                max_id = max(max_id, nuclear_charges2(i,a))
            endif
        enddo
    enddo
    
    max_two_body = max(maxval(n1), maxval(n2))**2

    allocate(two_body_distance1(max_id**2, max_two_body,nm1))
    allocate(two_body_distance2(max_id**2, max_two_body,nm2))
    allocate(two_body_scaling1(max_id**2, max_two_body,nm1))
    allocate(two_body_scaling2(max_id**2, max_two_body,nm2))
    allocate(two_body_counts1(max_id**2,nm1))
    allocate(two_body_counts2(max_id**2,nm2))

    two_body_distance1 = 0.0d0
    two_body_distance2 = 0.0d0
    two_body_scaling1 = 0.0d0
    two_body_scaling2 = 0.0d0
    two_body_counts1 = 0
    two_body_counts2 = 0
    
   
    if (three) then 
        max_three_body = max(maxval(n1), maxval(n2))**3

        allocate(three_body_distance1(max_id**3, max_three_body,nm1))
        allocate(three_body_distance2(max_id**3, max_three_body,nm2))
        allocate(three_body_scaling1(max_id**3,  max_three_body,nm1))
        allocate(three_body_scaling2(max_id**3,  max_three_body,nm2))
        allocate(three_body_counts1(max_id**3,nm1))
        allocate(three_body_counts2(max_id**3,nm2))
        
        three_body_distance1 = 0.0d0
        three_body_distance2 = 0.0d0
        three_body_scaling1 = 0.0d0
        three_body_scaling2 = 0.0d0
        three_body_counts1 = 0
        three_body_counts2 = 0
        allocate(acube1(maxval(n1(:)),maxval(n1(:)),maxval(n1(:))))
        allocate(acube2(maxval(n2(:)),maxval(n2(:)),maxval(n2(:))))

        acube1 = 0.0d0
        acube2 = 0.0d0
    
        ! write(*,*) max_id**3,  max_three_body,nm1

    endif

    allocate(dmat1(maxval(n1(:)),maxval(n1(:))))
    allocate(dmat2(maxval(n2(:)),maxval(n2(:))))
    

    dmat1 = 0.0d0
    dmat2 = 0.0d0

    ! write (*,*) size(x1, dim=1), size(x1, dim=2),size(x1, dim=3)
    ! write (*,*) size(x2, dim=1), size(x2, dim=2),size(x2, dim=3)
    
    ! write(*,*) 11
    !$OMP PARALLEL DO PRIVATE(dmat1,acube1)
    do a = 1, nm1
        ! write(*,*) 11, a

        dmat1(:,:) = distance_matrix(transpose(X1(a,2:4,:n1(a))))

        call initialize_two_body(transpose(X1(a,2:4,:n1(a))), nuclear_charges1(:n1(a),a), dmat1(:n1(a),:n1(a)), &
         & two_body_distance1(:,:,a), two_body_scaling1(:,:,a), two_body_counts1(:,a), max_id, d_power)

        s11(a) =  two_body_kernel( &
            & two_body_distance1(:,:,a), two_body_scaling1(:,:,a), two_body_counts1(:,a), &
            & two_body_distance1(:,:,a), two_body_scaling1(:,:,a), two_body_counts1(:,a), d_width)

        if (three) then
     
            call initialize_three_body(transpose(X1(a,2:4,:n1(a))), nuclear_charges1(:n1(a),a), dmat1(:n1(a),:n1(a)), acube1, &
                 & three_body_distance1(:,:,a), three_body_scaling1(:,:,a), three_body_counts1(:,a), max_id)
     
            s11(a) = s11(a) + t_weight * two_body_kernel( &
                & three_body_distance1(:,:,a), three_body_scaling1(:,:,a), three_body_counts1(:,a), &
                & three_body_distance1(:,:,a), three_body_scaling1(:,:,a), three_body_counts1(:,a), t_width)
        endif 

    enddo
    !$OMP END PARALLEL DO

    ! write(*,*) 22
    !$OMP PARALLEL DO private(dmat2, acube2)
    do a = 1, nm2
        ! write(*,*) 22, a
        dmat2(:,:) = distance_matrix(transpose(X2(a,2:4,:n2(a))))

        call initialize_two_body(transpose(X2(a,2:4,:n2(a))), nuclear_charges2(:n2(a),a), dmat2(:n2(a),:n2(a)), &
            & two_body_distance2(:,:,a), two_body_scaling2(:,:,a), two_body_counts2(:,a), max_id, d_power)
        
        s22(a) =  two_body_kernel( &
            & two_body_distance2(:,:,a), two_body_scaling2(:,:,a), two_body_counts2(:,a), &
            & two_body_distance2(:,:,a), two_body_scaling2(:,:,a), two_body_counts2(:,a), d_width)
        
        if (three) then
            call  ang_cube(transpose(X2(a,2:4,:n2(a))), acube2)

            call initialize_three_body(transpose(X2(a,2:4,:n2(a))), nuclear_charges2(:n2(a),a), dmat2(:n2(a),:n2(a)), acube2, &
                & three_body_distance2(:,:,a), three_body_scaling2(:,:,a), three_body_counts2(:,a), max_id)

            s22(a) = s22(a) + t_weight * two_body_kernel( &
                & three_body_distance2(:,:,a), three_body_scaling2(:,:,a), three_body_counts2(:,a), &
                & three_body_distance2(:,:,a), three_body_scaling2(:,:,a), three_body_counts2(:,a), t_width)
        endif


    enddo
    !$OMP END PARALLEL DO


    ! write(*,*) 12
    k = 0.0d0

    !$OMP PARALLEL DO PRIVATE(s12,l2) schedule(dynamic)
    do a = 1, nm1
        do b = 1, nm2
             
    ! write(*,*) 12, a, b
            s12 =  two_body_kernel( &
         & two_body_distance1(:,:,a), two_body_scaling1(:,:,a), two_body_counts1(:,a), &
         & two_body_distance2(:,:,b), two_body_scaling2(:,:,b), two_body_counts2(:,b), d_width)
        
        if (three) then
            s12 = s12 + t_weight * two_body_kernel( &
                & three_body_distance1(:,:,a), three_body_scaling1(:,:,a), three_body_counts1(:,a), &
                & three_body_distance2(:,:,b), three_body_scaling2(:,:,b), three_body_counts2(:,b), t_width)
        endif

            l2 = s11(a) + s22(b) - 2.0d0 * s12

            k(a,b) = exp( - l2 / (2.0d0 * sigma**2))

        enddo
    enddo
    !$OMP END PARALLEL DO

    ! write(*,*) "END"
    
    deallocate(nuclear_charges1)
    deallocate(nuclear_charges2)
    deallocate(n1)
    deallocate(n2)

    deallocate(two_body_distance1)
    deallocate(two_body_distance2)
    deallocate(two_body_scaling1)
    deallocate(two_body_scaling2)
    deallocate(two_body_counts1)
    deallocate(two_body_counts2)
    
    deallocate(dmat1)
    deallocate(dmat2)

    if (three) then
        deallocate(three_body_distance1)
        deallocate(three_body_distance2)
        deallocate(three_body_scaling1)
        deallocate(three_body_scaling2)
        deallocate(three_body_counts1)
        deallocate(three_body_counts2)
        
        deallocate(acube1)
        deallocate(acube2)

    endif
end subroutine boss_kernel



