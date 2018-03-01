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

function ang_cube(x) result(angmat)

    double precision, dimension(:,:), intent(in) :: x 

    double precision, allocatable, dimension(:,:,:) :: angmat

    double precision :: theta

    integer :: n, i, j, k

    n = size(x, dim=1)

    allocate(angmat(n,n,n))
    angmat = 0.0d0

    do i = 1, n
        do j = 1, n
            if (i == j) cycle

            do k = j+1, n
                if (i == k) cycle

                theta = calc_angle(x(j, :), x(i, :), x(k, :))

                ! write (*,*) i,j,k,theta
                angmat(i,j,k) = theta
                ! angmat(i,k,j) = theta

            enddo
        enddo
    enddo

end function ang_cube


subroutine initialize_two_body(coords, charges, dmat, distance, scaling, counts, max_id )

    implicit none
    
    double precision, dimension(:,:), intent(in) :: coords
    integer, dimension(:), intent(in) :: charges 
    double precision, dimension(:,:), intent(in) :: dmat

    double precision, dimension(:,:), intent(out) :: distance
    double precision, dimension(:,:), intent(out) :: scaling
    integer, dimension(:), intent(out) :: counts

    integer :: i, j, n, max_id, key
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
            scaling(key,counts(key)) = dmat(i,j)**(-6)


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

    integer i, j, a, b
        
    cut = 4.0d0*width
    k = 0.0d0

    inv_width = -1.0d0 / width**2

    j = 0

    do i = 1, size(two_body_counts1)

        if (two_body_counts1(i) < 1) cycle
        if (two_body_counts2(i) < 1) cycle

        do a = 1, two_body_counts1(i)
            do b = 1, two_body_counts2(i)

                mu = (two_body_distance1(i,a) - two_body_distance2(i,b))**2

                if (mu > cut) cycle
                scaling = two_body_scaling1(i,a) * two_body_scaling2(i,b)

                contrib = exp(mu * inv_width) * scaling
                
                k = k + contrib

                ! write (*,*), i, k, mu, scaling, contrib
                j = j + 1
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
    max_id = maxval(charges)

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

subroutine boss_kernel(coordinates1, nuclear_charges1, coordinates2, nuclear_charges2, k)

    use boss, only: distance_matrix, ang_cube, initialize_two_body, two_body_kernel, &
        & initialize_three_body

    implicit none

    double precision, dimension(:,:), intent(in) :: coordinates1
    double precision, dimension(:,:), intent(in) :: coordinates2

    integer, dimension(:),   intent(in) :: nuclear_charges1 
    integer, dimension(:),   intent(in) :: nuclear_charges2
    
    double precision, intent(out) :: k

    double precision, dimension(:,:), allocatable :: dmat1 
    double precision, dimension(:,:), allocatable :: dmat2
    
    double precision, dimension(:,:,:), allocatable :: acube1 
    double precision, dimension(:,:,:), allocatable :: acube2

    double precision, allocatable, dimension(:,:) :: two_body_distance1
    double precision, allocatable, dimension(:,:) :: two_body_scaling1
    integer, allocatable, dimension(:) :: two_body_counts1

    double precision, allocatable, dimension(:,:) :: two_body_distance2
    double precision, allocatable, dimension(:,:) :: two_body_scaling2
    integer, allocatable, dimension(:) :: two_body_counts2
    
    double precision, allocatable, dimension(:,:) :: three_body_distance1
    double precision, allocatable, dimension(:,:) :: three_body_scaling1
    integer, allocatable, dimension(:) :: three_body_counts1

    double precision, allocatable, dimension(:,:) :: three_body_distance2
    double precision, allocatable, dimension(:,:) :: three_body_scaling2
    integer, allocatable, dimension(:) :: three_body_counts2

    integer :: max_id
    integer :: max_two_body
    integer :: max_three_body

    double precision, parameter :: pi = 4.d0 * atan(1.d0)
    double precision, parameter :: d_width = 0.2d0
    double precision, parameter :: t_width = 0.2d0

    double precision :: k2, k3

    max_id = max(maxval(nuclear_charges1(:)), maxval(nuclear_charges2(:)))
    max_two_body = max(size(nuclear_charges1), size(nuclear_charges2))**2
    max_three_body = max(size(nuclear_charges1), size(nuclear_charges2))**3

    allocate(two_body_distance1(max_id**2, max_two_body))
    allocate(two_body_distance2(max_id**2, max_two_body))
    allocate(two_body_scaling1(max_id**2, max_two_body))
    allocate(two_body_scaling2(max_id**2, max_two_body))
    allocate(two_body_counts1(max_id**2))
    allocate(two_body_counts2(max_id**2))
    
    allocate(three_body_distance1(max_id**3, max_three_body))
    allocate(three_body_distance2(max_id**3, max_three_body))
    allocate(three_body_scaling1(max_id**3,  max_three_body))
    allocate(three_body_scaling2(max_id**3,  max_three_body))
    allocate(three_body_counts1(max_id**3))
    allocate(three_body_counts2(max_id**3))
    
    dmat1 = distance_matrix(coordinates1)
    dmat2 = distance_matrix(coordinates2)

    acube1 = ang_cube(coordinates1)
    acube2 = ang_cube(coordinates2)

    call initialize_two_body(coordinates1, nuclear_charges1, dmat1, &
        & two_body_distance1, two_body_scaling1, two_body_counts1, max_id)

    call initialize_two_body(coordinates2, nuclear_charges2, dmat2, &
        & two_body_distance2, two_body_scaling2, two_body_counts2, max_id)

    call initialize_three_body(coordinates1, nuclear_charges1, dmat1, acube1, &
        & three_body_distance1, three_body_scaling1, three_body_counts1, max_id)

    call initialize_three_body(coordinates2, nuclear_charges2, dmat2, acube2, &
        & three_body_distance2, three_body_scaling2, three_body_counts2, max_id)

    k2 =  1.0d0 * two_body_kernel(two_body_distance1, two_body_scaling1, two_body_counts1, &
                       & two_body_distance2, two_body_scaling2, two_body_counts2, d_width)

    k3 = 0.1d0 * two_body_kernel(three_body_distance1, three_body_scaling1, three_body_counts1, &
          & three_body_distance2, three_body_scaling2, three_body_counts2, t_width)

    ! write(*,*) max_id, max_two_body, k2, k3

    k = k2 + k3 

    ! k = exp( - (k)**2 / (2.0d0 * 100000.0d0))



    ! acube1 = ang_cube(coordinates1) 
    ! acube2 = ang_cube(coordinates2) 

    ! write(*,*) dmat1


end subroutine boss_kernel

