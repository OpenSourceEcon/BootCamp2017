module spl
! Spline module written by Viktor Tsyrennikov
! It is a version of Carl de Boor's code

    implicit none

contains


subroutine csint(x, y, breaks, coefs)
real(8) :: x(:), y(:), breaks(size(x)), coefs(4,size(x))
    breaks = x
    coefs(1,:) = y
    call cspl(breaks, coefs, size(x), 0, 0)
end subroutine


subroutine cspl ( tau, c, n, bcbeg, bcend )
! This code is based on C. de Boor's cubcpl.f
! Inputs:
!   n = number of data points. assumed to be > 3
!   (tau(i),c(1,i), i=1,...,n) = approximation data (tau str. increasing)
!   bcbeg, bcend = boundary condition indicators
!   c(2,1), c(2,n) = boundary condition information
!   bcbeg = 0 = not-a-knot condition x(2)
!           1 = first derivative at x(1) equals c(2,1)
!           2 = second derivative at x(1) equals c(2,1)
!   bcend = 0 = not-a-knot condition at x(n-1)
!           1 = first derivative at x(1) equals c(2,1)
!           2 = second derivative at x(1) equals c(2,1)
! Output:
!   c(j,i) = polynomial coefficient for x^{i-1} on interval j 
!           in the interval (tau(i), tau(i+1))
!           f(x) = c(1,i) + h*(c(2,i) + h*(c(3,i) + h*c(4,i)/3)/2), h = x-tau(i)

integer :: bcbeg, bcend, n, i, j, m
real(8) :: c(4,n), tau(n), divdf1, divdf3, dtau, g, h

    ! Compute first differences of tau sequence and store in c(3,:)
    ! Compute first divided difference of data and store in c(4,:)
    do m = 2, n
        c(3,m) = tau(m) - tau(m-1)
        c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
    end do

    ! Incorporate left boundary not-a-knot condition
    if (bcbeg == 0) then
        c(4,1) = c(3,3)
        c(3,1) = c(3,2) + c(3,3)
        c(2,1) = ((c(3,2) + 2d0*c(3,1))*c(4,2)*c(3,3) + c(3,2)**2*c(4,3))/c(3,1)
    end if

    ! Incorporate left boundary first derivative condition
    if (bcbeg == 1) then
        c(4,1) = 1d0
        c(3,1) = 0d0
    end if

    ! Incorporate left boundary second derivative condition
    if (bcbeg == 2) then
        c(4,1) = 2d0
        c(3,1) = 1d0
        c(2,1) = 3d0*c(4,2) - c(3,2)/2d0*c(2,1)
    end if

    ! Generate the equations for interior knots and carry out the forward
    ! pass of gauss elimination, after which the m-th equation reads
    !   c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m)
    do m = 2, n-1
        g = -c(3,m+1)/c(4,m-1)
        c(2,m) = g*c(2,m-1) + 3d0*(c(3,m)*c(4,m+1) + c(3,m+1)*c(4,m))
        c(4,m) = g*c(3,m-1) + 2d0*(c(3,m) + c(3,m+1))
    end do

    ! Incorporate right boundary not-a-knot condition
    if (bcend == 0) then
        g = c(3,n-1) + c(3,n)
        c(2,n) = ((c(3,n) + 2d0*g)*c(4,n)*c(3,n-1) &
            + c(3,n)**2*(c(1,n-1) - c(1,n-2))/c(3,n-1))/g
        g = -g/c(4,n-1)
        c(4,n) = c(3,n-1)
    end if

    ! Incorporate right boundary second derivative condition
    if (bcend == 2) then
        c(2,n) = 3d0*c(4,n) + c(3,n)/2d0*c(2,n)
        c(4,n) = 2d0

        g = -1./c(4,n-1)
    end if

    if (bcend .ne. 1) then
        c(4,n) = g*c(3,n-1) + c(4,n)
        c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
    end if

    ! Carry out back substitution
    j = n-1
    do while (j > 0)
        c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
        j = j - 1
    end do

    ! Generate cubic coefficients in each interval
    do i = 2, n
        dtau = c(3,i)
        divdf1 = (c(1,i) - c(1,i-1))/dtau
        divdf3 = c(2,i-1) + c(2,i) - 2d0*divdf1
        c(3,i-1) = 2d0*(divdf1 - c(2,i-1) - divdf3)/dtau
        c(4,i-1) = (divdf3/dtau)*(6d0/dtau)
    end do

    ! Update coefficients c(:,n)
    h = tau(n)-tau(n-1)
    c(1,n) = c(1,n-1)+(c(2,n-1)+(c(3,n-1)+c(4,n-1)*h/3)*h/2)*h
    c(2,n) = c(2,n-1)+(c(3,n-1)+c(4,n-1)*h/3)*h
    c(3,n) = c(3,n-1)/2+c(4,n-1)*h/3
    c(4,n) = c(4,n-1)
end subroutine

function csval(breaks, coefs, x) result(y)
    real(8) :: breaks(:), coefs(:,:), x, y
    integer :: j, k, m
    real(8) :: h

    ! Find index i of largest breakpoint to the left of x
    j = 1
    k = size(breaks)
    do while (k-j > 1)
        m = floor(0.5d0*(k+j))
        if (x > breaks(m)) then
            j = m
        else
            k = m
        end if
    end do

    h = x - breaks(j)
    y = coefs(4,j)
    do k = 3, 1, -1
        y = coefs(k,j) + y*h/k
    end do
end function

end module
