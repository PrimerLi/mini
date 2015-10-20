subroutine inverse(matrix, matrixInverse)
    implicit none
    complex(kind = 8), intent(in) :: matrix(:, :)
    complex(kind = 8), intent(out) :: matrixInverse(:, :)
    integer :: matrixDimension
    integer :: M, N
    integer :: LDA
    integer, dimension(:), allocatable :: ipiv
    integer :: info
    complex(kind = 8), dimension(:), allocatable :: work
    integer :: lwork
    complex(kind = 8), dimension(:, :), allocatable :: A

    matrixDimension = size(matrix, 1)
    if (size(matrix, 2) .ne. matrixDimension) then
    	write(*, *) "Only square matrix will be accepted. "
	stop
    endif
    if (size(matrixInverse, 1) .ne. size(matrixInverse, 2) .or. &
	    &size(matrixInverse, 1) .ne. matrixDimension) then
    	write(*, *) "Matrix dimension incompatible. "
	stop
    endif

    Allocate(A(matrixDimension, matrixDimension))
    Allocate(ipiv(matrixDimension))
    M = matrixDimension
    N = matrixDimension
    LDA = matrixDimension
    lwork = 2*N
    allocate(work(lwork))

    A = matrix

    call zgetrf(M, N, A, LDA, ipiv, info)

    if (info == 0) then
    	call zgetri(N, A, LDA, ipiv, work, lwork, info)
    else
    	write(*, *) "The matrix is singular. "
	stop
    endif

    matrixInverse = A

    deallocate(ipiv)
    deallocate(A)
    Deallocate(work)
end subroutine inverse

