# Documentation for `unify_eigvec.f90`

## Overview

`unify_eigvec.f90` is a Fortran module designed to ensure the uniqueness and orthogonality of eigenvectors, especially in cases where eigenvalues are degenerate (i.e., multiple eigenvectors correspond to the same eigenvalue). This is crucial in quantum mechanical calculations and other areas of physics and chemistry where degenerate states are common. The module processes a given set of eigenvalues and their corresponding eigenvectors, modifying the eigenvectors within each degenerate subspace to produce a unique, orthonormal set.

## Key Components

### `module m_unique_eigvec`

The main module containing all subroutines and functions for the eigenvector unification process.

#### `subroutine unique_eigvec(eigval, z, ierr)`
-   **Description:** This is the primary public subroutine. It takes an array of eigenvalues and a matrix of corresponding eigenvectors. It identifies groups of degenerate eigenvalues and then calls `unify_group_operator` for each group to ensure the eigenvectors in that degenerate subspace are unique and orthonormal. The input eigenvector matrix `z` is modified in place.
-   **Parameters:**
    -   `eigval (real(kind=8), intent(in))`: A 1D array containing the eigenvalues of the matrix.
    -   `z (real(kind=8), intent(inout))`: A 2D array (matrix) where each column represents an eigenvector corresponding to the eigenvalue at the same index in `eigval`. This matrix is modified in place.
    -   `ierr (integer, intent(out))`: An error flag. It returns 0 if the process is successful, and a non-zero value if an error occurs (e.g., failure to orthogonalize a degenerate subspace).

#### `subroutine unify_group_operator(beg_group, end_group, z, ierr)`
-   **Description:** This subroutine processes a specific group of eigenvectors that correspond to degenerate eigenvalues. It aims to make these eigenvectors orthonormal and unique. It first extracts the relevant eigenvectors, constructs an "RQ" matrix using `make_rq_mat`, performs a QR decomposition on this RQ matrix, and then uses the resulting Q matrix to transform the original degenerate eigenvectors into an orthonormal set.
-   **Parameters:**
    -   `beg_group (integer, intent(in))`: The starting column index (1-based) in the `z` matrix for the current group of degenerate eigenvectors.
    -   `end_group (integer, intent(in))`: The ending column index (1-based) in the `z` matrix for the current group of degenerate eigenvectors.
    -   `z (real(kind=8), intent(inout))`: The full matrix of eigenvectors. Only the columns from `beg_group` to `end_group` are modified by this subroutine.
    -   `ierr (integer, intent(out))`: An error flag. Returns 0 on success, non-zero on failure (e.g., if `make_rq_mat` fails or the subsequent QR decomposition fails).

#### `function make_groups(eigval) result(groups)`
-   **Description:** This function identifies groups of degenerate or near-degenerate eigenvalues. It iterates through the sorted eigenvalues and groups them if the absolute difference between consecutive eigenvalues is within a defined tolerance (e.g., `1e-8`).
-   **Parameters:**
    -   `eigval (real(kind=8), intent(in))`: A 1D array of eigenvalues, assumed to be sorted.
-   **Returns:**
    -   `groups (integer, allocatable)`: A 1D allocatable array. Each element `groups(i)` stores the size of the i-th degenerate group found. The sum of elements in `groups` should equal the total number of eigenvalues.

#### `subroutine make_rq_mat(eigvecs, rq_mat, ierr)`
-   **Description:** This subroutine constructs an "RQ" matrix from a given set of (degenerate) eigenvectors. The goal is to create a smaller square matrix (`rq_mat`) whose rows are selected components of the input eigenvectors. The selection strategy aims to find a set of rows that makes `rq_mat` well-conditioned and its rows linearly independent. This `rq_mat` is then used in a QR-like decomposition process (specifically, `dgeqrfp` followed by `dorgqr` on `rq_mat` transpose) to help orthogonalize the original `eigvecs`. The selection of rows involves picking rows with components that are above a certain `cutoff` value and iteratively checking for linear independency using `linear_independency_real`.
-   **Parameters:**
    -   `eigvecs (real(kind=8), intent(in))`: A 2D array where columns are the eigenvectors of a degenerate subspace.
    -   `rq_mat (real(kind=8), intent(inout), allocatable)`: The constructed RQ matrix. It will be allocated as an N-by-N matrix, where N is the number of eigenvectors (and thus the dimension of the degenerate subspace).
    -   `ierr (integer, intent(out))`: Error flag. Returns 0 on success, 1 if a linearly independent `rq_mat` cannot be formed.

#### `function linear_independency_real(mat) result(linear_independency_real)`
-   **Description:** This function assesses the linear independency of a real square matrix by computing its Singular Value Decomposition (SVD) using LAPACK's `dgesdd` routine. The smallest singular value is returned; a very small value (close to zero) indicates that the matrix is singular or nearly singular, meaning its rows/columns are linearly dependent.
-   **Parameters:**
    -   `mat (real(kind=8), intent(inout))`: The real square matrix to analyze. Note: The `dgesdd` routine can modify the input matrix `mat`.
-   **Returns:**
    -   `linear_independency_real (real(kind=8))`: The smallest singular value of the matrix `mat`.

#### `function linear_independency_cmplx(mat) result(linear_independency_cmplx)`
-   **Description:** This function assesses the linear independency of a complex square matrix by computing its Singular Value Decomposition (SVD) using LAPACK's `zgesdd` routine. Similar to the real version, it returns the smallest singular value.
-   **Parameters:**
    -   `mat (complex(kind=8), intent(inout))`: The complex square matrix to analyze. Note: The `zgesdd` routine can modify the input matrix `mat`.
-   **Returns:**
    -   `linear_independency_cmplx (real(kind=8))`: The smallest singular value of the matrix `mat` (singular values are always real).

#### `interface linear_independency`
-   **Description:** This generic interface allows the function `linear_independency` to be called with either a real matrix (dispatching to `linear_independency_real`) or a complex matrix (dispatching to `linear_independency_cmplx`). This provides flexibility, though `unify_eigvec.f90` primarily deals with real eigenvectors.

#### `subroutine print_mtx(mtx)`
-   **Description:** A utility subroutine used for debugging purposes. It prints the elements of a given real matrix to the standard output unit, formatted for readability.
-   **Parameters:**
    -   `mtx (real(kind=8), intent(in))`: The real matrix to be printed.

## Important Variables/Constants

-   **Tolerance in `make_groups` (e.g., `1e-8`):** This floating-point value is used as a threshold to determine if two eigenvalues are considered degenerate. If `abs(eigval(i) - eigval(i+1)) < tolerance`, they are grouped together.
-   **`cutoff` in `make_rq_mat` (e.g., `1e-6 * sqrt(1.0/size(eigvecs,1))`):** This threshold is used during the construction of the `rq_mat`. When selecting rows from the input `eigvecs` to form `rq_mat`, components (elements) of eigenvectors whose absolute value is below this `cutoff` are initially ignored or considered less significant. The scaling factor `sqrt(1.0/size(eigvecs,1))` (where `size(eigvecs,1)` is the dimension of the eigenvectors) helps to normalize this cutoff.
-   **Tolerance in `make_rq_mat` for `lindep < 1e-9`:** After attempting to construct `rq_mat` by selecting rows, its linear independency is checked. If the smallest singular value (`lindep`) of `rq_mat` is less than this tolerance (e.g., `1e-9`), the `rq_mat` is considered to be not sufficiently linearly independent, and the `make_rq_mat` subroutine may report an error or try a different strategy.

## Usage Examples

The `example.f90` file in the repository demonstrates how to use the `m_unique_eigvec` module.
A typical workflow involves:
1.  Diagonalizing a Hermitian matrix (e.g., using LAPACK's `dsyev` or a similar routine) to obtain eigenvalues and eigenvectors.
2.  Calling `call unique_eigvec(eigenvalues, eigenvectors, error_flag)` to process the eigenvectors, making them unique and orthonormal within degenerate subspaces.

```fortran
! Program to demonstrate the use of m_unique_eigvec
! Assume H is your Hermitian matrix, eig1 are eigenvalues, vecs1 are eigenvectors
!
! module m_example_diag
!   use m_unique_eigvec ! Import the module
!   implicit none
! contains
!   subroutine run_example()
!     real(kind=8), allocatable :: H(:,:), eig1(:), vecs1(:,:)
!     integer :: n, i, j, ierr
!
!     ! Example: Create a sample Hermitian matrix (e.g., one with degenerate eigenvalues)
!     n = 5 ! Size of the matrix
!     allocate(H(n,n), eig1(n), vecs1(n,n))
!
!     ! Populate H with some values - for a real example, this would be your physical system's matrix
!     ! For instance, a matrix that is known to produce degenerate eigenvalues:
!     H = 0.0
!     H(1,1) = 2.0; H(2,2) = 2.0; H(3,3) = 3.0; H(4,4) = 3.0; H(5,5) = 5.0
!     H(1,2) = 1.0e-9; H(2,1) = 1.0e-9 ! Small perturbation to break perfect numerical degeneracy if needed by solver
!     H(3,4) = 1.0e-9; H(4,3) = 1.0e-9
!     ! ... ensure H is Hermitian ...
!
!     ! Step 1: Diagonalize the matrix H to get eigenvalues (eig1) and eigenvectors (vecs1)
!     ! This typically involves calling a LAPACK routine like DSYEV.
!     ! For simplicity, assume eig1 and vecs1 are populated here.
!     ! e.g., call DSYEV('V', 'U', n, H, n, eig1, vecs1_lapack_output_wk, lwork, info)
!     ! vecs1 = H ! After DSYEV, eigenvectors are in H if jobz='V'
!
!     ! --- Placeholder for actual diagonalization ---
!     ! For this example, let's manually set some eigenvalues and mock eigenvectors
!     ! True eigenvalues: 2, 2, 3, 3, 5
!     eig1 = [2.0, 2.0, 3.0, 3.0, 5.0]
!     ! Mock eigenvectors (columns). For degenerate pairs, these might not be orthogonal.
!     vecs1 = reshape([ &
!         1.0, 0.0, 0.0, 0.0, 0.0, & ! vec for eig=2
!         0.9, 0.1, 0.0, 0.0, 0.0, & ! vec for eig=2 (potentially not orthogonal to first)
!         0.0, 0.0, 1.0, 0.0, 0.0, & ! vec for eig=3
!         0.0, 0.0, 0.8, 0.2, 0.0, & ! vec for eig=3 (potentially not orthogonal to third)
!         0.0, 0.0, 0.0, 0.0, 1.0  & ! vec for eig=5
!     ], [n, n])
!     ! --- End Placeholder ---
!
!     write(*,*) "Original Eigenvalues:"
!     write(*,'(5F10.6)') eig1
!     write(*,*) "Original Eigenvectors (columns):"
!     call print_mtx(vecs1) ! Assuming print_mtx is accessible or implemented
!
!     ! Step 2: Call unique_eigvec to process the eigenvectors
!     call unique_eigvec(eig1, vecs1, ierr)
!
!     if (ierr /= 0) then
!         write (*,*) "Error in unique_eigvec. Error code:", ierr
!     else
!         write (*,*) "Eigenvectors processed by unique_eigvec."
!         write (*,*) "Unified Eigenvectors (columns):"
!         call print_mtx(vecs1)
!         ! Here you would proceed with using the unified eigenvectors in vecs1
!         ! For example, check their orthogonality:
!         ! do i = 1, n
!         !   do j = i, n
!         !     write(*,*) "Dot product V(:,",i,") . V(:,",j,") = ", dot_product(vecs1(:,i), vecs1(:,j))
!         !   end do
!         ! end do
!     endif
!
!     deallocate(H, eig1, vecs1)
!   end subroutine run_example
! end module m_example_diag
!
! program main
!   use m_example_diag
!   call run_example()
! end program main
```

## Dependencies and Interactions

-   **OpenMP (`omp_lib`):** The module includes `use omp_lib`, suggesting that it may leverage OpenMP for parallelizing loops, particularly in `unify_group_operator` or other computationally intensive parts. If OpenMP directives (`!$OMP ...`) are used, compilation should include an OpenMP flag (e.g., `-fopenmp` for gfortran, `icc -qopenmp`).
-   **LAPACK (Linear Algebra PACKage):** The module relies heavily on LAPACK routines for performing numerical linear algebra operations. Key routines that are likely used (based on the described functionality) include:
    -   `dgeqrfp` (or `dgeqrf`): Computes a QR factorization of a general M-by-N real matrix. Used in `unify_group_operator` to orthogonalize the `rq_mat`.
    -   `dorgqr` (or `dormqr`): Generates an M-by-N real matrix Q with orthonormal columns from the reflectors returned by `dgeqrfp`, or applies Q. Used to construct the orthogonal transformation matrix from the QR factorization.
    -   `dgemm`: Performs matrix-matrix multiplication. Used in `unify_group_operator` to apply the orthogonal transformation (derived from Q) to the original degenerate eigenvectors.
    -   `dgesdd`: Computes the singular value decomposition (SVD) of a real general matrix A = U * SIGMA * V**T. Used in `linear_independency_real` to find the smallest singular value.
    -   `zgesdd`: Computes the SVD of a complex general matrix. Used in `linear_independency_cmplx`.
-   **Input Data:** The module expects eigenvalues and their corresponding eigenvectors as input. These are typically the result of a preceding matrix diagonalization step (e.g., from diagonalizing a Hamiltonian in physics). The module then modifies the input eigenvectors in place to ensure uniqueness and orthogonality within any identified degenerate subspaces.
```
