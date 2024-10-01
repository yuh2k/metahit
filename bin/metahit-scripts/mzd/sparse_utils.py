import logging
import numpy as np
import scipy.sparse as scisp
import sparse
from math import ceil

logger = logging.getLogger(__name__)


def is_hermitian(m, tol=1e-6):
    """
    Test that a sparse matrix is hermitian (also suffices for symmetric)

    :param m: square matrix
    :param tol: tolerance above zero for m - m.T < tol
    :return: True matrix is Hermitian
    """
    return np.all(~(np.abs(m - m.H) >= tol).todense())


def tensor_print(T):
    """
    Pretty print a dense (numpy) 4D matrix. Users should consider the size of the matrix before
    printing, as they can be large! More useful for smaller objects

    :param T: the tensor to print with dim: (N,M,n,m)
    """
    try:
        pw = int(np.ceil(np.log10(T.max())))
    except OverflowError:
        pw = 1
    for i in range(T.shape[0]):
        for k in range(T.shape[2]):
            print('|', end=' ')
            for j in range(T.shape[1]):
                print('[', end=' ')
                for l in range(T.shape[3]):
                    print(f'{T[i, j, k, l]:{pw}d}', end=' ')
                print(']', end=' ')
            print('|')
        if i < T.shape[1] - 1:
            print('+')
    print()


def downsample(m, block_size, method='mean'):
    """
    Perform a down-sampling of a 2D matrix (scipy.sparse or ndarray)
    by a factor of block_size in each dimension. Block size must be
    an integer larger than 1.

    :param m: a matrix (scipy.sparse or ndarray)
    :param block_size: an integer reduction factor
    :param method: per block mean or maximum.
    :return: scipy.sparse.csr_matrix
    """
    assert isinstance(m, (np.ndarray, scisp.spmatrix)), 'supplied array must be of type np.ndarray or scipy.spmatrix'
    assert block_size > 1 and isinstance(block_size, int), 'block_size must be an integer larger than 1'

    pad_size = lambda N, n: int(ceil(N / float(n)) * n) - N

    pad_row = pad_size(m.shape[0], block_size)
    pad_col = pad_size(m.shape[1], block_size)

    # only dok has resize method up until scipy 1.1.0
    if isinstance(m, np.ndarray):
        m = scisp.dok_matrix(m)
    else:
        m = m.todok()
    m.resize((m.shape[0] + pad_row, m.shape[1] + pad_col))

    # conversion to csr here appears necessary to properly preserve matrix
    m = sparse.COO(m.tocsr())

    if method == 'mean':
        m = m.reshape((m.shape[0] // block_size, block_size, m.shape[1] // block_size, block_size)).sum(axis=(1, 3)).tocsr()
        m *= 1.0 / block_size**2
    elif method == 'max':
        m = m.reshape((m.shape[0] // block_size, block_size, m.shape[1] // block_size, block_size)).max(axis=(1, 3)).tocsr()

    return m


def kr_biostochastic(m, tol=1e-6, x0=None, delta=0.1, Delta=3, max_iter=1000):
    """
    Normalize a matrix to be bistochastic using Knight-Ruiz algorithm.

    :param m: the input matrix (fully symmetric)
    :param tol: precision tolerance
    :param x0: an initial guess
    :param delta: how close balancing vector can get
    :param Delta: how far balancing vector can get
    :param max_iter: maximum number of iterations
    :return: tuple containing the bistochastic matrix and the scale factors
    """
    assert scisp.isspmatrix(m), 'input matrix must be sparse matrix from scipy.spmatrix'
    assert m.shape[0] == m.shape[1], 'input matrix must be square'

    _orig = m.copy()

    # replace 0 diagonals with 1
    m = m.tolil()
    is_zero = m.diagonal() == 0
    if np.any(is_zero):
        logger.warning(f'treating {is_zero.sum()} zeros on diagonal as ones')
        ix = np.where(is_zero)
        m[ix, ix] = 1

    if not scisp.isspmatrix_csr(m):
        m = m.tocsr()

    if not is_hermitian(m, tol):
        logger.warning('input matrix is expected to be fully symmetric')

    n = m.shape[0]
    e = np.ones(n)

    if x0 is None:
        x0 = e.copy()

    g = 0.9
    etamax = 0.1
    eta = etamax
    stop_tol = tol * 0.5

    x = x0.copy()
    rt = tol ** 2
    v = x * m.dot(x)

    rk = 1 - v
    rho_km1 = rk.T.dot(rk)
    rout = rho_km1
    rold = rout

    n_iter = 0
    i = 0
    y = np.empty_like(e)
    while rout > rt and n_iter < max_iter:
        i += 1
        k = 0
        y[:] = e

        inner_tol = max(rout * eta ** 2, rt)

        while rho_km1 > inner_tol:
            k += 1
            if k == 1:
                Z = rk / v
                p = Z
                rho_km1 = rk.T.dot(Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            w = x * m.dot(x * p) + v * p
            alpha = rho_km1 / p.T.dot(w)
            ap = alpha * p

            ynew = y + ap

            if np.amin(ynew) <= delta:
                if delta == 0:
                    break
                ind = np.where(ap < 0)[0]
                gamma = np.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            if np.amax(ynew) >= Delta:
                ind = np.where(ynew > Delta)[0]
                gamma = np.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            y = ynew
            rk = rk - alpha * w
            rho_km2 = rho_km1

            Z = rk * v
            rho_km1 = np.dot(rk.T, Z)

            if np.any(np.isnan(x)):
                raise RuntimeError('scale vector has developed invalid values (NANs)!')

        x *= y
        v = x * m.dot(x)

        rk = 1 - v
        rho_km1 = np.dot(rk.T, rk)
        rout = rho_km1
        n_iter += k + 1

        rat = rout / rold
        rold = rout
        res_norm = np.sqrt(rout)
        eta_o = eta
        eta = g * rat

        if g * eta_o ** 2 > 0.1:
            eta = max(eta, g * eta_o ** 2)
        eta = max(min(eta, etamax), stop_tol / res_norm)

    if n_iter > max_iter:
        raise RuntimeError(f'matrix balancing failed to converge in {n_iter} iterations')

    del m

    logger.debug(f'It took {n_iter} iterations to achieve bistochasticity')

    if n_iter >= max_iter:
        logger.warning(f'Warning: maximum number of iterations ({max_iter}) reached without convergence')

    X = scisp.spdiags(x, 0, n, n, 'csr')
    return X.T.dot(_orig.dot(X)), x



class Sparse2DAccumulator(object):
    def __init__(self, N):
        self.shape = (N, N)
        self.mat = {}
        self.dtype = np.uint32

    def setitem(self, index, value):
        assert len(index) == 2 and index[0] >= 0 and index[1] >= 0, 'invalid index: {}'.format(index)
        assert isinstance(value, (int, int)), 'values must be integers'
        self.mat[index] = value

    def getitem(self, index):
        return self.mat.get(index, 0)

    def get_coo(self, symm=True):
        """
        Create a COO format sparse representation of the accumulated values.

        :param symm: ensure matrix is symmetric on return
        :return: a scipy.coo_matrix sparse matrix
        """
        coords = [[], []]
        data = []
        for (i, j), v in self.mat.items():
            coords[0].append(i)
            coords[1].append(j)
            data.append(v)

        m = scisp.coo_matrix((data, coords), shape=self.shape, dtype=self.dtype)

        if symm:
            m += scisp.tril(m.T, k=-1)

        return m.tocoo()

def max_offdiag(_m):
    """
    Determine the maximum off-diagonal values of a given symmetric matrix. As this
    is assumed to be symmetric, we consider only the rows.

    :param _m: a scipy.sparse matrix
    :return: the off-diagonal maximum values
    """
    assert scisp.isspmatrix(_m), 'Input matrix is not a scipy.sparse object'
    _m = _m.tolil()
    _m.setdiag(0)
    return np.asarray(_m.tocsr().max(axis=0).todense()).ravel()


def compress(_m, _mask):
    """
    Remove rows and columns using a 1d boolean mask.

    :param _mask: True (keep), False (drop)
    :return: a coo_matrix of only the accepted rows/columns
    """
    assert scisp.isspmatrix(_m), 'Input matrix is not a scipy sparse matrix type'

    if not scisp.isspmatrix_coo(_m):
        _m = _m.tocoo()

    keep_row = []
    keep_col = []
    keep_data = []
    accept_index = set(np.where(_mask)[0])
    for i in range(_m.nnz):
        if _m.row[i] in accept_index and _m.col[i] in accept_index:
            keep_row.append(_m.row[i])
            keep_col.append(_m.col[i])
            keep_data.append(_m.data[i])

    shift = np.cumsum(~_mask)
    for i in range(len(keep_row)):
        keep_row[i] -= shift[keep_row[i]]
        keep_col[i] -= shift[keep_col[i]]

    return scisp.coo_matrix((keep_data, (keep_row, keep_col)), shape=_m.shape - shift[-1])


class Sparse4DAccumulator(object):
    """
    Simple square sparse tensor of dimension (N, N, 2, 2).
    """

    def __init__(self, N):
        self.shape = (N, N, 2, 2)
        self.mat = {}
        self.dtype = np.uint32

    def __setitem__(self, index, value):
        assert isinstance(index, tuple), 'index must be a list of indices'
        if len(index) == 4:
            assert 0 <= index[0] < self.shape[0] and \
                   0 <= index[1] < self.shape[1] and \
                   0 <= index[2] < 2 and \
                   0 <= index[3] < 2, f'invalid range {index} for dimension {self.shape}'
            if index[:2] not in self.mat and np.any(value != 0):
                self.mat.setdefault(index[:2], self._make_elem())[index[2:]] = value

        if len(index) == 2:
            assert 0 <= index[0] < self.shape[0] and \
                   0 <= index[1] < self.shape[1], f'invalid range {index} for dimension {self.shape}'
            if index not in self.mat:
                self.mat.setdefault(index, self._make_elem())[:] = value

    def __getitem__(self, index):
        return self.mat.setdefault(index, self._make_elem())

    def _make_elem(self):
        return np.zeros((2, 2), dtype=self.dtype)

    def get_coo(self, symm=True):
        """
        Create a COO format sparse representation of the accumulated values. NOTE: As scipy
        does not support multidimensional arrays, this object is from the "sparse" module.
        """
        _coords = [[], [], [], []]
        _data = []
        _m = self.mat
        _inner_indices = [[0, 0], [0, 1], [1, 0], [1, 1]]
        for i, j in _m.keys():
            for k, l in _inner_indices:
                v = _m[i, j][k, l]
                if v != 0:
                    _coords[0].append(i)
                    _coords[1].append(j)
                    _coords[2].append(k)
                    _coords[3].append(l)
                    _data.append(v)

        _m = sparse.COO(_coords, _data, self.shape, has_duplicates=False)

        if symm:
            _m = Sparse4DAccumulator.symm(_m)

        return _m

    @staticmethod
    def _flip(c_row):
        """
        Flip indices (coordinates) as pairs: (i,j), (k,l) -> (j,i), (l,k)
        """
        c_row = c_row.copy()
        c_row[0], c_row[1] = c_row[1], c_row[0]
        c_row[2], c_row[3] = c_row[3], c_row[2]
        return c_row

    @staticmethod
    def symm(_m):
        """
        Make a 4D COO matrix symmetric, all elements above and below the diagonal are included.
        Duplicate entries will be summed.
        """
        ix = np.where(~np.apply_along_axis(lambda x: x[0] == x[1], 0, _m.coords))[0]
        _coords = np.hstack((_m.coords, np.apply_along_axis(Sparse4DAccumulator._flip, 0, _m.coords[:, ix])))
        _data = np.hstack((_m.data, _m.data[ix]))
        return sparse.COO(_coords, _data, shape=_m.shape, has_duplicates=True)


def max_offdiag_4d(_m):
    """
    Determine the maximum off-diagonal summed signal, where "summed signal" refers to reducing the
    tensor to a 2D matrix by summing over the last two axes (2x2 submatrices).

    :param _m: a 4d sparse.COO or DOK matrix with dimension NxNx2x2.
    :return: a vector of length N containing off-diagonal maximums.
    """
    return max_offdiag(_m.sum(axis=(2, 3)).tocsr())


def flatten_tensor_4d(_m):
    """
    Flatten a 4D tensor into 2D by doubling the first two dimensions.
    """
    _coords = [[], []]
    _data = []
    for n in range(_m.nnz):
        i, j, k, l = _m.coords[:, n]
        ii = 2 * i
        jj = 2 * j
        _coords[0].append(ii + k)
        _coords[1].append(jj + l)
        _data.append(_m.data[n])

    _m = scisp.coo_matrix((_data, _coords), shape=(2 * _m.shape[0], 2 * _m.shape[1]))
    return _m


def compress_4d(_m, _mask):
    """
    Remove rows and columns of a sparse 4D matrix using a 1d boolean mask. Masking operates on
    only the first two primary axes (essentially a 2D matrix with 2x2 cells). The returned matrix is of type sparse.COO.
    """
    assert isinstance(_m, (sparse.COO, sparse.DOK)), 'Input matrix must be of sparse.COO or sparse.DOK type'
    if not isinstance(_m, sparse.COO):
        _m = _m.to_coo()

    keep_coords = []
    keep_data = []
    accept_index = set(np.where(_mask)[0])
    for i in range(_m.nnz):
        if _m.coords[0, i] in accept_index and _m.coords[1, i] in accept_index:
            keep_coords.append(_m.coords[:, i])
            keep_data.append(_m.data[i])
    keep_coords = np.array(keep_coords, dtype=int).T

    shift = np.cumsum(~_mask, dtype=int)
    keep_coords[:2, :] -= shift[keep_coords[:2, :]]

    new_shape = list(_m.shape)
    new_shape[:2] -= shift[-1]
    return sparse.COO(keep_coords, keep_data, shape=new_shape, has_duplicates=False)


def dotdot(_m, _a):
    """
    Assuming A is a vector representing the trace of a diagonal matrix, dotdot
    performs the transformation dot(A.T,dot(M,A)) on  a sparse matrix.
    """
    for n in range(_m.nnz):
        i, j = _m.coords[:2, n]
        _m.data[n] *= _a[i] * _a[j]
    return _m


def kr_biostochastic_4d(m4d, **kwargs):
    """
    Knight-Ruiz applied to a NxNx2x2 tensor.
    """
    m2d = m4d.astype(float).sum(axis=(2, 3)).tocsr()
    _, scl = kr_biostochastic(m2d, **kwargs)
    return dotdot(m4d.astype(float), scl), scl
