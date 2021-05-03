from numba import njit


def number_of_nodes(func, nmin, nmax):
    assert nmin >= 1, 'Error: number_of_nodes: nmin < 1'
    assert nmin+1 <= len(func), 'Error: number_of_nodes: nmin+1>len(func)'
    assert nmax-1 >= 1, 'Error: number_of_modes: nmax-1 < 1'
    assert nmax <= len(func), 'Error: number_of_nodes: nmax>size(func)'

    return number_of_nodes_calc(func, nmin, nmax)

@njit()
def number_of_nodes_calc(func, nmin, nmax):
    num_nodes = 0

    for i in range(nmin+1, nmax):
        if func[i-1] * func[i] < 0.0:
            num_nodes += 1

    return num_nodes