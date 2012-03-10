from mpi4py import MPI
comm = MPI.COMM_WORLD

def is_my_item(i, n, r):
    """i = the item ID in n total items. r is my MPI rank."""
    chunk_size = int( n / comm.Get_size() )
    start = chunk_size * r
    end = start + chunk_size
    remainder = n - comm.Get_size() * chunk_size
    if i < end and i >= start:
        return True
    elif (r == comm.Get_size()-1) and (i >= n-remainder):
        return True
    else:
        return False