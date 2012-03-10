from mpi4py import MPI
comm = MPI.COMM_WORLD

def is_my_item(i, l, r):
    """i = the item ID in a list (l) of all item IDs. r is my MPI rank."""
    chunk_size = int( l.__len__() / comm.Get_size() )
    start = chunk_size * r
    end = start + chunk_size
    remainder = l.__len__() - comm.Get_size() * chunk_size
    if i < end and i >= start:
        return True
    elif (r == comm.Get_size()-1) and (i >= l.__len__()-remainder):
        return True
    else:
        return False
    
def list_my_items(l, r):
    my_items = []
    for i in range(0, l.__len__()):
        if is_my_item(i, l, r):
            my_items.append(i)
    return my_items