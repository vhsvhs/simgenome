from mpi4py import MPI
comm = MPI.COMM_WORLD

def mpi_check():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    tally = 0
    data = {}
    
    if rank == 0 and comm.Get_size() > 1:
        data = {'a': 7, 'b': 3.14}
        for i in range(1, comm.Get_size()):
            comm.send(data, dest=i, tag=11)
        #print "Master sent it."
    elif rank > 0:
        data = comm.recv(source=0, tag=11)
        #print "Rank ", rank, "got it.", data

    nodename = MPI.Get_processor_name()

    print "\n. MPI process", rank, "of", comm.Get_size(), "is alive on", nodename

    #print "My rank is", rank
    if rank == 0:
        for i in range(1, comm.Get_size()):
            #print "Master is listening for", i
            recv_data = {}
            recv_data = comm.recv(source=i, tag=12)
            #print "Master got it from", i, recv_data
            if recv_data.__contains__('a'):
                if recv_data.__contains__('b'):
                    if recv_data['a'] == 7 and recv_data['b'] == 3.14:
                        tally += 1
    else:
        comm.send(data, dest=0, tag=12)
        #print "Rank ", rank, "sent it.", data
    
    comm.Barrier()
    
    if rank == 0:    
        if tally == comm.Get_size()-1:
            return
        else:
            print "\n. Hmmm, something is wrong with MPI.  Goodbye."
            exit(1)


def is_my_item_core(i, l, r, s):
    """i = the item ID in a list (l) of all item IDs. r is my MPI rank."""
    chunk_size = int( l.__len__() / (s-1 ) )
    start = chunk_size * (r-1)
    end = start + chunk_size
    remainder = l.__len__() - (s-1) * chunk_size
    if i < end and i >= start:
        return True
    elif (r == s-2) and (i >= l.__len__()-remainder):
        return True
    else:
        return False

def is_my_item(i, l, r):
    """i = the item ID, l = the list of all item IDs, and r = my MPI rank.
        This method returns True if i should be processed by r."""
    if comm.Get_size() == 1:
        return True
    else:
        return is_my_item_core(i, l, r, comm.Get_size())
    
def list_my_items(l, r):
    """l is a list of all item IDs, r = my MPI rank.
        This method returns a list of all the items within l that belong to r."""
    my_items = []
    for i in range(0, l.__len__()):
        if is_my_item(i, l, r):
            my_items.append(i)
    return my_items

