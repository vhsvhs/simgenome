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

    print "\n. MPI process", rank, "of", comm.Get_size()-1, "is alive."

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
        #print tally
        if tally == comm.Get_size()-1:
            return
        else:
            print "\n. Hmmm, something is wrong with MPI.  Goodbye."
            exit(1)



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