import pstats
import sys

p = pstats.Stats(sys.argv[1])
#help(p.sort_stats)
p.sort_stats('time', 'cum').print_stats(40)
#p.sort_stats('cumulative').print_stats(40)
