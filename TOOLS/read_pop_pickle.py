from test_include import *
ap = ArgParser(sys.argv)
p0 = pickle.load( open(sys.argv[1], "r"))
pop0 = Population()
pop0.uncollapse(p0)

