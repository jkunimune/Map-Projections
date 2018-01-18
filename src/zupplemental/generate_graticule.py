GRATICULE = 10
ADJUST_POLES = False #set to True to reduce the number of meridians as you approach the pole


def plot_meridian(lamd, cut=0):
	tag = '\t\t\t<path d=M{},{}"'.format(lamd, cut) + 'v1'*(180-2*cut) + '"/>'
	print(tag)


def plot_parallel(phid, cut=0):
	tag = '\t\t\t<path d=M{},{}"'.format(cut, phid) + 'h1'*(360-2*cut) + '"/>'
	print(tag)


if __name__ == '__main__':
	for x in range(0, 361, 10):
		plot_meridian(x)
	for y in range(-80, 90, 10):
		plot_parallel(y)