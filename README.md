# Map-Projections
A class to create custom maps of the Earth's surface. There are thousands of combinations of color-schemes, projections, and aspects. Includes Mercator, Gall-Peters, Orthographic, Peirce Quincuncial, and More!

## Installation
If you are a fancy Windows user, I recommend the convenient [fancy Windows binaries](https://github.com/jkunimune15/Map-Projections/releases). Double-click to install them and then keep pressing buttons until something good happens. If you see a map, you're in the right place.

If you are not on Windows or are otherwise not fancy enough to deserve such executables, simply double-click on the .jar files in the main directory and, if you have [Java](https://java.com/en/download/) installed (10/10 would recommend), it should just run without any set-up.

You could also compile and run the source code, but if you do, there are a few dependencies. The Java dependencies can be obtained as .jar files:

* [Apache Commons Mathematics Library](http://commons.apache.org/proper/commons-math/download_math.cgi)
* [Java Tools for Experimental Mathematics "ellipticFunctions" package](http://www3.math.tu-berlin.de/jtem/downloads.html), which requires their "mfc" package

If you want to use the Python, you'll need a couple of packages from [PyPI](https://pypi.python.org/pypi). Just install Python 3 and pip, and then call the following from a command line (or use Anaconda or something, I don't know. Up to you).
~~~~
pip3 install numpy pyshp
~~~~
`generate_coastlines.py` also takes input data from [naturalearthdata.com](http://www.naturalearthdata.com/downloads/), which should be placed in `src\\zupplemental\\data\\`.

## Usage
There are three executable Java files, three runnable Java scripts, and four runnable Python scripts. These are, in order:

* `MapDesignerRaster.jar` &ndash; The original program. Create custom oblique raster images of the Earth's surface using a variety of algorithms called _projections_.  
* `MapDesignerVector.jar` &ndash; The same idea, but working in vector images instead in case you want to cut a vinyl sticker or something.  
* `MapAnalyzer.jar` &ndash; See graphs and figures quantifying the amount of scale and angular distortion present in each map projection.  
* `MapPlotter.java` &ndash; Plot a large group of map projections by the amount of distortion they produce.  
* `MapOptimizer.java` &ndash; Run gradient descent on parametric projections to minimize their distortion.  
* `MapExplainer.java` &ndash; Generate an HTML blurb outlining and displaying every map projection.
* `generate_coastlines.py` &ndash; Generate an SVG string outlining the continents, islands, and major lakes of the world, to be used as vector input.
* `generate_graticule.py` &ndash; Generate an SVG string displaying a map graticule, to be used as vector input.
* `generate_indicatrices.py` &ndash; Generate an SVG string outlining an array of Tissot's indiatrices of distortion, to be used as vector input.
* `generate_orthodromes.py` &ndash; Generate a mesh of orthodromes in an Equirectangular projection, to be used as vector input.

The executable applications all have similar layouts that let you select an input equirectangular map, a projection, an aspect (where the North Pole is situated with respect to the projection), and parameters if applicable. Go crazy! There are a practically unlimited number of combinations.

The runnable scripts just kind of work on their own. Those ones aren't really meant for mass consumption.

## Wherefore?
I'll write a little blurb here later.

For some examples, check out the `output` folder. For more information, go to [jkunimune15.github.io/Map-Projections](https://jkunimune15.github.io/Map-Projections).

## Credits
While I wrote all of the code in this repository myself, and I created several of the simpler images from scratch, other people did help. Here's a comprehensive list.
* **The NASA** for [Basic.png](https://visibleearth.nasa.gov/view.php?id=57730), [Satellite.jpg](https://visibleearth.nasa.gov/view.php?id=57752), and [Altitude.png](https://asterweb.jpl.nasa.gov/gdem.asp),
* **Tom Patterson** for [Clouds.jpg](http://www.shadedrelief.com/natural3/pages/textures.html), [Rivers.png](http://www.shadedrelief.com/natural3/pages/extra.html),
* **Natural Earth** for [Pastel.png](http://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-2/) and their [detailed vector data](http://www.naturalearthdata.com/downloads/),
* **The Apache Commons** for their [complex mathematics code](https://commons.apache.org/proper/commons-math/),
* **Technische Universit&auml;t Berlin** for their [complex mathematics code](http://www3.math.tu-berlin.de/jtem/ellipticFunctions/),
* **Gene Keyes** for his [impressively in-depth documentation of his map projection](http://www.genekeyes.com/CKOG-OOo/7-CKOG-illus-&-coastline.html),
* **Robert Gray** for his [research on Buckminster Fuller's map projection](http://www.rwgrayprojects.com/rbfnotes/toc.html),
* **AuthaGraph Co., Ltd.** for their [vague information about their map projection](http://www.authagraph.com/projects/description/%E3%80%90%E4%BD%9C%E5%93%81%E8%A7%A3%E8%AA%AC%E3%80%91%E8%A8%98%E4%BA%8B01/?lang=en), and
* **Wikipedia** for their [substantial collection of map projection information and equations](https://en.wikipedia.org/wiki/List_of_map_projections).
