# Map-Projections
A class to create custom maps of the Earth's surface. There are thousands of combinations of color-schemes, projections, and aspects. Includes Mercator, Gall-Peters, Orthographic, Peirce Quincuncial, and More!

## Installation
If you are a fancy Windows user, I recommend the convenient [fancy Windows binaries](https://github.com/jkunimune15/Map-Projections/releases). Double-click to install them and then keep pressing buttons until something good happens. If you see a map, you're in the right place.

If you are not on Windows or are otherwise not fancy enough to deserve such executables, simply double-click on the .jar files in the main directory and, if you have [Java](https://java.com/en/download/) installed (10/10 would recommend), it should just run without any set-up.

You could also compile and run the source code, but if you do, there are a few dependencies. All of them can be obtained as .jar files:

* [Apache Commons Mathematics Library](http://commons.apache.org/proper/commons-math/download_math.cgi)
* [Java Tools for Experimental Mathematics "ellipticFunctions" package](http://www3.math.tu-berlin.de/jtem/downloads.html), which requires their "mfc" package

## Usage
There are three executable files and three other runnable Java scripts. These are, in order:

* `MapDesignerRaster.jar` &ndash; The original program. Create custom oblique raster images of the Earth's surface using a variety of algorithms called _projections_.  
* `MapDesignerVector.jar` &ndash; The same idea, but working in vector images instead in case you want to cut a vinyl sticker or something.  
* `MapAnalyzer.jar` &ndash; See graphs and figures quantifying the amount of scale and angular distortion present in each map projection.  
* `MapPlotter.java` &ndash; Plot a large group of map projections by the amount of distortion they produce.  
* `MapOptimizer.java` &ndash; Run gradient descent on parametric projections to minimize their distortion.  
* `MapExplainer.java` &ndash; Generate an HTML blurb outlining and displaying every map projection.

The executable applications all have similar layouts that let you select an input equirectangular map, a projection, an aspect (where the North Pole is situated with respect to the projection), and parameters if applicable. Go crazy! There are a practically unlimited number of combinations.

The runnable scripts just kind of work on their own. Those ones aren't really meant for mass consumption.

## Wherefore?
I'll write a little blurb here later.

For some examples, check out the `output` folder. For more information, go to [jkunimune15.github.io/Map-Projections](https://jkunimune15.github.io/Map-Projections).

## Credits
While I wrote all of the code in this repository myself, and I created several of the simpler images from scratch, other people did help. Here's a comprehensive list.
* **The NASA** for [Basic.png](https://visibleearth.nasa.gov/view.php?id=57730), [Satellite.jpg](https://visibleearth.nasa.gov/view.php?id=57752), and [Altitude.png](https://asterweb.jpl.nasa.gov/gdem.asp)
* **Tom Patterson** for [Pastel.jpg](http://www.shadedrelief.com/natural3/pages/textures.html) and [Rivers.png](http://www.shadedrelief.com/natural3/pages/extra.html)
* **RokerHRO** for [the indicatrix layer of Tissot.jpg](https://commons.wikimedia.org/wiki/File:Tissot_10deg.png)
* **Crates** for [Landmasses.svg and the landmass layer of Basic.svg](https://commons.wikimedia.org/wiki/File:World_map_blank_without_borders.svg)
* **The CIA** for [Political.svg and Political.png](https://commons.wikimedia.org/wiki/File:BlankMap-World6-Equirectangular.svg).
* **The Apache Commons** for their [complex mathematics code](https://commons.apache.org/proper/commons-math/)
* **Technische Universit&auml;t Berlin** for their [complex mathematics code](http://www3.math.tu-berlin.de/jtem/ellipticFunctions/)
* **AuthaGraph Co., Ltd.** for their [vague information about their map projection](http://www.authagraph.com/projects/description/%E3%80%90%E4%BD%9C%E5%93%81%E8%A7%A3%E8%AA%AC%E3%80%91%E8%A8%98%E4%BA%8B01/?lang=en), which inspired several of my own.
* **Wikipedia** for their [impressive collection of map projection information and equations](https://en.wikipedia.org/wiki/List_of_map_projections)
