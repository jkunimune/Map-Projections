# Map-Projections
A class to create custom maps of the Earth's surface. There are thousands of combinations of color-schemes, projections, and aspects. Includes Mercator, Gall-Peters, Orthographic, Peirce Quincuncial, and More!

## Installation
If you are a fancy Windows user, I recommend the convenient [fancy Windows binaries](https://example.org). Double-click to install them and then keep pressing buttons until something good happens. If you see a map, you're in the right place.

If you are not on Windows or are otherwise not fancy enough to deserve such executables, simply double-click on the .jar files in the main directory and, if you have [Java](https://java.com/en/download/) installed (10/10 would recommend), it should just run without any set-up.

## Features
There are three executable files and two other runnable Java scripts. These are:

* MapDesignerRaster.jar – The original program. Create custom oblique raster images of the Earth's surface using a variety of algorithms called _projections_.  
* MapDesignerVector.jar – The same idea, but working in vector images instead in case you want to cut a vinyl sticker or something.  
* MapAnalyzer.jar – See graphs and figures quantifying the amount of scale and angular distortion present in each map projection.  
* MapPlotter.java – Plot a large group of map projections by the amount of distortion they produce.  
* MapOptimizer.java – Run gradient descent on parametric projections to minimize their distortion.  

The executable applications all have similar layouts that let you select an input equirectangular map, a projection, an aspect (where the North Pole is situated with respect to the projection), and parameters if applicable. Go crazy! There are a practically unlimited number of combinations.

## Dependencies
While the excecutables are standalone, and the Jars require only Java, the source code makes use of several external libraries. These are

* Apache commons `math3`
* "mfc"
* "ellipticFunctions"

To be perfectly honest, I don't remember where I got most of these. Oops. It looks like they might be German. I would recommend looking the `de.jtem` package from (math.tu-berlin.de)[www3.math.tu.berlin.de/jtem/].

## Wherefore?
I'll write a little blurb here later.

For more information go to [jkunimune15.github.io/Map-Projections](https://jkunimune15.github.io/Map-Projections).
