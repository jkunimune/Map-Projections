# Map-Projections
A class to create custom maps of the Earth's surface. There are thousands of combinations of color-schemes, projections, and aspects. Includes Mercator, Gall-Peters, Orthographic, Peirce Quincuncial, and More!

<img src="output/Political equal-area.png" alt="Tobler hyperelliptical political map" width="300px"/>
<img src="output/AuthaGraph indicatrices.png" alt="AuthaGraph imitation with Tissot's indicatrices of distortion" width="300px"/>
<img src="output/Better Guyou.jpg" alt="Guyou physical map" width="300px"/>

## Running the programs

There are three main programs here:
* `MapDesignerRaster.jar` &ndash; The original program. Create custom oblique raster images of the Earth's surface using a variety of algorithms called _projections_.
* `MapDesignerVector.jar` &ndash; The same idea, but working in vector images instead in case you want to cut a vinyl sticker or something.
* `MapAnalyzer.jar` &ndash; See graphs and figures quantifying the amount of scale and angular distortion present in each map projection.

They all have similar layouts that let you select an input equirectangular map, a projection, an aspect (where the North Pole is situated with respect to the map), and parameters if applicable.
Go crazy! There are a practically unlimited number of combinations.

I will note that while I think the interface is mostly intuitive, there are a couple of things where I never got around to making the proper GUI elements, so you won't be able to figure out on your own.
The first is the fact that the graticule checkbox, available on Map Designer Raster, draws its formatting information from a file called `graticule.txt` in the input folder.
If you're fine with the default styling, or if you don't use that checkbox, don't worry about it.  But if you do use that checkbox, know that you can alter the color and width of the lines by editing that file.

The twoth is the existence of truncated inputs.
If you load an input map with the word "octant" in the filename (all lowercase), then the program will load it into the octant bounded by 0°N, 90°N, 0°E, and 90°E.
This is useful if you have very large inputs and/or memory constraints.
The output will still be sized as though the entire map were there, unless it's a projection that doesn't show the entire globe ("Cahill–Keyes (single octant)" does not show the entire globe and is in fact specifically designed to work with this feature.)

### Running by double-clicking

If you are a fancy Windows user, I recommend the convenient [fancy Windows binaries](https://github.com/jkunimune/Map-Projections/releases). Double-click to install them and then keep pressing buttons until something good happens. If you see a map, you're in the right place.

### Running from the command line

If you are not on Windows or are otherwise not fancy enough to deserve such binaries, there are also equivalent `.jar` files in the main directory.
For that, you'll need to download and install [Java](https://www.java.com/download/).
You can use any version 8 or higher.
Once you've installed it, you can run the programs like so *if* you're on Java 8 specifically:
~~~bash
java -jar MapDesignerRaster.jar
~~~
Change the filename at the end accordingly, of course.

If you have Java 9 or higher, that's fine too, but you'll need to also install [JavaFX](https://gluonhq.com/products/javafx/).
Make sure you get the SDK JavaFX, not the Jmods JavaFX; the version doesn't matter.
Once you've unzipped JavaFX into some directory – let’s say, for example, `/home/jkunimune/javafx-sdk-17` – then you can run the programs like so:

~~~bash
java --module-path '/home/jkunimune/javafx-sdk-17/lib' --add-modules javafx.controls,javafx.swing -jar MapDesignerRaster.jar
~~~

I think this syntax might be somewhat platform dependent, but I can’t really remember.
If you’re having problems, try forward slashes instead of backslashes or double quotes instead of single quotes.

## Building from source

If you want to edit the code, or use some of the deeper functionality not meant for mass consumption, you can also compile and run the Java source code.
In addition to three `.java` files corresponding to the three executables, there are also these runnable scripts:

* `MapPlotter.java` &ndash; Plot a large group of map projections by the amount of distortion they produce.  
* `src/app/MapOptimizer.java` &ndash; Run gradient descent on parametric projections to minimize their distortion.  
* `src/app/MapExplainer.java` &ndash; Generate an HTML blurb outlining and displaying every map projection.
* `src/app/MapConverter.java` &ndash; Generate a bunch of maps in one projection from a bunch of input images.

To run these, you’ll need to install some dependencies in addition to JavaFX as mentioned above; you can get them as `.jar` files:

* [Apache Commons Mathematics Library](http://commons.apache.org/proper/commons-math/download_math.cgi) (you can just download and unzip the whole thing)
* [Java Tools for Experimental Mathematics "ellipticFunctions" package](http://www3.math.tu-berlin.de/jtem/downloads.html) (you only need "ellipticFunctions.jar" and "mfc.jar")

Once you have those and put them in, for example, `/home/jkunimune/apache` and `/home/jkunimune/jtem`, you can compile and run with
~~~bash
javac --module-path '/home/jkunimune/javafx-sdk-17/lib:/home/jkunimune/apache/commons-math3-3.6.1.jar:/home/jkunimune/jtem' --add-modules javafx.controls,javafx.swing,ellipticFunctions --source-path=src src/apps/MapPlotter.java
java --class-path '/home/jkunimune/javafx-sdk-17/lib/javafx.controls.jar:/home/jkunimune/javafx-sdk-17/lib/javafx.swing.jar:/home/jkunimune/apache/commons-math3-3.6.1.jar:/home/jkunimune/jtem/ellipticFunctions.jar:src' apps.MapPlotter
~~~
As with running the .jar file, the syntax might be somewhat platform-dependent; the colons might need to be semicolons, for example.

To build the JAR and EXE files, you use the Ant file `build.xml`.
I don't know what commands are used for that; IntelliJ IDEA does it for me.
The main trick here is you *must* use the Java 8 JDK because the tool I use, Javapackager, was dropped in Java 9.
If you use Javapackager in conjunction with a newer JDK, it will seem to work at first, but the resulting EXEs won't run.
You may have to put `javapackager.exe` in your path and install the WiX Toolset 4 and Inno Setup 5 (not Inno Setup 6!).
Javapackager does often have problems but doesn't have useful error messages (every failure mode appears to be "builder did not produce a bundle"),
so I apologize about that.
Make sure you set the directories at the top of the file according to your file system, make sure the compiled Java is in the `bin/` directory, make sure none of the jar files have restricted read access or are being used by any other programs.

There are also some Python files used to generate SVG inputs for MapDesignerVector in the src/zupplemental directory.
To run those, you'll need a few packages from [PyPI](https://pypi.python.org/pypi).
Just install Python 3 and pip, and then call the following from a command line (or use Anaconda or something, I don't know. Up to you).
~~~~
pip3 install numpy pyshp shapely
~~~~

Note that `compose_maps.py` requires input data from [naturalearthdata.com](http://www.naturalearthdata.com/downloads/), which should be downloaded and placed in `src/zupplemental/shapefiles/`.

## Wherefore?
I'll write a little blurb here later.

For some examples, check out the `output` folder (all were created with this program but some also involved some postprocessing in paint.net). For more information, go to [jkunimune.github.io/Map-Projections](https://jkunimune15.github.io/Map-Projections).

## Credits
While I wrote all of the code in this repository myself, and I created several of the simpler images from scratch, other people did help. Here's a comprehensive list.
* **The NASA** for [Basic.png](https://visibleearth.nasa.gov/view.php?id=57730), [Large.jpg](https://visibleearth.nasa.gov/view.php?id=57752), and [Altitude.png](https://asterweb.jpl.nasa.gov/gdem.asp),
* **Tom Patterson** for [Clouds.jpg](http://www.shadedrelief.com/natural3/pages/textures.html), [Rivers.png](http://www.shadedrelief.com/natural3/pages/extra.html),
* **Natural Earth** for [Pastel.png](http://www.naturalearthdata.com/downloads/50m-raster-data/50m-natural-earth-2/) and their [detailed vector data](http://www.naturalearthdata.com/downloads/),
* **The Apache Commons** for their [complex mathematics code](https://commons.apache.org/proper/commons-math/),
* **Technische Universit&auml;t Berlin** for their [complex mathematics code](http://www3.math.tu-berlin.de/jtem/ellipticFunctions/),
* **Gene Keyes** for his [impressively in-depth documentation of his map projection](http://www.genekeyes.com/CKOG-OOo/7-CKOG-illus-&-coastline.html),
* **Robert Gray** for his [research on Buckminster Fuller's map projection](http://www.rwgrayprojects.com/rbfnotes/toc.html),
* **AuthaGraph Co., Ltd.** for their [vague information about their map projection](http://www.authagraph.com/projects/description/%E3%80%90%E4%BD%9C%E5%93%81%E8%A7%A3%E8%AA%AC%E3%80%91%E8%A8%98%E4%BA%8B01/?lang=en),
* **Wikipedia** for its [substantial collection of map projection information and equations](https://en.wikipedia.org/wiki/List_of_map_projections), and
* **John P. Snyder** for his [even more substantial collection of map projection information and equations](https://press.uchicago.edu/ucp/books/book/chicago/F/bo3632853.html).
