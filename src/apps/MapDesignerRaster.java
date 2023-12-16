/**
 * MIT License
 * 
 * Copyright (c) 2017 Justin Kunimune
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package apps;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Supplier;

import dialogs.MapConfigurationDialog;
import image.ImageUtils;
import image.PixelMap;
import image.SavableImage;
import image.TruncatedPixelMap;
import javafx.concurrent.Task;
import javafx.embed.swing.SwingFXUtils;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Separator;
import javafx.scene.image.ImageView;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import maps.Projection;
import utils.Flag;
import utils.MutableDouble;
import utils.Procedure;
import utils.Shape;

import static java.lang.Double.isNaN;
import static java.lang.Integer.parseInt;
import static java.lang.Math.PI;
import static java.lang.Math.toRadians;
import static utils.Math2.max;
import static utils.Math2.min;

/**
 * An application to make raster oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapDesignerRaster extends MapApplication {
	
	public static void main(String[] args) {
		launch(args);
	}
	
	
	
	private static final FileChooser.ExtensionFilter[] READABLE_TYPES = {
			new FileChooser.ExtensionFilter("All image types", "*.png","*.jpg","*.jpeg","*.jpe","*.jfif","*.gif"),
			new FileChooser.ExtensionFilter("PNG", "*.png"),
			new FileChooser.ExtensionFilter("JPG", "*.jpg","*.jpeg","*.jpe","*.jfif"),
			new FileChooser.ExtensionFilter("GIF", "*.gif") };
	private static final FileChooser.ExtensionFilter[] RASTER_TYPES = {
			new FileChooser.ExtensionFilter("PNG", "*.png"),
			new FileChooser.ExtensionFilter("JPG", "*.jpg"),
			new FileChooser.ExtensionFilter("GIF", "*.gif") };
	
	private static final double GRATICULE_PRECISION = 0.02;
	
	private Region aspectSelector;
	private double[] aspect;
	private Flag cropAtIDL;
	private MutableDouble graticuleSpacing;
	private PixelMap input;
	private ImageView display;
	private MapConfigurationDialog configDialog;
	
	
	
	public MapDesignerRaster() {
		super("Map Designer");
	}
	
	
	
	@Override
	public void start(Stage root) {
		super.start(root);
		Task<Void> task = setInputTask(new File("input/Basic.png"));
		task.setOnSucceeded((event) -> new Thread(calculateTaskForUpdate()).start());
		new Thread(task).start();
	}
	
	
	@Override
	protected Region makeWidgets() {
		this.aspect = new double[3];
		this.cropAtIDL = new Flag(true);
		this.graticuleSpacing = new MutableDouble();
		final Region inputSelector = buildInputSelector(READABLE_TYPES,
				RASTER_TYPES[0], this::setInputTask);
		final Region projectionSelector = buildProjectionSelector(this::updateAspect);
		this.aspectSelector = buildAspectSelector(this.aspect, Procedure.NONE);
		final Region parameterSelector = buildParameterSelector(Procedure.NONE);
		final Region optionPane = buildOptionPane(cropAtIDL, graticuleSpacing);
		final Region updateBtn = buildUpdateButton("Update Map", this::calculateTaskForUpdate);
		final Region saveMapBtn = buildSaveButton(
				"map", RASTER_TYPES, RASTER_TYPES[0],
				this::collectFinalSettings, this::calculateTaskForSaving);
		final HBox buttons = new HBox(H_SPACE, updateBtn, saveMapBtn);
		buttons.setAlignment(Pos.CENTER);
		
		final VBox layout = new VBox(V_SPACE,
				inputSelector, new Separator(), projectionSelector, new Separator(),
				aspectSelector, parameterSelector, new Separator(), optionPane, new Separator(),
				buttons);
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(GUI_WIDTH);
		
		this.display = new ImageView();
		this.display.setFitWidth(IMG_SIZE);
		this.display.setFitHeight(IMG_SIZE);
		this.display.setPreserveRatio(true);
		final StackPane pane = new StackPane(display);
		pane.setMinWidth(IMG_SIZE);
		pane.setMinHeight(IMG_SIZE);
		
		final HBox gui = new HBox(MARGIN, layout, pane);
		gui.setAlignment(Pos.CENTER);
		StackPane.setMargin(gui, new Insets(MARGIN));
		
		return gui;
	}
	
	
	private Task<Void> setInputTask(File file) {
		return new Task<Void>() {
			protected Void call() throws IOException {
				if (file.getName().contains("octant"))
					input = new TruncatedPixelMap(file);
				else
					input = new PixelMap(file);
				return null;
			}
			
			protected void failed() {
				showError("File not found!",
						"Couldn't find "+file.getAbsolutePath()+".");
			}
		};
	}
	
	
	private void updateAspect() {
		aspectSelector.setVisible(this.getProjection().hasAspect());
	}
	
	
	private boolean collectFinalSettings() {
		loadParameters();
		final double ratio = getProjection().getShape().aspectRatio;
		this.configDialog = new MapConfigurationDialog(ratio);
		this.configDialog.showAndWait();
		return this.configDialog.getResult();
	}
	
	
	private Task<SavableImage> calculateTaskForUpdate() {
		loadParameters();
		if (getProjection().getShape().aspectRatio >= 1) //fit it to an IMG_SIZE x IMG_SIZE box
			return calculateTask(
					IMG_SIZE, (int)max(1,IMG_SIZE/getProjection().getShape().aspectRatio), 1);
		else
			return calculateTask(
					(int)max(1,IMG_SIZE*getProjection().getShape().aspectRatio), IMG_SIZE, 1);
	}
	
	private Task<SavableImage> calculateTaskForSaving() {
		int[] outDims = configDialog.getDims();
		int step = configDialog.getSmoothing();
		return calculateTask(outDims[0], outDims[1], step);
	}
	
	private Task<SavableImage> calculateTask(int width, int height, int step) {
		return calculateTask(width, height, step,
				input, getProjection(), aspect.clone(), cropAtIDL.isSet(), graticuleSpacing.get(),
				display);
	}
	
	/**
	 * Prepare a task to create a new savable raster map.
	 * @param width - The desired map width.
	 * @param height - The desired map height.
	 * @param step - The desired amount of smoothing to apply.
	 * @param input - The input equirectangular image.
	 * @param proj - The Projection to do the mapping.
	 * @param aspect - The oblique axis of the map.
	 * @param crop - Should points with extreme longitudes be hidden?
	 * @param gratSpacing - The number of degrees between graticule lines, or 0 for no graticule.
	 * @param display - The ImageViewer in which to put the new image, or null if you don't want us
	 * 		to do that.
	 * @return the task that will call MapDesignerRaster.calculate()
	 */
	public static Task<SavableImage> calculateTask(int width, int height, int step,
			PixelMap input, Projection proj, double[] aspect, boolean crop, double gratSpacing,
			ImageView display) {
		return new Task<SavableImage>() {
			private BufferedImage map;

			protected SavableImage call() {
				map = MapDesignerRaster.calculate(
					  width, height, step, input, proj,
					  aspect, crop, gratSpacing,
					  this::updateProgress, this::updateMessage, this::isCancelled);
				return SavableImage.savable(map);
			}
			
			protected void failed() {
				getException().printStackTrace();
				if (this.getException() instanceof OutOfMemoryError)
					showError("Failure!", "Java memory constraints forbid that image size. I'm sorry; I never should have let you input a number that big. Please try something smaller.");
			}
			
			protected void succeeded() {
				if (display != null)
					display.setImage(SwingFXUtils.toFXImage(map, null));
			}
		};
	}

	/**
	 * Create a new savable raster map.
	 * @param width - The desired map width.
	 * @param height - The desired map height.
	 * @param step - The desired amount of smoothing to apply.
	 * @param input - The input equirectangular image.
	 * @param proj - The Projection to do the mapping.
	 * @param aspect - The oblique axis of the map.
	 * @param crop - Should points with extreme longitudes be hidden?
	 * @param gratSpacing - The number of degrees between graticule lines, or 0 for no graticule.
	 * @return the projected image
	 */
	public static BufferedImage calculate(int width, int height, int step,
										 PixelMap input, Projection proj,
										 double[] aspect, boolean crop,
										 double gratSpacing,
										 BiConsumer<Integer, Integer> updateProgress,
										 Consumer<String> updateMessage,
										 Supplier<Boolean> isCancelled) {
		if (updateProgress == null)
			updateProgress = (i, j) -> {};
		if (updateMessage == null)
			updateMessage = (s) -> {};
		if (isCancelled == null)
			isCancelled = () -> false;

		updateProgress.accept(-1, 1);
		updateMessage.accept("Generating map\u2026");

		Shape domain = proj.getShape();
		BufferedImage theMap = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB); //why is this a BufferedImage when the rest of this program uses JavaFX? Because the only JavaFX alternatives are WritableImage, which doesn't do anything but single-pixel-editing, and Canvas, which doesn't properly support transparency.
		for (int y = 0; y < theMap.getHeight(); y ++) { //iterate through the map, filling in pixels
			if (isCancelled.get()) 	return null;
			updateProgress.accept(y, theMap.getHeight());
			for (int x = 0; x < theMap.getWidth(); x ++) {
				int[] colors = new int[step*step];
				for (int dy = 0; dy < step; dy ++) {
					for (int dx = 0; dx < step; dx ++) {
						double X = domain.xMin + (domain.xMax - domain.xMin)*(x+(dx+.5)/step)/width;
						double Y = domain.yMax - (domain.yMax - domain.yMin)*(y+(dy+.5)/step)/height;
						double[] coords = proj.inverse(X, Y, aspect, crop);
						if (coords != null) { //if it is null, the default (0:transparent) is used
							if (isNaN(coords[0]) || isNaN(coords[1]))
								colors[step*dy+dx] = 0xffff00ff;
							else
								colors[step*dy+dx] = input.getArgb(coords[0], coords[1]);
						}
					}
				}
				theMap.setRGB(x, y, ImageUtils.blend(colors));
			}
		}

		if (gratSpacing != 0) { //draw the graticule, if desired
			if (isCancelled.get()) 	return null;
			updateProgress.accept(-1, 1);
			updateMessage.accept("Drawing graticule\u2026");

			int r = 255, g = 255, b = 255, a = 255;
			float lineWidth = (float)(min(width, height)/300);
			BufferedReader fileReader = null;
			try {
				fileReader = new BufferedReader(new FileReader("input/graticule.txt"));
				r = parseInt(fileReader.readLine().split(":")[1].trim());
				g = parseInt(fileReader.readLine().split(":")[1].trim());
				b = parseInt(fileReader.readLine().split(":")[1].trim());
				a = parseInt(fileReader.readLine().split(":")[1].trim());
				lineWidth = Float.parseFloat(fileReader.readLine().split(":")[1]);
			} catch (NumberFormatException | IOException e) {
				e.printStackTrace();
			} finally {
				if (fileReader != null)
					try {
						fileReader.close();
					} catch (IOException ignored) { }
			}

			ImageUtils.drawSVGPath(
				  proj.drawGraticule(toRadians(gratSpacing), GRATICULE_PRECISION,
									 width, height, PI/2, PI, aspect),
				  new Color(r, g, b, a), lineWidth,
				  true, (Graphics2D)theMap.getGraphics());
		}

		return theMap;
	}
}
