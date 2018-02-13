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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Path2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import dialogs.MapConfigurationDialog;
import image.ImageUtils;
import image.PixelMap;
import image.SavableImage;
import image.SVGMap.Command;
import image.SVGMap.Path;
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

/**
 * An application to make raster oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapDesignerRaster extends MapApplication {
	
	public static final void main(String[] args) {
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
	
	private static final float GRATICULE_WIDTH = 0.5f;
	private static final Color GRATICULE_COLOR = Color.WHITE;
	
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
		this.cropAtIDL = new Flag(false);
		this.graticuleSpacing = new MutableDouble();
		final Region inputSelector = buildInputSelector(READABLE_TYPES,
				RASTER_TYPES[0], this::setInputTask);
		final Region projectionSelector = buildProjectionSelector(this::updateAspect);
		this.aspectSelector = buildAspectSelector(this.aspect, Procedure.NONE);
		final Region parameterSelector = buildParameterSelector(Procedure.NONE);
		final Region optionPane = buildOptionPane(cropAtIDL, graticuleSpacing);
		final Region updateBtn = buildUpdateButton("Update Map", this::calculateTaskForUpdate);
		final Region saveMapBtn = buildSaveButton(true, "map", RASTER_TYPES,
				RASTER_TYPES[0], this::collectFinalSettings, this::calculateTaskForSaving);
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
		final double ratio = getProjection().getAspectRatio();
		this.configDialog = new MapConfigurationDialog(ratio);
		this.configDialog.showAndWait();
		return this.configDialog.getResult();
	}
	
	
	private Task<SavableImage> calculateTaskForUpdate() {
		if (getProjection().isLandscape()) //either fit it to an IMG_SIZE x IMG_SIZE box
			return calculateTask(
					IMG_SIZE, (int)Math.max(IMG_SIZE/getProjection().getAspectRatio(),1), 1);
		else
			return calculateTask(
					(int)Math.max(IMG_SIZE/getProjection().getAspectRatio(),1), IMG_SIZE, 1);
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
	 * @return
	 */
	public static Task<SavableImage> calculateTask(int width, int height, int step,
			PixelMap input, Projection proj, double[] aspect, boolean crop, double gratSpacing,
			ImageView display) {
		System.out.println("Let's make a task!");
		return new Task<SavableImage>() {
			private BufferedImage theMap;
			
			protected SavableImage call() {
				System.out.println("Here we go!");
				updateProgress(-1, 1);
				updateMessage("Generating map\u2026");
				
				theMap = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB); //why is this a BufferedImage when the rest of this program uses JavaFX? Because the only JavaFX alternatives are WritableImage, which doesn't do anything but single-pixel-editing, and Canvas, which doesn't properly support transparency.
				for (int y = 0; y < theMap.getHeight(); y ++) { //iterate through the map, filling in pixels
					if (isCancelled()) 	return null;
					updateProgress(y, theMap.getHeight());
					for (int x = 0; x < theMap.getWidth(); x ++) {
						int[] colors = new int[step*step];
						for (int dy = 0; dy < step; dy ++) {
							for (int dx = 0; dx < step; dx ++) {
								double X = ((x+(dx+.5)/step)/width - 1/2.) *proj.getWidth();
								double Y = (1/2. - (y+(dy+.5)/step)/height) *proj.getHeight();
								double[] coords = proj.inverse(X, Y, aspect, crop);
								if (coords != null) //if it is null, the default (0:transparent) is used
									colors[step*dy+dx] = input.getArgb(coords[0], coords[1]);
							}
						}
						theMap.setRGB(x, y, ImageUtils.blend(colors));
					}
				}
				
				if (gratSpacing != 0) { //draw the graticule, if desired
					if (isCancelled()) 	return null;
					updateProgress(-1, 1);
					updateMessage("Drawing graticule\u2026");
					
					Graphics2D g = (Graphics2D)theMap.getGraphics();
					g.setStroke(new BasicStroke(GRATICULE_WIDTH));
					g.setColor(GRATICULE_COLOR);
					g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
					Path svgPath = proj.drawGraticule(
							Math.toRadians(gratSpacing), .02, width, height, Math.PI/2, Math.PI, aspect);
					Path2D awtPath = new Path2D.Double(Path2D.WIND_NON_ZERO, svgPath.size());
					for (Command svgCmd: svgPath) {
						if (isCancelled()) 	return null;
						switch (svgCmd.type) {
						case 'M':
							awtPath.moveTo(svgCmd.args[0], svgCmd.args[1]);
							break;
						case 'L':
							awtPath.lineTo(svgCmd.args[0], svgCmd.args[1]);
							break;
						case 'Z':
							awtPath.closePath();
						}
					}
					g.draw(awtPath);
				}
				
				return SavableImage.savable(theMap);
			}
			
			protected void failed() {
				getException().printStackTrace();
			}
			
			protected void succeeded() {
				if (display != null)
					display.setImage(SwingFXUtils.toFXImage(theMap, null));
			}
		};
	}
}
