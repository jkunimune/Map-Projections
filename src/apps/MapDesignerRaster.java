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

import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;

import dialogs.MapConfigurationDialog;
import dialogs.ProgressBarDialog;
import javafx.application.Platform;
import javafx.embed.swing.SwingFXUtils;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.Separator;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.image.PixelReader;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import maps.Azimuthal;
import maps.Conic;
import maps.Cylindrical;
import maps.Misc;
import maps.MyProjections;
import maps.Projection;
import maps.Pseudocylindrical;
import maps.Robinson;
import maps.Tetrahedral;
import maps.Tobler;
import maps.WinkelTripel;
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
			new FileChooser.ExtensionFilter("JPG", "*.jpg","*.jpeg","*.jpe","*.jfif"),
			new FileChooser.ExtensionFilter("GIF", "*.gif") };
	
	private static final Projection[] PROJ_ARR = { Cylindrical.MERCATOR,
			Cylindrical.EQUIRECTANGULAR, Cylindrical.EQUAL_AREA, Cylindrical.GALL,
			Azimuthal.STEREOGRAPHIC, Azimuthal.POLAR, Azimuthal.EQUAL_AREA, Azimuthal.ORTHOGRAPHIC,
			Azimuthal.GNOMONIC, Conic.LAMBERT, Conic.EQUIDISTANT, Conic.ALBERS, Tetrahedral.LEE,
			Tetrahedral.TETRAGRAPH, Tetrahedral.AUTHAGRAPH, Pseudocylindrical.SINUSOIDAL,
			Pseudocylindrical.MOLLWEIDE, Tobler.TOBLER, Misc.AITOFF, Misc.VAN_DER_GRINTEN,
			Robinson.ROBINSON, WinkelTripel.WINKEL_TRIPEL, Misc.PEIRCE_QUINCUNCIAL, Misc.GUYOU,
			Misc.TWO_POINT_EQUIDISTANT, Misc.HAMMER_RETROAZIMUTHAL, Pseudocylindrical.LEMONS,
			MyProjections.MAGNIFIER, MyProjections.EXPERIMENT, MyProjections.PSEUDOSTEREOGRAPHIC,
			MyProjections.HYPERELLIPOWER, Tetrahedral.TETRAPOWER, Tetrahedral.TETRAFILLET };
	
	
	private Node aspectSelector;
	private Button updateBtn, saveMapBtn;
	private double[] aspect;
	private Image input;
	private ImageView display;
	private MapConfigurationDialog configDialog;
	
	
	
	public MapDesignerRaster() {
		super("Map Designer");
	}
	
	
	
	@Override
	public void start(Stage root) {
		super.start(root);
		new Thread(() -> {
			setInput(new File("input/basic.jpg")); //TODO: this should cause the buttons to grey out
			updateMap();
		}).start();
	}
	
	
	@Override
	protected Node makeWidgets() {
		this.aspect = new double[3];
		final Node inputSelector = buildInputSelector(READABLE_TYPES,
				RASTER_TYPES[0], this::setInput);
		final Node projectionSelector = buildProjectionSelector(PROJ_ARR,
				this::hideAspect);
		this.aspectSelector = buildAspectSelector(this.aspect, Procedure.NONE);
		final Node parameterSelector = buildParameterSelector(Procedure.NONE);
		this.updateBtn = buildUpdateButton(this::updateMap);
		this.saveMapBtn = buildSaveButton(true, "map", RASTER_TYPES,
				RASTER_TYPES[0], this::collectFinalSettings, this::calculateAndSaveMap);
		final HBox buttons = new HBox(5, updateBtn, saveMapBtn);
		buttons.setAlignment(Pos.CENTER);
		aspectSelector.managedProperty().bind(aspectSelector.visibleProperty());
		
		final VBox layout = new VBox(5,
				inputSelector, new Separator(), projectionSelector,
				new Separator(), aspectSelector, parameterSelector,
				new Separator(), buttons);
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(GUI_WIDTH);
		
		this.display = new ImageView();
		this.display.setFitWidth(IMG_WIDTH);
		this.display.setFitHeight(IMG_WIDTH);
		this.display.setPreserveRatio(true);
		final StackPane pane = new StackPane(display);
		pane.setMinWidth(IMG_WIDTH);
		
		final HBox gui = new HBox(10, layout, pane);
		gui.setAlignment(Pos.CENTER);
		StackPane.setMargin(gui, new Insets(10));
		
		return gui;
	}
	
	
	private void setInput(File file) {
		updateBtn.setDisable(true);
		saveMapBtn.setDisable(true);
		
		try {
			input = new Image(file.toURI().toString());
		} catch (IllegalArgumentException e) {
			final Alert alert = new Alert(Alert.AlertType.ERROR);
			alert.setHeaderText("File not found!");
			alert.setContentText("Couldn't find "+file.getAbsolutePath()+".");
			Platform.runLater(alert::showAndWait);
		} finally {
			updateBtn.setDisable(false);
			saveMapBtn.setDisable(false);
		}
	}
	
	
	private void hideAspect() {
		aspectSelector.setVisible(this.getProjection().hasAspect());
	}
	
	
	private void updateMap() {
		loadParameters();
		display.setImage(makeImage(this.getProjection().map(IMG_WIDTH, aspect.clone())));
	}
	
	
	private boolean collectFinalSettings() {
		loadParameters();
		final double ratio = getProjection().getAspectRatio();
		this.configDialog = new MapConfigurationDialog(ratio);
		this.configDialog.showAndWait();
		return this.configDialog.getResult();
	}
	
	
	private void calculateAndSaveMap(File file, ProgressBarDialog pBar) {
		final int[] outDims = configDialog.getDims();
		final int smoothing = configDialog.getSmoothing();
		double[][][] points =this.getProjection().map(
				outDims[0]*smoothing, outDims[1]*smoothing,
				aspect.clone(), pBar::setProgress);
		pBar.setProgress(-1);
		Image theMap = makeImage(points, smoothing); //calculate
		
		final String filename = file.getName();
		final String extension = filename.substring(filename.lastIndexOf('.')+1);
		try {
			ImageIO.write(SwingFXUtils.fromFXImage(theMap,null), extension, file); //save
		} catch (IOException e) {
			showError("Failure!",
					"Could not access "+file.getAbsolutePath()+". It's possible that another program has it open.");
		}
	}
	
	
	private Image makeImage(double[][][] points) {
		return makeImage(points, 1);
	}
	
	private Image makeImage(double[][][] points, int step) {
		final PixelReader in = input.getPixelReader();
		final WritableImage out = new WritableImage(points[0].length/step, points.length/step);
		for (int y = 0; y < out.getHeight(); y ++) {
			for (int x = 0; x < out.getWidth(); x ++) {
				int[] colors = new int[step*step];
				for (int dy = 0; dy < step; dy ++) {
					for (int dx = 0; dx < step; dx ++) {
						colors[step*dy+dx] = getArgb(points[step*y+dy][step*x+dx], in, input.getWidth(), input.getHeight());
					}
				}
				out.getPixelWriter().setArgb(x, y, blend(colors));
			}
		}
		return out;
	}
	
	
	private static int getArgb(double[] coords, PixelReader in,
			double inWidth, double inHeight) { // returns the color of any coordinate on earth
		if (coords == null) 	return 0;
		
		double x = 1/2.0 + coords[1]/(2*Math.PI);
		x = (x - Math.floor(x)) * inWidth;
		
		double y = inHeight/2.0 - coords[0]*inHeight/(Math.PI);
		if (y < 0)
			y = 0;
		else if (y >= inHeight)
			y = inHeight - 1;
		
		return in.getArgb((int) x, (int) y);
	}
	
	
	private static final int blend(int[] colors) {
		int a_tot = 0;
		int r_tot = 0;
		int g_tot = 0;
		int b_tot = 0;
		for (int argb: colors) {
			double a = ((argb >> 24) & 0xFF);
			a_tot += a;
			r_tot += a*((argb >> 16) & 0xFF);
			g_tot += a*((argb >> 8) & 0xFF);
			b_tot += a*((argb >> 0) & 0xFF);
		}
		if (a_tot == 0)	return 0;
		else
			return (a_tot/colors.length<<24) +
					(r_tot/a_tot<<16) + (g_tot/a_tot<<8) + (b_tot/a_tot);
	}

}