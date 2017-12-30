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

import java.awt.image.BufferedImage;
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
import javafx.scene.image.ImageView;
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
import maps.Snyder;
import maps.Tetrahedral;
import maps.Tobler;
import maps.WinkelTripel;
import utils.ImageUtils;
import utils.PixelMap;
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
	
	private static final Projection[] PROJ_ARR = { Cylindrical.MERCATOR,
			Cylindrical.EQUIRECTANGULAR, Cylindrical.EQUAL_AREA, Cylindrical.GALL,
			Azimuthal.STEREOGRAPHIC, Azimuthal.POLAR, Azimuthal.EQUAL_AREA, Azimuthal.GNOMONIC,
			Azimuthal.PERSPECTIVE, Conic.LAMBERT, Conic.EQUIDISTANT, Conic.ALBERS, Tetrahedral.LEE,
			Tetrahedral.ACTUAUTHAGRAPH, Tetrahedral.AUTHAPOWER, Tetrahedral.AUTHAGRAPH,
			Pseudocylindrical.SINUSOIDAL, Pseudocylindrical.MOLLWEIDE, Tobler.TOBLER, Misc.AITOFF,
			Misc.VAN_DER_GRINTEN, Robinson.ROBINSON, WinkelTripel.WINKEL_TRIPEL,
			Misc.PEIRCE_QUINCUNCIAL, Misc.TWO_POINT_EQUIDISTANT, Misc.HAMMER_RETROAZIMUTHAL, Snyder.GS50,
			Pseudocylindrical.LEMONS, MyProjections.EXPERIMENT, MyProjections.PSEUDOSTEREOGRAPHIC,
			MyProjections.TWO_POINT_EQUALIZED };
	
	
	private Node aspectSelector;
	private Button updateBtn, saveMapBtn;
	private double[] aspect;
	private PixelMap input;
	private ImageView display;
	private MapConfigurationDialog configDialog;
	
	
	
	public MapDesignerRaster() {
		super("Map Designer");
	}
	
	
	
	@Override
	public void start(Stage root) {
		super.start(root);
		new Thread(() -> {
			setInput(new File("input/Basic.jpg")); //TODO: this should cause the buttons to grey out
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
		pane.setMinHeight(IMG_WIDTH);
		
		final HBox gui = new HBox(10, layout, pane);
		gui.setAlignment(Pos.CENTER);
		StackPane.setMargin(gui, new Insets(10));
		
		return gui;
	}
	
	
	private void setInput(File file) {
		updateBtn.setDisable(true);
		saveMapBtn.setDisable(true);
		
		try {
			input = new PixelMap(file);
		} catch (IOException e) {
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
		display.setImage(SwingFXUtils.toFXImage(makeImage(), null));
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
		BufferedImage theMap = makeImage(outDims[0], outDims[1], smoothing, pBar); //calculate
		
		pBar.setProgress(-1);
		
		final String filename = file.getName();
		final String extension = filename.substring(filename.lastIndexOf('.')+1);
		try {
			ImageIO.write(theMap, extension, file); //save
		} catch (IOException e) {
			showError("Failure!",
					"Could not access "+file.getAbsolutePath()+". It's possible that another program has it open.");
		}
	}
	
	
	private BufferedImage makeImage() {
		final double aspRat = this.getProjection().getAspectRatio();
		if (aspRat >= 1)
			return makeImage(IMG_WIDTH, (int)Math.max(IMG_WIDTH/aspRat,1), 1, null);
		else
			return makeImage((int)Math.max(IMG_WIDTH*aspRat,1), IMG_WIDTH, 1, null);
	}
	
	private BufferedImage makeImage(int width, int height, int step, ProgressBarDialog pBar) {
		final double[] pole = aspect.clone();
		final Projection pjc = this.getProjection();
		final BufferedImage out = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		for (int y = 0; y < out.getHeight(); y ++) {
			if (pBar != null)
				pBar.setProgress((double) y/out.getHeight());
			for (int x = 0; x < out.getWidth(); x ++) {
				int[] colors = new int[step*step];
				for (int dy = 0; dy < step; dy ++) {
					for (int dx = 0; dx < step; dx ++) {
						final double X = ((x+(dx+.5)/step)/out.getWidth() - 1/2.) *pjc.getWidth();
						final double Y = (1/2. - (y+(dy+.5)/step)/out.getHeight()) *pjc.getHeight();
						final double[] coords = this.getProjection().inverse(X, Y, pole);
						if (coords != null) //if it is null, the default (0:transparent) is used
							colors[step*dy+dx] = input.getArgb(coords[0], coords[1]);
					}
				}
				out.setRGB(x, y, ImageUtils.blend(colors));
			}
		}
		return out;
	}
}
