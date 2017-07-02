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
import javafx.scene.control.Alert.AlertType;
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
import maps.Projection;
import util.Procedure;

/**
 * An application to make raster oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapDesignerRaster extends MapApplication {

	public static final void main(String[] args) {
		launch(args);
	}
	
	
	
	private static final FileChooser.ExtensionFilter[] RASTER_TYPES = {
			new FileChooser.ExtensionFilter("JPG", "*.jpg; *.jpeg; *.jpe; *.jfif"),
			new FileChooser.ExtensionFilter("PNG", "*.png") };
	
	private static final Projection[] PROJ_ARR = { Projection.MERCATOR, Projection.PLATE_CARREE, Projection.HOBO_DYER,
			Projection.GALL, Projection.STEREOGRAPHIC, Projection.POLAR, Projection.E_A_AZIMUTH,
			Projection.ORTHOGRAPHIC, Projection.GNOMONIC, Projection.LAMBERT_CONIC, Projection.E_D_CONIC,
			Projection.ALBERS, Projection.LEE, Projection.TETRAGRAPH, Projection.AUTHAGRAPH, Projection.SINUSOIDAL,
			Projection.MOLLWEIDE, Projection.TOBLER, Projection.AITOFF, Projection.VAN_DER_GRINTEN, Projection.ROBINSON,
			Projection.WINKEL_TRIPEL, Projection.PEIRCE_QUINCUNCIAL, Projection.GUYOU, Projection.LEMONS,
			Projection.MAGNIFIER, Projection.EXPERIMENT };
	
	
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
		final Node inputSelector = buildInputSelector(RASTER_TYPES, this::setInput);
		final Node projectionSelector = buildProjectionSelector(PROJ_ARR, Projection.MERCATOR, Procedure.NONE);
		final Node aspectSelector = buildAspectSelector(this.aspect, Procedure.NONE);
		this.updateBtn = buildUpdateButton(this::updateMap);
		this.saveMapBtn = buildSaveButton(true, "map", RASTER_TYPES,
				this::collectFinalParameters, this::calculateAndSaveMap);
		final HBox buttons = new HBox(5, updateBtn, saveMapBtn);
		buttons.setAlignment(Pos.CENTER);
		
		final VBox layout = new VBox(5,
				inputSelector, new Separator(), projectionSelector,
				new Separator(), aspectSelector, new Separator(), buttons);
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(CONT_WIDTH);
		
		this.display = new ImageView();
		this.display.setFitWidth(IMG_WIDTH);
		this.display.setFitHeight(IMG_WIDTH);
		this.display.setPreserveRatio(true);
		
		final HBox gui = new HBox(layout, this.display);
		gui.setAlignment(Pos.CENTER);
		gui.setSpacing(10);
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
			alert.showAndWait();
		} finally {
			updateBtn.setDisable(false);
			saveMapBtn.setDisable(false);
		}
	}
	
	
	private void updateMap() {
		display.setImage(map());
	}
	
	
	protected void collectFinalParameters() {
		final double ratio = getProjection().getAspectRatio(getProjection().getDefaultParameters()); //TODO: get real parameters
		this.configDialog = new MapConfigurationDialog(ratio);
		this.configDialog.showAndWait();
	}
	
	
	protected void calculateAndSaveMap(File file, ProgressBarDialog pBar) {
		pBar.setContentText("Finalizing map...");
		Image theMap = map(
				configDialog.getDims(), configDialog.getSmoothing(), pBar); //calculate
		Platform.runLater(() -> pBar.setProgress(-1));
		pBar.setContentText("Saving...");
		try {
			ImageIO.write(SwingFXUtils.fromFXImage(theMap,null), "jpg", file); //save
		} catch (IOException e) {
			final Alert alert = new Alert(AlertType.ERROR);
			alert.setHeaderText("Failure!");
			alert.setContentText("Could not access "+file.getAbsolutePath()+". It's possible that another program has it open.");
			alert.showAndWait();
		}
	}
	
	
	public Image map() {
		final double a = this.getProjection().getAspectRatio(getProjection().getDefaultParameters()); //TODO: get real parameters
		return map(new int[] { IMG_WIDTH, (int)(IMG_WIDTH/a) }, 1);
	}
	
	public Image map(int[] outputDims, int smoothing) {
		return map(outputDims,smoothing, null);
	}
	
	public Image map(int[] outDims, int smoothing,
			ProgressBarDialog pbar) {
		final Projection proj = this.getProjection();
		final double[] params = proj.getDefaultParameters(); //TODO: get real parameters
		final PixelReader ref = input.getPixelReader();
		final int[] refDims = {(int)input.getWidth(), (int)input.getHeight()};
		
		WritableImage img = new WritableImage(outDims[0], outDims[1]);
		
		for (int x = 0; x < outDims[0]; x ++) {
			for (int y = 0; y < outDims[1]; y ++) {
				int[] colors = new int[smoothing*smoothing];
				int i = 0;
				for (double dx = 0.5/smoothing; dx < 1; dx += 1.0/smoothing) {
					for (double dy = .5/smoothing; dy < 1; dy += 1.0/smoothing) {
						colors[i] = getArgb(x+dx, y+dy,
								proj,params,aspect,ref,refDims,outDims);
						i ++;
					}
				}
				img.getPixelWriter().setArgb(x, y, blend(colors));
			}
			if (pbar != null)	pbar.setProgress((double)(x+1)/outDims[0]);
		}
		
		return img;
	}
	
	
	public static int getArgb(double x, double y, Projection proj, double[] params, double[] pole,
			PixelReader ref, int[] refDims, int[] outDims) {
		final double[] coords = proj.inverse(
				2.*x/outDims[0]-1, 1-2.*y/outDims[1], params, pole);
		if (coords != null)
			return getColorAt(coords, ref, refDims);
		else
			return 0;
	}
	
	
	public static int getColorAt(double[] coords, PixelReader ref, int[] refDims) { // returns the color of any coordinate on earth
		double x = 1/2.0 + coords[1]/(2*Math.PI);
		x = (x - Math.floor(x)) * refDims[0];
		
		double y = refDims[1]/2.0 - coords[0]*refDims[1]/(Math.PI);
		if (y < 0)
			y = 0;
		else if (y >= refDims[1])
			y = refDims[1] - 1;
		
		return ref.getArgb((int)x, (int)y);
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