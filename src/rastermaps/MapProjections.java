package rastermaps;
import java.io.File;
import java.io.IOException;
import java.util.Optional;
import javax.imageio.ImageIO;

import ellipticFunctions.Jacobi;
import javafx.application.Application;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.concurrent.Task;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Dialog;
import javafx.scene.control.Label;
import javafx.scene.control.MenuButton;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Separator;
import javafx.scene.control.Slider;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.image.PixelReader;
import javafx.scene.image.WritableImage;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.scene.input.KeyEvent;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import mfc.field.Complex;
import util.Dixon;
import util.ProgressBarDialog;
import util.Robinson;
import util.WinkelTripel;

/**
 * An application to make raster oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapProjections extends Application {

	private static final int CONT_WIDTH = 300;
	private static final int IMG_WIDTH = 500;
	
	
	private static final KeyCombination ctrlO = new KeyCodeCombination(KeyCode.O, KeyCodeCombination.CONTROL_DOWN);
	private static final KeyCombination ctrlS = new KeyCodeCombination(KeyCode.S, KeyCodeCombination.CONTROL_DOWN);
	
	
	private static final String[] PROJ_ARR = { "Mercator", "Equirectangular", "Hobo-Dyer", "Gall Stereographic",
			"Stereographic", "Polar", "Azimuthal Equal-Area", "Orthographic", "Gnomonic", "Conformal Conic",
			"Equidistant Conic", "Albers", "Lee", "TetraGraph", "AuthaGraph", "Mollweide", "Hammer", "Sinusoidal", "Van der Grinten", "Robinson",
			"Winkel Tripel", "Lemons", "Pierce Quincuncial", "Guyou", "Magnifier", "Experimental" };
	private static final double[] DEFA = { 1, 2, 1.977, 4/3., 1, 1, 1, 1, 1, 2,
			2, 2, Math.sqrt(3), Math.sqrt(3), 4/Math.sqrt(3), 2, 2, 2, 1, 1.9716, Math.PI/2, 2, 1, 2, 1, 1 };
	private static final String[] DESC = { "A conformal cylindrical map", "An equidistant cylindrical map",
			"An equal-area cylindrical map", "A compromising cylindrical map", "A conformal azimuthal map",
			"An equidistant azimuthal map", "An equal-area azimuthal map",
			"Represents earth viewed from an infinite distance",
			"Every straight line on the map is a straight line on the sphere", "A conformal conic map",
			"An equidistant conic map", "An equal-area conic map", "A conformal tetrahedral map that really deserves more attention",
			"An equidistant tetrahedral map that I invented", "An almost-equal-area tetrahedral map", "An equal-area map shaped like an ellipse",
			"An equal-area map shaped like an ellipse", "An equal-area map shaped like a sinusoid","A circular compromise map",
			"A visually pleasing piecewise compromise map",
			"The compromise map used by National Geographic", 
			"BURN LIFE'S HOUSE DOWN!", "A conformal square map that uses complex math",
			"A rearranged Pierce Quincuncial map",
			"A novelty map that swells the center to disproportionate scale",
			"What happens when you apply a complex differentiable function to a stereographic projection?" };
	
	private static final String[] AXES = { "Standard", "Transverse", "Center of Mass", "Jerusalem", "Point Nemo",
			"Longest Line", "Longest Line Transverse", "Cylindrical", "Conic", "Tetrahedral", "Quincuncial", "Antipode", "Random" };
	private static final double[] DEF_LATS = { 90, 0, 29.9792, 31.7833, 48.8767, -28.5217, -46.4883, -35, -10, 47, 60 };
	private static final double[] DEF_LONS = { 0, 0, 31.1344, 35.216, 56.6067, 141.451, 16.5305, -13.6064, 65, -173, -6 };
	private static final double[] DEF_THTS = { 0, 0, -32, -35, -45, 161.5, 137, 145, -150, 138, -10 };
	
	
	private Stage stage;
	private FileChooser inputChooser, saver;
	private Text inputLabel;
	private Button changeInput;
	private ComboBox<String> projectionChooser;
	private Text projectionDesc;
	private Slider latSlider, lonSlider, thtSlider;
	private Button update, saveMap;
	private Image input;
	private ImageView output;
	
	
	
	public static final void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage primaryStage) {
		stage = primaryStage;
		stage.setTitle("Map Designer");
		
		final VBox layout = new VBox();
		layout.setSpacing(5);
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(CONT_WIDTH);
		
		Label lbl = new Label("Current input:");
		inputLabel = new Text("None");
		layout.getChildren().add(new HBox(3, lbl, inputLabel));
		
		inputChooser = new FileChooser();
		inputChooser.setInitialDirectory(new File("input"));
		inputChooser.setTitle("Choose an input map");
		inputChooser.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("JPG", "*.jpg; *.jpeg; *.jpe; *.jfif"),
				new FileChooser.ExtensionFilter("PNG", "*.png"));
		
		changeInput = new Button("Choose input...");
		changeInput.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				chooseInput();
			}
		});
		changeInput.setTooltip(new Tooltip(
				"Choose the image to determine your map's color scheme"));
		stage.addEventHandler(KeyEvent.KEY_RELEASED, new EventHandler<KeyEvent>() {	// ctrl-O opens
			public void handle(KeyEvent event) {
				if (ctrlO.match(event))	changeInput.fire();
			}
		});
		layout.getChildren().add(changeInput);
		
		layout.getChildren().add(new Separator());
		
		lbl = new Label("Projection:");
		ObservableList<String> items = FXCollections.observableArrayList(PROJ_ARR);
		projectionChooser = new ComboBox<String>(items);
		projectionChooser.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				for (int i = 0; i < PROJ_ARR.length; i ++) {
					if (PROJ_ARR[i].equals(projectionChooser.getValue())) {
						projectionDesc.setText(DESC[i]);
						break;
					}
				}
			}
		});
		projectionChooser.setValue("Mercator");
		layout.getChildren().add(new HBox(3, lbl, projectionChooser));
		
		projectionDesc = new Text(DESC[1]);
		projectionDesc.setWrappingWidth(CONT_WIDTH);
		layout.getChildren().add(projectionDesc);
		
		layout.getChildren().add(new Separator());
		
		final MenuButton defAxes = new MenuButton("Aspect Presets");
		for (String preset: AXES) {
			MenuItem m = new MenuItem(preset);
			m.setOnAction(new EventHandler<ActionEvent>() {
				public void handle(ActionEvent event) {
					setAxisByPreset(((MenuItem) event.getSource()).getText());
				}
			});
			defAxes.getItems().add(m);
		}
		defAxes.setTooltip(new Tooltip(
				"Set the aspect sliders based on a preset"));
		layout.getChildren().add(defAxes);
		
		latSlider = new Slider(-90, 90, 90);
		lonSlider = new Slider(-180,180,0);
		thtSlider = new Slider(-180,180,0);
		Tooltip aspTlTp = new Tooltip("Change the aspect of the map");
		latSlider.setTooltip(aspTlTp);
		lonSlider.setTooltip(aspTlTp);
		thtSlider.setTooltip(aspTlTp);
		
		GridPane grid = new GridPane();
		grid.addRow(0, new Text("Latitude:"), latSlider);
		grid.addRow(1, new Text("Longitude:"), lonSlider);
		grid.addRow(2, new Text("Orientation:"), thtSlider);
		GridPane.setHgrow(latSlider, Priority.ALWAYS);
		GridPane.setHgrow(lonSlider, Priority.ALWAYS);
		GridPane.setHgrow(thtSlider, Priority.ALWAYS);
		//GridPane.setFillWidth(latSlider, Boolean.valueOf(true));
		layout.getChildren().add(grid);
		
		layout.getChildren().add(new Separator());
		
		update = new Button("Update Map");
		update.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				updateMap();
			}
		});
		update.setTooltip(new Tooltip(
				"Update the current map with your parameters."));
		update.setDefaultButton(true);
		
		saver = new FileChooser();
		saver.setInitialDirectory(new File("output"));
		saver.setInitialFileName("myMap.jpg");
		saver.setTitle("Save Map");
		saver.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("JPG", "*.jpg"),
				new FileChooser.ExtensionFilter("PNG", "*.png"));
		
		saveMap = new Button("Save Map...");
		saveMap.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				startFinalizingMap();
			}
		});
		saveMap.setTooltip(new Tooltip("Save the map with current settings."));
		stage.addEventHandler(KeyEvent.KEY_RELEASED, new EventHandler<KeyEvent>() {	// ctrl-S saves
			public void handle(KeyEvent event) {
				if (ctrlS.match(event))	saveMap.fire();
			}
		});
		
		HBox box = new HBox(5, update, saveMap);
		box.setAlignment(Pos.CENTER);
		layout.getChildren().add(box);
		
		output = new ImageView();
		output.setFitWidth(IMG_WIDTH);
		output.setFitHeight(IMG_WIDTH);
		output.setPreserveRatio(true);
		
		new Thread(() -> {
			setInput("basic.jpg", "input/basic.jpg");
			update.fire();
		}).start();
		
		final HBox gui = new HBox(layout, output);
		gui.setAlignment(Pos.CENTER);
		gui.setSpacing(10);
		StackPane.setMargin(gui, new Insets(10));
		stage.setScene(new Scene(new StackPane(gui)));
		stage.show();
	}
	
	
	private void chooseInput() {
		final File f = inputChooser.showOpenDialog(stage);
		if (f != null) {
			new Thread(() -> {
				setInput(f.getName(), f.getAbsolutePath());
			}).start();
		}
	}
	
	
	private void setInput(String name, String filename) {
		changeInput.setDisable(true);
		update.setDisable(true);
		saveMap.setDisable(true);
		
		try {
			input = new Image("file:"+filename);
			inputLabel.setText(name);
		} catch (IllegalArgumentException e) {
			final Alert alert = new Alert(Alert.AlertType.ERROR);
			alert.setHeaderText("File not found!");
			alert.setContentText("Couldn't find "+filename+".");
		}
		
		changeInput.setDisable(false);
		update.setDisable(false);
		saveMap.setDisable(false);
	}
	
	
	private void setAxisByPreset(String preset) {
		if (preset.equals("Antipode")) {
			latSlider.setValue(-latSlider.getValue());
			lonSlider.setValue((lonSlider.getValue()+360)%360-180);
			thtSlider.setValue(-thtSlider.getValue());
			return;
		}
		if (preset.equals("Random")) {
			latSlider.setValue(Math.toDegrees(Math.asin(Math.random()*2-1)));
			lonSlider.setValue(Math.random()*360-180);
			thtSlider.setValue(Math.random()*360-180);
			return;
		}
		for (int i = 0; i < AXES.length; i ++) {
			if (AXES[i].equals(preset)) {
				latSlider.setValue(DEF_LATS[i]);
				lonSlider.setValue(DEF_LONS[i]);
				thtSlider.setValue(DEF_THTS[i]);
				break;
			}
		}
	}
	
	
	private void updateMap() {
		update.setDisable(true);
		new Thread(new Task<Void>() {
			protected Void call() {
				output.setImage(map());
				update.setDisable(false);
				return null;
			}
		}).start();
	}
	
	
	private void startFinalizingMap() {
		int p = 0;
		while (p < PROJ_ARR.length &&
				!PROJ_ARR[p].equals(projectionChooser.getValue()))
			p ++;
		
		Dialog<Thread> dialog = new MapConfigurationDialog(DEFA[p], this);
		Optional<Thread> mapMaker = dialog.showAndWait();
		if (mapMaker.isPresent())	mapMaker.get().start();
		else						return;
	}
	
	
	public void saveImage(Image img, ProgressBarDialog pBar) {	// call from the main thread!
		pBar.close();
		
		final File f = saver.showSaveDialog(stage);
		if (f != null) {
			new Thread(() -> {
				try {
					saveMap.setDisable(true);
					ImageIO.write(SwingFXUtils.fromFXImage(img,null), "png", f);
					saveMap.setDisable(false);
				} catch (IOException e) {}
			}).start();
		}
	}
	
	
	public Image map() {
		int p = 0;
		while (p < PROJ_ARR.length &&
				!PROJ_ARR[p].equals(projectionChooser.getValue()))
			p ++;
		return map(IMG_WIDTH, (int)(IMG_WIDTH/DEFA[p]), 1);
	}
	
	
	public Image map(int outputWidth, int outputHeight, int smoothing) {
		return map(outputWidth,outputHeight,smoothing, null);
	}
	
	
	public Image map(int outputWidth, int outputHeight, int smoothing,
			ProgressBarDialog pbar) {
		final String proj = projectionChooser.getValue();
		final double[] pole = {Math.toRadians(latSlider.getValue()),
							Math.toRadians(lonSlider.getValue()),
							Math.toRadians(thtSlider.getValue())};
		final PixelReader ref = input.getPixelReader();
		final int[] refDims = {(int)input.getWidth(), (int)input.getHeight()};
		final int[] outDims = {outputWidth, outputHeight};
		
		WritableImage img = new WritableImage(outputWidth, outputHeight);
		
		for (int x = 0; x < outputWidth; x ++) {
			for (int y = 0; y < outputHeight; y ++) {
				int[] colors = new int[smoothing*smoothing];
				int i = 0;
				for (double dx = 0; dx < 1; dx += 1.0/smoothing) {
					for (double dy = 0; dy < 1; dy += 1.0/smoothing) {
						colors[i] = getArgb(x+dx, y+dy,
								proj,pole,ref,refDims,outDims);
						i ++;
					}
				}
				img.getPixelWriter().setArgb(x, y, blend(colors));
			}
			if (pbar != null)	pbar.setProgress((double)(x+1)/outputWidth);
		}
		
		return img;
	}
	
	
	public static int getArgb(double x, double y, String proj, double[] pole,
			PixelReader ref, int[] refDims, int[] outDims) {
		double[] coords = project(x, y, proj, outDims);
		
		if (coords != null)
			return getColorAt(obliquify(pole, coords), ref, refDims);
		else
			return 0;
	}
	
	
	public static double[] project(double x, double y, String p, int[] dims) { //apply a map projection
		final double X = 2.0*x/dims[0]-1;
		final double Y = 1-2.0*y/dims[1];
		
		if (p.equals("Pierce Quincuncial"))
			return quincuncial(X, Y);
		else if (p.equals("Equirectangular"))
			return equirectangular(X, Y);
		else if (p.equals("Mercator"))
			return mercator(X, Y);
		else if (p.equals("Polar"))
			return polar(X, Y);
		else if (p.equals("Gall Stereographic"))
			return gall(X, Y);
		else if (p.equals("Sinusoidal"))
			return sinusoidal(X, Y);
		else if (p.equals("Stereographic"))
			return stereographic(X, Y);
		else if (p.equals("Gnomonic"))
			return gnomonic(X, Y);
		else if (p.equals("Orthographic"))
			return orthographic(X, Y);
		else if (p.equals("Hobo-Dyer"))
			return eaCylindrical(X, Y);
		else if (p.equals("Conformal Conic"))
			return lambert(X, Y);
		else if (p.equals("Lemons"))
			return lemons(X, Y);
		else if (p.equals("Azimuthal Equal-Area"))
			return eaAzimuth(X, Y);
		else if (p.equals("Guyou"))
			return quinshift(X, Y);
		else if (p.equals("Mollweide"))
			return mollweide(X, Y);
		else if (p.equals("Winkel Tripel"))
			return winkel_tripel(X, Y);
		else if (p.equals("Van der Grinten"))
			return grinten(X, Y);
		else if (p.equals("Magnifier"))
			return magnus(X, Y);
		else if (p.equals("Hammer"))
			return hammer(X, Y);
		else if (p.equals("AuthaGraph"))
			return authagraph(X, Y);
		else if (p.equals("TetraGraph"))
			return tetragraph(X, Y);
		else if (p.equals("Experimental"))
			return experiment(X, Y);
		else if (p.equals("Equidistant Conic"))
			return edConic(X, Y);
		else if (p.equals("Albers"))
			return albers(X, Y);
		else if (p.equals("Robinson"))
			return robinson(X, Y);
		else if (p.equals("Lee"))
			return lee(X, Y);
		else
			throw new IllegalArgumentException(p);
	}
	
	
	public static int getColorAt(double[] coords,
			PixelReader ref, int[] refDims) { // returns the color of any coordinate on earth
		double x = 1/2.0 + coords[1]/(2*Math.PI);
		x = (x - Math.floor(x)) * refDims[0];
		
		double y = refDims[1]/2.0 - coords[0]*refDims[1]/(Math.PI);
		if (y < 0)
			y = 0;
		else if (y >= refDims[1])
			y = refDims[1] - 1;
		
		return ref.getArgb((int)x, (int)y);
	}
	
	
	public static final double[] obliquify(double[] pole, double[] coords) {
		final double lat0 = pole[0];
		final double lon0 = pole[1];
		final double tht0 = pole[2];
		double lat1 = coords[0];
		double lon1 = coords[1];
		lon1 += tht0;
		double latf = Math.asin(Math.sin(lat0)*Math.sin(lat1) - Math.cos(lat0)*Math.cos(lon1)*Math.cos(lat1));
		double lonf;
		double innerFunc = Math.sin(lat1)/Math.cos(lat0)/Math.cos(latf) - Math.tan(lat0)*Math.tan(latf);
		if (lat0 == Math.PI / 2) // accounts for special case when lat0 = pi/2
			lonf = lon1+lon0;
		else if (lat0 == -Math.PI / 2) // accounts for special case when lat0 = -pi/2
			lonf = -lon1+lon0 + Math.PI;
		else if (Math.abs(innerFunc) > 1) { // accounts for special case when cos(lat1) = --> 0
			if ((lon1 == 0 && lat1 < -lat0) || (lon1 != 0 && lat1 < lat0))
				lonf = lon0 + Math.PI;
			else
				lonf = lon0;
		}
		else if (Math.sin(lon1) > 0)
			lonf = lon0 +
					Math.acos(innerFunc);
		else
			lonf = lon0 -
					Math.acos(innerFunc);
		
		double thtf = 0;
		if (pole.length >= 3)
			thtf += pole[2];
		
		double[] output = {latf, lonf, thtf};
		return output;
	}
	
	
	private static double[] quincuncial(double x, double y) { // a tessalatable square map
		Complex u = new Complex(1.854 * (x+1), 1.854 * y); // 1.854 is approx K(sqrt(1/2)
		Complex k = new Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus
		Complex ans = Jacobi.cn(u, k);
		double p = 2 * Math.atan(ans.abs());
		double theta = ans.arg() - Math.PI/2;
		double lambda = Math.PI/2 - p;
		return new double[] {lambda, theta};
	}
	
	private static double[] experiment(double x, double y) { // just some random complex plane stuff
		Complex z = new Complex(x*3, y*3);
		Complex ans = z.sin();
		double p = 2 * Math.atan(ans.abs());
		double theta = ans.arg();
		double lambda = Math.PI/2 - p;
		return new double[] {lambda, theta};
	}
	
	private static double[] equirectangular(double x, double y) { // a linear scale
		return new double[] {y*Math.PI/2, x*Math.PI};
	}
	
	private static double[] mercator(double x, double y) { // a popular shape-preserving map
		double phi = Math.atan(Math.sinh(y*Math.PI));
		return new double[] {phi, x*Math.PI};
	}
	
	private static double[] polar(double x, double y) { // the projection used on the UN flag
		double phi = Math.PI/2 - Math.PI * Math.hypot(x, y);
		if (phi > -Math.PI/2)
			return new double[] {phi, Math.atan2(y, x) + Math.PI/2};
		else
			return null;
	}
	
	private static double[] gall(double x, double y) { // a compromise map, similar to mercator
		return new double[] {2*Math.atan(y), x*Math.PI};
	}
	
	private static double[] sinusoidal(double x, double y) { // a map shaped like a sinusoid
		return new double[] {y*Math.PI/2, x*Math.PI / Math.cos(y*Math.PI/2)};
	}
	
	private static double[] stereographic(double x, double y) { // a shape-preserving infinite map
		return new double[] {Math.PI/2 - 2*Math.atan(2*Math.hypot(x, y)),
				Math.atan2(y, x) + Math.PI/2};
	}
	
	private static double[] gnomonic(double x, double y) { // map where straight lines are straight
		return new double[] {Math.PI/2 - Math.atan(2*Math.hypot(x, y)),
				Math.atan2(y, x) + Math.PI/2};
	}
	
	private static double[] orthographic(double x, double y) { // a map that mimics the view from space
		double R = Math.hypot(x, y);
		if (R <= 1)
			return new double[] {Math.acos(R), Math.atan2(y, x) + Math.PI/2};
		else
			return null;
	}
	
	private static double[] eaCylindrical(double x, double y) { // an equal-area cylindrical map
		return new double[] {Math.asin(y), x*Math.PI};
	}
	
	private static double[] lambert(double x, double y) { // a conical projection
		y = (y-1)/2;
		return new double[] {
				Math.PI/2 - 2*Math.atan(Math.pow(1.5*Math.hypot(x, y), 2)),
				2*(Math.atan2(y, x) + Math.PI/2)};
	}
	
	private static double[] lemons(double x, double y) { // a simple map that is shaped like lemons
		x = x+2;
		final double lemWdt = 1/6.0;
		
		if (Math.abs(x % lemWdt - lemWdt / 2.0) <= Math.cos(y*Math.PI/2) * lemWdt/2.0) // if it is in
			return new double[] {y*Math.PI/2,	// a sine curve
					Math.PI * (x%lemWdt - lemWdt/2.0) / (Math.cos(y*Math.PI/2))
							+ (int)(x/lemWdt) * Math.PI/6};
		else
			return null;
	}
	
	private static double[] eaAzimuth(double x, double y) { // the lambert azimuthal equal area projection
		double R = Math.hypot(x, y);
		if (R <= 1)
			return new double[] {Math.asin(1-2*R*R), Math.atan2(y,x)+Math.PI/2};
		else
			return null;
	}
	
	private static double[] quinshift(double x, double y) { // a tessalatable rectangle map
		Complex u = new Complex(1.8558*(x - y/2 - 0.5), 1.8558*(x + y/2 + 0.5)); // don't ask me where 3.7116 comes from
		Complex k = new Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus
		Complex ans = Jacobi.cn(u, k);
		double p = 2 * Math.atan(ans.abs());
		double theta = ans.arg();
		double lambda = Math.PI/2 - p;
		return new double[] {lambda, theta};
	}
	
	private static double[] lee(double x, double y) { //a tessalatable rectangle map
		final double x1, y1;
		if (y > x+1) {
			x1 = -x - 4/3.;
			y1 = y - 1;
		}
		else if (y < -x-1) {
			x1 = x + 2/3.;
			y1 = -y - 1;
		}
		else if (y > 1-x) {
			x1 = -4/3. + x;
			y1 = 1 - y;
		}
		else if (y < x-1) {
			x1 = 2/3. - x;
			y1 = 1 + y;
		}
		else if (x < 0) {
			x1 = x + 2/3.;
			y1 = -y - 1;
		}
		else {
			x1 = 2/3. - x;
			y1 = y + 1;
		}
		Complex w = new Complex(x1, y1/Math.sqrt(3)).times(2.655);
		Complex ans = Dixon.leeFunc(w).times(Math.pow(2, -5/6.));
		double p = 2 * Math.atan(ans.abs());
		double theta = ans.arg();
		double lambda = p - Math.PI/2;
		return new double[] {lambda, theta};
	}
	
	private static double[] mollweide(double x, double y) {
		double tht = Math.asin(y);
		return new double[] {
				Math.asin((2*tht + Math.sin(2*tht)) / Math.PI),
				Math.PI * x / Math.cos(tht)};
	}
	
	private static double[] winkel_tripel(double x, double y) {
		final double tolerance = 10e-2;
		
		double X = x * (1 + Math.PI/2);
		double Y = y * Math.PI/2;
		if ((x<0) == (y<0))	Y = -Y; // this makes it converge faster for some reason
		double phi = Y;
		double lam = X; // I used equirectangular for my initial guess
		double error = Math.PI;
		
		for (int i = 0; i < 6 && error > tolerance; i++) {
			final double f1 = WinkelTripel.f1pX(phi, lam) - X;
			final double f2 = WinkelTripel.f2pY(phi, lam) - Y;
			final double df1dP = WinkelTripel.df2dphi(phi, lam);
			final double df1dL = WinkelTripel.df1dlam(phi, lam);
			final double df2dP = WinkelTripel.df2dphi(phi, lam);
			final double df2dL = WinkelTripel.df2dlam(phi, lam);
			
			phi -= (f1*df2dL - f2*df1dL) / (df1dP*df2dL - df2dP*df1dL);
			lam -= (f2*df1dP - f1*df2dP) / (df1dP*df2dL - df2dP*df1dL);
			
			error = Math.hypot(f1, f2);
		}
		if ((x<0) == (y<0))	phi = -phi;
		
		if (error > tolerance) // if it aborted due to timeout
			return null;
		else // if it aborted due to convergence
			return new double[] {phi, lam};
	}
	
	private static double[] grinten(double x, double y) {
		if (y == 0) // special case 1: equator
			return new double[] {0, x*Math.PI};
		
		if (x == 0) // special case 3: meridian
			return new double[] {Math.PI/2 * Math.sin(2*Math.atan(y)), 0};
		
		double c1 = -Math.abs(y) * (1 + x*x + y*y);
		double c2 = c1 - 2*y*y + x*x;
		double c3 = -2 * c1 + 1 + 2*y*y + Math.pow(x*x + y*y, 2);
		double d = y*y / c3 + 1 / 27.0 * (2*Math.pow(c2 / c3, 3) - 9*c1*c2 / (c3*c3));
		double a1 = 1 / c3*(c1 - c2*c2 / (3*c3));
		double m1 = 2 * Math.sqrt(-a1 / 3);
		double tht1 = Math.acos(3*d / (a1 * m1)) / 3;
		
		return new double[] {
				Math.signum(y) * Math.PI * (-m1 * Math.cos(tht1 + Math.PI/3) - c2 / (3*c3)),
				Math.PI*(x*x + y*y - 1 + Math.sqrt(1 + 2*(x*x - y*y) + Math.pow(x*x + y*y, 2)))
						/ (2*x)};
	}
	
	private static double[] magnus(double x, double y) { // a novelty map that magnifies the center profusely
		double R = Math.hypot(x, y);
		if (R <= 1)
			return new double[] {
					Math.PI/2 * (1 - R*.2 - R*R*R*1.8),
					Math.atan2(y, x) + Math.PI/2};
		else
			return null;
	}
	
	private static double[] hammer(double x, double y) { // similar to Mollweide, but moves distortion from the poles to the edges
		final double X = x * Math.sqrt(8);
		final double Y = y * Math.sqrt(2);
		final double z = Math.sqrt(1 - Math.pow(X/4, 2) - Math.pow(Y/2, 2));
		return new double[] {
				Math.asin(z * Y),
				2*Math.atan(0.5*z*X / (2*z*z - 1))};
	}
	
	private static double[] authagraph(double x, double y) { // a modern Japanese almost-equal-area map
		final double[] faceCenter = new double[3];
		final double rot, localX, localY;
		if (y-1 < 4*x && y-1 < -4*x) {
			faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
			faceCenter[1] = 0;
			rot = 0;
			localX = 4/Math.sqrt(3)*x;
			localY = y+1/3.0;
		}
		else if (y-1 < -4*(x+1)) {
			faceCenter[0] = -Math.PI/2;
			faceCenter[1] = Math.PI;
			rot = 0;
			localX = 4/Math.sqrt(3)*(x+1);
			localY = y+1/3.0;
		}
		else if (y-1 < 4*(x-1)) {
			faceCenter[0] = -Math.PI/2;
			faceCenter[1] = Math.PI;
			rot = 0;
			localX = 4/Math.sqrt(3)*(x-1);
			localY = y+1/3.0;
		}
		else if (x < 0) {
			faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
			faceCenter[1] = 4*Math.PI/3;
			rot = Math.PI/3;
			localX = 4/Math.sqrt(3)*(x+0.5);
			localY = y-1/3.0;
		}
		else {
			faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
			faceCenter[1] = 2*Math.PI/3;
			rot = -Math.PI/3;
			localX = 4/Math.sqrt(3)*(x-0.5);
			localY = y-1/3.0;
		}
		faceCenter[2] = 0;
		
		final double t = Math.atan2(localY, localX) + rot;
		final double t0 = Math.floor((t+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
		final double dt = t-t0;
		final double z = 2.49*Math.hypot(localX, localY)*Math.cos(dt);
		
		final double g = 0.03575*z*z*z + 0.0219*z*z + 0.4441*z;
		
		double[] triCoords = {
				Math.PI/2 - Math.atan(Math.tan(g)/Math.cos(dt)),
				Math.PI/2 + t0 + dt};
		return obliquify(faceCenter, triCoords);
	}
	
	private static double[] tetragraph(double x, double y) { // a tetrahedral compromise
		final double[] faceCenter = new double[3];
		final double rot, localX, localY;
		if (y < x-1) {
			faceCenter[0] = -Math.PI/2;
			faceCenter[1] = 0;
			rot = -Math.PI/2;
			localX = Math.sqrt(3)*(x-2/3.0);
			localY = y+1;
		}
		else if (y < -x-1) {
			faceCenter[0] = -Math.PI/2;
			faceCenter[1] = 0;
			rot = Math.PI/2;
			localX = Math.sqrt(3)*(x+2/3.0);
			localY = y+1;
		}
		else if (y > -x+1) {
			faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
			faceCenter[1] = Math.PI;
			rot = -Math.PI/2;
			localX = Math.sqrt(3)*(x-2/3.0);
			localY = y-1;
		}
		else if (y > x+1) {
			faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
			faceCenter[1] = Math.PI;
			rot = Math.PI/2;
			localX = Math.sqrt(3)*(x+2/3.0);
			localY = y-1;
		}
		else if (x < 0) {
			faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
			faceCenter[1] = -Math.PI/3;
			rot = Math.PI/6;
			localX = Math.sqrt(3)*(x+1/3.0);
			localY = y;
		}
		else {
			faceCenter[0] = Math.PI/2-Math.asin(Math.sqrt(8)/3);
			faceCenter[1] = Math.PI/3;
			rot = -Math.PI/6;
			localX = Math.sqrt(3)*(x-1/3.0);
			localY = y;
		}
		faceCenter[2] = 0;
		
		final double t = Math.atan2(localY, localX) + rot;
		final double t0 = Math.floor((t+Math.PI/2)/(2*Math.PI/3)+0.5)*(2*Math.PI/3) - Math.PI/2;
		final double dt = t-t0;
		
		double[] triCoords = {
				Math.PI/2 - Math.atan(Math.tan(1.654*Math.hypot(localX, localY)*Math.cos(dt))/Math.cos(dt)),
				Math.PI/2 + t0 + dt};
		return obliquify(faceCenter, triCoords);
	}
	
	
	private static double[] edConic(double x, double y) {
		y = (y-1)/2;
		final double r = Math.hypot(x, y);
		if (r < 0.2 || r > 1)	return null;
		return new double[] {
				3*Math.PI/4 - 5*Math.PI/4*r,
				2*Math.atan2(y, x)+Math.PI};
	}
	
	
	private static double[] albers(double x, double y) {
		y = (y-1)/2;
		final double r = Math.hypot(x, y);
		if (r < Math.sqrt(1/11.) || r > 1)	return null;
		return new double[] {
				Math.asin(1.2-2.2*Math.pow(r, 2)),
				2*Math.atan2(y, x)};
	}
	
	
	private static double[] robinson(double x, double y) {
		return new double[] {
				Robinson.latFromPdfe(Math.abs(y))*Math.signum(y),
				x/Robinson.plenFromPdfe(Math.abs(y))*Math.PI };
	}
	
	
	
	public static final int blend(int[] colors) {
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
	
	
	public static final double asinh(double x) {
		return Math.log(x+Math.sqrt(1+x*x));
	}
}