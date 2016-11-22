import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;

import ellipticFunctions.Jacobi;
import javafx.application.Application;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Label;
import javafx.scene.control.MenuButton;
import javafx.scene.control.MenuItem;
import javafx.scene.control.Separator;
import javafx.scene.control.Slider;
import javafx.scene.control.Tooltip;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.text.Text;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import mfc.field.Complex;

/**
 * 
 */

/**
 * An application to make raster oblique aspects of map projections
 * @author Justin Kunimune
 */
public class MapProjections extends Application {

	private static final int CONT_WIDTH = 300;
	private static final int IMG_WIDTH = 500;
	
	
	private static final String[] PROJ_ARR = { "Equirectangular", "Mercator", "Gall Stereographic",
			"Cylindrical Equal-Area", "Polar", "Stereographic", "Azimuthal Equal-Area", "Orthographic", "Gnomonic",
			"Lambert Conical", "Winkel Tripel", "Van der Grinten", "Mollweide", "Hammer", "Sinusoidal", "Lemons",
			"Pierce Quincuncial", "Guyou Hemisphere-in-a-Square", "AuthaGraph", "Magnifier", "Experimental" };
	private static final double[] DEFA = { 2, 1, 4/3.0, 2, 1, 1, 1, 1, 1, 2, Math.PI/2, 1, 2, 2,
			2, 2, 1, 2, 4.0/Math.sqrt(3), 1, 1 };
	private static final String[] DESC = { "An equidistant cylindrical map", "A conformal cylindrical map",
			"A compromising cylindrical map", "An equal-area cylindrical map", "An equidistant azimuthal map",
			"A conformal azimuthal map", "An equal-area azimuthal map",
			"Represents earth viewed from an infinite distance",
			"Every straight line on the map is a straight line on the sphere", "A conformal conical map",
			"The compromise map used by National Geographic (caution: very slow)", "A circular compromise map",
			"An equal-area map shaped like an ellipse", "An equal-area map shaped like an elipse",
			"An equal-area map shaped like a sinusoid", "BURN LIFE'S HOUSE DOWN!",
			"A conformal square map that uses complex math",
			"A reorganized version of Pierce Quincuncial and actually the best map ever",
			"An almost-equal-area map based on a tetrahedron.",
			"A novelty map that swells the center to disproportionate scale",
			"What happens when you apply a complex differentiable function to a stereographic projection?" };
	
	private static final String[] AXES = { "Standard", "Transverse", "Center of Mass", "Jerusalem", "Point Nemo",
			"Longest Line", "Longest Line Transverse", "Cylindrical", "Conical", "Quincuncial", "Antipode", "Random" };
	private static final double[] DEF_LATS = { 90, 0, 29.9792, 31.7833, 48.8767, -28.5217, -46.4883, -35, -10, 60 };
	private static final double[] DEF_LONS = { 0, 0, 31.1344, 35.216, 56.6067, 141.451, 16.5305, -13.6064, 65, -6 };
	private static final double[] DEF_THTS = { 0, 0, -32, -35, -45, 161.5, 137, 145, -150, -10 };
	
	
	private FileChooser inputChooser, saver;
	private Text inputLabel;
	private ComboBox<String> projectionChooser;
	private Text projectionDesc;
	private Slider latSlider, lonSlider, thtSlider;
	private Button update;
	private Image input;
	private int outputWidth, outputHeight;
	private ImageView output;
	
	
	
	public static void main(String[] args) {
		launch(args);
	}
	
	
	@Override
	public void start(Stage stage) {
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
				new FileChooser.ExtensionFilter("All Images", "*.*"),
				new FileChooser.ExtensionFilter("JPG", "*.jpg"),
				new FileChooser.ExtensionFilter("PNG", "*.png"));
		
		final Button changeInput = new Button("Choose input...");
		changeInput.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				final File f = inputChooser.showOpenDialog(stage);
				if (f == null)	return;
				input = new Image("file:"+f.getAbsolutePath());
				inputLabel.setText(f.getName());
				update.setDisable(false);
			}
		});
		changeInput.setTooltip(new Tooltip(
				"Choose the image to determine your map's color scheme"));
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
		projectionChooser.setPrefWidth(200);
		projectionChooser.setValue(PROJ_ARR[1]);
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
		layout.getChildren().add(grid);
		
		layout.getChildren().add(new Separator());
		
		update = new Button("Update Map");
		update.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				output.setImage(map(projectionChooser.getValue(),
						latSlider.getValue(),
						lonSlider.getValue(),
						thtSlider.getValue()));
			}
		});
		update.setTooltip(new Tooltip(
				"Update the current map with your parameters."));
		update.setDefaultButton(true);
		layout.getChildren().add(update);
		
		saver = new FileChooser();
		saver.setInitialDirectory(new File("output"));
		saver.setInitialFileName("myMap.png");
		saver.setTitle("Save Map");
		saver.getExtensionFilters().addAll(
				new FileChooser.ExtensionFilter("PNG", "*.png"));
		
		final Button saveMap = new Button("Save Map...");
		saveMap.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				final File f = saver.showSaveDialog(stage);
				if (f == null)	return;
				try {
					ImageIO.write(
							SwingFXUtils.fromFXImage(output.getImage(),null),
							"png", f);
				} catch (IOException e) {}
			}
		});
		saveMap.setTooltip(new Tooltip("Save your new custom map!"));
		layout.getChildren().add(saveMap);
		
		output = new ImageView();
		output.setFitWidth(IMG_WIDTH);
		output.setFitHeight(IMG_WIDTH);
		output.setPreserveRatio(true);
		
		try {
			input = new Image("file:input/basic.jpg");
			inputLabel.setText("basic.jpg");
			update.getOnAction().handle(null);
		} catch (IllegalArgumentException e) {
			update.setDisable(true);
		}
		
		final HBox gui = new HBox(layout, output);
		gui.setAlignment(Pos.CENTER);
		gui.setSpacing(10);
		StackPane.setMargin(gui, new Insets(10));
		stage.setScene(new Scene(new StackPane(gui)));
		stage.show();
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
	
	
	public Image map(String projName,
			double latD, double lonD, double thtD) {
		int p = 0;
		for (int i = 0; i < PROJ_ARR.length; i ++)
			if (PROJ_ARR[i].equals(projName))
				p = i;
		
		outputWidth = (int)(600*Math.sqrt(DEFA[p]));
		outputHeight = (int)(600/Math.sqrt(DEFA[p]));
		
		WritableImage img = new WritableImage(outputWidth, outputHeight);
		
		for (int x = 0; x < outputWidth; x ++)
			for (int y = 0; y < outputHeight; y ++)
				img.getPixelWriter().setArgb(x, y, getArgb(x, y));
		return img;
	}
	
	
	public int getArgb(int x, int y) {
		final double lat0 = Math.toRadians(latSlider.getValue());
		final double lon0 = Math.toRadians(lonSlider.getValue());
		final double tht0 = Math.toRadians(thtSlider.getValue());
		final double pole[] = {lat0, lon0, tht0};
		final double width = input.getWidth();
		final double height = input.getHeight();
		final int[] refDims = {(int)width, (int)height};
		final String p = projectionChooser.getValue();
		final double X = 2.0*x/outputWidth-1;
		final double Y = 1-2.0*y/outputHeight;
		
		double[] coords;
		if (p.equals("Pierce Quincuncial"))
			coords = quincuncial(X, Y);
		else if (p.equals("Equirectangular"))
			coords = equirectangular(X, Y);
		else if (p.equals("Mercator"))
			coords = mercator(X, Y);
		else if (p.equals("Polar"))
			coords = polar(X, Y);
		else if (p.equals("Gall Stereographic"))
			coords = gall(X, Y);
		else if (p.equals("Sinusoidal"))
			coords = sinusoidal(X, Y);
		else if (p.equals("Stereographic"))
			coords = stereographic(X, Y);
		else if (p.equals("Gnomonic"))
			coords = gnomonic(X, Y);
		else if (p.equals("Orthographic"))
			coords = orthographic(X, Y);
		else if (p.equals("Cylindrical Equal-Area"))
			coords = eaCylindrical(X, Y);
		else if (p.equals("Lambert Conical"))
			coords = lambert(X, Y);
		else if (p.equals("Lemons"))
			coords = lemons(X, Y);
		else if (p.equals("Azimuthal Equal-Area"))
			coords = eaAzimuth(X, Y);
		else if (p.equals("Guyou Hemisphere-in-a-Square"))
			coords = quinshift(X, Y);
		else if (p.equals("Mollweide"))
			coords = mollweide(X, Y);
		else if (p.equals("Winkel Tripel"))
			coords = winkel_tripel(X, Y);
		else if (p.equals("Van der Grinten"))
			coords = grinten(X, Y);
		else if (p.equals("Magnifier"))
			coords = magnus(X, Y);
		else if (p.equals("Hammer"))
			coords = hammer(X, Y);
		else if (p.equals("AuthaGraph"))
			coords = authagraph(X, Y);
		else if (p.equals("Experimental"))
			coords = experiment(X, Y);
		else
			throw new IllegalArgumentException(p);
		
		if (coords != null)
			return getColor(pole, coords, refDims, input);
		else
			return 0;
	}
	
	
	public static int getColor(final double[] pole, double[] coords,
			int[] refDims, Image input) { // returns the color of any coordinate on earth
		final double[] convCoords = obliquify(pole, coords);
		
		double x = 1/2.0 + convCoords[1]/(2*Math.PI);
		x = (x - Math.floor(x)) * refDims[0];
		
		double y = refDims[1]/2.0 - convCoords[0]*refDims[1]/(Math.PI);
		if (y < 0)
			y = 0;
		else if (y >= refDims[1])
			y = refDims[1] - 1;
		
		return input.getPixelReader().getArgb((int)x, (int)y);
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
		
		double P = Math.sin(lat0)*Math.cos(latf)-Math.cos(lat0)*Math.sin(latf)*Math.cos(lonf-lon0);
		double thtf = Math.acos(P/Math.cos(lat1));
		thtf = 0;
		if (coords.length >= 3)
			thtf += coords[2];
		
		double[] output = {latf, lonf, thtf};
		return output;
	}
	
	
	private static double[] quincuncial(double x, double y) { // a tessalatable square map
		Complex u = new Complex(1.854 * (x+1), 1.854 * y); // 1.854 is approx K(sqrt(1/2)
		Complex k = new Complex(Math.sqrt(0.5)); // the rest comes from some fancy complex calculus
		Complex ans = Jacobi.cn(u, k);
		double p = 2 * Math.atan(ans.abs());
		double theta = Math.atan2(ans.getIm(), ans.getRe()) - Math.PI/2;
		double lambda = Math.PI/2 - p;
		return new double[] {lambda, theta};
	}
	
	private static double[] experiment(double x, double y) { // just some random complex plane stuff
		Complex z = new Complex(x*3, y*3);
		Complex ans = z;
		double p = 2 * Math.atan(ans.abs());
		double theta = Math.atan2(ans.getIm(), ans.getRe()) + Math.PI/2;
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
		
		if (Math.abs(x % lemWdt - lemWdt / 2.0) < Math.cos(y*Math.PI/2) * lemWdt/2.0) // if it is in
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
		double theta = Math.atan2(ans.getIm(), ans.getRe());
		double lambda = Math.PI/2 - p;
		return new double[] {lambda, theta};
	}
	
	private static double[] mollweide(double x, double y) {
		double tht = Math.asin(y);
		return new double[] {
				Math.asin((2*tht + Math.sin(2*tht)) / Math.PI),
				Math.PI * x / Math.cos(tht)};
	}
	
	private static double[] winkel_tripel(double x, double y) {
		final double tolerance = 0.001;
		
		double phi = y * Math.PI/2;
		double lam = x * Math.PI; // I used equirectangular for my initial guess
		double xf = x * (2 + Math.PI);
		double yf = y * Math.PI;
		double error = Math.PI;
		
		for (int i = 0; i < 100 && error > tolerance; i++) {
			final double X = WinkelTripel.X(phi, lam);
			final double Y = WinkelTripel.Y(phi, lam);
			final double dXdP = WinkelTripel.dXdphi(phi, lam);
			final double dYdP = WinkelTripel.dYdphi(phi, lam);
			final double dXdL = WinkelTripel.dXdlam(phi, lam);
			final double dYdL = WinkelTripel.dYdlam(phi, lam);
			
			phi -= (dYdL*(X - xf) - dXdL*(Y - yf)) / (dXdP*dYdL - dXdL*dYdP);
			lam -= (dXdP*(Y - yf) - dYdP*(X - xf)) / (dXdP*dYdL - dXdL*dYdP);
			
			error = Math.hypot(X - xf, Y - yf);
		}
		if (error >= tolerance) // if it aborted due to timeout
			return null;
		else // if it aborted due to convergence
			return new double[] {phi, lam + Math.PI};
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
		final double localX, localY;
		if (y-1 < 4*x && y-1 < -4*x) {
			//faceCenter[0] = Math.asin(Math.sqrt(8)/3)-Math.PI/2;
			faceCenter[0] = 0.7;
			faceCenter[1] = 0;
			faceCenter[2] = 0;
			localX = 2*x;
			localY = y-1/3.0;
		}
		else if (y-1 < -4*(x+1)) {
			faceCenter[0] = Math.PI/2;
			faceCenter[1] = 0;
			faceCenter[2] = 0;
			localX = 2*(x+1);
			localY = y-1/3.0;
		}
		else if (y-1 < 4*(x-1)) {
			faceCenter[0] = Math.PI/2;
			faceCenter[1] = 0;
			faceCenter[2] = 0;
			localX = 2*(x-1);
			localY = y-1/3.0;
		}
		else if (x < 0) {
			faceCenter[0] = Math.asin(Math.sqrt(8)/3)-Math.PI/2;
			faceCenter[1] = 4*Math.PI/3;
			faceCenter[2] = Math.PI/3;
			localX = 2*(x+0.5);
			localY = y+1/3.0;
		}
		else {
			faceCenter[0] = Math.asin(Math.sqrt(8)/3)-Math.PI/2;
			faceCenter[1] =2*Math.PI/3;
			faceCenter[2] = -Math.PI/3;
			localX = 2*(x-0.5);
			localY = y+1/3.0;
		}
		
		double[] triCoords = {
				Math.PI/2 - Math.atan(2.35*Math.hypot(localX, localY)),
				Math.atan2(localY, localX) + Math.PI/2};
		return obliquify(faceCenter, triCoords);
	}
}