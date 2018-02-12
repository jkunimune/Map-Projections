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
import java.util.LinkedList;
import java.util.List;
import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import image.SVGMap;
import image.SavableImage;
import image.SVGMap.Command;
import image.SVGMap.Path;
import javafx.concurrent.Task;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Separator;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import maps.Projection;
import utils.Math2;

/**
 * An application to make vector oblique aspects of map projections
 * 
 * @author Justin Kunimune
 */
public class MapDesignerVector extends MapApplication {

	public static final void main(String[] args) {
		launch(args);
	}
	
	
	
	private static final int DEF_MAX_VTX = 5000;
	private static final int FAST_MAX_VTX = 2000;
	
	private static final FileChooser.ExtensionFilter[] VECTOR_TYPES = {
			new FileChooser.ExtensionFilter("SVG", "*.svg") };
	
	private Region aspectSelector;
	private double[] aspect;
	private SVGMap input;
	private StackPane viewer;
	
	
	
	public MapDesignerVector() {
		super("Map Designer");
	}
	
	
	
	@Override
	public void start(Stage root) {
		super.start(root);
		new Thread(setInputTask(new File("input/Basic.svg"))).start(); //this automatically updates the map
	}
	
	
	@Override
	protected Region makeWidgets() {
		this.aspect = new double[3];
		final Region inputSelector = buildInputSelector(VECTOR_TYPES,
				VECTOR_TYPES[0], this::setInputTask);
		final Region projectionSelector = buildProjectionSelector(this::updateAspect);
		this.aspectSelector = buildAspectSelector(this.aspect, this::updateMap);
		final Region parameterSelector = buildParameterSelector(this::updateMap);
		final Region saveBtn = buildSaveButton(true, "map", VECTOR_TYPES,
				VECTOR_TYPES[0], ()->true, this::calculateTaskForSaving);
		
		final VBox layout = new VBox(V_SPACE,
				inputSelector, new Separator(), projectionSelector,
				new Separator(), aspectSelector, parameterSelector,
				new Separator(), saveBtn);
		
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(GUI_WIDTH);
		
		viewer = new StackPane();
		viewer.setMinWidth(IMG_SIZE);
		viewer.setMinHeight(IMG_SIZE);
		
		final HBox gui = new HBox(MARGIN, layout, viewer);
		gui.setAlignment(Pos.CENTER);
		StackPane.setMargin(gui, new Insets(MARGIN));
		
		return gui;
	}
	
	
	private Task<Void> setInputTask(File file) {
		return new Task<Void>() {
			protected Void call() throws IOException, SAXException, ParserConfigurationException {
				input = new SVGMap(file);
				return null;
			}
			
			protected void failed() {
				if (getException() instanceof IOException)
					showError("File not found!",
							"We couldn't find "+file.getAbsolutePath()+".");
				else if (getException() instanceof SAXException)
					showError("Unreadable file!",
							"We couldn't read "+file.getAbsolutePath()+". It may be corrupt or an unreadable format.");
				else if (getException() instanceof ParserConfigurationException)
					showError("Parser Configuration Error!",
							"My parser configured incorrectly. I blame you.");
				else
					showError("Unexpected error!", getException().getMessage());
			}
			
			protected void succeeded() {
				updateMap();
			}
		};
		
	}
	
	
	private void updateAspect() { //check the visibility of the aspect sliders
		aspectSelector.setVisible(this.getProjection().hasAspect());
		updateMap();
	}
	
	
	private void updateMap() { //execute a new calculation Task immediately
		new Thread(calculateTaskForUpdate()).start();
	}
	
	
	private Task<SavableImage> calculateTaskForUpdate() {
		int maxVtx = this.getParamsChanging() ? FAST_MAX_VTX : DEF_MAX_VTX;
		int step = input.length()/maxVtx+1;
		return calculateTask(step, true);
	}
	
	private Task<SavableImage> calculateTaskForSaving() {
		return calculateTask(1, false);
	}
	
	private Task<SavableImage> calculateTask(int step, boolean render) {
		loadParameters();
		final Projection proj = this.getProjection();
		final double[] aspect = this.aspect.clone();
		
		return new Task<SavableImage>() {
			private Canvas rendered;
			
			protected SavableImage call() {
				updateProgress(-1, 1);
				updateMessage("Generating map\u2026");
				
				List<Path> theMap = new LinkedList<Path>();
				int i = 0;
				for (Path path0: input) {
					updateProgress(i, input.numCurves());
					if (path0.size() <= step) 	continue; //don't bother drawing singular points
					Path path1 = new Path();
					int counter = 0;
					for (Command cmdS: path0) {
						counter --; //skip the requisite number of 'L's
						if (counter > 0 && cmdS.type != 'M' && cmdS.type != 'Z') 	continue;
						counter = step;
						
						Command cmdP = new Command(cmdS.type, new double[cmdS.args.length]);
						for (int k = 0; k < cmdS.args.length; k += 2) {
							double[] coords = proj.project(cmdS.args[k+1], cmdS.args[k], aspect);
							cmdP.args[k] =
									Math.max(Math.min(coords[0], proj.getWidth()), -proj.getWidth());
							cmdP.args[k+1] =
									Math.max(Math.min(coords[1], proj.getHeight()), -proj.getHeight());
							if (Double.isNaN(cmdP.args[k]) || Double.isNaN(cmdP.args[k+1]))
								System.err.println(getProjection()+" returns "+cmdP.args[k]+","+cmdP.args[k+1]+" at "+cmdS.args[k+1]+","+cmdS.args[k]+"!");
						}
						path1.add(cmdP); //TODO: if I was smart, I would divide landmasses that hit an interruption so that I didn't get those annoying lines that cross the map, and then run adaptive resampling to make sure the cuts look clean and not polygonal (e.g. so Antarctica extends all the way to the bottom), but that sounds really hard.
					}
					theMap.add(path1);
					
					i ++;
				}
				
				if (render) {
					updateProgress(-1, 1);
					updateMessage("Rendering map\u2026");
					
					int width, height;
					if (getProjection().isLandscape()) { //fit it to an IMG_SIZE x IMG_SIZE box
						width = IMG_SIZE;
						height = (int)Math.max(IMG_SIZE/getProjection().getAspectRatio(), 1);
					}
					else {
						width = (int)Math.max(IMG_SIZE/getProjection().getAspectRatio(), 1);
						height = IMG_SIZE;
					}
					rendered = drawImage(theMap, width, height);
				}
				
				return new SavableImage() {
					public void save(File file) throws IOException {
						Projection proj = getProjection();
						SVGMap altered = input.replace("Equirectangular", proj.getName());
						altered.save(theMap, file, -proj.getWidth()/2, proj.getHeight()/2,
								proj.getWidth(), proj.getHeight()); //save
					}
				};
			}
			
			protected void succeeded() { //draw the image once successful
				if (render) {
					viewer.getChildren().clear();
					viewer.getChildren().add(rendered);
				}
			}
		};
	}
	
	
	private Canvas drawImage(Iterable<Path> paths, int width, int height) { //parse the SVG path, with a few modifications
		final double mX = this.getProjection().getWidth()/2;
		final double mY = this.getProjection().getHeight()/2;
		Canvas c = new Canvas(width, height);
		GraphicsContext g = c.getGraphicsContext2D();
		g.clearRect(0, 0, c.getWidth(), c.getHeight());
		g.beginPath();
		for (Path path: paths) {
			double startX = 0, startY = 0, lastX = 0, lastY = 0;
			for (Command cmd: path) {
				final double[] args = new double[cmd.args.length];
				for (int i = 0; i < args.length; i ++)
					if (i%2 == 0)
						args[i] = Math2.linInterp(cmd.args[i], -mX, mX, 0, c.getWidth());
					else
						args[i] = Math2.linInterp(cmd.args[i], -mY, mY, c.getHeight(), 0);
				
				switch (cmd.type) {
				case 'M':
					startX = args[0];
					startY = args[1];
					g.moveTo(args[0], args[1]);
					break;
				case 'L':
				case 'T':
					for (int i = 0; i < args.length; i += 2)
						if (Math.hypot(args[i+0]-lastX, args[i+1]-lastY) < IMG_SIZE/4) // break lines that are too long
							g.lineTo(args[i+0], args[i+1]); //TODO: I really need to actually look for interruptions or something
						else
							g.moveTo(args[i+0], args[i+1]);
					break;
				case 'Q':
				case 'S':
					for (int i = 0; i < args.length; i += 4)
						g.quadraticCurveTo(args[i+0], args[i+1], args[i+2], args[i+3]);
					break;
				case 'C':
					for (int i = 0; i < args.length; i += 6)
						g.bezierCurveTo(
								args[i+0], args[i+1], args[i+2], args[i+3], args[i+4], args[i+5]);
					break;
				case 'Z':
					if (Math.hypot(startX-lastX, startY-lastY) < IMG_SIZE/4)
						g.lineTo(startX, startY);
					break;
				default:
					System.err.println("Unsupported movement type: "+cmd.type); //I don't do arcs; they just don't work well with projection
				}
				
				if (args.length > 0) {
					lastX = args[0];
					lastY = args[1];
				}
			}
		}
		g.stroke();
		return c;
	}
}
