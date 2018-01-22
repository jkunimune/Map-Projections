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
import java.util.function.DoubleConsumer;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import dialogs.ProgressBarDialog;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.control.Button;
import javafx.scene.control.Separator;
import javafx.scene.layout.HBox;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import maps.Projection;
import utils.Math2;
import utils.Procedure;
import utils.SVGMap;
import utils.SVGMap.Command;
import utils.SVGMap.Path;

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
	
	private Node aspectSelector;
	private Button saveBtn;
	private double[] aspect;
	private SVGMap input;
	private Canvas viewer;
	
	
	
	public MapDesignerVector() {
		super("Map Designer");
	}
	
	
	
	@Override
	public void start(Stage root) {
		super.start(root);
		new Thread(() -> {
			setInput(new File("input/Basic.svg")); //this automatically updates the map//TODO: this should cause the buttons to grey out
		}).start();
	}
	
	
	@Override
	protected Node makeWidgets() {
		this.aspect = new double[3];
		final Node inputSelector = buildInputSelector(VECTOR_TYPES,
				VECTOR_TYPES[0], this::setInput);
		final Node projectionSelector = buildProjectionSelector(
				Procedure.concat(this::updateMap, this::hideAspect));
		this.aspectSelector = buildAspectSelector(this.aspect, this::updateMap);
		final Node parameterSelector = buildParameterSelector(this::updateMap);
		this.saveBtn = buildSaveButton(true, "map", VECTOR_TYPES,
				VECTOR_TYPES[0], ()->true, this::calculateAndSaveMap);
		
		final VBox layout = new VBox(V_SPACE,
				inputSelector, new Separator(), projectionSelector,
				new Separator(), aspectSelector, parameterSelector,
				new Separator(), saveBtn);
		
		layout.setAlignment(Pos.CENTER);
		layout.setPrefWidth(GUI_WIDTH);
		
		viewer = new Canvas(1,1);
		final StackPane pane = new StackPane(viewer);
		pane.setMinWidth(IMG_WIDTH);
		pane.setMinHeight(IMG_WIDTH);
		
		final HBox gui = new HBox(MARGIN, layout, pane);
		gui.setAlignment(Pos.CENTER);
		StackPane.setMargin(gui, new Insets(MARGIN));
		
		return gui;
	}
	
	
	private void setInput(File file) {
		saveBtn.setDisable(true);
		
		try {
			input = new SVGMap(file);
		} catch (IOException e) {
			showError("File not found!",
					"We couldn't find "+file.getAbsolutePath()+".");
		} catch (SAXException e) {
			showError("Unreadable file!",
					"We couldn't read "+file.getAbsolutePath()+". It may be corrupt or an unreadable format.");
		} catch (ParserConfigurationException e) {
			// TODO: Handle this
			e.printStackTrace();
		} finally {
			saveBtn.setDisable(false);
		}
		
		updateMap();
	}
	
	
	
	private void hideAspect() {
		aspectSelector.setVisible(this.getProjection().hasAspect());
	}
	
	
	private void updateMap() {
		loadParameters();
		int maxVtx = this.getParamsChanging() ? FAST_MAX_VTX : DEF_MAX_VTX;
		final Iterable<Path> transformed = map(input, input.length()/maxVtx+1, aspect.clone(), null);
		
		if (this.getProjection().getWidth() > this.getProjection().getHeight()) {
			viewer.setWidth(IMG_WIDTH);
			viewer.setHeight(IMG_WIDTH/this.getProjection().getAspectRatio());
		}
		else {
			viewer.setWidth(IMG_WIDTH*this.getProjection().getAspectRatio());
			viewer.setHeight(IMG_WIDTH);
		}
		drawImage(transformed, viewer);
	}
	
	
	private void calculateAndSaveMap(File file, ProgressBarDialog pBar) {
		loadParameters();
		final List<Path> transformed = map(input, 1, aspect.clone(), pBar::setProgress); //calculate
		try {
			Projection proj = this.getProjection();
			SVGMap altered = input.replace("Equirectangular", proj.getName());
			altered.save(transformed, file, -proj.getWidth()/2, proj.getHeight()/2,
					proj.getWidth(), proj.getHeight(), pBar::setProgress); //save
		} catch (IOException e) {
			showError("Failure!",
					"Could not access "+file.getAbsolutePath()+". It's possible that another program has it open.");
		}
	}
	
	
	private List<Path> map(SVGMap input, int step, double[] pole, DoubleConsumer tracker) {
		Projection proj = this.getProjection();
		List<Path> output = new LinkedList<Path>();
		int i = 0;
		for (Path path0: input) {
			Path path1 = new Path();
			int counter = 0;
			for (Command cmdS: path0) {
				counter --;
				if (counter > 0 && cmdS.type != 'M' && cmdS.type != 'Z') 	continue;
				counter = step;
				
				Command cmdP = new Command(cmdS.type, new double[cmdS.args.length]);
				for (int k = 0; k < cmdS.args.length; k += 2) {
					double[] coords = proj.project(cmdS.args[k+1], cmdS.args[k], pole);
					cmdP.args[k] =
							Math.max(Math.min(coords[0], proj.getWidth()), -proj.getWidth());
					cmdP.args[k+1] =
							Math.max(Math.min(coords[1], proj.getHeight()), -proj.getHeight());
					if (Double.isNaN(cmdP.args[k]) || Double.isNaN(cmdP.args[k+1]))
						System.err.println(this.getProjection()+" returns "+cmdP.args[k]+","+cmdP.args[k+1]+" at "+cmdS.args[k+1]+","+cmdS.args[k]+"!");
				}
				path1.add(cmdP); //TODO: if I was smart, I would divide landmasses that hit an interruption so that I didn't get those annoying lines that cross the map, and then run adaptive resampling to make sure the cuts look clean and not polygonal (e.g. so Antarctica extends all the way to the bottom), but that sounds really hard.
			}
			output.add(path1);
			
			if (tracker != null) {
				i ++;
				tracker.accept((double)i/input.numCurves());
			}
		}
		return output;
	}
	
	
	private void drawImage(
			Iterable<Path> paths, Canvas c) { //parse the SVG path, with a few modifications
		final double mX = this.getProjection().getWidth()/2;
		final double mY = this.getProjection().getHeight()/2;
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
						if (Math.hypot(args[i+0]-lastX, args[i+1]-lastY) < IMG_WIDTH/4) // break lines that are too long
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
					if (Math.hypot(startX-lastX, startY-lastY) < IMG_WIDTH/4)
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
	}

}
