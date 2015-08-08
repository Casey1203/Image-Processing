import java.awt.AlphaComposite;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.awt.geom.Ellipse2D;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Cursor;
import java.awt.Frame;
import java.awt.Paint;
import java.awt.GridLayout;
import java.awt.Composite;
import java.awt.Point;
import javax.swing.JWindow;
import javax.swing.UIManager;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import javax.swing.SwingUtilities;
import java.awt.BorderLayout;
import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import javax.imageio.ImageIO;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

public class Assignment2UI extends JFrame {
	JPopupMenu viewportPopup;
	Assignment2 imgProcessor;
	BufferedImage img;
	JWindow toolTip;
	JLabel label;
	ArrayList<double[]> corners = new ArrayList<double[]>();
	ArrayList<double[]> pointList2D = new ArrayList<double[]>();
	ArrayList<double[]> pointList3D = new ArrayList<double[]>();
	private JLabel infoLabel = new JLabel(" ");
	public Assignment2UI() {
		super("COMP 7502 - Assignment 2");
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		JScrollPane scroller = new JScrollPane(new ImagePanel());
		this.add(scroller);
		this.add(infoLabel, BorderLayout.SOUTH);
		this.setSize(750, 600);
		initToolTip();
		this.setVisible(true);
	}

	private void initToolTip() {
		label = new JLabel(" ");
		label.setOpaque(true);
		label.setBackground(UIManager.getColor("ToolTip.background"));
		toolTip = new JWindow(new Frame());
		toolTip.getContentPane().add(label);
		toolTip.addMouseMotionListener(new ToolTipMouseListener(toolTip));
	}
	
	class ToolTipMouseListener implements MouseMotionListener {
		JWindow toolTip;
		public ToolTipMouseListener(JWindow toolTip) {
			this.toolTip = toolTip;
		}
		public void mouseDragged(MouseEvent e) {
		}
		public void mouseMoved(MouseEvent e) {
			toolTip.setVisible(false);
		}
		
	}


	public static void main(String args[]) {
		javax.swing.SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				new Assignment2UI();
			}
		});
	}

	class ImagePanel extends JPanel implements MouseListener, ActionListener, MouseMotionListener {
		int row;
		int column;
		public ImagePanel() {
			imgProcessor = new Assignment2();
			this.addMouseListener(this);
			this.addMouseMotionListener(this);
		}

		public void mouseDragged(MouseEvent e) {}

		public void mouseMoved(MouseEvent e) {
			column = e.getX();
			row = e.getY();
			infoLabel.setText("(x = column = "+column+", y = row = "+row+")");
			if (hitTest(e)) {
				int SELFUZZ = 5;
				double[] selectedPoint = hitTest(e.getX(), e.getY(), SELFUZZ);
				int indexOfSelectedPoint = pointList2D.indexOf(selectedPoint);
				double[] point3d = pointList3D.get(indexOfSelectedPoint);
				double[] point2d = pointList2D.get(indexOfSelectedPoint);
				NumberFormat formatter = new DecimalFormat("#000.00");   
				label.setText("(x="+formatter.format(point2d[0])+", y="+formatter.format(point2d[1])+")->(X="+point3d[0]+", Y="+point3d[1]+", Z="+point3d[2]+")");
				label.repaint();
				toolTip.pack();
				if (!toolTip.isVisible()) {
					toolTip.setVisible(true);
				}
				Point p = e.getPoint();
				SwingUtilities.convertPointToScreen(p, this);
				toolTip.setLocation(p.x + 2, p.y - toolTip.getHeight() - 2);
				//
			} else {
				toolTip.setVisible(false);
				//selectedPoint = null;
			}
			this.invalidate();
			this.repaint();
		}

		private boolean hitTest(MouseEvent e) {
			int SELFUZZ = 5;
			if (hitTest(e.getX(), e.getY(), SELFUZZ) != null) return true;
			else return false;
		}

		private double[] hitTest(double x, double y, double fuzz) {
			if (img != null) {
				for (double[] currentPoint : pointList2D) {
					if (Math.abs(x - currentPoint[0]) < fuzz && Math.abs(y - currentPoint[1]) < fuzz)
						return currentPoint;
				}
			}
			return null;
		}

		public Dimension getPreferredSize() {
			if (img != null) {
				return (new Dimension(img.getWidth(), img.getHeight()));
			} else {
				return (new Dimension(0, 0));
			}
		}

		public void paintComponent(Graphics g) {
			super.paintComponent(g);
			if (img != null) {
				g.drawImage(img, 0, 0, this);
				renderCorners((Graphics2D)g);
				renderCorrespondencePoints((Graphics2D)g, 8.0);
			}

		}
		
		public void renderCorners(Graphics2D g2d) {
			double crossLength = 2.0;
			g2d.setColor(Color.RED);
			if (corners!=null) {
				for (double[] p : corners) {
					Line2D l = new Line2D.Double(p[0]-crossLength, p[1]-crossLength, p[0]+crossLength, p[1]+crossLength);
					g2d.draw(l);
					l = new Line2D.Double(p[0]-crossLength, p[1]+crossLength, p[0]+crossLength, p[1]-crossLength);
					g2d.draw(l);
				}
			}
		}

		public void renderCorrespondencePoints(Graphics2D g, double radius) {
			Composite originalComposite = g.getComposite();
			Paint oldPaint = g.getPaint();

			g.setComposite(makeComposite(0.75f));
			g.setPaint(Color.MAGENTA);
			for (double[] current2DPoint : pointList2D) {
				g.fill(new Ellipse2D.Double(current2DPoint[0] - (radius / 2), current2DPoint[1] - (radius / 2), (radius), (radius)));
			}

			g.setPaint(oldPaint);
			g.setComposite(originalComposite);
		}

		private AlphaComposite makeComposite(float alpha) {
			int type = AlphaComposite.SRC_OVER;
			return (AlphaComposite.getInstance(type, alpha));
		}
		

		private void showPopup(MouseEvent e) {
			JPopupMenu.setDefaultLightWeightPopupEnabled(false);
			viewportPopup = new JPopupMenu();

			JMenuItem openImageMenuItem = new JMenuItem("open image ...");
			openImageMenuItem.addActionListener(this);
			openImageMenuItem.setActionCommand("open image");
			viewportPopup.add(openImageMenuItem);

			JMenuItem loadTestImage1ImageMenuItem = new JMenuItem("load test image 1 from web");
			loadTestImage1ImageMenuItem.addActionListener(this);
			loadTestImage1ImageMenuItem.setActionCommand("load test image 1");
			viewportPopup.add(loadTestImage1ImageMenuItem);

			JMenuItem loadTestImage2ImageMenuItem = new JMenuItem("load test image 2 from web");
			loadTestImage2ImageMenuItem.addActionListener(this);
			loadTestImage2ImageMenuItem.setActionCommand("load test image 2");
			viewportPopup.add(loadTestImage2ImageMenuItem);

			viewportPopup.addSeparator();

			JMenuItem gaussianSmoothMenuItem = new JMenuItem("Task 1: gaussianSmooth");
			gaussianSmoothMenuItem.addActionListener(this);
			gaussianSmoothMenuItem.setActionCommand("gaussianSmooth");
			viewportPopup.add(gaussianSmoothMenuItem);


			JMenuItem detectCornorsMenuItem = new JMenuItem("Task 2: detectCorners");
			detectCornorsMenuItem.addActionListener(this);
			detectCornorsMenuItem.setActionCommand("detect corners");
			viewportPopup.add(detectCornorsMenuItem);

			JMenuItem determineProjectionMatrixMenuItem = new JMenuItem("Task 3: determineProjectionMatrix");
			determineProjectionMatrixMenuItem.addActionListener(this);
			determineProjectionMatrixMenuItem.setActionCommand("determineProjectionMatrix");
			viewportPopup.add(determineProjectionMatrixMenuItem);

			viewportPopup.addSeparator();

			JMenuItem exitMenuItem = new JMenuItem("exit");
			exitMenuItem.addActionListener(this);
			exitMenuItem.setActionCommand("exit");
			viewportPopup.add(exitMenuItem);

			viewportPopup.show(e.getComponent(), e.getX(), e.getY());
		}

		public void mouseClicked(MouseEvent e) {}
		public void mouseEntered(MouseEvent e) {}
		public void mouseExited(MouseEvent e) {}
		public void mouseReleased(MouseEvent e) {}

		public void mousePressed(MouseEvent e) {
			if (viewportPopup != null) {
				viewportPopup.setVisible(false);
				viewportPopup = null;
			} else {
				showPopup(e);
			}
		}

		public void actionPerformed(ActionEvent e) {
			if (e.getActionCommand().equals("open image")) {
				final JFileChooser fc = new JFileChooser();
				FileFilter imageFilter = new FileNameExtensionFilter("Image files", "bmp", "gif", "jpg");
				fc.addChoosableFileFilter(imageFilter);
				fc.setDragEnabled(true);
				fc.setMultiSelectionEnabled(false);
				fc.showOpenDialog(this);
				File file = fc.getSelectedFile();
				try {

					corners.clear();

					pointList2D = new ArrayList<double[]>();
					pointList2D.add(new double[]{145, 128});
					pointList2D.add(new double[]{147, 296});
					pointList2D.add(new double[]{316, 254});
					pointList2D.add(new double[]{318, 103});
					pointList2D.add(new double[]{336, 103});
					pointList2D.add(new double[]{332, 256});
					pointList2D.add(new double[]{484, 317});
					pointList2D.add(new double[]{494, 145});
				
					pointList3D = new ArrayList<double[]>();
					pointList3D.add(new double[]{9.5,0.0,7.5});
					pointList3D.add(new double[]{9.5, 0.0, 0.5});
					pointList3D.add(new double[]{0.5, 0.0, 0.5});
					pointList3D.add(new double[]{0.5, 0.0, 7.5});
					pointList3D.add(new double[]{0.0, 0.5, 7.5});
					pointList3D.add(new double[]{0.0, 0.5, 0.5});
					pointList3D.add(new double[]{0.0, 9.5, 0.5});
					pointList3D.add(new double[]{0.0, 9.5, 7.5});

					img = ImageIO.read(file);
					img = colorToGray(img);
				} catch (Exception ee) {
				}
			} else if (e.getActionCommand().equals("load test image 1")) {
				try {
					img = ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/samplescene1.jpg"));
					img = colorToGray(img);
					///
					corners.clear();

					pointList2D = new ArrayList<double[]>();
					pointList2D.add(new double[]{145, 128});
					pointList2D.add(new double[]{147, 296});
					pointList2D.add(new double[]{316, 254});
					pointList2D.add(new double[]{318, 103});
					pointList2D.add(new double[]{336, 103});
					pointList2D.add(new double[]{332, 256});
					pointList2D.add(new double[]{484, 317});
					pointList2D.add(new double[]{494, 145});
				
					pointList3D = new ArrayList<double[]>();
					pointList3D.add(new double[]{9.5,0.0,7.5});
					pointList3D.add(new double[]{9.5, 0.0, 0.5});
					pointList3D.add(new double[]{0.5, 0.0, 0.5});
					pointList3D.add(new double[]{0.5, 0.0, 7.5});
					pointList3D.add(new double[]{0.0, 0.5, 7.5});
					pointList3D.add(new double[]{0.0, 0.5, 0.5});
					pointList3D.add(new double[]{0.0, 9.5, 0.5});
					pointList3D.add(new double[]{0.0, 9.5, 7.5});
				
					///
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("load test image 2")) {
				try {
					img = ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/samplescene2.jpg"));
					img = colorToGray(img);


					corners.clear();

					pointList2D = new ArrayList<double[]>();
					pointList2D.add(new double[]{133.7037037037037, 120.3982985305491});
					pointList2D.add(new double[]{263.3333333333333, 83.35266821345706});
					pointList2D.add(new double[]{279.6296296296296, 82.98221191028615});
					pointList2D.add(new double[]{461.85185185185185, 111.50734725444701});
					pointList2D.add(new double[]{451.85185185185185, 270.4331013147718});
					pointList2D.add(new double[]{277.4074074074074, 230.79427687548335});
					pointList2D.add(new double[]{262.96296296296293, 231.16473317865427});
					pointList2D.add(new double[]{138.5185185185185, 287.84454756380507});

					pointList3D = new ArrayList<double[]>();
					pointList3D.add(new double[]{9.5, 0.0, 7.5});
					pointList3D.add(new double[]{0.5, 0.0, 7.5});
					pointList3D.add(new double[]{0.0, 0.5, 7.5});
					pointList3D.add(new double[]{0.0, 9.5, 7.5});
					pointList3D.add(new double[]{0.0, 9.5, 0.5});
					pointList3D.add(new double[]{0.0, 0.5, 0.5});
					pointList3D.add(new double[]{0.5, 0.0, 0.5});
					pointList3D.add(new double[]{9.5, 0.0, 0.5});
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("detect corners")) {
				if (img != null) {
					new DetectCornersGUI().createAndShowGUI();
				}
			} else if (e.getActionCommand().equals("gaussianSmooth")) {
				if (img!=null) {
					double sigma = 1.0;
					boolean notOk = true;
					String s = (String)JOptionPane.showInputDialog(this, "Please enter sigma!", ""+sigma);
					if (s==null) {
						return;
					}
					try {
						sigma = Double.parseDouble(s);
					} catch (Exception ee){
					}
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.gaussianSmooth(imgData, img.getWidth(), img.getHeight(), sigma);
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				}
			} else if (e.getActionCommand().equals("determineProjectionMatrix")) {
				if (img != null && corners.size() > 0) {
					Matrix P = imgProcessor.determineProjectionMatrix(corners, pointList2D, pointList3D);
					Matrix K = new Matrix(3, 3);
					Matrix RT = new Matrix(3, 4);
					decomposeProjectionMatrix(P, K, RT);					
					Matrix R = RT.subMat(0, 2, 0, 2);
					Matrix T = RT.subMat(0, 2, 3, 3);
					Matrix C = R.transpose().mul(T).mul(-1);
					this.updateUI();
					JOptionPane.showMessageDialog(this, "Image size = ("+img.getWidth()+" x "+img.getHeight()+")\n\nP = \n" + P + "\n\nK = \n" + K + "\n\n[R|T] = \n" + RT +"\n \n C = \n"+C, "Calibration Result", JOptionPane.INFORMATION_MESSAGE);

				}
			}else if (e.getActionCommand().equals("exit")) {
				System.exit(0);
			}
			viewportPopup = null;
			this.updateUI();
		}
		public BufferedImage colorToGray(BufferedImage source) {
	        BufferedImage returnValue = new BufferedImage(source.getWidth(), source.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
	        Graphics g = returnValue.getGraphics();
	        g.drawImage(source, 0, 0, null);
	        return returnValue;
	    }
	}

	public void decomposeProjectionMatrix(final Matrix P_in, final Matrix K_out, final Matrix RT_out) {
		Matrix Q = new Matrix(3, 3);
		Matrix R = new Matrix(3, 3);
		P_in.subMat(0, 2, 0, 2).transpose().QR2(Q, R);

		K_out.set(R.transpose());
		Matrix tempR = Q.transpose();

		// focal length in x-direction should not be negative
		if (K_out.get(0, 0) < 0) {
			for (int i = 0; i < 3; i++) {
				K_out.set(i, 0, -K_out.get(i, 0));
				tempR.set(0, i, -tempR.get(0, i));
			}
		}

		// focal length in y-direction should not be negative
		if (K_out.get(1, 1) < 0) {
			for (int i = 0; i < 3; i++) {
				K_out.set(i, 1, -K_out.get(i, 1));
				tempR.set(1, i, -tempR.get(1, i));
			}
		}

		// principal point should not be negative
		if (K_out.get(0, 2) < 0) {
			for (int i = 0; i < 3; i++) {
				K_out.set(i, 2, -K_out.get(i, 2));
				tempR.set(2, i, -tempR.get(2, i));
			}
		}

		Matrix t = K_out.inverse().mul(P_in.subMat(0, 2, 3, 3));

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				RT_out.set(i, j, tempR.get(i, j));
		for (int i = 0; i < 3; i++)
			RT_out.set(i, 3, t.get(i, 0));

		K_out.set(K_out.div(K_out.get(2, 2)));
	}
	
	class DetectCornersGUI extends JPanel implements ActionListener {
	    JButton okButton;
	    JButton applyButton;
		JDialog dialog = null;
		JSpinner sigmaSpinner;
		JSpinner thresholdSpinner;
		public DetectCornersGUI() {
			super(new GridLayout(0,1));
	        okButton = new JButton("Ok");
	        okButton.setActionCommand("ok");
	        okButton.addActionListener(this);
	        applyButton = new JButton("Apply");
	        applyButton.setActionCommand("apply");
	        applyButton.addActionListener(this);
	        JPanel rootPanel = new JPanel();
	        GridLayout gridLayout = new GridLayout(2,1, 3, 3);
	        JPanel detectCornersPanel = new JPanel(gridLayout);
	        JLabel sigmaLabel = new JLabel("Sigma: ");
	        double valueSigma = 1.5;
	        double minSigma = 0.01;
	        double maxSigma = 5.0;
	        double stepSigma = 0.25;
	        SpinnerNumberModel sigmaSpinnerNumberModel = new SpinnerNumberModel( valueSigma, minSigma, maxSigma, stepSigma);
	        sigmaSpinner = new JSpinner( sigmaSpinnerNumberModel );
	        detectCornersPanel.add(sigmaLabel);
	        detectCornersPanel.add(sigmaSpinner);
	        JLabel thresholdLabel = new JLabel("Threshold: ");
	        double valueThreshold = 10000000.0;
	        double minThreshold = -100000000;
	        double maxThreshold = 100000000.0;
	        double stepThreshold = 100000;
	        SpinnerNumberModel thresholdSpinnerNumberModel = new SpinnerNumberModel( valueThreshold, minThreshold, maxThreshold, stepThreshold );
	        thresholdSpinner = new JSpinner( thresholdSpinnerNumberModel );
	        detectCornersPanel.add(thresholdLabel);
	        detectCornersPanel.add(thresholdSpinner);
			rootPanel.add(detectCornersPanel);
	        JPanel panel = new JPanel();
	        panel.add(applyButton);
	        rootPanel.add(panel);
	        add(rootPanel);
	        GUIValuesToProgram();
		}
		
	    public void createAndShowGUI() {
	        dialog = new JDialog(Assignment2UI.this, "Corner Detection", false);
	        dialog.setLocation(new Point(Assignment2UI.this.getLocation().x+100, Assignment2UI.this.getLocation().y+100));
	        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
	        dialog.setModal(true);
	        dialog.setContentPane(this);
	        dialog.pack();
	        dialog.setResizable(false);
	        dialog.setVisible(true);
	    }
	    public void GUIValuesToProgram() {
			double sigma = (Double)sigmaSpinner.getValue();
			double threshold = (Double)thresholdSpinner.getValue();
			byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
			corners = imgProcessor.detectCorners(imgData, img.getWidth(), img.getHeight(), sigma, threshold);
			Assignment2UI.this.repaint();
	    }
		public void actionPerformed(ActionEvent e) {
			if (e.getActionCommand() == "apply") {
				GUIValuesToProgram();
				repaint();
			}
		}
	}
}