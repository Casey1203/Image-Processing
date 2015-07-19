import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.File;
import java.net.URL;

import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

public class Workshop2UI extends JFrame {

	private JPopupMenu viewportPopup;
	private JLabel infoLabel = new JLabel("");

	public Workshop2UI() {
		super("COMP 7502 - Workshop 2");
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		JScrollPane scroller = new JScrollPane(new ImagePanel());
		this.add(scroller);
		this.add(infoLabel, BorderLayout.SOUTH);
		this.setSize(500, 500);
		this.setVisible(true);
	}

	public static void main(String args[]) {
		new Workshop2UI();
	}

	private class ImagePanel extends JPanel implements MouseListener, ActionListener {
		private BufferedImage img;
		private Workshop2 imgProcessor;

		public ImagePanel() {
			imgProcessor = new Workshop2();
			this.addMouseListener(this);
		}

		public Dimension getPreferredSize() {
			if (img != null) return (new Dimension(img.getWidth(), img.getHeight()));
			else return (new Dimension(0, 0));
		}

		public void paintComponent(Graphics g) {
			super.paintComponent(g);
			if (img != null)
				g.drawImage(img, 0, 0, this);
		}

		private void showPopup(MouseEvent e) {
			JPopupMenu.setDefaultLightWeightPopupEnabled(false);
			viewportPopup = new JPopupMenu();
			JMenuItem openImageMenuItem = new JMenuItem("load any image ...");
			openImageMenuItem.addActionListener(this);
			openImageMenuItem.setActionCommand("open image");
			viewportPopup.add(openImageMenuItem);

			JMenuItem loadDefaultImageMenuItem = new JMenuItem("load camel image from web");
			loadDefaultImageMenuItem.addActionListener(this);
			loadDefaultImageMenuItem.setActionCommand("load default image");
			viewportPopup.add(loadDefaultImageMenuItem);

			JMenuItem loadIllusionImageMenuItem = new JMenuItem("load illusion image from web");
			loadIllusionImageMenuItem.addActionListener(this);
			loadIllusionImageMenuItem.setActionCommand("load illusion image");
			viewportPopup.add(loadIllusionImageMenuItem);			
			 
			viewportPopup.addSeparator();

			JMenuItem thresholdingMenuItem = new JMenuItem("Task 1: thresholdTransformation");
			thresholdingMenuItem.addActionListener(this);
			thresholdingMenuItem.setActionCommand("thresholding");
			viewportPopup.add(thresholdingMenuItem);

			JMenuItem negativeTransformationMenuItem = new JMenuItem("Task 2: negativeTransformation");
			negativeTransformationMenuItem.addActionListener(this);
			negativeTransformationMenuItem.setActionCommand("negativeTransformation");
			viewportPopup.add(negativeTransformationMenuItem);			

			JMenuItem logTransformationMenuItem = new JMenuItem("Task 3: logTransformation");
			logTransformationMenuItem.addActionListener(this);
			logTransformationMenuItem.setActionCommand("logTransformation");
			viewportPopup.add(logTransformationMenuItem);

			JMenuItem bitPlaneSlicingMenuItem = new JMenuItem("Task 4: bitPlaneSlicing");
			bitPlaneSlicingMenuItem.addActionListener(this);
			bitPlaneSlicingMenuItem.setActionCommand("bitPlaneSlicing");
			viewportPopup.add(bitPlaneSlicingMenuItem);

			viewportPopup.addSeparator();
			
			JMenuItem showHistogramMenuItem = new JMenuItem("Task 5: histogram");
			showHistogramMenuItem.addActionListener(this);
			showHistogramMenuItem.setActionCommand("show histogram");
			viewportPopup.add(showHistogramMenuItem);

			JMenuItem histogramEqualizationMenuItem = new JMenuItem("Task 6: histogramEqualization");
			histogramEqualizationMenuItem.addActionListener(this);
			histogramEqualizationMenuItem.setActionCommand("histogram equalization");
			viewportPopup.add(histogramEqualizationMenuItem);
			
			viewportPopup.addSeparator();
			
			JMenuItem boxSmoothMenuItem = new JMenuItem("Homework 1: boxSmoothFilter");
			boxSmoothMenuItem.addActionListener(this);
			boxSmoothMenuItem.setActionCommand("box smooth");
			viewportPopup.add(boxSmoothMenuItem);
			
			JMenuItem medianFilterMenuItem = new JMenuItem("Homework 2: medianFilter");
			medianFilterMenuItem.addActionListener(this);
			medianFilterMenuItem.setActionCommand("median filter");
			viewportPopup.add(medianFilterMenuItem);
			
			JMenuItem laplacianFilterMenuItem = new JMenuItem("Homework 3: laplacianFilter");
			laplacianFilterMenuItem.addActionListener(this);
			laplacianFilterMenuItem.setActionCommand("laplacian filter");
			viewportPopup.add(laplacianFilterMenuItem);
			
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
			} else
				showPopup(e);
		}

		@Override
		public void actionPerformed(ActionEvent e) {
			if (e.getActionCommand().equals("open image")) {
				final JFileChooser fc = new JFileChooser();
				FileFilter imageFilter = new FileNameExtensionFilter("Image files", "bmp", "gif", "jpg");
				fc.addChoosableFileFilter(imageFilter);
				fc.setDragEnabled(true);
				fc.setMultiSelectionEnabled(false);
				int result =  fc.showOpenDialog(this);
				if (result == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					try {
						long start = System.nanoTime();
						img = colorToGray((ImageIO.read(file)));
						double seconds = (System.nanoTime() - start) / 1000000000.0;
						infoLabel.setText(seconds+"");
					} catch (Exception ee) {
						ee.printStackTrace();
					}
				}
			} else if (e.getActionCommand().equals("load default image")) {
				try {
					long start = System.nanoTime();
					img = colorToGray(ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/georgesteinmetz.jpg")));
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("load illusion image")) {
				try {
					long start = System.nanoTime();
					img = colorToGray(ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/illusion.png")));
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("show histogram")) {
				if (img!=null) {
					JFrame frame = new JFrame();
					frame.setTitle("Histogram");					
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					frame.add(new HistogramPanel(imgProcessor.histogram(imgData)));
					frame.pack();
					frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
					frame.setResizable(false);
					frame.setVisible(true);
				}
			}  else if (e.getActionCommand().equals("histogram equalization")) {
				if (img!=null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.histogramEqualization(imgData);
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");	
				}
			} else if (e.getActionCommand().equals("bitPlaneSlicing")) {
				if (img!=null) {
					int mask = 126;
					boolean notOk = true;
					while (notOk) {
						String s = (String)JOptionPane.showInputDialog(this, "Please enter a mask (in decimal)! \n(0 <= mask <= 255)", ""+mask);
						if (s==null) {
							return;
						}
						try {
							int i = Integer.parseInt(s);
							if (0<=i && i<=255) {
								notOk = false;
								mask = i;
							}
						} catch (Exception ee){
						}
					}
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.bitPlaneSlicing(imgData, mask);
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				}
			} else if (e.getActionCommand().equals("thresholding")) {
				if (img!=null) {
					int threshold = 126;
					boolean notOk = true;
					while (notOk) {
						String s = (String)JOptionPane.showInputDialog(this, "Please enter a threshold! \n(0 <= threshold <= 255)", ""+threshold);
						if (s==null) {
							return;
						}
						try {
							int i = Integer.parseInt(s);
							if (0<=i && i<=255) {
								notOk = false;
								threshold = i;
							}
						} catch (Exception ee){
						}
					}
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.thresholdTransformation(imgData, threshold);
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				}
			} else if (e.getActionCommand().equals("box smooth")) {
				if (img!=null) {
					int filterSize = 7;
					boolean notOk = true;
					while (notOk) {
						String s = (String)JOptionPane.showInputDialog(this, "Please enter a filter size! \n(Must be odd and smaller than image width and height)", ""+filterSize);
						if (s==null) {
							return;
						}
						try {
							int i = Integer.parseInt(s);
							if (i%2==1 && i<img.getWidth() && i<img.getHeight()) {
								notOk = false;
								filterSize = i;
							}
						} catch (Exception ee){
						}
					}
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.boxSmoothFilter(imgData, img.getWidth(), img.getHeight(), filterSize);
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				}
			} else if (e.getActionCommand().equals("median filter")) {
				if (img!=null) {
					int filterSize = 7;
					boolean notOk = true;
					while (notOk) {
						String s = (String)JOptionPane.showInputDialog(this, "Please enter a filter size! \n(Must be odd and smaller than image width and height)", ""+filterSize);
						if (s==null) {
							return;
						}
						try {
							int i = Integer.parseInt(s);
							if (i%2==1 && i<img.getWidth() && i<img.getHeight()) {
								notOk = false;
								filterSize = i;
							}
						} catch (Exception ee){
						}
					}
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.medianFilter(imgData, img.getWidth(), img.getHeight(), filterSize);
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				}
			} else if (e.getActionCommand().equals("laplacian filter")) {
				if (img!=null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.laplacianFilter(imgData, img.getWidth(), img.getHeight());
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");	
				}
			} else if (e.getActionCommand().equals("logTransformation")) {
				if (img!=null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.logTransformation(imgData);
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");	
				}
			} else if (e.getActionCommand().equals("negativeTransformation")) {
				if (img!=null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.negativeTransformation(imgData);
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");	
				}
			} else if (e.getActionCommand().equals("exit")) {
				System.exit(0);
			}
			viewportPopup = null;
			this.updateUI();

		}
		
		public BufferedImage colorToGray(BufferedImage source) {
	        BufferedImage returnValue = new BufferedImage(source.getWidth(), source.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
	        Graphics g = returnValue.getGraphics();
	        g.drawImage(source, 0, 0, null);
	        //g.dispose();
	        return returnValue;
	    }
	}
	
	private class HistogramPanel extends JPanel {
		
		private int[] histogram;
		int width = 512;
		int height = 256;
		
		public HistogramPanel(int[] histogram) {
			this.histogram = histogram;
		}
		
		public Dimension getPreferredSize() {
			return (new Dimension(width, height));
		}
		
		public int rescale(int in, int max) {
			return (int)(in*(255.0)/(max));
		}
		
		public void paintComponent(Graphics g) {
			super.paintComponent(g);
			Graphics2D g2 = (Graphics2D) g;
			
			int max = 0;
			for (int i=0;i<histogram.length;i++) {
				if (histogram[i]>max)
					max = histogram[i];
			}
			
			for (int i=0;i<histogram.length;i++) {
				g2.draw(new Line2D.Double(2*i, height, 2*i, height-rescale(histogram[i], max)));
				g2.draw(new Line2D.Double(2*i+1, height, 2*i+1, height-rescale(histogram[i], max)));
			}	
		}	
	}
}

