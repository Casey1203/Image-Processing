


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

public class Workshop3UI extends JFrame {

	private JPopupMenu viewportPopup;
	private JLabel infoLabel = new JLabel("");

	public Workshop3UI() {
		super("COMP 7502 - Workshop 3");
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		JScrollPane scroller = new JScrollPane(new ImagePanel());
		this.add(scroller);
		this.add(infoLabel, BorderLayout.SOUTH);
		this.setSize(500, 500);
		this.setVisible(true);
	}

	public static void main(String args[]) {
		new Workshop3UI();
	}

	private class ImagePanel extends JPanel implements MouseListener, ActionListener {
		private BufferedImage img;
		private Workshop3 imgProcessor;

		public ImagePanel() {
			imgProcessor = new Workshop3();
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

			JMenuItem loadLennaMenuItem = new JMenuItem("load lenna from web");
			loadLennaMenuItem.addActionListener(this);
			loadLennaMenuItem.setActionCommand("load lenna from web");
			viewportPopup.add(loadLennaMenuItem);

			JMenuItem loadLennaSmallMenuItem = new JMenuItem("load lenna small from web");
			loadLennaSmallMenuItem.addActionListener(this);
			loadLennaSmallMenuItem.setActionCommand("load lenna small from web");
			viewportPopup.add(loadLennaSmallMenuItem);

			JMenuItem loadSEMMenuItem = new JMenuItem("load SEM from web");
			loadSEMMenuItem.addActionListener(this);
			loadSEMMenuItem.setActionCommand("load SEM from web");
			viewportPopup.add(loadSEMMenuItem);

			JMenuItem loadSEMSmallMenuItem = new JMenuItem("load SEM small from web");
			loadSEMSmallMenuItem.addActionListener(this);
			loadSEMSmallMenuItem.setActionCommand("load SEM small from web");
			viewportPopup.add(loadSEMSmallMenuItem);
			 
			viewportPopup.addSeparator();

			JMenuItem fourierTransformMenuItem = new JMenuItem("Task 1: fourierTransform");
			fourierTransformMenuItem.addActionListener(this);
			fourierTransformMenuItem.setActionCommand("fourierTransform");
			viewportPopup.add(fourierTransformMenuItem);			

			JMenuItem changeImageMenuItem = new JMenuItem("Task 2: changeImage");
			changeImageMenuItem.addActionListener(this);
			changeImageMenuItem.setActionCommand("change image");
			viewportPopup.add(changeImageMenuItem);

			JMenuItem filteringMenuItem = new JMenuItem("Task 3: filtering");
			filteringMenuItem.addActionListener(this);
			filteringMenuItem.setActionCommand("filtering");
			viewportPopup.add(filteringMenuItem);

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
			} else if (e.getActionCommand().equals("load lenna from web")) {
				try {
					long start = System.nanoTime();
					img = colorToGray(ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/lenna.png")));
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("load lenna small from web")) {
				try {
					long start = System.nanoTime();
					img = colorToGray(ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/lenna_small.png")));
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("load SEM from web")) {
				try {
					long start = System.nanoTime();
					img = colorToGray(ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/SEM.png")));
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("load SEM small from web")) {
				try {
					long start = System.nanoTime();
					img = colorToGray(ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/SEM_small.png")));
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("change image")) {
				if (img!=null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.changeImage(imgData, img.getWidth(), img.getHeight());
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");	
				}
			} else if (e.getActionCommand().equals("fourierTransform")) {
				if (img!=null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.fourierTransform(imgData, img.getWidth(), img.getHeight());
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");	
				}
			} else if (e.getActionCommand().equals("filtering")) {
				if (img!=null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					long start = System.nanoTime();
					imgProcessor.filtering(imgData, img.getWidth(), img.getHeight());
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
	
}

