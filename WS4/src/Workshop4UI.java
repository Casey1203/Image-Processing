import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
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
import javax.swing.UIManager;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

public class Workshop4UI extends JFrame {

	private JPopupMenu viewportPopup;
	private JLabel infoLabel = new JLabel("");

	public Workshop4UI() {
		super("COMP 7502 - Workshop 4");
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		JScrollPane scroller = new JScrollPane(new ImagePanel());
		this.add(scroller);
		this.add(infoLabel, BorderLayout.SOUTH);
		this.setSize(500, 500);
		this.setVisible(true);
	}

	public static void main(String args[]) {
		try {
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		} catch (Exception e) {}
		new Workshop4UI();
	}

	private class ImagePanel extends JPanel implements MouseListener, ActionListener, MouseMotionListener {
		private BufferedImage img;
		private Workshop4 imgProcessor;
		int row;
		int column;
		
		int otherRow;
		int otherColumn;

		public ImagePanel() {
			imgProcessor = new Workshop4();
			this.addMouseListener(this);
			this.addMouseMotionListener(this);
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

			JMenuItem loadDefaultImageMenuItem = new JMenuItem("load office image from web");
			loadDefaultImageMenuItem.addActionListener(this);
			loadDefaultImageMenuItem.setActionCommand("load office image");
			viewportPopup.add(loadDefaultImageMenuItem);

			JMenuItem loadIllusionImageMenuItem = new JMenuItem("load checker image from web");
			loadIllusionImageMenuItem.addActionListener(this);
			loadIllusionImageMenuItem.setActionCommand("load checker image");
			viewportPopup.add(loadIllusionImageMenuItem);			
			
			viewportPopup.addSeparator();

			JMenuItem gradientImageMenuItem = new JMenuItem("Task 1: gradient image");
			gradientImageMenuItem.addActionListener(this);
			gradientImageMenuItem.setActionCommand("gradientImage");
			viewportPopup.add(gradientImageMenuItem);

			JMenuItem cornerResponseImageMenuItem = new JMenuItem("Task 2: corner response image");
			cornerResponseImageMenuItem.addActionListener(this);
			cornerResponseImageMenuItem.setActionCommand("cornerResponseImage");
			viewportPopup.add(cornerResponseImageMenuItem);
			
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
			} else if (e.getActionCommand().equals("load office image")) {
				try {
					long start = System.nanoTime();
					img = colorToGray(ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/office.png")));
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("load checker image")) {
				try {
					long start = System.nanoTime();
					img = colorToGray(ImageIO.read(new URL("http://www.cs.hku.hk/~sdirk/checker.jpg")));
					double seconds = (System.nanoTime() - start) / 1000000000.0;
					infoLabel.setText(seconds+"");
				} catch (Exception ee) {
					JOptionPane.showMessageDialog(this, "Unable to fetch image from URL", "Error",
							JOptionPane.ERROR_MESSAGE);
					ee.printStackTrace();
				}
			} else if (e.getActionCommand().equals("gradientImage")) {
				if (img!= null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					imgProcessor.gradientImage(imgData, img.getWidth(), img.getHeight());
				}
			} else if (e.getActionCommand().equals("cornerResponseImage")) {
				if (img!= null) {
					byte[] imgData = ((DataBufferByte)img.getRaster().getDataBuffer()).getData();
					imgProcessor.cornerResponseImage(imgData, img.getWidth(), img.getHeight());
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
	        return returnValue;
	    }

		public void mouseDragged(MouseEvent e) {}

		public void mouseMoved(MouseEvent e) {
			column = e.getX();
			row = e.getY();
			infoLabel.setText("("+row+","+column+")");
		}
	}
}

