package chrisdaw.samFastqIO;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

public class CConsole extends JDialog implements ActionListener {

	private static final long serialVersionUID = 10010L;
	public JTextArea cConsole;
	private JButton cancel; //, test;
	
	public CConsole(){ //default constructor to start gui
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); //dispose when jframe is closed
		this.setTitle("Console"); //Title bar
		this.setResizable(false); //don't let user screw with the size of the window
		
		startGUI();
	}

	private void startGUI() { //add all gui components etc
		this.setLayout(new GridBagLayout()); //layout manager
		this.setPreferredSize(new Dimension(450, 500)); //size of the CConsole window
		this.getContentPane().setBackground(Color.white); //set background to white
		GridBagConstraints c = new GridBagConstraints();  //layout manager constraint for layout
		
		//Label the output
		JLabel title = new JLabel("Output: ");
		title.setFont(title.getFont().deriveFont(20.0f));
		c.gridy=0; c.gridwidth = 2;
		c.insets = new Insets(0, 0, 10, 0);
		this.add(title, c);
		
		//area to print to with information from both C and Java
		cConsole = new JTextArea();
    	cConsole.setEditable(false);
    	cConsole.setText("Status: ");
    	cConsole.setLineWrap(true);
    	cConsole.setWrapStyleWord(true);
    	cConsole.setFont(cConsole.getFont().deriveFont(12.5f)); 
    	cConsole.setBorder(javax.swing.BorderFactory.createEmptyBorder());
    	c.ipadx = 20;
    	c.gridx=0; c.gridy=1; 
    	c.gridwidth = 2;
    	JScrollPane scrolly = new JScrollPane(cConsole); //to make the text area scroll
    	scrolly.setPreferredSize(new Dimension(350, 325));
    	scrolly.setBackground(Color.white);
    	this.add(scrolly, c);
    	
    	//button to cancel or close the window (depending on if it is working or not
    	c = new GridBagConstraints(); //reset constraints
    	cancel = new JButton("cancel");
    	cancel.addActionListener(this);
    	c.gridy= 2; c.insets = new Insets(10, 0, 0, 0);
    	c.ipadx = 10; c.ipady = 10;
    	this.add(cancel, c);
		
	}
	
	public void updateConsole(String update, boolean fromJava){ //output to console gui
		if(Display.cv.getState() || fromJava){ //if -v option is checked or if update comes from java ("initialized" etc.)
			cConsole.setText(cConsole.getText() + "\n" + update);
		}
		if(!Display.cworker.isWorking()){ //when it's done working, cancel no longer cancels, so close.
			cancel.setText("close");
		}
	}

	public void actionPerformed(ActionEvent e) { //handle actions
		if(e.getSource() == cancel){ //handle cancel
			if(Display.cworker.isWorking()){ //if it's working
				if (JOptionPane.showConfirmDialog(null, "Are you sure?\nThis will close the program.", "Cancellation", 
						JOptionPane.YES_NO_OPTION, JOptionPane.PLAIN_MESSAGE, null)
						== JOptionPane.YES_OPTION){
					System.out.println("User Exit.");
					System.exit(0);
					//when a user cancels, console get's "User Exit." and program terminates.
				}else{
					//Do nothing, (continue)
				}
			}else{ //if it's not working, it'll say close. If it's close, then it will dispose the panel.
				this.dispose();
			}
		}
	}
}
