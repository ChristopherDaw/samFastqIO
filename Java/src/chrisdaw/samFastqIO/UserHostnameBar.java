package chrisdaw.samFastqIO;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;

public class UserHostnameBar extends JPanel{
	private static final long serialVersionUID = -2553851197679275766L;
	
	private JTextField user, hostname, pathToFile;
	private String userDefault = "user", hostnameDefault = "hostname", pathDefault = "path to file";

	UserHostnameBar(){ //default constructor does mostly everything
		this.setBackground(Color.white); //white background
		this.setLayout(new GridBagLayout()); //layout manager
		
		GridBagConstraints c = new GridBagConstraints();
		
		//first text field, takes user name
		user = new JTextField();
		user.setFont(user.getFont().deriveFont(12.5f));
		user.setText(userDefault);
		user.setHorizontalAlignment(SwingConstants.CENTER);
		user.addFocusListener(new FocusListener(){ //to make the hint text go away and come back if it's empty 
			public void focusGained(FocusEvent e) { //when it gains and loses focus respectively
				if(user.getText().equals(userDefault)){
					user.setText("");
				}else if(user.getText().length()>0){
					user.setText(user.getText());
				}else{
					user.setText("");
				}
			}
			public void focusLost(FocusEvent e) {
				if(user.getText().length()==0){
					user.setText(userDefault);
				}
			}
		});
		c.gridx = 0;
		c.ipadx = 60;
		c.ipady = 5;
		c.gridy = 0;
		c.anchor = GridBagConstraints.LINE_START;
		this.add(user, c);
		
		//Font object to make our lives easier (uses the font from the user text area.
		Font f = user.getFont();
		
		//Just a label to denote there will be an '@' in the address
		JLabel at = new JLabel("@");
		at.setFont(f);
		c.anchor = GridBagConstraints.LINE_START;
		c.gridx = 1;
		this.add(at, c);
		
		//text area to collect the hostname
		hostname = new JTextField();
		hostname.setText(hostnameDefault);
		hostname.setHorizontalAlignment(SwingConstants.CENTER);
		hostname.addFocusListener(new FocusListener(){ //exactly the same as user's focus gained/lost functions
			public void focusGained(FocusEvent e) {
				if(hostname.getText().equals(hostnameDefault)){
					hostname.setText("");
				}else if(hostname.getText().length()>0){
					hostname.setText(hostname.getText());
				}else{
					hostname.setText("");
				}
			}
			public void focusLost(FocusEvent e) {
				if(hostname.getText().length()==0){
					hostname.setText(hostnameDefault);
				}
			}
		});
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0, 10, 0, 0);
		c.ipadx = 50;
		c.gridx = 1;
		this.add(hostname, c);
		
		//label to denote a ':' will be added to the address here
		JLabel colon = new JLabel(":");
		at.setFont(f);
		c.anchor = GridBagConstraints.LINE_START;
		c.insets = new Insets(0, 0, 0, 0);
		c.gridx=2;
		this.add(colon, c);
		
		//text area to collect the path to the specified file
		pathToFile = new JTextField();
		pathToFile.setText(pathDefault);
		pathToFile.setHorizontalAlignment(SwingConstants.CENTER);
		pathToFile.addFocusListener(new FocusListener(){ //exactly the same as user's and hostname's functions
			public void focusGained(FocusEvent e) {
				if(pathToFile.getText().equals(pathDefault)){
					pathToFile.setText("");
				}else if(pathToFile.getText().length()>0){
					pathToFile.setText(pathToFile.getText());
				}else{
					pathToFile.setText("");
				}
			}
			public void focusLost(FocusEvent e) {
				if(pathToFile.getText().length()==0){
					pathToFile.setText(pathDefault);
				}
			}
		});
		c.anchor = GridBagConstraints.CENTER;
		c.ipadx = 20;
		c.gridx = 1;
		c.gridy = 1;
		c.gridwidth = 4;
		c.fill = GridBagConstraints.BOTH;
		this.add(pathToFile, c);
	}
	
	public String getAbsoluteFilePath(){ //function to collect all given data from the UserHostnameBar object.
		return user.getText() + "@" + hostname.getText() + ":" + pathToFile.getText();
	}
	
}
