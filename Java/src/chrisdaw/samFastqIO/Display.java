package chrisdaw.samFastqIO;

import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.Scanner;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;

public class Display extends JPanel implements ActionListener,ItemListener{
	
	public static char[][] chars = new char[25][10]; //chars that will go to CWorker then to C
	public static int varCount = 0;  //amount of the variables in the chars array
	public static CWorker cworker;   //local cworker thread thing
	public CConsole console;         //it's the output console window
	
	private static final long serialVersionUID = -8547241010206408606L; //don't worry about it
	
	//JFrames and JPanels and components to fill said JFrames and JPanels
	private static JFrame mainFrame;
	private JPanel fileSelectPanel, optionsPanel, goPanel;
	private UserHostnameBar inputuhb, outputuhb;
	private JFileChooser fc, fcoutput, fcreference;
	private JLabel selInput, selOutput, selReference;
	private JButton inputButton, outputButton, referenceButton;
	public JButton go;
	private boolean go1 = false, go2 = false, go3 = false;
	public Checkbox cx, cc, cd, cu, cs, cD;
	public static Checkbox cv;
	public static Font f;
	public JTextField tc, td, tu, ts;
	public JRadioButton m, l, a;
	private ButtonGroup bg = new ButtonGroup();
	
	//Strings to help the Display (text that is either default or a file location)
	public String cText = "ratio", dText = "rate",                  //default string
				  uText = "rate", sText = "rate"; 					//default string
	private String ratioDefault = "ratio", rateDefault = "rate";    //default string
	public static String inputFile, outputLocation, referenceFile;  //file locations
	
	private void addCommponentsToFrame(Container frame){ //add panels to main frame
		frame.setLayout(new GridBagLayout());
		frame.setBackground(Color.white);
		addFileSelection();
		addOptions();
		addGo();
		
		GridBagConstraints c = new GridBagConstraints();
		
		c.gridx = 0; c.weightx = 0.1;
		c.weighty=1; c.gridy = 1;
		c.anchor = GridBagConstraints.LINE_START;
		c.insets = new Insets(0, 30, 0, 0);
		frame.add(fileSelectPanel, c);
		
		c.gridx = 1; c.weightx = 0.5;
		c.anchor = GridBagConstraints.CENTER;
		c.insets = new Insets(0, 0, 0, 0);
		frame.add(optionsPanel, c);
		
		c.gridx = 2; c.weightx = 0.1;
		c.anchor = GridBagConstraints.LINE_START;
		c.insets = new Insets(0, 0, 0, 30);
		frame.add(goPanel, c);
		
	}
	
	private void addFileSelection(){  //define file selection components
		GridBagConstraints cFile = new GridBagConstraints();
		
		fileSelectPanel = new JPanel();
		fileSelectPanel.setBackground(Color.white);
		fileSelectPanel.setLayout(new GridBagLayout());
		
		
		////Input
		JLabel input = new JLabel("<html>input file<font color='red'>*</font></html>");
		input.setFont(input.getFont().deriveFont(20.0f));
		cFile.gridy=0; cFile.insets = new Insets(0, 30, 0, 0);
		cFile.anchor = GridBagConstraints.FIRST_LINE_START; 
		cFile.fill = GridBagConstraints.HORIZONTAL;
		fileSelectPanel.add(input, cFile);
		
		fc = new JFileChooser();
		fc.setAcceptAllFileFilterUsed(false);
		fc.addChoosableFileFilter(new SAMFilter());
		inputButton = new JButton("Choose input file...");
		inputButton.addActionListener(this);
		cFile.gridy=1; cFile.insets = new Insets(0, 23, 0, 0);
		fileSelectPanel.add(inputButton, cFile);
		
		selInput = new JLabel("selected input file...");
		selInput.setForeground(Color.gray);
		cFile.gridy=2; cFile.insets = new Insets(0, 50, 0, 0);
		cFile.fill = GridBagConstraints.HORIZONTAL;
		fileSelectPanel.add(selInput, cFile);
		////
		
		inputuhb = new UserHostnameBar();
		inputuhb.setVisible(false);
		cFile.gridy=1; cFile.insets = new Insets(0, 0, 0, 0);
		fileSelectPanel.add(inputuhb, cFile);
		
		outputuhb = new UserHostnameBar();
		outputuhb.setVisible(false);
		cFile.gridy = 4;
		fileSelectPanel.add(outputuhb, cFile);

		////Output
		JLabel output = new JLabel("<html>output location<font color='red'>*</font></html>");
		output.setFont(output.getFont().deriveFont(20.0f));
		cFile.gridy=3; cFile.insets = new Insets(15, 30, 0, 0);
		cFile.fill = GridBagConstraints.HORIZONTAL;
		fileSelectPanel.add(output, cFile);
		
		fcoutput = new JFileChooser();
		fcoutput.setAcceptAllFileFilterUsed(false);
		fcoutput.addChoosableFileFilter(new IdoFilter());
		outputButton = new JButton("Choose output location...");
		outputButton.addActionListener(this);
		cFile.gridy=4; cFile.insets = new Insets(0, 23, 0, 0);
		fileSelectPanel.add(outputButton, cFile);
		
		selOutput = new JLabel("selected output file...");
		selOutput.setForeground(Color.gray);
		cFile.gridy=5; cFile.insets = new Insets(0, 50, 0, 0);
		cFile.fill = GridBagConstraints.HORIZONTAL;
		fileSelectPanel.add(selOutput, cFile);
		////
		
		////Reference
		JLabel reference = new JLabel("<html>reference file<font color='red'>*</font></html>");
		reference.setFont(reference.getFont().deriveFont(20.0f));
		cFile.gridy=6; cFile.insets = new Insets(15, 30, 0, 0);
		cFile.fill = GridBagConstraints.NONE;
		fileSelectPanel.add(reference, cFile);
		
		fcreference = new JFileChooser();
		fcreference.setAcceptAllFileFilterUsed(false);
		fcreference.addChoosableFileFilter(new FastaFilter());
		referenceButton = new JButton("Choose reference file...");
		referenceButton.addActionListener(this);
		cFile.gridy=7; cFile.insets = new Insets(0, 23, 0, 0);
		fileSelectPanel.add(referenceButton, cFile);
		
		selReference = new JLabel("selected reference file...");
		selReference.setForeground(Color.gray);
		selReference.setVisible(true);
		cFile.gridy=8; cFile.insets = new Insets(0, 50, 0, 0);
		cFile.fill = GridBagConstraints.HORIZONTAL;
		fileSelectPanel.add(selReference, cFile);
		////
	}
	 
	private void addOptions(){  //define file selection components
		GridBagConstraints oc = new GridBagConstraints();
		oc.weightx = 1.0;
		oc.weighty = 1.0;
		
		optionsPanel = new JPanel();
		optionsPanel.setBackground(Color.white);
		optionsPanel.setLayout(new GridBagLayout());
		
		JLabel optionsLabel = new JLabel("options: ");
		optionsLabel.setFont(optionsLabel.getFont().deriveFont(17.5f));
		oc.gridx=0; oc.gridy=0; oc.anchor = GridBagConstraints.LINE_START;
		optionsPanel.add(optionsLabel, oc);
		
		JLabel helpOptions = new JLabel("(hover over options for more info)");
		helpOptions.setToolTipText("like that, but over the letters");
		helpOptions.setForeground(Color.darkGray);
		helpOptions.setFont(helpOptions.getFont().deriveFont(12.5f));
		oc.gridx=0; oc.gridy=1; 
		oc.fill = GridBagConstraints.BOTH;
		optionsPanel.add(helpOptions, oc);
		
		f = optionsLabel.getFont();
		
		oc = new GridBagConstraints();
		JLabel v = new JLabel("-v");
		v.setFont(f);
		v.setToolTipText("View extended output.");
		oc.anchor = GridBagConstraints.LINE_START;
		oc.gridx = 0; oc.gridy = 2;
		oc.insets = new Insets(20, 0, 0, 0);
		optionsPanel.add(v, oc);
		
		cv = new Checkbox();
		cv.setState(true);
		oc.insets = new Insets (19, 25, 0, 0);
		optionsPanel.add(cv, oc);
		
		oc.fill = GridBagConstraints.NONE;
		JLabel x = new JLabel("-x");
		x.setFont(f);
		x.setToolTipText("Regenerate file from compressed file");
		oc.gridy = 3; oc.insets = new Insets(5, 0, 0, 0);
		oc.gridx=0; oc.gridwidth = 3; 
		optionsPanel.add(x, oc);
		
		cx = new Checkbox();
		oc.gridx=0; 
		oc.insets = new Insets(4, 25, 0, 0);
		optionsPanel.add(cx, oc);
		
		JLabel c = new JLabel("-c");
		c.setFont(f);
		c.setToolTipText("Compress using [ratio] bits per bit of input entropy per symbol");
		oc.gridy=4; oc.insets = new Insets(0, 0, 0, 0);
		oc.gridx=0;
		optionsPanel.add(c, oc);
		
		cc = new Checkbox();
		cc.addItemListener(this);
		oc.gridx=0;
		oc.insets = new Insets(0, 25, 0, 0);
		optionsPanel.add(cc, oc);
		
		tc = new JTextField(cText);
		tc.setEnabled(false);
		tc.setFont(f);
		tc.setColumns(5);
		tc.setText(ratioDefault);
		tc.setHorizontalAlignment(SwingConstants.CENTER);
		tc.addFocusListener(new FocusListener(){
			public void focusGained(FocusEvent e) {
				if(tc.getText().equals(ratioDefault)){
					tc.setText("");
				}else if(tc.getText().length()>0){
					tc.setText(tc.getText());
				}else{
					tc.setText("");
				}
			}
			public void focusLost(FocusEvent e) {
				if(tc.getText().length()==0){
					tc.setText(ratioDefault);
				}
			}
		});
		oc.anchor = GridBagConstraints.CENTER;
		oc.insets = new Insets(0, 40, 0, 0);
		optionsPanel.add(tc, oc);
		
		oc.anchor = GridBagConstraints.LINE_START;
		JLabel d = new JLabel("-d");
		d.setFont(f);
		d.setToolTipText("Decompress and Download file from remote [output]");
		oc.gridy = 5; oc.gridx=0;
		oc.insets = new Insets(0, 0, 0, 0);
		optionsPanel.add(d, oc);
		
		cd = new Checkbox();
		cd.addItemListener(this);
		oc.gridx=0;
		oc.insets = new Insets(0, 25, 0, 0);
		optionsPanel.add(cd, oc);
		
		td = new JTextField(cText);
		td.setEnabled(false);
		td.setFont(f);
		td.setColumns(5);
		td.setText(rateDefault);
		td.setHorizontalAlignment(SwingConstants.CENTER);
		td.addFocusListener(new FocusListener(){
			public void focusGained(FocusEvent e) {
				if(td.getText().equals(rateDefault)){
					td.setText("");
				}else if(td.getText().length()>0){
					td.setText(td.getText());
				}else{
					td.setText("");
				}
			}
			public void focusLost(FocusEvent e) {
				if(td.getText().length()==0){
					td.setText(rateDefault);
				}
			}
		});
		oc.anchor = GridBagConstraints.CENTER;
		oc.insets = new Insets(0, 40, 0, 0);
		optionsPanel.add(td, oc);
		
		oc.anchor = GridBagConstraints.LINE_START;
		JLabel u = new JLabel("-u");
		u.setFont(f);
		u.setToolTipText("Compress and Upload file to remote [output]");
		oc.gridy=6; oc.gridx=0;
		oc.insets = new Insets(0, 0, 0, 0);
		optionsPanel.add(u, oc);
		
		cu = new Checkbox();
		cu.addItemListener(this);
		oc.gridx=0;
		oc.insets = new Insets(0, 25, 0, 0);
		optionsPanel.add(cu, oc);
		
		tu = new JTextField(cText);
		tu.setEnabled(false);
		tu.setFont(f);
		tu.setColumns(5);
		tu.setText(rateDefault);
		tu.setHorizontalAlignment(SwingConstants.CENTER);
		tu.addFocusListener(new FocusListener(){
			public void focusGained(FocusEvent e) {
				if(tu.getText().equals(rateDefault)){
					tu.setText("");
				}else if(tu.getText().length()>0){
					tu.setText(tu.getText());
				}else{
					tu.setText("");
				}
			}
			public void focusLost(FocusEvent e) {
				if(tu.getText().length()==0){
					tu.setText(rateDefault);
				}
			}
		});
		oc.anchor = GridBagConstraints.CENTER;
		oc.insets = new Insets(0, 40, 0, 0);
		optionsPanel.add(tu, oc);
		
		oc.anchor = GridBagConstraints.LINE_START;
		JLabel s = new JLabel("-s");
		s.setFont(f);
		s.setToolTipText("Stream file to remote [output]");
		oc.gridy =7; oc.gridx=0;
		oc.insets = new Insets(0, 0, 0, 0);
		optionsPanel.add(s, oc);
		
		cs = new Checkbox();
		cs.addItemListener(this);
		oc.gridx=0;
		oc.insets = new Insets(0, 25, 0, 0);
		optionsPanel.add(cs, oc);
		
		ts = new JTextField(cText);
		ts.setEnabled(false);
		ts.setFont(f);
		ts.setColumns(5);
		ts.setText(rateDefault);
		ts.setHorizontalAlignment(SwingConstants.CENTER);
		ts.addFocusListener(new FocusListener(){
			public void focusGained(FocusEvent e) {
				if(ts.getText().equals(rateDefault)){
					ts.setText("");
				}else if(ts.getText().length()>0){
					ts.setText(ts.getText());
				}else{
					ts.setText("");
				}
			}
			public void focusLost(FocusEvent e) {
				if(ts.getText().length()==0){
					ts.setText(rateDefault);
				}
			}
		});
		oc.anchor = GridBagConstraints.CENTER;
		oc.insets = new Insets(0, 40, 0, 0);
		optionsPanel.add(ts, oc);
		
		oc.anchor = GridBagConstraints.LINE_START;
		JLabel D = new JLabel("-D");
		D.setFont(f);
		D.setToolTipText("Optimize for MSE, Lot(0+L0), L0 distortions, respectively (default: MSE)");
		oc.gridy=8; oc.gridx=0;
		oc.insets = new Insets(0, 0, 0, 0);
		optionsPanel.add(D, oc);
		
		cD = new Checkbox();
		cD.addItemListener(this);
		oc.gridx=0;
		oc.insets = new Insets(0, 25, 0, 0);
		optionsPanel.add(cD, oc);
		
		JPanel buttonPanel = new JPanel();
		buttonPanel.setBackground(Color.white);
		buttonPanel.setLayout(new FlowLayout());
		m = new JRadioButton("m");
		m.setEnabled(false);
		buttonPanel.add(m);
		bg.add(m);
		l = new JRadioButton("l");
		l.setEnabled(false);
		buttonPanel.add(l);
		bg.add(l);
		a = new JRadioButton("a");
		a.setEnabled(false);
		buttonPanel.add(a);
		bg.add(a);
		oc.anchor = GridBagConstraints.CENTER;
		oc.insets = new Insets(0, 40, 0, 0);
		optionsPanel.add(buttonPanel, oc);
	}
	
	private void addGo(){  //define go button component
		goPanel = new JPanel();
		goPanel.setBackground(Color.white);
		
		go = new JButton("Go");
		go.addActionListener(this);
		go.setEnabled(false);
		
		goPanel.add(go);
	}
	
	public void actionPerformed(ActionEvent e) {  //define what happens when things are pressed
		if (e.getSource() == inputButton) { //input file selection
			boolean useSam = !cx.getState();

			if(useSam){ //input defaults to SAM so if it's not that it must be ido.
				fc.setFileFilter(new SAMFilter());
			}else{
				fc.setFileFilter(new IdoFilter());
			}

			int returnVal = fc.showOpenDialog(Display.this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				inputFile = fc.getSelectedFile().getAbsolutePath();

				go1 = true;

				if(inputFile.length() > 30){//if the label is too long to display
					int length = inputFile.length();
					String shortInputFile = "..." + inputFile.subSequence(length-30, length);
					selInput.setText(shortInputFile);
				}else{
					selInput.setText(inputFile);
				}
				selInput.setToolTipText(inputFile);
			}

		}else if (e.getSource() == outputButton) { //output location selection
			boolean useIdo = !cx.getState();

			if(useIdo){ //output defaults to ido so if it's not that it must be SAM.
				fcoutput.setFileFilter(new IdoFilter());
			}else{
				fcoutput.setFileFilter(new SAMFilter());
			}

			int returnVal = fcoutput.showSaveDialog(Display.this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				if(fcoutput.getSelectedFile().getName().endsWith((useIdo ? ".ido" : ".sam"))){
					outputLocation = fcoutput.getSelectedFile().getAbsolutePath();
				}else{
					outputLocation = fcoutput.getSelectedFile().getAbsolutePath() + (useIdo ? ".ido" : ".sam");
				}
			
				go2 = true;
	
				if(outputLocation.length() > 30){//if the label is too long to display
					int length = outputLocation.length();
					String shortOutput = "..." + outputLocation.subSequence(length-20, length); //display the last bit
					selOutput.setText(shortOutput);
				}else{
					selOutput.setText(outputLocation);
				}
				selOutput.setToolTipText(outputLocation);
			}
        }else if (e.getSource() == referenceButton) { //reference file selection
            int returnVal = fcreference.showOpenDialog(Display.this);
            
            if (returnVal == JFileChooser.APPROVE_OPTION) {
            	referenceFile = fcreference.getSelectedFile().getAbsolutePath();
                go3 = true;
                
                if(referenceFile.length() > 25){//if the label is too long to display
            		int length = referenceFile.length(); //get the length
            		String shortReference = "..." + referenceFile.subSequence(length-25, length); //put the end on an ellipse
            		selReference.setText(shortReference); //show the shorter name
                }else{
                	selReference.setText(referenceFile); //otherwise just display it normally
                }
                selReference.setToolTipText(referenceFile); //show the full file path if user hovers over it
            }
        }else if(e.getSource() == go){ //go button actions
        	if(cu.getState()){
        		outputLocation = outputuhb.getAbsoluteFilePath();
        	}else if(cs.getState()){
        		outputLocation = outputuhb.getAbsoluteFilePath();
        		inputFile = inputuhb.getAbsoluteFilePath();
        	}else if(cd.getState()){
        		inputFile = inputuhb.getAbsoluteFilePath();
        	}
        	
        	boolean go = false; //do not go yet
        	if(!cu.getState() && !cs.getState() && !cd.getState()){
	        	if(fcoutput.getSelectedFile().getAbsoluteFile().exists()){ //check if output already exists
	        		if (JOptionPane.showConfirmDialog(null, "The file " + outputLocation + " already exists.\n"
	        				+ "Do you wish to overwrite it?", "Overwrite", 
							JOptionPane.YES_NO_OPTION, JOptionPane.PLAIN_MESSAGE, null)
							== JOptionPane.YES_OPTION){
	        			go = true; //yes overwrite
	        		}else{
	        			//don't overwrite, choose a new output
	        			selOutput.setText("");
	        		}
	        	}else{
	        		go = true; //file doesn't exist, go ahead
	        	}
        	}else{
        		go = true;
        	}
        	
        	if(go){ //open the console and run the CWorker 
        		console = new CConsole(); //make a new CConsole object
            	console.pack(); //put everything in order
        		console.setLocationRelativeTo(null); //put in center of screen
        		console.setVisible(true); //make it visible
        		
            	chars = buildForPush(); //see buildForPush()
            	new Thread(new Runnable() {
                    public void run() {
                        SwingUtilities.invokeLater(new Runnable() {
    					    public void run() {
    					    	
    					    	cworker = new CWorker(varCount, chars, console);
    					    	cworker.execute();
    					    }
    					});
                    }
                }).start();
        	}
        	
        }
		if(go1 && go2 && go3){
			go.setEnabled(true); //if all files are selected, user can select go
		}
	}
	
	private char[][] buildForPush(){ //resets variables and reinits. Returns settings in char[][] from getSettings()
		varCount = 0;
		return getSettings();
	}
	
	private char[][] getSettings(){ //return a 2d char array of settings (to send to C)
		StringBuilder settings = new StringBuilder();
		settings.append((cx.getState() ? "-x " : "") + (cc.getState() ? ("-c " + tc.getText() + " ") : "") + 
				(cd.getState() ? ("-d " + td.getText() + " ") : "") + (cu.getState() ? "-u " + tu.getText() + " ": "") + 
				(cs.getState() ? "-s " + ts.getText() + " ": "") + 
				(cD.getState() ? "-D " + (m.isSelected() ? "m " : (l.isSelected() ? "l " : "a ")): ""));
		
		chars[0] = new String("ChrisDaw").toCharArray();
		Scanner delimSpace = new Scanner(settings.toString());
		delimSpace.useDelimiter(" ");
		int currentLeft = 1;
		while(delimSpace.hasNext()){
			chars[currentLeft]=delimSpace.next().toCharArray();
			currentLeft++;
		}
		//put input, output, and reference files on the end. pick if from UHBs or not
		chars[currentLeft] = inputFile.toCharArray(); currentLeft++;
		chars[currentLeft] = outputLocation.toCharArray(); currentLeft++;
		chars[currentLeft] = referenceFile.toCharArray(); 
		
		//to get rid of whitespace and the empty latter half of chars array
		int longest = 0;
		for(char[] arr : chars){
			if(arr.length > longest){
				longest = arr.length;
			}
			if(arr[0] > 0){
				varCount++;
			}
		}
		char[][] temp = new char[varCount][longest];
		for(int i = 0; i < varCount; i++){
			for(int j = 0; j < chars[i].length; j++){
				temp[i][j] = chars[i][j];
			}
		}
		delimSpace.close();
		return temp;
	}
	
	public static void createAndShowGUI(){
		mainFrame = new JFrame("samFastqIO");
		mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		Display display = new Display();
		display.addCommponentsToFrame(mainFrame.getContentPane());
		
		mainFrame.setPreferredSize(new Dimension(750, 500));
		mainFrame.setResizable(false);
		mainFrame.pack();
		mainFrame.setLocationRelativeTo(null);
		mainFrame.setVisible(true);
	}

	public void itemStateChanged(ItemEvent e) {
		if(e.getSource() == cs){ //if checkbox -s is changing
			go2 = cs.getState() ? true : false;
			ts.setEnabled(cs.getState());
			if(!cs.getState()){ //if it is turning off
				ts.setText(rateDefault);
				
				cu.setEnabled(true);
				cd.setEnabled(true);
				
				inputuhb.setVisible(false);
				outputuhb.setVisible(false);
				inputButton.setVisible(true);
				selInput.setVisible(true);
				outputButton.setVisible(true);
				selOutput.setVisible(true);
			}else{  //if it is turning on
				cu.setEnabled(false);
				cd.setEnabled(false);
				
				inputuhb.setVisible(true);
				outputuhb.setVisible(true);
				inputButton.setVisible(false);
				selInput.setVisible(false);
				outputButton.setVisible(false);
				selOutput.setVisible(false);
			}
		}
		if(e.getSource() == cD){ //if checkbox -D is changing
			m.setEnabled(cD.getState());
			l.setEnabled(cD.getState());
			a.setEnabled(cD.getState());
		}
		if(e.getSource() == cu){ //if checkbox -u is changing
			go2 = cu.getState() ? true : false;
			tu.setEnabled(cu.getState());
			if(!cu.getState()){ //if it is turning off
				tu.setText(rateDefault);
				
				cs.setEnabled(true);
				cd.setEnabled(true);
				
				outputuhb.setVisible(false);
				outputButton.setVisible(true);
				selOutput.setVisible(true);
			}else{ //if it is turning on
				cs.setEnabled(false);
				cd.setEnabled(false);
				
				outputuhb.setVisible(true);
				outputButton.setVisible(false);
				selOutput.setVisible(false);
			}
		}
		if(e.getSource() == cd){ //if checkbox -d is changing
			go2 = cd.getState() ? true : false;
			td.setEnabled(cd.getState());
			if(!cd.getState()){ //if it is turning off
				td.setText(rateDefault);
				
				cs.setEnabled(true);
				cu.setEnabled(true);
				
				inputuhb.setVisible(false);
				inputButton.setVisible(true);
				selInput.setVisible(true);
			}else{ //if it is turning on
				cs.setEnabled(false);
				cu.setEnabled(false);
				
				inputuhb.setVisible(true);
				inputButton.setVisible(false);
				selInput.setVisible(false);
			}
		}
		if(e.getSource() == cc){ //if checkbox -c is changing
			tc.setEnabled(cc.getState()); 
			if(!cc.getState()){ //if it is turning off
				tc.setText(ratioDefault);
			}
		}
		if(go1 && go2 && go3){
			go.setEnabled(true); //if all files are selected, user can select go
		}
	}
}

//File filters to restrict file selection to acceptable files

class SAMFilter extends FileFilter{
	
	public boolean accept(File f){
		return f.isDirectory() || f.getAbsolutePath().endsWith(".sam");
	}

	public String getDescription() {
		return "SAM file (*.sam)";
	}
}

class IdoFilter extends FileFilter{
	public boolean accept(File f){
		return f.isDirectory() || f.getAbsolutePath().endsWith(".ido");
	}
	public String getDescription(){
		return "Ido file (*.ido)";
	}
}

class FastaFilter extends FileFilter{
	
	public boolean accept(File f){
		return f.isDirectory() || f.getAbsolutePath().endsWith(".fa") || 
				f.getAbsolutePath().endsWith(".fasta");
	}

	public String getDescription() {
		return "Fasta file (*.fa, *.fasta)";
	}
}
