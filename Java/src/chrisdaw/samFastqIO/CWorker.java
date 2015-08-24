package chrisdaw.samFastqIO;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.SwingWorker;

public class CWorker extends SwingWorker<Integer, String>{

	public static native void runmain(int varCount, byte[][] var, int sizeX, int sizeY);
	
	private char[][] arguments; //to be changed to bytes
	private int argCount, sizeX, sizeY;
	private byte[][] argToC;  //to be sent to C
	private static CConsole console;
	private static boolean isWorking = false;
	
	public CWorker(int argc, char[][] argv, CConsole cc){
		arguments = argv;
		argCount = argc; 
		sizeX = getSizeX(arguments);
		sizeY = getSizeY(arguments);
		argToC = charToByteArray(arguments);
		console = cc;
		isWorking = true;
	}
	
	private int getSizeX(char[][] arr){ //returns size of x length of 2D array
		int size = 0;
		for(char[] array : arr){
			if(array[0] > 0){
				size++;
			}
		}
		return size;
	}
	
	private int getSizeY(char[][] arr){ //returns size of y length of 2D array
		int longest = 0;
		for(char[] array : arr){
			if(array.length > longest){
				longest = array.length;
			}
		}
		return longest;
	}
	
	private byte[][] charToByteArray(char[][] chars){  //transform char[][] to byte[][]
		byte[][] bytes = new byte[sizeX][sizeY];
		for(int i = 0; i < getSizeX(chars); i++){
			for(int j = 0; j < chars[i].length; j++){
				bytes[i][j] = (byte) chars[i][j];
				System.out.print(chars[i][j]);
			}System.out.print(" ");
		}System.out.println();
		return bytes;
	}
	
	private static byte[] charToByteArray(char[] chars){
		byte[] bytes = new byte[chars.length];
		for(int i = 0; i < chars.length; i++){
			bytes[i] = (byte)chars[i];
		}
		return bytes;
	}

	
	public static void updateCConsole(String[] strings){ //used in the C code
		for(String string : strings){
			console.updateConsole(string, false);
		}
	}
	public static void updateCConsole(String string){ //used in below method
		console.updateConsole(string, true);
	}
	
	public boolean isWorking(){ //return true if is still working
		return isWorking;
	}
	
	public static void finishWorking(){ //called from C, finishes Console output
		isWorking = false;
		updateCConsole("File output to " + Display.outputLocation +"\n");
		updateCConsole("Done Working");
	}

	public static byte[] askForPassword(){ //called from C if it needs an ssh password
		JPanel panel = new JPanel();
		JLabel label = new JLabel("Password: ");
		JPasswordField pass = new JPasswordField(15);
		panel.add(label);
		panel.add(pass);
		String[] options = new String[]{"Enter", "Cancel"};
		int option = JOptionPane.showOptionDialog(null, panel, "Password",
				JOptionPane.NO_OPTION, JOptionPane.PLAIN_MESSAGE,
				null, options, options[1]);
		if(option == 0){ // pressing OK button
			System.out.println(pass.getPassword());
			return charToByteArray(pass.getPassword());

		}else if(JOptionPane.showConfirmDialog(null, "Are you sure?\nThis will close the program.", "Cancellation", 
				JOptionPane.YES_NO_OPTION, JOptionPane.PLAIN_MESSAGE, null)
				== JOptionPane.YES_OPTION){
			System.out.println("User Exit.");
			System.exit(0);
		}
		return charToByteArray(pass.getPassword());
	}
	
	protected Integer doInBackground() throws Exception {  //main part of CWorker
		updateCConsole("Initialized\n"); //to inform user
		
		runmain(argCount, argToC, sizeX, sizeY);  //call to C code
		
		return argCount;
	}
}
