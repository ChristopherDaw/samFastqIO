package chrisdaw.samFastqIO;

import java.lang.reflect.InvocationTargetException;

import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;

public class Main {
	
	static{
		System.loadLibrary("SamFastqIO");
	}
	
	public static void main(String[] args) throws InvocationTargetException, InterruptedException {
		//Dislay
		javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
			public void run() {
				Display.createAndShowGUI();
			}
		});

		try {
			for (LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) { //set look and feel
				if ("Windows".equals(info.getName())) {
					UIManager.setLookAndFeel(info.getClassName());
					break;
				} else if ("Mac OS X".equals(info.getName())) {
					UIManager.setLookAndFeel(info.getClassName());
					break;
				}
			}
		} catch (Exception e) { //not likely
			System.out.println("THERE WAS A PROBLEM ask Chris for help");
		}
	}
}
