Sam Fastq IO JNI Program Documentation

by Chris Daw
Overview:
This program is a Graphical User Interface (GUI) designed in Java that communicates with the native C code.  The Java program is nothing other than a helpful tool for anybody who is not familiar with the terminal.  This program is meant to make compressing/decompressing with the C algorithm easier to do.
I used Java Native Interface (JNI) to interface between the Java and C codes.  Since this interface is built into Java, this program can run on whatever computer has Java (and the OS that the C code was compiled in). 
Thank you for contributing to the code in whichever way you will do, are doing, or have done!


Java Classes and Functions:
1.  Main
   1. main(String[] args)
      1. Initializes the GUI and sets the look and feel of it.
   1. loadLibrary call (System function)
      1. Loads the dynamic library (.dylib) that tells Java about the C code.
1. Display
   1. addComponentsToFrame(Container frame)
      1. Adds every JPanel to the main JFrame as well as laying them out.
   1. addFileSelection()
      1. Adds and lays out all components necessary to file selection to fileSelectPanel.
   1. addOptions()
      1. Adds and lays out all components necessary to the options to its optionsPanel.
   1. addGo()
      1. Adds the go button to the goPanel.
   1. actionPerformed(ActionEvent e)
      1. Whenever a button is pushed (an action), this method is called through the ActionListener interface. 
      2. This function defines what all the buttons do.
   1. buildForPush()
      1. returns: char[][]
      2. Resets the variable varCount back to 0 for multiple runs in one program loop.
      3. calls getSettings() [2.g.] for the char[][] return value
   1. getSettings()
      1. returns: char[][]
      2. Creates a String using a StringBuilder and appending all the options selected as the time of the call of getSettings.
      3. Creates a char[][] using a delimiter of spaces to do so.
      4. Removes whitespace characters (‘\0’) from the original char(25)(10)
   1. createAndShowGUI()
      1. Creates, lays out, positions, and shows the GUI
   1. itemStateChanged(ItemEvent e)
      1. Same as actionPerformed() [2.e.] but for checkboxes and radio buttons.
      2. This function defines what to do when the different option checkboxes are changed (e.g. -v, -c, -u).
1. <Type>Filter
   1. File filters to restrict the user from choosing incorrect files.  Nothing other than that.
1. CWorker
   1. runmain(int varCount, byte[][] var, int sizeX, int sizeY)
      1. Native function that inits the C algorithm.
      2. varCount: number of arguments passed to C
      3. var: the list of arguments (argv) passed to C
      4. sizeX: the size of the byte[][] array
      5. sizeY: the longest byte[] element in the byte[][] array
   1. CWorker(int argc, char[][] argv, CConsole cc)
      1. Inits the CWorker object with the given arguments.
      2. argc: number of arguments in argv
      3. argv: lists of arguments to pass to C
      4. cc: the output console to be tied to this CWorker
   1. getSizeX(char[][] arr)
      1. returns: int
      2. Returns size of char array arr.
   1. getSizeY(char[][] arr)
      1. returns: int
      2. Returns length of longest char[] element in arr.
   1. charToByteArray(char[][] chars)
      1. returns: byte[][]
      2. Converts a 2d char array to a 2d byte array.
         1. This is because bytes are one bit (they convert to the char that C has) while java chars are two.
   1. charToByteArray(char[] chars)
      1. returns: byte[]
      2. Converts a char array to a byte array.
         1. See [4.e.II.1] for reasoning.
   1. updateCConsole(String[] strings)
      1. Sends strings to the linked output console.
      2. This method is solely called from the C code.
   1. updateCConsole(String string)
      1. Sends a string to the linked output console.
      2. This method is solely called from the Java code.
   1. isWorking()
      1. returns: boolean
      2. Returns state of CWorker (running = true and not = false)
   1. finishWorking()
      1. Sets boolean isWorking to false
      2. Uses updateCConsole(String string) to let the user know that it has finished working.
   1. askForPassword()
      1. returns: byte[]
      2. This method is used solely from C when ssh asks for a password.
      3. It opens a popup window that contains a password field (to hide their password from any sneaky onlookers).  The contents of this password field are then sent to the C code to enter the ssh server password.
   1. doInBackground()
      1. Required function of SwingWorker()s
      2. Lets user know (through CConsole) that the process has started)
      3. Runs runmain() with arguments compiled with the above methods.
1. CConsole
   1. CConsole()
      1. Inits CConsole object and runs startGUI() see [5.b.]
   1. startGUI()
      1. Lays out the entire GUI for the initial CConsole (before text)
   1. updateConsole(String update, boolean fromJava)
      1. Sets the text of the JTextArea to what it previously was plus the string update.
      2. The boolean fromJava is true when the text comes from Java, thus meaning it will be printed no matter what.
      3. The only time text wouldn’t be printed is if the option -v is not chosen from the options menu.
   1. actionPerformed(ActionEvent e)
      1. Similar to every other actionPerformed() function See [2.e.].
1. UsernameHostBar
   1. The UsernameHostBar object is for options that need to connect to an ssh server.  The object is made to avoid redundancies in the Display class (and because I’m lazy).
   2. UsernameHostBar()
      1. Inits the UsernameHostBar object and lays out the bars and text
   1. getAbsoluteFilePath()
      1. returns: String
      2. Concatenates everything into a pretty string.