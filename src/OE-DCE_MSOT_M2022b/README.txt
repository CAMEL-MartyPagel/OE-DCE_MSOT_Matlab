+-----------------------------+
| msotlib.itheramedical.com   |
+-----------------------------+



MATLAB Interfacing Library to the ViewMSOT Proprietary File Format





1. INSTALL
--------------------------------------------------------------------------------------

In order to use the functions within the “msotlib” library, please follow those two steps: 

1.	Add the following lines to your MATLAB startup.m * file (please create this file if it does not exist yet): 

javaaddpath <DIRECTORY>\MSOTBeans\xbean.jar
javaaddpath <DIRECTORY>\MSOTBeans\msotbeans.jar 

Within those command lines, replace <DIRECTORY> with the path to the location where you saved the “msotlib” library.

* Your startup.m file should be located in your MATLAB startup folder. 
In order to identify your MATLAB startup folder, you can type pwd in the MATLAB command window immediately after starting MATLAB and before typing any other command. 

2.	Add the “msotlib” library to your MATLAB search path using the dialog box “Set Path”.






2. CONTENTS
--------------------------------------------------------------------------------------

Listing of functions:
- listMSOT	    Lists contents of a study folder
- loadMSOT          Loads MSOT META information
- loadMSOTRecon     Loads MSOT Reconstructions
- loadMSOTMsp       Loads MSOT MSP (multispectrally processed) data
- loadMSOTSignals   Loads optoacoustic Signals as acquired by the transducers






3. USAGE
--------------------------------------------------------------------------------------

For usage please refer to the function documentation using help <function> and to our instruction manual