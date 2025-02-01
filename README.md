# IOFT User guide
Inclination-only fold test

1. Copy the ***ioft.m*** file to your MATLAB working directory.

2. Import your data into an Excel file following the format of the example file (***test.xlsx***) and save it in ***.xlsx*** format. Ensure the column headers match the case sensitivity of the example. Then, copy the data file into your MATLAB working directory.

3. In the MATLAB command window, enter

   ```
   [Tilt_Correction,Stepwise_Unfold] = ioft(filename)
   ```

   and press Enter to run the function. Here, filename is the name of the data spreadsheet in string format, excluding the file extension. For example, if the data file is named ***data.xlsx***, you should enter 

   ```
   [Tilt_Correction,Stepwise_Unfold] = ioft(‘data’)
   ```

4. Upon successful execution, the program will generate three figures and two tables as outputs. Figure 1 displays the equal-area projection before and after tilt correction. Figure 2 shows the K values obtained from the stepwise unfolding of the entire dataset. Figure 3 presents the results of stepwise unfolding for each limb of the fold separately. In the MATLAB command window, two tables, **Tilt_Correction**, and **Stepwise_Unfold** will be output. Both tables will be saved in CSV format in the MATLAB working directory. The three figures can be saved manually by selecting **File** > **Export** in the figure window and choosing the desired format and location. 

5. Additionally, you can view the built-in documentation for the program by entering 

   ```
   help ioft
   ```

    in the MATLAB command window. This documentation provides detailed explanations of the function, its input parameters, output formats, and usage instructions.
