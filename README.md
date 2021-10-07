# ERTProtocolsCreator
In this repository, you will find different Python scripts to generate ERT protocols for given geometries or process data originating from the ABME Terrameter (LS1 and LS2).

# Create 2D arrays:
They are two functions to create 2D arrays:
- ```CreateDNN(nMax, nElec, parralelizedMin) -> np.ndarray```
- ```CreateGradient(s, nElec) -> np.ndarray```

Both can be used with their default parameters creating a simple 64 electrodes profile.

# Create 3D arrays:
There, we base ourself on the creation of 2D protocols (gradient and/or dipole). 2 functions perform the basics 3D protocols (not for roll-alongs):
- ```CreateMultiLineGradient(nLines, nElecLine, s) -> np.ndarray```
- ```CreateMultiLineDDN(nLines, nElecLine, nMax) -> np.ndarray```

You can build anything you want from this (see ```CreateTom()```).

# Utilities:
All created protocols can be saved:
- ```SaveArrayTXT(filename, array)```
- ```ArrayToXML(array)```

You can also sort an array using the function ```SortArray(array, type)```

Creating the reciprocals array is done using:
```python
CreateReciprocals(array) # Full reciprocals
SampleReciprocals(array, sampling) # Sampled injections for reciprocals
```

# ExtractLines:
Extract in-line measurements from a 3D array that has been measured using the ABEM Terrameter LS1 or LS2.

Change the parameters atop the file and run the code to process the files.

