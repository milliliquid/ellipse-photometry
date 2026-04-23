# ellipse-photometry
A python tool that takes a set of inputs and performs photometry using a fitted elliptical aperture and returns the data.

The inputs can be inputted through the command line or through input.txt. Modify the example parameters and set Active to 1 to use input.txt. This will stop the code from asking you for the parameters and will use the values from the input file. Data will be outputted in three different files. _data is for data, _params is for the ellipse parameters and _bg is for the background count data. They will use the same name as your FITS file input. This only works with a FITS file.
