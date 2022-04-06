# HSTDownload
This code helps finding and downloading HST images for a given object or around given set of coordinates.

## License and Referencing
HSTDownload is released under a BSD 2-clause license. HSTDownload is freely available and you can freely modify, extend, and improve the HSTDownload source code. However, if you use it for published research, you are requested to cite del Pino et al. 2022 (ApJ submitted), where the method is described.

## Execution
Just download execute the script. Some examples are:

   ./hstdownload.py --name "Sculptor dSph"

   ./hstdownload.py --name "NGC 5053" --quiet

   ./hstdownload.py --ra 201.405 --dec -47.667 --search_radius 2 --hst_filters "F625W"

For more information:

   ./hstdownload.py --help
