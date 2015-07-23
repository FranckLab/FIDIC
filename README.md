The Fast Iterative Digital Image Correlation Algorithm (FIDIC) is a 2D version of FIDVC algorithm ((please see [Bar-Kochba, Toyjanova et al., Exp. Mechanics, 2014](http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst) for more details)) to find dispalcements fields in a 2D image. 

* [Download latest version v1.0!](https://github.com/FranckLab/FIDIC/releases)
* [FAQ](https://github.com/FranckLab/FIDIC/blob/master/README.md#faq)
* [Questions/Issues](https://github.com/FranckLab/FIDIC/issues)
* [Bug Fixes/history](https://github.com/FranckLab/FIDIC/wiki/Bug-Fixes!)
* [Franck Lab](http://franck.engin.brown.edu)
 
## Purpose
The following implementation contains the Matlab m-files for our FIDIC algorithm along with example images. The FIDIC algorithm determines the 2D displacement fields between consecutive images. 

## Running FIDIC

### Sofware Requirement
Matlab 2011b and the associated Image Processing Toolbox are the minimum requirement to run this code.  

### Input Image Requirements
* To check if the images have the required speckle pattern and intensity values for correlation please use our [DIC simulator](https://github.com/FranckLab/DIC-Simulator).
* We recommend that the input image stack at each dimension should have at least 3 times of the subset size as the number of pixels. The default subset size is 64x64, so we recommend that the minimum input image size should be 192x192.
* The size of the input image stack should be divisible by 0.5 times the size of the subset. 

### Running including example case
1. Make sure that the main files and the supplemental m files (from file exchange) are in the working directory on Matlab. 
2. Run the `exampleRunFile.m` file to get 2D displacement fields between the image. Note that the displacement output is in the form of a generic uniaxial tension test. 


## Files
* Main files
 - addDisplacements_2D.m
 - checkConvergenceSSD_2D.m
 - DIC.m
 - filterDisplacements_2D.m
 - funIDIC.m
 - IDIC.m
 - removeOutliers_2D.m
 - areaMapping_2D.m

* Supplement m files from the MATLAB file exchange:
 - inpaint_nans.m
 - mirt2D_mexinterp.m  (Optional, not currently in use)

* Example Run files
 - exampleRunFile.m
 - tiff2mat.m
 - Example test images

## FAQ
**What are the requirements for the input images?**

Please refer to [input image requirement](https://github.com/FranckLab/FIDIC#input-image-requirements).

**Can I use FIDIC for finding displacement fields in 3D images?**

No. But you can use [FIDVC](https://github.com/FranckLab/FIDVC), this finds 3D displacements in 3D image stack.


## Cite
If used please cite:
[Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast iterative digital volume correlation algorithm for large deformations. Experimental Mechanics. doi: 10.1007/s11340-014-9874-2](http://link.springer.com/article/10.1007/s11340-014-9874-2?sa_campaign=email/event/articleAuthor/onlineFirst)

```bibtex
@article{bar2014fast,
  title={A fast iterative digital volume correlation algorithm for large deformations},
  author={Bar-Kochba, E and Toyjanova, J and Andrews, E and Kim, K-S and Franck, C},
  journal={Experimental Mechanics},
  pages={1--14},
  year={2014},
  publisher={Springer}
}
```

## Contact and support
For questions, please first refer to [FAQ](https://github.com/FranckLab/FIDIC#faq) and [Questions/Issues](https://github.com/FranckLab/FIDIC/issues). Add a new question if similar issue hasn't been reported. We shall help you at the earliest. The author's contact information can be found at [Franck Lab](http://franck.engin.brown.edu).
