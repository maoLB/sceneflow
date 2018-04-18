####################################################################################
# This is free software; you can redistribute it and/or modify it under the		     #
# terms of the GNU General Public License as published by the Free Software        #
# Foundation; either version 2 of the License, or any later version.               #
#                                                                                  #
# This code is distributed in the hope that it will be useful, but WITHOUT ANY     #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  #
# PARTICULAR PURPOSE. See the GNU General Public License for more details.         #
#                                                                                  #
# You should have received a copy of the GNU General Public License along with     #
# this code; if not, write to the Free Software Foundation, Inc., 51 Franklin      #
# Street, Fifth Floor, Boston, MA 02110-1301, USA                                  #
####################################################################################

This is the implementation of 'Object Scene Flow for Autonomous Vehicles' as 
described in 

@INPROCEEDINGS{Menze2015CVPR,
  author = {Moritz Menze and Andreas Geiger},
  title = {Object Scene Flow for Autonomous Vehicles},
  booktitle = {Conference on Computer Vision and Pattern Recognition (CVPR)},
  year = {2015}
} 

The only difference to the paper is that we have replaced the proprietary SGM
implementation with the freely available spsstereo implementation (see references
below) which improves the numbers of our algorithm a little bit on the KITTI
scene flow test set. However, for consistency with the paper, the results submitted
to the benchmark are calculated with SGM as input.

Building the mex files
======================

1) run 'make_mex.m'
2) check ./build for compiled mex files

The code was developed and tested under Matlab 2013b and Ubuntu 14.04.
We do not provide any support for compilation issues or other OS.
Note that you need a processor that supports SSE 4.2 instructions (for
instance an Intel i7, otherwise the mex wrappers might crash Matlab due
to the missing popcnt implementation).


Running the Demo
================

To run our algorithm on the provided demo data (in ./data) simply execute demo.m.
Results will be stored in the newly created directory ./results


Processing the complete dataset
===============================

You can use the script run_osf.m to run our algorithm on the KITTI 2015 scene flow 
dataset: http://www.cvlibs.net/download.php?file=data_scene_flow.zip
Please copy the folders training and testing to ./data or adjust the paths in
run_osf.m accordingly.

There are two different parameter settings. You can toggle the comments in lines 14, 15
to execute the full method or the fast variant respectively.


Dependencies
============

To make the zip file self-contained, it includes copies of 

a) spsstereo: http://ttic.uchicago.edu/~dmcallester/SPS/ 
   as described in
   @InProceedings{Yamaguchi2014ECCV,
			Title = {Efficient joint segmentation, occlusion labeling, stereo and flow estimation},
			Author = {Yamaguchi, Koichiro and McAllester, David and Urtasun, Raquel},
			Booktitle = ECCV,
			Year = {2014}
	 }
   
b) libviso2: http://www.cvlibs.net/software/libviso/
   as described in
   @INPROCEEDINGS{Geiger2011IV,
			author = {Andreas Geiger and Julius Ziegler and Christoph Stiller},
			title = {StereoScan: Dense 3D Reconstruction in Real-time},
			booktitle = {Intelligent Vehicles Symposium (IV)},
			year = {2011}
	 } 
		
c) TRW-S: http://pub.ist.ac.at/~vnk/papers/TRW-S.html
	 as described in 
	 @Article{Kolmogorov2006PAMI,
			Title                    = {Convergent Tree-Reweighted Message Passing for Energy Minimization},
			Author                   = {Kolmogorov, Vladimir},
			Journal                  = PAMI,
			Year                     = {2006},
			Number                   = {10},
			Pages                    = {1568-1583},
			Volume                   = {28}
	 }
	 
d) Other dependencies (installable via the package manager):
	 png++:	 http://savannah.nongnu.org/projects/pngpp/
   libpng: http://www.libpng.org/pub/png/libpng.html
	 zlib: 	 http://www.zlib.net/

Citation
========

If you find this project, data or code useful, we would be happy if you cite us:

@INPROCEEDINGS{Menze2015CVPR,
  author = {Moritz Menze and Andreas Geiger},
  title = {Object Scene Flow for Autonomous Vehicles},
  booktitle = {Conference on Computer Vision and Pattern Recognition (CVPR)},
  year = {2015}
}

