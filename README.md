# CellFLEX-FM
Cell Fiber Lattice Experiment for Force Measurements

CellFLEX-FM is a Matlab code that can be used to extract 3D force measurements and 3D force maps from cell experiments on deformable photopolymerized fibers.
It was used in our publication that describes the calibration and testing of this new method for mechanobiology. 

You will find instructions on how to preprocess your image stacks and execute the code. You can also apply it on the provided test sample (in this case, skip directly to section B).

## A. Preprocessing your image stack (ImageJ)

1. Merge the 2 image stacks together (Cell and Fibers).

2. Rotate the image to correct the tilt.

3. Crop roughly around your cell.  
*NB: If you have an important drift, make sure that the cell and fibers are contained within the cropped area for all timepoints.  
NB: Do not include the black area from the image rotation in the crop.*

4. Rescale the image to a target x-y pixel size of 0.3 μm (Fiji >Image > Scale).  
Apply a scaling factor of (initial pixel size)/(final pixel size).   
*NB: Do not rescale in the z direction. Isometric voxel size will be obtained from the Matlab code.  
NB: It is possible to work with a different target pixel size. Rescaling is performed to reduce the computation time for the subsequent image analysis steps.*

5. Correct the drift using the ImageJ plugin Correct 3D drift (Fiji > Plugins > Registration > Correct 3D drift). Use the fiber channel as the reference.  
*NB: The plugin is functional on ImageJ 1.53c  
NB: The drift is corrected on the cropped resized image and not the original image. This is in order to considerably speed up the registration process.*

6. Crop finely the image stack so that it stops at the fibers ends. 

7. Split the 2 channels.

8. Subtract the background from the Fibers stack (Fiji > Process > Subtract background). For an x-y pixel size of 0.3 μm, the size of the structuring element can be set to around 10 pixels.

![/assets/pretreatment.png](https://github.com/uclapierre/CellFLEX-FM/blob/main/pretreatment.png)

9.	Save the Fibers stack as an image sequence in your working repertory (e.g. “Work”).

10.	Save the Cell stack as an image sequence in a sub-repertory (e.g. Work>cell).

## B.	Segmentation, Tracking, Intersection, Deflection (Matlab)

### 1.	In the "main", fill the different parameters of “Section 1: Initialization”:

**Z** = number of z stacks  
**xy_pix_size** = xy pixel size (e.g. 0.3e-6 m)  
**z_pix_size** = z pixel size (e.g. 1e-6 m)  
**L** = Fiber length (e.g. 80e-6 m)  
**Nb_nodes** = number of nodes in the Finite Element model (e.g. 41)  
**Param_fiber** = exposure time used for fabrication (e.g. 600 μs)  

Execute “Section 1: Initialization”.   
*NB: This step will load the 3D+T Fibers image stack and resample it to a final isometric voxel size. This step will also create folders in the “Work” repertory to host the outputs.*

### 2.	Execute “Section 2: Segmentation and tracking”.  
a.	Select a point on the image to define the threshold used for fiber segmentation.   
*NB: Choose the point on the fibers that has the weakest intensity, generally close to the ends. If you select an area with an intensity too close to the background signal, the segmentation will go wrong.*

b.	Select the areas that need to be erased if needed. Press Enter.  
*NB: It can for instance be used to erase a dust or, as in the example, a fiber with very weak intensity that will not be segmented properly.*

c.	Delineate the anchor points from the two extremities of the fibers from layer 1. Wait for processing. Delineate the anchor points from the two extremities of the fibers from layer 2. Wait for processing.  
*NB: The anchor points will be used to perform the subsequent tracking step, as they are a fixed part of the fibers.  They will also be used to handle fiber segmentation in the case where several fibers are contacting each other. The delineated area should typically be around 1/10 of the fiber length.*

d.	Wait for the segmentation step (in the tested example, the 245x258x21 x 61 time points stack is segmented in ~162 seconds).  
*NB: The segmented image is stored in Work>Seg and can be opened on Fiji to check the quality of the segmentation.*

e.	Wait for the tracking step (in the tested example, ~133 seconds).  
*NB: The segmented image is stored in Work>Track and can be opened on Fiji to check the quality of the segmentation.*
![/assets/segmentation.png](https://github.com/uclapierre/CellFLEX-FM/blob/main/segmentation.png)



f.	Wait for the computation of Lsmooth, the 3D+T tracked image where fibers are reconstructed as smooth cylinders (in the example, ~170 seconds).

### 3.	Execute “Section 3: Intersection and deflection”.  
a.	Select the repertory where the Cell stack image sequence is stored (e.g. Work>cell).

b.	If Mask_mode is set to 1, the binary mask of the cell will be generated by setting manually a threshold, by clicking on the cell image similar to step B.2.
Alternatively, (and preferably), if Mask_mode is set to 0, the binary mask of the cell needs to be generated in Fiji and saved as an image sequence prior to the execution of Section 3. In this case, open the Cell stack (from Work>cell) in Fiji. Apply some filtering to smooth the image (e.g. Fiji>Process>Filters>Median 3D with size 2). Apply a threshold manually (Fiji>Image>Adjust>Threshold). If needed, clear the borders of the binary images by applying a rectangular selection on the areas to be cleared (Fiji>Edit>Clear). Save the binary image as an image sequence in Work>cell mask.  
*NB: It does not matter if the thresholding leaves some small artefacts. The binary image will undergo morphological opening in Matlab which will remove them.*

![/assets/CellMask.png](https://github.com/uclapierre/CellFLEX-FM/blob/main/CellMask.png)
  

### 4.	Wait for the calculation of the Intersection (< 1 min).  

Wait for the calculation of the deflection vectors (here, for 12 fibers x 61 time points, 
~14 min). The deflection map will be stored in Work>Def3D  
![/assets/Deflection.png](https://github.com/uclapierre/CellFLEX-FM/blob/main/Deflection.png)
 



## C.	Traction force recovery (Matlab)

### 1.	Execute “Section 4: Traction force recovery”.
Wait for the computation of traction forces (in the tested example, ~30 min).
The traction force map will be stored in Work>Force3D.
The traction force vectors are stored in the variable pres (Nb_fibers x Nb_timepoints cell).    
![/assets/Force.png](https://github.com/uclapierre/CellFLEX-FM/blob/main/Force.png)




