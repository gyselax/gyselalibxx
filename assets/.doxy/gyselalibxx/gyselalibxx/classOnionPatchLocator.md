

# Class OnionPatchLocator

**template &lt;class MultipatchIdxRanges, class LogicalToPhysicalMapping, class PhysicalToLogicalMapping, class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**OnionPatchLocator**](classOnionPatchLocator.md)



[_**Patch**_](structPatch.md) _locator specialised for "onion" geometry._[More...](#detailed-description)


































































## Detailed Description


We define an "onion" geometry a set of patches mapping to the physical domain in a shape of concentrical rings. The first patch is supposed to be a disk containing the O-point. The other patches are ordered as concentrical rings drawing away from the O-point. The order of patches is made by MultipatchIdxRanges and it is important for the dichotomy method.


We also suppose that we can define a global logical grid that we can split into the different logical grids of the patches.


This operator locates on which patch a given physical coordinate is.




**Warning:**

The operator can works on GPU or CPU according to the given execution space ExecSpace. The ExecSpace by default is device. The constructor will still need to be called from CPU, but the operator() needs to be called from the given ExecSpace.


 

**Template parameters:**


* `MultipatchIdxRanges` A [**MultipatchType**](classMultipatchType.md) type containing the 2D index ranges on each patch. 
* `LogicalToPhysicalMapping` A mapping type for all the patches. 
* `ExecSpace` The space (CPU/GPU) where the calculations are carried out. By default it is on device. 




    

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/onion_patch_locator.hpp`

