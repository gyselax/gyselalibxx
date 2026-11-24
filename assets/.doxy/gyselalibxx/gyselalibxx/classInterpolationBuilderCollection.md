

# Class InterpolationBuilderCollection

**template &lt;ddc::concepts::type\_seq PatchSeq, ddc::concepts::type\_seq BuilderSeq&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationBuilderCollection**](classInterpolationBuilderCollection.md)



_A class to call all the builders of all the patches once._ [More...](#detailed-description)


































































## Detailed Description


We need to instantiate all the builders for all the patches in the main code. We process the same way for the Fields containing the interpolation coefficients and the values of the function on each patch. The Fields are stored in MultipatchField objects. The builders are stored in this class. This class is instantiated with all the builders. The operator() allows to call all the builders stored in the member of this class in one single line.


This function is useful to avoid calling all the builders individually, especially in a multipatch geometry with several patches.


Any 1D interpolation builder (see concepts::InterpolationBuilder1D) can be used. Different builder types can be used on different patches (e.g. to use different boundary conditions).




**Template parameters:**


* `PatchSeq` A ddc::TypeSeq of the patches. 
* `BuilderSeq` A ddc::TypeSeq of the 1D interpolation builders on each patch, in the same order as the patches in PatchSeq. 




    

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interpolation/interpolation_builder_collection.hpp`

