

# Class InterpolationEvaluator2DCollection

**template &lt;ddc::concepts::type\_seq PatchSeq, ddc::concepts::type\_seq EvaluatorSeq, class PatchLocator, class [**ExtrapolationRule**](namespaceExtrapolationRule.md)&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationEvaluator2DCollection**](classInterpolationEvaluator2DCollection.md)



_A class to evaluate all the interpolations of all the patches at once._ [More...](#detailed-description)


































































## Detailed Description


This class computes the evaluations of the interpolations defined on each patch at given coordinates, using the interpolation evaluators given for each patch. The coordinates are stored in fields on each patch. On a given storing patch, we do not enforce that a coordinate is physically located on this storing patch. These fields of coordinates can be the result of characteristics equations solving.


This function is useful to avoid calling all the evaluators individually, especially in a multipatch geometry with several patches. It also manages to evaluate the right interpolation at a coordinate physically located on another patch than the storing patch.


Additionally, methods to compute the first derivatives and cross-derivatives are implemented.


Any 2D interpolation evaluator (see concepts::InterpolationEvaluator) can be used. Different evaluator types can be used on different patches.




**Warning:**

This operator does not work on batched domain.




**Warning:**

This operator only works for patches defined on the same continuous dimensions!




**Template parameters:**


* `PatchSeq` A ddc::TypeSeq of the patches. 
* `EvaluatorSeq` A ddc::TypeSeq of the 2D interpolation evaluators on each patch, in the same order as the patches in PatchSeq. 
* `PatchLocator` A operator that finds the patch where a given coordinate is physically located. 
* [**ExtrapolationRule**](namespaceExtrapolationRule.md) The extrapolation rule type for outside of the global domain. 




    

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interpolation/interpolation_evaluator_2d_collection.hpp`

