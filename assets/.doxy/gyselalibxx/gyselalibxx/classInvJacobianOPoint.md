

# Class InvJacobianOPoint

**template &lt;class Mapping, class CoordRTheta&gt;**



[**ClassList**](annotated.md) **>** [**InvJacobianOPoint**](classInvJacobianOPoint.md)



_An operator for calculating the inverse of the Jacobian at an O-point. This class is used in_ [_**CombinedMapping**_](classCombinedMapping.md) _to calculate the inverse of the Jacobian at an O-point when one of the mappings does not allow the evaluation of its Jacobian/inverse Jacobian at the O-point._[More...](#detailed-description)


































































## Detailed Description


Specialisations of this class must implement:
* A constructor taking the Mapping as an argument
* An operator() returning the Jacobian matrix at the O-point.






**Template parameters:**


* `Mapping` The mapping for which the inverse of the Jacobian is calculated. 
* `CoordRTheta` The coordinate system in which the inverse of the Jacobian is calculated. 




    

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/inv_jacobian_o_point.hpp`

