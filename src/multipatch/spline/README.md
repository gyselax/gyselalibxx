# Spline on multipatch geometry

## Spline builder 

The MultipatchSplineBuilder allows all the given spline builders ot be called in one single line. 

### Example of use 
We define all the builders for each patch, 
```cpp
SplineBuilderType1 builder_1(domain_1); 
SplineBuilderType2 builder_2(domain_2); 
SplineBuilderType3 builder_3(domain_3); 
...
```

We can instantiate a MultipatchSplineBuilder from these builders, 
```cpp
MultipatchSplineBuilder builder (builder_1, builder_2, builder_3); 
```

We allocate memory for the spline representations on each patch,
```cpp
SplineType1 spline_1(spline_domain_1); 
SplineType2 spline_2(spline_domain_2); 
SplineType3 spline_3(spline_domain_3); 
...
std::tuple splines_tuple = {spline_1, spline_2, spline_3}; 
```

We do the same with the values of the function, 
```cpp
ValuesType1 values_1(domain_1); 
ValuesType2 values_2(domain_2); 
ValuesType3 values_3(domain_3); 
...
// initialisation 
...
std::tuple values_tuple = {values_1, values_2, values_3}; 
```

And we build the spline representation on every patch using, 
```cpp
builder(splines_tuple, values_tuple); 
```
